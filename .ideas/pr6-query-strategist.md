# PR #6: Query Strategist Agent - Gap Analysis

**Sprint:** Week 2, Days 6-7  
**Goal:** Intelligent query generation based on KG gaps

## What to Build

1. KG gap analysis (weak evidence, missing links)
2. Entity extraction from user questions
3. Smart query generation
4. Integration with Fetcher Agent

## Deliverable

```bash
# User question → Analyzes KG → Generates targeted queries
# Queries sent to Fetcher automatically
```

## agents/query_strategist/agent.py

```python
#!/usr/bin/env python3
"""
Query Strategist Agent - Analyzes gaps and generates smart queries
"""

import sys
import os
import time
import json
import re
import redis

sys.path.insert(0, '/app/shared')
from neo4j_client import KnowledgeGraph
from llm_service import get_llm_provider


class QueryStrategistAgent:
    """Agent that analyzes knowledge gaps and generates smart queries"""
    
    def __init__(self):
        self.redis = redis.from_url(os.getenv("REDIS_URL"))
        self.kg = KnowledgeGraph(
            uri=os.getenv("NEO4J_URI"),
            user=os.getenv("NEO4J_USER"),
            password=os.getenv("NEO4J_PASSWORD")
        )
        self.llm = get_llm_provider()
    
    def extract_entities_from_question(self, question: str) -> list:
        """Extract key entities from user question"""
        prompt = f"""Extract key biological/scientific entities from this question.

Question: {question}

Return ONLY a JSON array of entity names:
["entity1", "entity2", "entity3"]

Focus on: organisms, diseases, molecules, biological processes.
Maximum 5 entities."""

        try:
            response = self.llm.generate(prompt)
            
            # Extract JSON array
            json_match = re.search(r'\[.*?\]', response, re.DOTALL)
            if json_match:
                entities = json.loads(json_match.group())
                return entities[:5]  # Max 5
        
        except Exception as e:
            print(f"   ⚠️  Entity extraction error: {e}")
        
        # Fallback: extract capitalized words
        entities = re.findall(r'\b[A-Z][a-z]+(?:\s+[a-z]+)?\b', question)
        return entities[:5]
    
    def analyze_kg_coverage(self, entities: list) -> tuple:
        """Analyze existing knowledge in graph"""
        if not entities:
            return [], []
        
        print(f"   🔍 Checking KG for: {', '.join(entities)}")
        
        # Query for existing relationships
        coverage = self.kg.graph.run("""
            MATCH (start)-[r]->(end)
            WHERE start.name IN $entities OR end.name IN $entities
            RETURN start.name as subject,
                   type(r) as predicate,
                   end.name as object,
                   avg(r.confidence) as avg_confidence,
                   count(*) as paper_count
            ORDER BY paper_count DESC
        """, entities=entities).data()
        
        # Identify gaps
        gaps = []
        
        for row in coverage:
            # Weak evidence
            if row['avg_confidence'] < 0.5:
                gaps.append({
                    'type': 'weak_evidence',
                    'subject': row['subject'],
                    'predicate': row['predicate'],
                    'object': row['object'],
                    'confidence': round(row['avg_confidence'], 2),
                    'priority': 'high'
                })
            
            # Understudied relationships
            if row['paper_count'] < 3:
                gaps.append({
                    'type': 'understudied',
                    'subject': row['subject'],
                    'predicate': row['predicate'],
                    'object': row['object'],
                    'paper_count': row['paper_count'],
                    'priority': 'medium'
                })
        
        # Check for completely missing entities
        found_entities = set()
        for row in coverage:
            found_entities.add(row['subject'])
            found_entities.add(row['object'])
        
        missing_entities = set(entities) - found_entities
        
        for entity in missing_entities:
            gaps.append({
                'type': 'missing_entity',
                'entity': entity,
                'priority': 'high'
            })
        
        print(f"      Found {len(coverage)} existing relationships")
        print(f"      Identified {len(gaps)} gaps")
        
        return coverage, gaps
    
    def generate_queries(self, question: str, coverage: list, gaps: list) -> list:
        """Generate targeted search queries to fill gaps"""
        
        # Format coverage for prompt
        coverage_summary = "\n".join([
            f"- {c['subject']} {c['predicate']} {c['object']} "
            f"(confidence: {c['avg_confidence']:.2f}, papers: {c['paper_count']})"
            for c in coverage[:5]
        ])
        
        # Format gaps for prompt
        gaps_summary = "\n".join([
            f"- {g['type']}: {g.get('subject', g.get('entity', 'unknown'))} "
            f"({g.get('priority', 'unknown')} priority)"
            for g in gaps[:5]
        ])
        
        prompt = f"""You are a scientific literature search expert.

User question: {question}

Current knowledge graph coverage:
{coverage_summary if coverage_summary else 'No existing knowledge'}

Knowledge gaps identified:
{gaps_summary if gaps_summary else 'No gaps identified'}

Generate 3 targeted PubMed search queries to answer the question and fill gaps.

Rules:
1. Make queries specific and scientific
2. Use Boolean operators (AND, OR) where helpful
3. Target the biggest gaps first
4. Each query should have a clear purpose

Return ONLY JSON:
[
  {{"query": "specific search terms", "reason": "fills X gap", "priority": "high"}},
  {{"query": "specific search terms", "reason": "fills Y gap", "priority": "medium"}},
  {{"query": "specific search terms", "reason": "verifies Z", "priority": "low"}}
]"""

        try:
            response = self.llm.generate(prompt)
            
            # Extract JSON
            json_match = re.search(r'\[.*\]', response, re.DOTALL)
            if json_match:
                queries = json.loads(json_match.group())
                return queries[:3]  # Max 3 queries
        
        except Exception as e:
            print(f"   ⚠️  Query generation error: {e}")
        
        # Fallback: create basic query from question
        return [{
            'query': question,
            'reason': 'direct search',
            'priority': 'medium'
        }]
    
    def send_to_fetcher(self, queries: list, question_id: str):
        """Send queries to Fetcher Agent"""
        print(f"   📤 Sending {len(queries)} queries to Fetcher:")
        
        for i, query_data in enumerate(queries):
            query_id = f"{question_id}_q{i}"
            
            fetch_payload = {
                'id': query_id,
                'query': query_data['query'],
                'question_id': question_id,
                'max_results': 20
            }
            
            self.redis.rpush("fetch_queue", json.dumps(fetch_payload))
            
            print(f"      {i+1}. {query_data['query']}")
            print(f"         Reason: {query_data['reason']}")
    
    def process_question(self, question_id: str, question: str):
        """Process a user question"""
        print(f"\n❓ Processing question {question_id}")
        print(f"   Question: {question}")
        
        # Extract entities
        print(f"   🔎 Extracting entities...")
        entities = self.extract_entities_from_question(question)
        print(f"      Entities: {', '.join(entities)}")
        
        # Analyze KG
        print(f"   📊 Analyzing knowledge graph...")
        coverage, gaps = self.analyze_kg_coverage(entities)
        
        # Generate queries
        print(f"   💡 Generating search queries...")
        queries = self.generate_queries(question, coverage, gaps)
        
        # Send to fetcher
        self.send_to_fetcher(queries, question_id)
        
        # Store analysis results
        analysis = {
            'question': question,
            'entities': entities,
            'coverage_count': len(coverage),
            'gaps_count': len(gaps),
            'gaps': gaps[:10],  # Store first 10 gaps
            'queries': queries
        }
        
        self.redis.set(
            f"question:{question_id}:analysis",
            json.dumps(analysis),
            ex=3600  # 1 hour expiry
        )
        
        # Store gaps separately for easy access
        self.redis.set(
            f"question:{question_id}:gaps",
            json.dumps(gaps),
            ex=3600
        )
        
        print(f"   ✅ Question {question_id} processed\n")
    
    def run(self):
        """Main agent loop"""
        print("="*60)
        print("🧭 Query Strategist Agent Started")
        print("="*60)
        print("Analyzing knowledge gaps and generating smart queries...")
        print("Waiting for questions...\n")
        
        while True:
            try:
                # Block and wait for questions
                result = self.redis.blpop("question_queue", timeout=5)
                
                if result:
                    _, question_json = result
                    question_data = json.loads(question_json)
                    
                    self.process_question(
                        question_data['id'],
                        question_data['question']
                    )
            
            except KeyboardInterrupt:
                print("\n\n👋 Shutting down gracefully...")
                break
            
            except Exception as e:
                print(f"❌ Error in main loop: {e}")
                time.sleep(5)


if __name__ == "__main__":
    agent = QueryStrategistAgent()
    agent.run()
```

## scripts/test_query_strategist.py

```python
#!/usr/bin/env python3
"""
Test script for Query Strategist Agent
"""

import redis
import json
import time
import sys


def test_strategist():
    """Test the Query Strategist agent"""
    r = redis.from_url("redis://localhost:6379")
    
    print("Testing Query Strategist Agent...")
    
    # Create test question
    question = {
        "id": "test_question_1",
        "question": "How do gut bacteria influence obesity?"
    }
    
    print(f"\n1. Sending test question: {question['question']}")
    r.rpush("question_queue", json.dumps(question))
    print("   ✅ Question queued")
    
    # Wait for processing
    print("\n2. Waiting for strategist to analyze (20 seconds)...")
    time.sleep(20)
    
    # Check analysis
    print("\n3. Checking analysis results...")
    
    analysis = r.get(f"question:{question['id']}:analysis")
    if analysis:
        analysis_data = json.loads(analysis)
        
        print(f"   ✅ Analysis complete:")
        print(f"      Entities: {', '.join(analysis_data['entities'])}")
        print(f"      Coverage: {analysis_data['coverage_count']} relationships")
        print(f"      Gaps: {analysis_data['gaps_count']} identified")
        
        print(f"\n   Generated queries:")
        for q in analysis_data['queries']:
            print(f"      - {q['query']}")
            print(f"        ({q['reason']})")
    else:
        print("   ⚠️  No analysis found")
    
    # Check gaps
    print("\n4. Checking identified gaps...")
    gaps = r.get(f"question:{question['id']}:gaps")
    if gaps:
        gaps_data = json.loads(gaps)
        print(f"   Found {len(gaps_data)} gaps:")
        for gap in gaps_data[:3]:
            print(f"      - {gap['type']}: {gap}")
    
    # Check if queries were sent to fetcher
    print("\n5. Checking Fetcher queue...")
    queue_length = r.llen("fetch_queue")
    print(f"   Fetch queue: {queue_length} queries")


if __name__ == "__main__":
    print("="*60)
    print("Query Strategist Test")
    print("="*60)
    print("\nMake sure services are running:")
    print("  docker-compose up -d")
    print("\nPress Enter to start test...")
    input()
    
    try:
        test_strategist()
        print("\n" + "="*60)
        print("✅ Test complete!")
        print("="*60)
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)
```

## Manual Testing

```bash
# Start strategist
docker-compose up -d query_strategist

# Watch logs
docker-compose logs -f query_strategist

# Send question via Redis
redis-cli RPUSH question_queue '{
  "id": "manual_q1",
  "question": "What is the role of Akkermansia in metabolic health?"
}'

# Check analysis
redis-cli GET question:manual_q1:analysis | jq

# Check gaps
redis-cli GET question:manual_q1:gaps | jq

# Check Fetcher queue
redis-cli LRANGE fetch_queue 0 -1
```

## Integration Test

Test the full pipeline:

```bash
# 1. Load a seed (from PR #9)
curl -X POST http://localhost:8000/load_seed \
  -H "Content-Type: application/json" \
  -d '{"seed": "Gut Microbiome Obesity"}'

# 2. Send question
redis-cli RPUSH question_queue '{
  "id": "integration_test",
  "question": "How does diet affect gut bacteria and obesity?"
}'

# 3. Watch the pipeline
docker-compose logs -f query_strategist fetcher kg_builder

# 4. Check results after 60 seconds
redis-cli GET question:integration_test:analysis
```

## Acceptance Criteria

- [ ] Agent starts and connects to Redis/Neo4j
- [ ] Receives questions from Redis queue
- [ ] Extracts entities from questions
- [ ] Analyzes KG for coverage and gaps
- [ ] Generates targeted search queries
- [ ] Sends queries to Fetcher queue
- [ ] Stores analysis results in Redis
- [ ] Test script runs successfully
- [ ] Can handle questions with no KG coverage
- [ ] Can handle questions with full coverage
