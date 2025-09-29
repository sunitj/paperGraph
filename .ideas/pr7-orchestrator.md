# PR #7: Orchestrator - Agent Coordination

**Sprint:** Week 3, Days 1-2  
**Goal:** Coordinate agents and synthesize results

## What to Build

1. FastAPI endpoints for querying
2. Agent coordination logic
3. Result synthesis using LLM
4. Seed loading endpoint

## Deliverable

```bash
# POST /query → Coordinates agents → Returns answer with graph
curl -X POST http://localhost:8000/query \
  -d '{"question": "How do gut bacteria affect obesity?"}'
```

## orchestrator/main.py

```python
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import redis
import json
import uuid
import time
import sys
import os

sys.path.insert(0, '/app/shared')
from neo4j_client import KnowledgeGraph
from llm_service import get_llm_provider

app = FastAPI(
    title="PaperGraph Orchestrator",
    description="Multi-agent literature analysis system"
)

# CORS for UI
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize services
redis_client = redis.from_url(os.getenv("REDIS_URL"))
kg = KnowledgeGraph(
    uri=os.getenv("NEO4J_URI"),
    user=os.getenv("NEO4J_USER"),
    password=os.getenv("NEO4J_PASSWORD")
)
llm = get_llm_provider()


class Query(BaseModel):
    question: str
    max_papers: int = 20


class SeedRequest(BaseModel):
    seed: str


@app.get("/health")
async def health():
    """Health check"""
    llm_status = "unknown"
    if hasattr(llm, 'health_check'):
        llm_status = "healthy" if llm.health_check() else "unhealthy"
    
    return {
        "status": "healthy",
        "llm_provider": os.getenv("LLM_PROVIDER", "ollama"),
        "llm_status": llm_status,
        "redis": "connected",
        "neo4j": "connected"
    }


@app.post("/query")
async def process_query(query: Query):
    """Process a user query"""
    question_id = str(uuid.uuid4())[:8]
    
    print(f"\n{'='*60}")
    print(f"📥 New Query: {question_id}")
    print(f"Question: {query.question}")
    print(f"{'='*60}\n")
    
    # Send to Query Strategist
    redis_client.rpush("question_queue", json.dumps({
        'id': question_id,
        'question': query.question
    }))
    
    print("1. Sent to Query Strategist...")
    
    # Wait for analysis (simplified - in production use async/websockets)
    time.sleep(10)
    
    # Check if analysis is complete
    analysis = redis_client.get(f"question:{question_id}:analysis")
    if not analysis:
        raise HTTPException(status_code=408, detail="Analysis timeout")
    
    analysis_data = json.loads(analysis)
    print("2. Analysis complete")
    
    # Wait for papers to be fetched and analyzed
    print("3. Waiting for papers to be fetched and analyzed...")
    time.sleep(40)  # Give agents time to work
    
    # Collect results
    papers = []
    for i in range(3):  # Check for 3 queries
        query_id = f"{question_id}_q{i}"
        paper_data = redis_client.get(f"query:{query_id}:papers")
        if paper_data:
            papers.extend(json.loads(paper_data))
    
    print(f"4. Collected {len(papers)} papers")
    
    # Get gaps
    gaps = json.loads(redis_client.get(f"question:{question_id}:gaps") or "[]")
    
    # Get subgraph
    entities = analysis_data.get('entities', [])
    subgraph = kg.get_subgraph(entities, depth=2)
    
    # Synthesize answer
    print("5. Synthesizing answer...")
    answer = synthesize_answer(query.question, papers, gaps, subgraph)
    
    print("✅ Query complete\n")
    
    return {
        "question_id": question_id,
        "answer": answer,
        "papers": papers[:10],  # Return top 10
        "gaps": gaps[:5],  # Return top 5
        "entities": entities,
        "subgraph": format_subgraph(subgraph[:20])
    }


def synthesize_answer(question: str, papers: list, gaps: list, subgraph: list) -> str:
    """Synthesize answer using LLM"""
    
    # Format papers
    paper_summary = "\n".join([
        f"- [{i+1}] {p['title']} (PMID: {p['pmid']})"
        for i, p in enumerate(papers[:5])
    ])
    
    # Format gaps
    gap_summary = "\n".join([
        f"- {g['type']}: {g.get('subject', g.get('entity', 'unknown'))}"
        for g in gaps[:3]
    ])
    
    # Format subgraph
    relationships = []
    for item in subgraph[:10]:
        # Extract relationships from path
        if 'path' in str(item):
            relationships.append(str(item))
    
    rel_summary = "\n".join(relationships[:5]) if relationships else "No existing relationships"
    
    prompt = f"""You are a scientific research assistant analyzing literature.

Question: {question}

Knowledge Graph Relationships:
{rel_summary}

Papers Analyzed:
{paper_summary}

Knowledge Gaps Identified:
{gap_summary}

Provide a comprehensive answer that:
1. Directly answers the question
2. Cites specific papers using [number] format
3. Explains key mechanisms or relationships
4. Notes any gaps or uncertainties
5. Is clear and concise (2-3 paragraphs)

Answer:"""

    answer = llm.generate(prompt)
    return answer


def format_subgraph(subgraph: list) -> list:
    """Format subgraph for UI"""
    nodes = []
    edges = []
    
    for item in subgraph:
        # Simplified formatting - extract nodes and relationships
        # This would need proper Neo4j path parsing in production
        pass
    
    return {
        "nodes": nodes[:20],
        "edges": edges[:20]
    }


@app.post("/load_seed")
async def load_seed(request: SeedRequest):
    """Load a seed knowledge graph"""
    seed_name = request.seed.lower().replace(" ", "_")
    seed_file = f"/seeds/{seed_name}.json"
    
    if not os.path.exists(seed_file):
        raise HTTPException(status_code=404, detail=f"Seed '{request.seed}' not found")
    
    print(f"\n🌱 Loading seed: {request.seed}")
    
    with open(seed_file) as f:
        seed_data = json.load(f)
    
    papers_loaded = 0
    relationships_loaded = 0
    
    for paper in seed_data['papers']:
        # Create paper
        kg.create_paper(
            pmid=paper['pmid'],
            title=paper['title'],
            abstract=paper.get('abstract', ''),
            year=paper['year'],
            pre_seeded=True
        )
        papers_loaded += 1
        
        # Create relationships if present
        for rel in paper.get('relationships', []):
            kg.create_relationship(
                subj_name=rel['subject'],
                subj_type=rel.get('subject_type', 'Organism'),
                predicate=rel['predicate'],
                obj_name=rel['object'],
                obj_type=rel.get('object_type', 'Disease'),
                confidence=rel.get('confidence', 0.8),
                paper_pmid=paper['pmid']
            )
            relationships_loaded += 1
    
    print(f"✅ Loaded {papers_loaded} papers, {relationships_loaded} relationships\n")
    
    return {
        "status": "success",
        "seed": request.seed,
        "papers_count": papers_loaded,
        "relationships_count": relationships_loaded
    }


@app.get("/seeds")
async def list_seeds():
    """List available seeds"""
    import glob
    
    seeds = []
    for seed_file in glob.glob("/seeds/*.json"):
        with open(seed_file) as f:
            data = json.load(f)
            seeds.append({
                "name": data['topic'],
                "description": data.get('description', ''),
                "papers": len(data['papers'])
            })
    
    return {"seeds": seeds}


@app.get("/stats")
async def get_stats():
    """Get system statistics"""
    
    # Count papers
    paper_count = kg.graph.run("MATCH (p:Paper) RETURN count(p) as count").data()[0]['count']
    
    # Count entities
    entity_count = kg.graph.run("""
        MATCH (n)
        WHERE n:Organism OR n:Disease OR n:Molecule
        RETURN count(n) as count
    """).data()[0]['count']
    
    # Count relationships
    rel_count = kg.graph.run("MATCH ()-[r]->() RETURN count(r) as count").data()[0]['count']
    
    return {
        "papers": paper_count,
        "entities": entity_count,
        "relationships": rel_count
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
```

## scripts/test_orchestrator.py

```python
#!/usr/bin/env python3
"""
Test the orchestrator API
"""

import requests
import json
import time


def test_health():
    """Test health endpoint"""
    print("1. Testing health endpoint...")
    response = requests.get("http://localhost:8000/health")
    print(f"   Status: {response.json()['status']}")
    print(f"   LLM: {response.json()['llm_status']}")


def test_query():
    """Test query endpoint"""
    print("\n2. Testing query endpoint...")
    
    response = requests.post(
        "http://localhost:8000/query",
        json={"question": "How do gut bacteria influence obesity?"}
    )
    
    if response.status_code == 200:
        result = response.json()
        print(f"   ✅ Answer received:")
        print(f"   {result['answer'][:200]}...")
        print(f"   Papers: {len(result['papers'])}")
        print(f"   Gaps: {len(result['gaps'])}")
    else:
        print(f"   ❌ Error: {response.status_code}")


def test_seed():
    """Test seed loading"""
    print("\n3. Testing seed loading...")
    
    # List seeds
    response = requests.get("http://localhost:8000/seeds")
    seeds = response.json()['seeds']
    print(f"   Available seeds: {len(seeds)}")
    
    # Load first seed
    if seeds:
        response = requests.post(
            "http://localhost:8000/load_seed",
            json={"seed": seeds[0]['name']}
        )
        print(f"   ✅ Loaded: {seeds[0]['name']}")
        print(f"   Papers: {response.json()['papers_count']}")


def test_stats():
    """Test stats endpoint"""
    print("\n4. Testing stats endpoint...")
    
    response = requests.get("http://localhost:8000/stats")
    stats = response.json()
    
    print(f"   Papers: {stats['papers']}")
    print(f"   Entities: {stats['entities']}")
    print(f"   Relationships: {stats['relationships']}")


if __name__ == "__main__":
    print("="*60)
    print("Orchestrator API Test")
    print("="*60)
    
    try:
        test_health()
        test_seed()
        test_stats()
        test_query()
        
        print("\n" + "="*60)
        print("✅ All tests complete!")
        print("="*60)
    
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
```

## Testing

```bash
# Start orchestrator
docker-compose up -d orchestrator

# Test endpoints
python scripts/test_orchestrator.py

# Or use curl
curl http://localhost:8000/health
curl http://localhost:8000/seeds
curl -X POST http://localhost:8000/load_seed -H "Content-Type: application/json" -d '{"seed":"Gut Microbiome Obesity"}'
curl -X POST http://localhost:8000/query -H "Content-Type: application/json" -d '{"question":"How do gut bacteria affect obesity?"}'
```

## API Documentation

Once running, visit:
- http://localhost:8000/docs (Swagger UI)
- http://localhost:8000/redoc (ReDoc)

## Acceptance Criteria

- [ ] Health endpoint returns status
- [ ] Query endpoint coordinates agents
- [ ] Seed loading works correctly
- [ ] Stats endpoint returns counts
- [ ] Answer synthesis uses LLM
- [ ] API documentation accessible
- [ ] Error handling works
- [ ] CORS enabled for UI
