# PR #5: KG Builder Agent - Entity Extraction

**Sprint:** Week 2, Days 4-5  
**Goal:** Extract entities and build knowledge graph from papers

## What to Build

1. spaCy-based entity extraction
2. LLM-based relationship extraction
3. Neo4j graph updates
4. Confidence scoring

## Deliverable

```bash
# Reads papers from Neo4j → Extracts entities → Updates graph
docker-compose logs -f kg_builder
# Shows: "Analyzed paper PMID123: 5 relationships"
```

## agents/kg_builder/agent.py

```python
#!/usr/bin/env python3
"""
KG Builder Agent - Extracts entities and builds knowledge graph
"""

import sys
import os
import time
import json
import re
import spacy
import redis

sys.path.insert(0, '/app/shared')
from neo4j_client import KnowledgeGraph
from llm_service import get_llm_provider

# Load spaCy model
nlp = spacy.load("en_core_web_sm")


class KGBuilderAgent:
    """Agent that builds knowledge graph from papers"""
    
    def __init__(self):
        self.redis = redis.from_url(os.getenv("REDIS_URL"))
        self.kg = KnowledgeGraph(
            uri=os.getenv("NEO4J_URI"),
            user=os.getenv("NEO4J_USER"),
            password=os.getenv("NEO4J_PASSWORD")
        )
        self.llm = get_llm_provider()
        
        # Biological entity keywords
        self.organism_keywords = [
            'bacteri', 'firmicutes', 'bacteroidetes', 'akkermansia',
            'coli', 'aureus', 'lactobacillus', 'bifidobacterium',
            'microbiome', 'microbiota', 'probiotic'
        ]
        
        self.disease_keywords = [
            'obesity', 'diabetes', 'cancer', 'disease', 'disorder',
            'syndrome', 'infection', 'inflammation', 'resistance'
        ]
    
    def extract_entities(self, text: str) -> dict:
        """Extract biological entities using spaCy + heuristics"""
        doc = nlp(text[:5000])  # Limit text length for performance
        
        organisms = []
        diseases = []
        molecules = []
        
        for ent in doc.ents:
            entity_text = ent.text
            entity_lower = entity_text.lower()
            
            # Classify entities
            if self._is_organism(entity_lower):
                organisms.append(entity_text)
            
            elif self._is_disease(entity_lower) or ent.label_ == "DISEASE":
                diseases.append(entity_text)
            
            elif ent.label_ in ["PRODUCT", "SUBSTANCE"]:
                molecules.append(entity_text)
        
        # Also extract biological terms from text patterns
        bio_patterns = self._extract_bio_patterns(text)
        organisms.extend(bio_patterns['organisms'])
        molecules.extend(bio_patterns['molecules'])
        
        return {
            'organisms': list(set(organisms))[:10],  # Limit to top 10
            'diseases': list(set(diseases))[:10],
            'molecules': list(set(molecules))[:10]
        }
    
    def _is_organism(self, text: str) -> bool:
        """Check if text is likely an organism name"""
        return any(kw in text for kw in self.organism_keywords)
    
    def _is_disease(self, text: str) -> bool:
        """Check if text is likely a disease name"""
        return any(kw in text for kw in self.disease_keywords)
    
    def _extract_bio_patterns(self, text: str) -> dict:
        """Extract biological entities using regex patterns"""
        organisms = []
        molecules = []
        
        # Pattern: Genus species (e.g., "Escherichia coli")
        genus_species = re.findall(r'\b[A-Z][a-z]+ [a-z]+\b', text)
        for match in genus_species:
            if any(kw in match.lower() for kw in ['bacteri', 'coccus', 'ella']):
                organisms.append(match)
        
        # Pattern: Chemical compounds
        chemical_patterns = [
            r'\b[A-Z]+-\d+\b',  # e.g., PD-1, IL-6
            r'\b[a-z]+-[a-z]+\b',  # e.g., short-chain fatty acids
        ]
        
        for pattern in chemical_patterns:
            matches = re.findall(pattern, text)
            molecules.extend(matches)
        
        return {
            'organisms': organisms[:5],
            'molecules': molecules[:5]
        }
    
    def extract_relationships(self, text: str, entities: dict) -> list:
        """Extract relationships using LLM"""
        all_entities = (
            entities['organisms'] + 
            entities['diseases'] + 
            entities['molecules']
        )
        
        if len(all_entities) < 2:
            return []
        
        # Limit text length for LLM
        text_snippet = text[:1000]
        
        prompt = f"""Extract biological relationships from this scientific text.

Text: {text_snippet}

Entities found: {', '.join(all_entities[:10])}

Output a JSON array of relationships. Use ONLY these predicates:
- PRODUCES (A creates/synthesizes B)
- AFFECTS (A influences/modulates B)
- INHIBITS (A blocks/prevents B)
- INCREASES (A upregulates B)
- DECREASES (A downregulates B)

Format:
[
  {{"subject": "entity1", "predicate": "PRODUCES", "object": "entity2", "confidence": 0.8}},
  ...
]

Rules:
1. Only include relationships explicitly stated in the text
2. Confidence should reflect how clearly stated (0.3-1.0)
3. Return ONLY the JSON array, nothing else
4. Maximum 5 relationships"""

        try:
            response = self.llm.generate(prompt)
            
            # Extract JSON from response
            json_match = re.search(r'\[.*\]', response, re.DOTALL)
            if json_match:
                relationships = json.loads(json_match.group())
                
                # Validate relationships
                valid_relationships = []
                for rel in relationships[:5]:  # Max 5
                    if all(k in rel for k in ['subject', 'predicate', 'object', 'confidence']):
                        valid_relationships.append(rel)
                
                return valid_relationships
        
        except Exception as e:
            print(f"      ⚠️  Relationship extraction error: {e}")
        
        return []
    
    def get_entity_type(self, entity: str, entities: dict) -> str:
        """Determine entity type"""
        entity_lower = entity.lower()
        
        if entity in entities['organisms']:
            return "Organism"
        elif entity in entities['diseases']:
            return "Disease"
        elif entity in entities['molecules']:
            return "Molecule"
        elif self._is_organism(entity_lower):
            return "Organism"
        elif self._is_disease(entity_lower):
            return "Disease"
        else:
            return "Molecule"  # Default
    
    def process_paper(self, pmid: str):
        """Process a single paper"""
        # Get paper from Neo4j
        result = self.kg.graph.run("""
            MATCH (p:Paper {pmid: $pmid})
            RETURN p.title as title, p.abstract as abstract
        """, pmid=pmid).data()
        
        if not result:
            print(f"   ⚠️  Paper {pmid} not found in Neo4j")
            return
        
        paper = result[0]
        
        if not paper['abstract'] or len(paper['abstract']) < 50:
            print(f"   ⚠️  Paper {pmid} has no abstract")
            self.kg.mark_paper_analyzed(pmid)
            return
        
        print(f"\n📄 Analyzing paper {pmid}")
        print(f"   Title: {paper['title'][:60]}...")
        
        # Combine title and abstract
        text = f"{paper['title']}. {paper['abstract']}"
        
        # Extract entities
        print(f"   🔍 Extracting entities...")
        entities = self.extract_entities(text)
        
        total_entities = sum(len(v) for v in entities.values())
        print(f"      Found {total_entities} entities:")
        print(f"      - Organisms: {len(entities['organisms'])}")
        print(f"      - Diseases: {len(entities['diseases'])}")
        print(f"      - Molecules: {len(entities['molecules'])}")
        
        # Extract relationships
        print(f"   🔗 Extracting relationships...")
        relationships = self.extract_relationships(text, entities)
        print(f"      Found {len(relationships)} relationships")
        
        # Update graph
        for rel in relationships:
            try:
                subj_type = self.get_entity_type(rel['subject'], entities)
                obj_type = self.get_entity_type(rel['object'], entities)
                
                self.kg.create_relationship(
                    subj_name=rel['subject'],
                    subj_type=subj_type,
                    predicate=rel['predicate'],
                    obj_name=rel['object'],
                    obj_type=obj_type,
                    confidence=float(rel['confidence']),
                    paper_pmid=pmid
                )
                
                print(f"      ✅ {rel['subject']} -[{rel['predicate']}]-> {rel['object']}")
            
            except Exception as e:
                print(f"      ⚠️  Error creating relationship: {e}")
        
        # Mark as analyzed
        self.kg.mark_paper_analyzed(pmid)
        self.redis.delete(f"paper:{pmid}:needs_analysis")
        
        print(f"   ✅ Analysis complete")
    
    def run(self):
        """Main agent loop"""
        print("="*60)
        print("🧠 KG Builder Agent Started")
        print("="*60)
        print("Building knowledge graph from papers...")
        print("Waiting for papers to analyze...\n")
        
        while True:
            try:
                # Find papers needing analysis
                keys = self.redis.keys("paper:*:needs_analysis")
                
                if keys:
                    for key in keys[:5]:  # Process up to 5 at a time
                        try:
                            pmid = key.decode().split(':')[1]
                            self.process_paper(pmid)
                            time.sleep(2)  # Rate limiting
                        
                        except Exception as e:
                            print(f"   ❌ Error processing paper: {e}")
                
                else:
                    # No papers to process, wait
                    time.sleep(5)
            
            except KeyboardInterrupt:
                print("\n\n👋 Shutting down gracefully...")
                break
            
            except Exception as e:
                print(f"❌ Error in main loop: {e}")
                time.sleep(5)


if __name__ == "__main__":
    agent = KGBuilderAgent()
    agent.run()
```

## scripts/test_kg_builder.py

```python
#!/usr/bin/env python3
"""
Test script for KG Builder Agent
"""

import redis
import sys
import time

sys.path.insert(0, '.')
from agents.shared.neo4j_client import KnowledgeGraph


def create_test_paper():
    """Create a test paper for analysis"""
    kg = KnowledgeGraph()
    r = redis.from_url("redis://localhost:6379")
    
    print("Creating test paper...")
    
    test_paper = {
        'pmid': 'TEST_KG_001',
        'title': 'Firmicutes Bacteria Produce Short-Chain Fatty Acids Affecting Obesity',
        'abstract': '''
        Firmicutes bacteria in the gut microbiome produce short-chain fatty acids (SCFAs) 
        that significantly affect obesity and metabolic health. Our study shows that 
        Firmicutes increases SCFA production, which decreases inflammation and affects 
        insulin resistance. These findings suggest that Bacteroidetes to Firmicutes ratio 
        influences obesity outcomes through metabolic pathways.
        ''',
        'year': 2024
    }
    
    kg.create_paper(
        pmid=test_paper['pmid'],
        title=test_paper['title'],
        abstract=test_paper['abstract'],
        year=test_paper['year']
    )
    
    # Mark for analysis
    r.set(f"paper:{test_paper['pmid']}:needs_analysis", "true")
    
    print(f"✅ Created test paper: {test_paper['pmid']}")
    print("\nWaiting for KG Builder to process (30 seconds)...")
    time.sleep(30)
    
    # Check results
    print("\nChecking results in Neo4j...")
    
    result = kg.graph.run("""
        MATCH (p:Paper {pmid: $pmid})
        RETURN p.analyzed as analyzed
    """, pmid=test_paper['pmid']).data()
    
    if result and result[0]['analyzed']:
        print("✅ Paper marked as analyzed")
    else:
        print("⚠️  Paper not yet analyzed")
    
    # Check entities
    entities = kg.graph.run("""
        MATCH (e)
        WHERE e.name IN ['Firmicutes', 'obesity', 'short-chain fatty acids']
        RETURN e.name as name, labels(e) as type
    """).data()
    
    print(f"\n✅ Found {len(entities)} entities:")
    for ent in entities:
        print(f"   - {ent['name']} ({ent['type'][0]})")
    
    # Check relationships
    rels = kg.graph.run("""
        MATCH (a)-[r]->(b)
        WHERE r.paper_pmid = $pmid
        RETURN a.name as subject, type(r) as predicate, b.name as object, r.confidence as confidence
    """, pmid=test_paper['pmid']).data()
    
    print(f"\n✅ Found {len(rels)} relationships:")
    for rel in rels:
        print(f"   - {rel['subject']} -[{rel['predicate']}]-> {rel['object']} (confidence: {rel['confidence']})")


if __name__ == "__main__":
    print("="*60)
    print("KG Builder Agent Test")
    print("="*60)
    print("\nMake sure KG Builder is running:")
    print("  docker-compose up -d kg_builder")
    print("\nPress Enter to start test...")
    input()
    
    try:
        create_test_paper()
        print("\n" + "="*60)
        print("✅ Test complete!")
        print("\nView in Neo4j Browser:")
        print("  MATCH (n)-[r]->(m) WHERE r.paper_pmid = 'TEST_KG_001' RETURN n,r,m")
        print("="*60)
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)
```

## Visualization in Neo4j

Run this query in Neo4j Browser to see the knowledge graph:

```cypher
// View all relationships for test paper
MATCH (n)-[r]->(m) 
WHERE r.paper_pmid = 'TEST_KG_001' 
RETURN n, r, m

// View all organisms and their relationships
MATCH (o:Organism)-[r]->(target)
RETURN o, r, target
LIMIT 50

// Find weak evidence (potential gaps)
MATCH (a)-[r]->(b)
WHERE r.confidence < 0.5
RETURN a.name, type(r), b.name, r.confidence
ORDER BY r.confidence
```

## Testing

```bash
# Start KG Builder
docker-compose up -d kg_builder

# Watch logs
docker-compose logs -f kg_builder

# Run test
python scripts/test_kg_builder.py

# Check Neo4j
open http://localhost:7474
```

## Performance Tuning

For better entity extraction:

```bash
# Install better spaCy model (larger, more accurate)
pip install https://github.com/explosion/spacy-models/releases/download/en_core_web_md-3.7.0/en_core_web_md-3.7.0-py3-none-any.whl

# Update in Dockerfile:
RUN python -m spacy download en_core_web_md
```

## Acceptance Criteria

- [ ] Agent starts and polls for papers
- [ ] Extracts biological entities from abstracts
- [ ] Generates relationships using LLM
- [ ] Creates nodes and relationships in Neo4j
- [ ] Marks papers as analyzed
- [ ] Handles errors gracefully
- [ ] Test script runs successfully
- [ ] Knowledge graph visible in Neo4j browser
