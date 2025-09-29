# PR #3: Neo4j Setup & Basic Schema

**Sprint:** Week 1, Days 5-7  
**Goal:** Setup knowledge graph database with schema

## What to Build

1. Neo4j connection helper class
2. Schema initialization (constraints, indexes)
3. Basic CRUD operations
4. Test script to verify setup

## Deliverable

```bash
python scripts/init_schema.py
# Creates constraints and indexes in Neo4j
```

## agents/shared/neo4j_client.py

```python
"""
Neo4j Knowledge Graph Client
Handles all graph database operations
"""

from py2neo import Graph, Node, Relationship
from typing import List, Dict, Any
import os


class KnowledgeGraph:
    """Knowledge graph client for PaperGraph"""
    
    def __init__(self, uri: str = None, user: str = None, password: str = None):
        self.uri = uri or os.getenv("NEO4J_URI", "bolt://localhost:7687")
        self.user = user or os.getenv("NEO4J_USER", "neo4j")
        self.password = password or os.getenv("NEO4J_PASSWORD", "papergraph123")
        
        self.graph = Graph(self.uri, auth=(self.user, self.password))
    
    def init_schema(self):
        """Initialize database schema with constraints and indexes"""
        print("Initializing Neo4j schema...")
        
        # Constraints (ensure uniqueness)
        constraints = [
            "CREATE CONSTRAINT organism_name IF NOT EXISTS FOR (o:Organism) REQUIRE o.name IS UNIQUE",
            "CREATE CONSTRAINT disease_name IF NOT EXISTS FOR (d:Disease) REQUIRE d.name IS UNIQUE",
            "CREATE CONSTRAINT molecule_name IF NOT EXISTS FOR (m:Molecule) REQUIRE m.name IS UNIQUE",
            "CREATE CONSTRAINT paper_pmid IF NOT EXISTS FOR (p:Paper) REQUIRE p.pmid IS UNIQUE",
        ]
        
        for constraint in constraints:
            try:
                self.graph.run(constraint)
                print(f"✅ Created: {constraint.split('FOR')[1].split('REQUIRE')[0].strip()}")
            except Exception as e:
                print(f"⚠️  Constraint may already exist: {e}")
        
        # Indexes (improve query performance)
        indexes = [
            "CREATE INDEX paper_year IF NOT EXISTS FOR (p:Paper) ON (p.year)",
            "CREATE INDEX paper_title IF NOT EXISTS FOR (p:Paper) ON (p.title)",
            "CREATE INDEX organism_ncbi IF NOT EXISTS FOR (o:Organism) ON (o.ncbi_id)",
        ]
        
        for index in indexes:
            try:
                self.graph.run(index)
                print(f"✅ Created: {index.split('FOR')[1].strip()}")
            except Exception as e:
                print(f"⚠️  Index may already exist: {e}")
        
        print("✅ Schema initialization complete!")
    
    def create_paper(self, pmid: str, title: str, abstract: str, year: int, 
                     pre_seeded: bool = False) -> Node:
        """Create or update a paper node"""
        result = self.graph.run("""
            MERGE (p:Paper {pmid: $pmid})
            SET p.title = $title,
                p.abstract = $abstract,
                p.year = $year,
                p.pre_seeded = $pre_seeded,
                p.analyzed = false
            RETURN p
        """, pmid=pmid, title=title, abstract=abstract, year=year, 
             pre_seeded=pre_seeded).data()
        
        return result[0]['p'] if result else None
    
    def create_entity(self, name: str, entity_type: str, **properties) -> Node:
        """Create an entity node (Organism, Disease, Molecule)"""
        # Sanitize entity type
        valid_types = ['Organism', 'Disease', 'Molecule']
        if entity_type not in valid_types:
            entity_type = 'Entity'
        
        query = f"""
            MERGE (e:{entity_type} {{name: $name}})
            SET e += $properties
            RETURN e
        """
        
        result = self.graph.run(query, name=name, properties=properties).data()
        return result[0]['e'] if result else None
    
    def create_relationship(self, subj_name: str, subj_type: str, 
                          predicate: str, obj_name: str, obj_type: str,
                          confidence: float = 1.0, paper_pmid: str = None,
                          **properties) -> Relationship:
        """Create a relationship between entities"""
        # Sanitize relationship type (Neo4j doesn't allow spaces or special chars)
        predicate = predicate.upper().replace(' ', '_').replace('-', '_')
        
        query = f"""
            MERGE (subj:{subj_type} {{name: $subj_name}})
            MERGE (obj:{obj_type} {{name: $obj_name}})
            MERGE (subj)-[r:{predicate}]->(obj)
            SET r.confidence = $confidence,
                r.paper_pmid = $paper_pmid,
                r += $properties
            RETURN r
        """
        
        result = self.graph.run(
            query,
            subj_name=subj_name,
            obj_name=obj_name,
            confidence=confidence,
            paper_pmid=paper_pmid,
            properties=properties
        ).data()
        
        return result[0]['r'] if result else None
    
    def get_entity_coverage(self, entities: List[str]) -> List[Dict[str, Any]]:
        """Get existing knowledge about entities"""
        query = """
            MATCH (start)-[r]->(end)
            WHERE start.name IN $entities OR end.name IN $entities
            RETURN start.name as subject,
                   type(r) as predicate,
                   end.name as object,
                   avg(r.confidence) as avg_confidence,
                   count(*) as paper_count
        """
        
        return self.graph.run(query, entities=entities).data()
    
    def get_subgraph(self, entities: List[str], depth: int = 2) -> List[Dict[str, Any]]:
        """Get subgraph around entities"""
        query = f"""
            MATCH path = (start)-[r*1..{depth}]-(end)
            WHERE start.name IN $entities
            RETURN path
            LIMIT 50
        """
        
        return self.graph.run(query, entities=entities).data()
    
    def mark_paper_analyzed(self, pmid: str):
        """Mark a paper as analyzed"""
        self.graph.run("""
            MATCH (p:Paper {pmid: $pmid})
            SET p.analyzed = true
        """, pmid=pmid)
    
    def get_unanalyzed_papers(self, limit: int = 10) -> List[Dict[str, str]]:
        """Get papers that need analysis"""
        result = self.graph.run("""
            MATCH (p:Paper)
            WHERE p.analyzed = false
            RETURN p.pmid as pmid, p.title as title, p.abstract as abstract
            LIMIT $limit
        """, limit=limit).data()
        
        return result
    
    def clear_all(self):
        """Clear all nodes and relationships (USE WITH CAUTION!)"""
        self.graph.run("MATCH (n) DETACH DELETE n")
        print("⚠️  All data cleared!")
```

## scripts/init_schema.py

```python
#!/usr/bin/env python3
"""
Initialize Neo4j schema and test operations
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from agents.shared.neo4j_client import KnowledgeGraph


def test_basic_operations(kg: KnowledgeGraph):
    """Test basic CRUD operations"""
    print("\nTesting basic operations...")
    
    # Create a test paper
    print("\n1. Creating test paper...")
    paper = kg.create_paper(
        pmid="TEST123",
        title="Test Paper: Gut Bacteria and Health",
        abstract="This is a test abstract about gut bacteria.",
        year=2024
    )
    print(f"✅ Created paper: {paper['title']}")
    
    # Create entities
    print("\n2. Creating entities...")
    organisms = kg.create_entity("Firmicutes", "Organism", ncbi_id="1239")
    disease = kg.create_entity("obesity", "Disease", mesh_id="D009765")
    molecule = kg.create_entity("short-chain fatty acids", "Molecule")
    
    print(f"✅ Created organism: {organisms['name']}")
    print(f"✅ Created disease: {disease['name']}")
    print(f"✅ Created molecule: {molecule['name']}")
    
    # Create relationships
    print("\n3. Creating relationships...")
    rel1 = kg.create_relationship(
        subj_name="Firmicutes",
        subj_type="Organism",
        predicate="PRODUCES",
        obj_name="short-chain fatty acids",
        obj_type="Molecule",
        confidence=0.85,
        paper_pmid="TEST123"
    )
    print(f"✅ Created: Firmicutes -[PRODUCES]-> short-chain fatty acids")
    
    rel2 = kg.create_relationship(
        subj_name="short-chain fatty acids",
        subj_type="Molecule",
        predicate="AFFECTS",
        obj_name="obesity",
        obj_type="Disease",
        confidence=0.78,
        paper_pmid="TEST123"
    )
    print(f"✅ Created: short-chain fatty acids -[AFFECTS]-> obesity")
    
    # Query coverage
    print("\n4. Testing coverage query...")
    coverage = kg.get_entity_coverage(["Firmicutes", "obesity"])
    print(f"✅ Found {len(coverage)} relationships")
    for row in coverage:
        print(f"   {row['subject']} -[{row['predicate']}]-> {row['object']} (confidence: {row['avg_confidence']:.2f})")
    
    # Get subgraph
    print("\n5. Testing subgraph query...")
    subgraph = kg.get_subgraph(["Firmicutes"], depth=2)
    print(f"✅ Found {len(subgraph)} paths in subgraph")
    
    print("\n✅ All operations successful!")


def visualize_schema():
    """Print schema visualization"""
    schema = """
    
    Knowledge Graph Schema:
    
    Nodes:
    ├── Paper {pmid, title, abstract, year, pre_seeded, analyzed}
    ├── Organism {name, ncbi_id}
    ├── Disease {name, mesh_id}
    └── Molecule {name, chebi_id}
    
    Relationships:
    ├── -[:PRODUCES {confidence, paper_pmid}]->
    ├── -[:AFFECTS {confidence, paper_pmid}]->
    ├── -[:INHIBITS {confidence, paper_pmid}]->
    ├── -[:CORRELATES_WITH {effect_size, p_value, paper_pmid}]->
    └── -[:MENTIONED_IN {context}]->
    
    Indexes:
    ├── Paper.pmid (unique)
    ├── Paper.year
    ├── Organism.name (unique)
    ├── Disease.name (unique)
    └── Molecule.name (unique)
    """
    print(schema)


if __name__ == "__main__":
    print("="*60)
    print("PaperGraph Neo4j Schema Initialization")
    print("="*60)
    
    try:
        # Initialize knowledge graph
        kg = KnowledgeGraph()
        
        # Create schema
        kg.init_schema()
        
        # Run tests
        test_basic_operations(kg)
        
        # Show schema
        visualize_schema()
        
        print("\n" + "="*60)
        print("✅ Schema initialization complete!")
        print("Neo4j Browser: http://localhost:7474")
        print("Username: neo4j")
        print("Password: papergraph123")
        print("="*60)
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        print("\nMake sure Neo4j is running:")
        print("  docker-compose up -d neo4j")
        sys.exit(1)
```

## scripts/test_neo4j.py

```python
#!/usr/bin/env python3
"""
Quick test script for Neo4j connection
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from agents.shared.neo4j_client import KnowledgeGraph


def test_connection():
    """Test Neo4j connection"""
    try:
        kg = KnowledgeGraph()
        
        # Simple query
        result = kg.graph.run("RETURN 'Hello from Neo4j!' as message").data()
        
        if result:
            print(f"✅ {result[0]['message']}")
            return True
        else:
            print("❌ No response from Neo4j")
            return False
            
    except Exception as e:
        print(f"❌ Connection failed: {e}")
        return False


if __name__ == "__main__":
    print("Testing Neo4j connection...")
    success = test_connection()
    sys.exit(0 if success else 1)
```

## Testing

```bash
# Start Neo4j
docker-compose up -d neo4j

# Wait for it to be ready (takes ~30 seconds)
sleep 30

# Test connection
python scripts/test_neo4j.py

# Initialize schema
python scripts/init_schema.py

# View in browser
open http://localhost:7474
# Login: neo4j / papergraph123

# Query test data
# Run in Neo4j browser:
MATCH (n) RETURN n LIMIT 25
```

## Acceptance Criteria

- [ ] Neo4j starts successfully
- [ ] Schema initialization runs without errors
- [ ] Test operations create nodes and relationships
- [ ] Can query data via Neo4j browser
- [ ] Constraints and indexes are created
- [ ] Knowledge graph client is reusable across agents
