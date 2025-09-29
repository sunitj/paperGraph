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
    organism = kg.create_entity("Firmicutes", "Organism", ncbi_id="1239")
    disease = kg.create_entity("obesity", "Disease", mesh_id="D009765")
    molecule = kg.create_entity("short-chain fatty acids", "Molecule")

    print(f"✅ Created organism: {organism['name']}")
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