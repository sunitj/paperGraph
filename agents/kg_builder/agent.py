#!/usr/bin/env python3
"""
KG Builder Agent - Extracts entities and relationships (stub implementation)
"""

import sys
import os
import time
import json
import redis

sys.path.insert(0, '/app/shared')
from neo4j_client import KnowledgeGraph
from llm_service import get_llm_provider


class KGBuilderAgent:
    """Agent that extracts entities and relationships from papers"""

    def __init__(self):
        self.redis = redis.from_url(os.getenv("REDIS_URL", "redis://localhost:6379"))
        self.kg = KnowledgeGraph(
            uri=os.getenv("NEO4J_URI"),
            user=os.getenv("NEO4J_USER"),
            password=os.getenv("NEO4J_PASSWORD")
        )
        self.llm = get_llm_provider()

    def extract_entities(self, text: str) -> dict:
        """Extract entities from text (stub)"""
        print(f"   🧠 Extracting entities (stub) from: {text[:50]}...")

        # Stub implementation - return fake entities
        fake_entities = {
            "organisms": ["Firmicutes", "Bacteroidetes"],
            "diseases": ["obesity", "diabetes"],
            "molecules": ["short-chain fatty acids", "glucose"]
        }

        print(f"   ✅ Extracted {sum(len(v) for v in fake_entities.values())} entities (stub)")
        return fake_entities

    def extract_relationships(self, text: str, entities: dict) -> list:
        """Extract relationships between entities (stub)"""
        print(f"   🔗 Extracting relationships (stub)...")

        # Stub implementation - return fake relationships
        fake_relationships = [
            {
                "subject": "Firmicutes",
                "subject_type": "Organism",
                "predicate": "PRODUCES",
                "object": "short-chain fatty acids",
                "object_type": "Molecule",
                "confidence": 0.85
            },
            {
                "subject": "short-chain fatty acids",
                "subject_type": "Molecule",
                "predicate": "AFFECTS",
                "object": "obesity",
                "object_type": "Disease",
                "confidence": 0.72
            }
        ]

        print(f"   ✅ Extracted {len(fake_relationships)} relationships (stub)")
        return fake_relationships

    def process_paper(self, pmid: str):
        """Process a single paper for entity/relationship extraction"""
        print(f"\n📄 Processing paper {pmid}")

        # Get paper from Neo4j (stub - would query actual database)
        paper_text = f"Stub abstract for paper {pmid}"

        # Extract entities
        entities = self.extract_entities(paper_text)

        # Create entity nodes in graph
        for entity_type, entity_list in entities.items():
            for entity_name in entity_list:
                self.kg.create_entity(entity_name, entity_type.rstrip('s').capitalize())

        # Extract relationships
        relationships = self.extract_relationships(paper_text, entities)

        # Create relationships in graph
        for rel in relationships:
            self.kg.create_relationship(
                subj_name=rel["subject"],
                subj_type=rel["subject_type"],
                predicate=rel["predicate"],
                obj_name=rel["object"],
                obj_type=rel["object_type"],
                confidence=rel["confidence"],
                paper_pmid=pmid
            )

        # Mark paper as analyzed
        self.kg.mark_paper_analyzed(pmid)

        # Remove from Redis analysis queue
        self.redis.delete(f"paper:{pmid}:needs_analysis")

        print(f"✅ Paper {pmid} analysis complete\n")

    def run(self):
        """Main agent loop"""
        print("="*60)
        print("🧠 KG Builder Agent Started (stub)")
        print("="*60)
        print("Waiting for papers to analyze...\n")

        while True:
            try:
                # Look for papers that need analysis
                keys = self.redis.keys("paper:*:needs_analysis")

                if keys:
                    for key in keys[:5]:  # Process up to 5 papers at a time
                        pmid = key.decode().split(':')[1]
                        self.process_paper(pmid)
                else:
                    time.sleep(5)  # Wait 5 seconds before checking again

            except KeyboardInterrupt:
                print("\n\n👋 Shutting down gracefully...")
                break

            except Exception as e:
                print(f"❌ Error in main loop: {e}")
                time.sleep(5)


if __name__ == "__main__":
    agent = KGBuilderAgent()
    agent.run()