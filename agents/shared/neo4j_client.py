"""
Neo4j Knowledge Graph Client
Handles all graph database operations
"""

from py2neo import Graph, Node, Relationship
from typing import List, Dict, Any, Optional
import os


class KnowledgeGraph:
    """Knowledge graph client for PaperGraph"""

    def __init__(self, uri: Optional[str] = None, user: Optional[str] = None, password: Optional[str] = None):
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
                     pre_seeded: bool = False) -> Optional[Node]:
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

    def create_entity(self, name: str, entity_type: str, **properties) -> Optional[Node]:
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
                          confidence: float = 1.0, paper_pmid: Optional[str] = None,
                          **properties) -> Optional[Relationship]:
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

    def health_check(self) -> bool:
        """Check if Neo4j is accessible"""
        try:
            self.graph.run("RETURN 1").data()
            return True
        except Exception:
            return False