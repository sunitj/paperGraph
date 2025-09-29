#!/usr/bin/env python3
"""
Fetcher Agent - Retrieves papers from PubMed (stub implementation)
"""

import sys
import os
import time
import json
import redis

sys.path.insert(0, '/app/shared')
from neo4j_client import KnowledgeGraph


class FetcherAgent:
    """Agent that fetches papers from PubMed"""

    def __init__(self):
        self.redis = redis.from_url(os.getenv("REDIS_URL", "redis://localhost:6379"))
        self.kg = KnowledgeGraph(
            uri=os.getenv("NEO4J_URI"),
            user=os.getenv("NEO4J_USER"),
            password=os.getenv("NEO4J_PASSWORD")
        )
        self.max_results = int(os.getenv("MAX_PAPERS_PER_QUERY", "20"))

    def fetch_papers(self, query: str, max_results: int = None) -> list:
        """Fetch papers from PubMed (stub)"""
        max_results = max_results or self.max_results
        print(f"🔍 Searching PubMed (stub): {query}")

        # Stub implementation - return fake papers
        fake_papers = [
            {
                'pmid': f'STUB{i:06d}',
                'title': f'Stub Paper {i}: {query}',
                'abstract': f'This is a stub abstract for paper {i} about {query}.',
                'year': 2024
            }
            for i in range(1, min(max_results + 1, 6))  # Return up to 5 fake papers
        ]

        print(f"   ✅ Found {len(fake_papers)} papers (stub)")
        return fake_papers

    def store_papers(self, papers: list, query_id: str):
        """Store papers in Neo4j and mark for analysis"""
        stored = []

        for paper in papers:
            try:
                # Store in Neo4j (stub)
                self.kg.create_paper(
                    pmid=paper['pmid'],
                    title=paper['title'],
                    abstract=paper['abstract'],
                    year=paper['year'],
                    pre_seeded=False
                )

                # Mark for analysis in Redis
                self.redis.set(
                    f"paper:{paper['pmid']}:needs_analysis",
                    "true",
                    ex=86400  # Expire after 24 hours
                )

                stored.append({
                    'pmid': paper['pmid'],
                    'title': paper['title']
                })

            except Exception as e:
                print(f"   ⚠️  Error storing paper {paper['pmid']}: {e}")

        # Store results for query
        self.redis.set(
            f"query:{query_id}:papers",
            json.dumps(stored),
            ex=3600  # Expire after 1 hour
        )

        print(f"   ✅ Stored {len(stored)} papers")
        return stored

    def process_query(self, query_data: dict):
        """Process a single query from the queue"""
        query_id = query_data.get('id')
        query = query_data.get('query')
        max_results = query_data.get('max_results', self.max_results)

        print(f"\n📥 Processing query {query_id}: {query}")

        # Fetch papers
        papers = self.fetch_papers(query, max_results)

        if papers:
            # Store papers
            stored = self.store_papers(papers, query_id)

            # Mark query as complete
            self.redis.set(
                f"query:{query_id}:status",
                "complete",
                ex=3600
            )

            print(f"✅ Query {query_id} complete: {len(stored)} papers\n")
        else:
            print(f"⚠️  Query {query_id}: No papers found\n")

            self.redis.set(
                f"query:{query_id}:status",
                "no_results",
                ex=3600
            )

    def run(self):
        """Main agent loop"""
        print("="*60)
        print("🔬 Fetcher Agent Started (stub)")
        print("="*60)
        print(f"Max papers per query: {self.max_results}")
        print(f"Listening on: fetch_queue")
        print("Waiting for queries...\n")

        while True:
            try:
                # Block and wait for queries (with timeout)
                result = self.redis.blpop("fetch_queue", timeout=5)

                if result:
                    _, query_json = result
                    query_data = json.loads(query_json)
                    self.process_query(query_data)

            except KeyboardInterrupt:
                print("\n\n👋 Shutting down gracefully...")
                break

            except Exception as e:
                print(f"❌ Error in main loop: {e}")
                time.sleep(5)  # Brief pause before retrying


if __name__ == "__main__":
    agent = FetcherAgent()
    agent.run()