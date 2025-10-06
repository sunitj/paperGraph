#!/usr/bin/env python3
"""
Fetcher Agent - Retrieves papers from PubMed using Biopython
"""

import sys
import os
import time
import json
import redis
from Bio import Entrez

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

        # Configure Entrez for PubMed API
        Entrez.email = os.getenv("PUBMED_EMAIL", "papergraph@example.com")
        api_key = os.getenv("PUBMED_API_KEY")
        if api_key:
            Entrez.api_key = api_key
            self.rate_limit_delay = 0.11  # 10 requests/sec with API key
        else:
            self.rate_limit_delay = 0.34  # 3 requests/sec without API key

        self.last_request_time = 0

    def _rate_limit(self):
        """Enforce NCBI rate limiting"""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit_delay:
            time.sleep(self.rate_limit_delay - elapsed)
        self.last_request_time = time.time()

    def fetch_papers(self, query: str, max_results: int = None) -> list:
        """Fetch papers from PubMed using Biopython"""
        max_results = max_results or self.max_results

        try:
            print(f"🔍 Searching PubMed: {query}")

            # Step 1: Search PubMed for PMIDs
            self._rate_limit()
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()

            pmids = record.get("IdList", [])
            print(f"   📄 Found {len(pmids)} PMIDs")

            if not pmids:
                print(f"   ⚠️  No papers found")
                return []

            # Step 2: Fetch paper details
            self._rate_limit()
            handle = Entrez.efetch(
                db="pubmed",
                id=pmids,
                rettype="abstract",
                retmode="xml"
            )
            articles = Entrez.read(handle)
            handle.close()

            # Step 3: Parse articles
            papers = []
            for article in articles.get('PubmedArticle', []):
                try:
                    paper_data = self._parse_article(article)
                    if paper_data:
                        papers.append(paper_data)
                except Exception as e:
                    print(f"   ⚠️  Error parsing article: {e}")
                    continue

            print(f"   ✅ Successfully parsed {len(papers)} papers")
            return papers

        except Exception as e:
            print(f"   ❌ PubMed error: {e}")
            return []

    def _parse_article(self, article):
        """Parse PubMed article XML"""
        medline = article.get('MedlineCitation', {})

        # Extract PMID
        pmid = str(medline.get('PMID', ''))
        if not pmid:
            return None

        # Extract article data
        article_data = medline.get('Article', {})

        # Title
        title = article_data.get('ArticleTitle', '')
        if isinstance(title, list):
            title = ' '.join(title)

        # Abstract
        abstract_data = article_data.get('Abstract', {})
        abstract_text = abstract_data.get('AbstractText', [])

        if isinstance(abstract_text, list):
            # Handle structured abstracts
            abstract = ' '.join([str(text) for text in abstract_text])
        else:
            abstract = str(abstract_text)

        # Year
        pub_date = article_data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
        year = pub_date.get('Year', pub_date.get('MedlineDate', '2024'))

        # Extract just the year if it's a range like "2023-2024"
        if isinstance(year, str) and '-' in year:
            year = year.split('-')[0]

        try:
            year = int(year)
        except:
            year = 2024

        return {
            'pmid': pmid,
            'title': title,
            'abstract': abstract,
            'year': year
        }

    def store_papers(self, papers: list, query_id: str):
        """Store papers in Neo4j and mark for analysis"""
        stored = []

        for paper in papers:
            try:
                # Store in Neo4j
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
        print("🔬 Fetcher Agent Started")
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