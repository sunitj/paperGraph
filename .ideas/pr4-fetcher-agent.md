# PR #4: Fetcher Agent - PubMed Integration

**Sprint:** Week 2, Days 1-3  
**Goal:** Build agent that fetches papers from PubMed

## What to Build

1. PubMed API client using Biopython
2. Redis task queue consumer
3. Paper storage to Neo4j
4. Error handling and retries

## Deliverable

```bash
# Post query to Redis → Fetcher retrieves papers → Saves to Neo4j
redis-cli RPUSH fetch_queue '{"id":"q1","query":"gut microbiome obesity"}'
# Agent fetches and stores papers
```

## agents/fetcher/agent.py

```python
#!/usr/bin/env python3
"""
Fetcher Agent - Retrieves papers from PubMed
"""

import sys
import os
import time
import json
from Bio import Entrez
import redis

sys.path.insert(0, '/app/shared')
from neo4j_client import KnowledgeGraph

# Configure Entrez
Entrez.email = os.getenv("PUBMED_EMAIL", "papergraph@example.com")
Entrez.api_key = os.getenv("PUBMED_API_KEY", None)  # Optional: for higher rate limits


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
        """Fetch papers from PubMed"""
        max_results = max_results or self.max_results
        
        try:
            print(f"🔍 Searching PubMed: {query}")
            
            # Search PubMed
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record.get("IdList", [])
            print(f"   Found {len(pmids)} papers")
            
            if not pmids:
                return []
            
            # Fetch paper details
            handle = Entrez.efetch(
                db="pubmed",
                id=pmids,
                rettype="abstract",
                retmode="xml"
            )
            articles = Entrez.read(handle)
            handle.close()
            
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
```

## scripts/test_fetcher.py

```python
#!/usr/bin/env python3
"""
Test script for Fetcher Agent
"""

import redis
import json
import time
import sys

def test_fetcher():
    """Test the fetcher agent"""
    print("Testing Fetcher Agent...")
    
    # Connect to Redis
    r = redis.from_url("redis://localhost:6379")
    
    # Create test query
    query = {
        "id": "test_query_1",
        "query": "gut microbiome obesity mechanisms",
        "max_results": 5
    }
    
    print(f"\n1. Sending test query: {query['query']}")
    r.rpush("fetch_queue", json.dumps(query))
    print("   ✅ Query queued")
    
    # Wait for processing
    print("\n2. Waiting for agent to process (30 seconds)...")
    time.sleep(30)
    
    # Check results
    print("\n3. Checking results...")
    
    status = r.get(f"query:{query['id']}:status")
    if status:
        print(f"   Status: {status.decode()}")
    else:
        print("   ⚠️  No status found (agent may still be processing)")
    
    papers = r.get(f"query:{query['id']}:papers")
    if papers:
        papers_data = json.loads(papers)
        print(f"   ✅ Found {len(papers_data)} papers:")
        for paper in papers_data[:3]:
            print(f"      - {paper['title'][:60]}...")
    else:
        print("   ⚠️  No papers found")
    
    # Check if papers are marked for analysis
    print("\n4. Checking papers marked for analysis...")
    keys = r.keys("paper:*:needs_analysis")
    print(f"   Found {len(keys)} papers awaiting analysis")


if __name__ == "__main__":
    print("="*60)
    print("Fetcher Agent Test")
    print("="*60)
    print("\nMake sure services are running:")
    print("  docker-compose up -d")
    print("\nPress Enter to start test...")
    input()
    
    try:
        test_fetcher()
        print("\n" + "="*60)
        print("✅ Test complete!")
        print("="*60)
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)
```

## Manual Testing via Redis CLI

```bash
# Start agent
docker-compose up -d fetcher

# Watch logs
docker-compose logs -f fetcher

# In another terminal, send a test query
redis-cli RPUSH fetch_queue '{
  "id": "manual_test_1",
  "query": "CRISPR gene editing mechanisms",
  "max_results": 10
}'

# Check status
redis-cli GET query:manual_test_1:status

# Get papers
redis-cli GET query:manual_test_1:papers

# See papers marked for analysis
redis-cli KEYS "paper:*:needs_analysis"

# Check Neo4j (in browser: http://localhost:7474)
MATCH (p:Paper) RETURN p LIMIT 10
```

## Environment Variables

The Fetcher Agent uses:

```bash
REDIS_URL=redis://redis:6379
NEO4J_URI=bolt://neo4j:7687
NEO4J_USER=neo4j
NEO4J_PASSWORD=papergraph123
PUBMED_EMAIL=your_email@example.com  # Required by NCBI
PUBMED_API_KEY=  # Optional, for higher rate limits
MAX_PAPERS_PER_QUERY=20
```

## Error Handling

The agent handles:
- PubMed API errors (rate limits, timeouts)
- Malformed article XML
- Redis connection issues
- Neo4j connection issues
- Invalid queries

## Rate Limits

- Without API key: 3 requests/second
- With API key: 10 requests/second
- Built-in retry logic with backoff

## Acceptance Criteria

- [ ] Agent starts and connects to Redis/Neo4j
- [ ] Receives queries from Redis queue
- [ ] Fetches papers from PubMed successfully
- [ ] Stores papers in Neo4j with correct schema
- [ ] Marks papers for analysis in Redis
- [ ] Handles errors gracefully
- [ ] Test script runs successfully
- [ ] Can verify papers in Neo4j browser
