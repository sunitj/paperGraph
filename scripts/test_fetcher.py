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
