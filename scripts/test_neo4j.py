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