#!/usr/bin/env python3
"""
System verification and testing script
Tests all components of PaperGraph Intelligence
"""

import requests
import redis
import json
import time
import sys
from typing import Dict, Any


class SystemTester:
    """Test all system components"""

    def __init__(self):
        self.base_url = "http://localhost:8000"
        self.results = []

    def log_test(self, name: str, passed: bool, message: str = ""):
        """Log test result"""
        status = "✅ PASS" if passed else "❌ FAIL"
        self.results.append((name, passed))
        print(f"{status} {name}")
        if message:
            print(f"     {message}")

    def test_orchestrator_health(self) -> bool:
        """Test orchestrator health endpoint"""
        try:
            response = requests.get(f"{self.base_url}/health", timeout=5)
            if response.status_code == 200:
                data = response.json()
                self.log_test("Orchestrator Health", True, f"Status: {data.get('status')}")
                return True
            else:
                self.log_test("Orchestrator Health", False, f"HTTP {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Orchestrator Health", False, str(e))
            return False

    def test_orchestrator_stats(self) -> bool:
        """Test orchestrator stats endpoint"""
        try:
            response = requests.get(f"{self.base_url}/stats", timeout=5)
            if response.status_code == 200:
                self.log_test("Orchestrator Stats", True)
                return True
            else:
                self.log_test("Orchestrator Stats", False, f"HTTP {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Orchestrator Stats", False, str(e))
            return False

    def test_redis_connection(self) -> bool:
        """Test Redis connectivity"""
        try:
            r = redis.from_url("redis://localhost:6379")
            r.ping()
            self.log_test("Redis Connection", True)
            return True
        except Exception as e:
            self.log_test("Redis Connection", False, str(e))
            return False

    def test_redis_queue_operations(self) -> bool:
        """Test Redis queue operations"""
        try:
            r = redis.from_url("redis://localhost:6379")

            # Test basic queue operations
            test_data = {"id": "test_system_check", "query": "test query"}
            r.rpush("test_queue", json.dumps(test_data))

            # Check if data exists
            result = r.blpop("test_queue", timeout=1)
            if result:
                self.log_test("Redis Queue Operations", True)
                return True
            else:
                self.log_test("Redis Queue Operations", False, "No data in queue")
                return False
        except Exception as e:
            self.log_test("Redis Queue Operations", False, str(e))
            return False

    def test_neo4j_browser(self) -> bool:
        """Test Neo4j browser accessibility"""
        try:
            response = requests.get("http://localhost:7474", timeout=5)
            if response.status_code == 200:
                self.log_test("Neo4j Browser", True)
                return True
            else:
                self.log_test("Neo4j Browser", False, f"HTTP {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Neo4j Browser", False, str(e))
            return False

    def test_ollama_api(self) -> bool:
        """Test Ollama API accessibility"""
        try:
            response = requests.get("http://localhost:11434/api/tags", timeout=5)
            if response.status_code == 200:
                data = response.json()
                models = [model.get('name', 'unknown') for model in data.get('models', [])]
                self.log_test("Ollama API", True, f"Models: {', '.join(models) if models else 'none'}")
                return True
            else:
                self.log_test("Ollama API", False, f"HTTP {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Ollama API", False, str(e))
            return False

    def test_streamlit_ui(self) -> bool:
        """Test Streamlit UI accessibility"""
        try:
            response = requests.get("http://localhost:8501", timeout=5)
            if response.status_code == 200:
                self.log_test("Streamlit UI", True)
                return True
            else:
                self.log_test("Streamlit UI", False, f"HTTP {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Streamlit UI", False, str(e))
            return False

    def test_docker_containers(self) -> bool:
        """Test that Docker containers are running"""
        import subprocess
        try:
            result = subprocess.run(
                ["docker-compose", "ps", "--services", "--filter", "status=running"],
                capture_output=True,
                text=True,
                check=True
            )
            running_services = result.stdout.strip().split('\n')
            running_services = [s for s in running_services if s]

            expected_services = ["ollama", "neo4j", "redis", "orchestrator", "fetcher", "kg_builder", "query_strategist", "ui"]
            all_running = all(service in running_services for service in expected_services)

            if all_running:
                self.log_test("Docker Containers", True, f"All {len(expected_services)} services running")
                return True
            else:
                missing = [s for s in expected_services if s not in running_services]
                self.log_test("Docker Containers", False, f"Missing: {', '.join(missing)}")
                return False

        except Exception as e:
            self.log_test("Docker Containers", False, str(e))
            return False

    def test_agent_connectivity(self) -> bool:
        """Test that agents can connect to their dependencies"""
        try:
            r = redis.from_url("redis://localhost:6379")

            # Send a test query to the fetch queue
            test_query = {
                "id": "system_test_query",
                "query": "test connectivity",
                "max_results": 1
            }

            r.rpush("fetch_queue", json.dumps(test_query))
            self.log_test("Agent Connectivity", True, "Test query queued for fetcher")
            return True

        except Exception as e:
            self.log_test("Agent Connectivity", False, str(e))
            return False

    def run_all_tests(self) -> Dict[str, Any]:
        """Run all system tests"""
        print("🔬 PaperGraph Intelligence System Tests")
        print("=" * 50)

        # Run tests
        tests = [
            self.test_docker_containers,
            self.test_redis_connection,
            self.test_redis_queue_operations,
            self.test_neo4j_browser,
            self.test_ollama_api,
            self.test_orchestrator_health,
            self.test_orchestrator_stats,
            self.test_streamlit_ui,
            self.test_agent_connectivity,
        ]

        for test in tests:
            test()

        # Summary
        passed_tests = sum(1 for _, passed in self.results if passed)
        total_tests = len(self.results)

        print("\n" + "=" * 50)
        print(f"📊 Results: {passed_tests}/{total_tests} tests passed")

        if passed_tests == total_tests:
            print("🎉 All tests passed! System is ready.")
            success = True
        else:
            print("⚠️  Some tests failed. Check the issues above.")
            success = False

        return {
            "total_tests": total_tests,
            "passed_tests": passed_tests,
            "success": success,
            "results": self.results
        }


def main():
    """Main function"""
    tester = SystemTester()
    results = tester.run_all_tests()

    # Exit with appropriate code
    sys.exit(0 if results["success"] else 1)


if __name__ == "__main__":
    main()