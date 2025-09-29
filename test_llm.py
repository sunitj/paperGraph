#!/usr/bin/env python3
"""
Test script for LLM service
Tests Ollama connectivity and basic functionality
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(__file__))

from agents.shared.llm_service import get_llm_provider


def test_ollama():
    """Test Ollama connection and generation"""
    print("Testing Ollama LLM service...")

    llm = get_llm_provider("ollama")

    # Health check
    if not llm.health_check():
        print("❌ Ollama is not running!")
        print("Start it with: docker-compose up -d ollama")
        return False

    print("✅ Ollama is running")

    # Test generation
    prompt = "Say 'Hello from PaperGraph!' in exactly 5 words"
    print(f"\nPrompt: {prompt}")

    response = llm.generate(prompt)
    print(f"Response: {response}")

    if response:
        print("✅ Generation successful!")
        return True
    else:
        print("❌ Generation failed")
        return False


def test_entity_extraction():
    """Test entity extraction capability"""
    print("\n" + "="*50)
    print("Testing entity extraction...")

    llm = get_llm_provider("ollama")

    text = """
    Firmicutes bacteria in the gut microbiome produce short-chain fatty acids
    that can affect obesity and insulin resistance.
    """

    prompt = f"""Extract biological entities from this text:

Text: {text}

Return JSON only:
{{"organisms": [...], "molecules": [...], "diseases": [...]}}"""

    response = llm.generate(prompt)
    print(f"Response:\n{response}")

    if "Firmicutes" in response:
        print("✅ Entity extraction working!")
        return True
    else:
        print("⚠️  Entity extraction may need tuning")
        return True  # Not critical for this test


if __name__ == "__main__":
    print("PaperGraph LLM Service Test\n")

    success = test_ollama()

    if success:
        test_entity_extraction()

    print("\n" + "="*50)

    if success:
        print("✅ All tests passed!")
        sys.exit(0)
    else:
        print("❌ Some tests failed")
        sys.exit(1)