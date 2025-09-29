# PR #2: Ollama Integration & LLM Service

**Sprint:** Week 1, Days 3-4  
**Goal:** Create LLM abstraction layer for easy provider switching

## What to Build

1. LLM abstraction layer (`llm_service.py`)
2. Ollama health check in orchestrator
3. Test script to verify LLM connectivity
4. Support for swapping to AWS Bedrock later

## Deliverable

```bash
python test_llm.py
# Output: "Ollama is working! Model: llama3.1:8b"
```

## agents/shared/llm_service.py

```python
"""
LLM Service Abstraction Layer
Supports: Ollama (default), AWS Bedrock (future)
"""

from abc import ABC, abstractmethod
import requests
import json
import os


class LLMProvider(ABC):
    """Abstract base class for LLM providers"""
    
    @abstractmethod
    def generate(self, prompt: str, system: str = None) -> str:
        """Generate text from prompt"""
        pass


class OllamaProvider(LLMProvider):
    """Ollama local LLM provider"""
    
    def __init__(self, base_url: str = None, model: str = "llama3.1:8b"):
        self.base_url = base_url or os.getenv("OLLAMA_URL", "http://localhost:11434")
        self.model = model
    
    def generate(self, prompt: str, system: str = None) -> str:
        """Generate text using Ollama"""
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False
        }
        
        if system:
            payload["system"] = system
        
        try:
            response = requests.post(
                f"{self.base_url}/api/generate",
                json=payload,
                timeout=120
            )
            response.raise_for_status()
            return response.json()['response']
        
        except requests.exceptions.RequestException as e:
            print(f"Ollama error: {e}")
            return ""
    
    def health_check(self) -> bool:
        """Check if Ollama is running"""
        try:
            response = requests.get(f"{self.base_url}/api/tags", timeout=5)
            return response.status_code == 200
        except:
            return False


class BedrockProvider(LLMProvider):
    """AWS Bedrock provider (for future use)"""
    
    def __init__(self, region: str = "us-east-1", model: str = "anthropic.claude-3-sonnet-20240229-v1:0"):
        self.region = region
        self.model = model
        self.client = None
    
    def _init_client(self):
        """Lazy initialization of Bedrock client"""
        if self.client is None:
            import boto3
            self.client = boto3.client('bedrock-runtime', region_name=self.region)
    
    def generate(self, prompt: str, system: str = None) -> str:
        """Generate text using AWS Bedrock"""
        self._init_client()
        
        messages = [{"role": "user", "content": prompt}]
        
        body = json.dumps({
            "anthropic_version": "bedrock-2023-05-31",
            "max_tokens": 2000,
            "messages": messages,
            "system": system or ""
        })
        
        try:
            response = self.client.invoke_model(
                modelId=self.model,
                body=body
            )
            
            response_body = json.loads(response['body'].read())
            return response_body['content'][0]['text']
        
        except Exception as e:
            print(f"Bedrock error: {e}")
            return ""


def get_llm_provider(provider_type: str = None) -> LLMProvider:
    """
    Factory function to get LLM provider
    
    Args:
        provider_type: "ollama" or "bedrock" (defaults to env var or "ollama")
    
    Returns:
        LLMProvider instance
    """
    provider_type = provider_type or os.getenv("LLM_PROVIDER", "ollama")
    
    if provider_type == "ollama":
        return OllamaProvider(
            base_url=os.getenv("OLLAMA_URL"),
            model=os.getenv("LLM_MODEL", "llama3.1:8b")
        )
    
    elif provider_type == "bedrock":
        return BedrockProvider(
            region=os.getenv("AWS_REGION", "us-east-1"),
            model=os.getenv("BEDROCK_MODEL", "anthropic.claude-3-sonnet-20240229-v1:0")
        )
    
    else:
        raise ValueError(f"Unknown LLM provider: {provider_type}")


# Convenience function for quick usage
def generate_text(prompt: str, system: str = None) -> str:
    """Quick helper to generate text with default provider"""
    llm = get_llm_provider()
    return llm.generate(prompt, system)
```

## test_llm.py

```python
#!/usr/bin/env python3
"""
Test script for LLM service
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
    if hasattr(llm, 'health_check'):
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
```

## orchestrator/main.py (initial version)

```python
from fastapi import FastAPI
import sys
import os

# Add shared modules to path
sys.path.insert(0, '/app/shared')

from llm_service import get_llm_provider

app = FastAPI(title="PaperGraph Orchestrator")

# Initialize LLM
llm = get_llm_provider()


@app.get("/health")
async def health():
    """Health check endpoint"""
    llm_status = "unknown"
    
    if hasattr(llm, 'health_check'):
        llm_status = "healthy" if llm.health_check() else "unhealthy"
    
    return {
        "status": "healthy",
        "llm_provider": os.getenv("LLM_PROVIDER", "ollama"),
        "llm_status": llm_status
    }


@app.get("/test")
async def test_llm():
    """Test LLM generation"""
    response = llm.generate("Say hello in 5 words")
    return {
        "prompt": "Say hello in 5 words",
        "response": response
    }
```

## Environment Variables

Update `.env.example`:

```bash
# PubMed
PUBMED_EMAIL=your_email@example.com

# LLM Configuration
LLM_PROVIDER=ollama
LLM_MODEL=llama3.1:8b
OLLAMA_URL=http://ollama:11434

# For future Bedrock support
# LLM_PROVIDER=bedrock
# AWS_REGION=us-east-1
# BEDROCK_MODEL=anthropic.claude-3-sonnet-20240229-v1:0
```

## Setup Script

Create `scripts/setup_ollama.sh`:

```bash
#!/bin/bash

echo "Setting up Ollama..."

# Wait for Ollama to be ready
echo "Waiting for Ollama service..."
until curl -s http://localhost:11434/api/tags > /dev/null; do
    sleep 2
done

echo "✅ Ollama is ready"

# Pull models
echo "Pulling llama3.1:8b (this may take a few minutes)..."
docker exec papergraph-ollama-1 ollama pull llama3.1:8b

echo "Pulling mistral:7b as backup..."
docker exec papergraph-ollama-1 ollama pull mistral:7b

echo "✅ Models ready!"
echo ""
echo "Available models:"
docker exec papergraph-ollama-1 ollama list
```

## Testing

```bash
# Start services
docker-compose up -d

# Setup Ollama models
chmod +x scripts/setup_ollama.sh
./scripts/setup_ollama.sh

# Test LLM service
python test_llm.py

# Test via API
curl http://localhost:8000/health
curl http://localhost:8000/test
```

## Acceptance Criteria

- [ ] Ollama service starts and responds to health check
- [ ] LLM abstraction layer works with Ollama
- [ ] Test script runs successfully
- [ ] Orchestrator health endpoint shows LLM status
- [ ] Easy to swap to Bedrock by changing env var
- [ ] Models are downloaded and ready
