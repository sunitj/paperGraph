"""
LLM Service Abstraction Layer
Supports: Ollama (default), AWS Bedrock (future)
"""

from abc import ABC, abstractmethod
import requests
import json
import os
from typing import Optional


class LLMProvider(ABC):
    """Abstract base class for LLM providers"""

    @abstractmethod
    def generate(self, prompt: str, system: Optional[str] = None) -> str:
        """Generate text from prompt"""
        pass

    @abstractmethod
    def health_check(self) -> bool:
        """Check if the LLM service is available"""
        pass


class OllamaProvider(LLMProvider):
    """Ollama local LLM provider"""

    def __init__(self, base_url: Optional[str] = None, model: str = "llama3.1:8b"):
        self.base_url = base_url or os.getenv("OLLAMA_URL", "http://localhost:11434")
        self.model = model

    def generate(self, prompt: str, system: Optional[str] = None) -> str:
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
            return response.json().get('response', '')

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

    def generate(self, prompt: str, system: Optional[str] = None) -> str:
        """Generate text using AWS Bedrock"""
        # Stub implementation for now
        return "Bedrock provider not yet implemented"

    def health_check(self) -> bool:
        """Check if Bedrock is available"""
        # Stub implementation for now
        return False


def get_llm_provider(provider_type: Optional[str] = None) -> LLMProvider:
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
def generate_text(prompt: str, system: Optional[str] = None) -> str:
    """Quick helper to generate text with default provider"""
    llm = get_llm_provider()
    return llm.generate(prompt, system)