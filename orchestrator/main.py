#!/usr/bin/env python3
"""
PaperGraph Intelligence Orchestrator
FastAPI coordination service for multi-agent system
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import sys
import os
import redis
from typing import Dict, Any

# Add shared modules to path
sys.path.insert(0, '/app/shared')

try:
    from llm_service import get_llm_provider
    llm = get_llm_provider()
except ImportError:
    print("Warning: LLM service not available")
    llm = None

app = FastAPI(
    title="PaperGraph Orchestrator",
    description="Multi-agent literature analysis coordination service",
    version="0.1.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize Redis connection
try:
    redis_client = redis.from_url(os.getenv("REDIS_URL", "redis://localhost:6379"))
    redis_client.ping()
except Exception as e:
    print(f"Warning: Redis connection failed: {e}")
    redis_client = None


@app.get("/health")
async def health() -> Dict[str, Any]:
    """Health check endpoint"""

    # Check Redis
    redis_status = "disconnected"
    try:
        if redis_client:
            redis_client.ping()
            redis_status = "connected"
    except:
        pass

    # Check Neo4j (stub for now)
    neo4j_status = "unknown"

    # Check LLM
    llm_status = "unknown"
    llm_provider = os.getenv("LLM_PROVIDER", "ollama")

    if llm and hasattr(llm, 'health_check'):
        try:
            llm_status = "healthy" if llm.health_check() else "unhealthy"
        except Exception:
            llm_status = "error"

    return {
        "status": "healthy",
        "version": "0.1.0",
        "llm_provider": llm_provider,
        "llm_status": llm_status,
        "redis": redis_status,
        "neo4j": neo4j_status,
        "timestamp": "2024-01-01T00:00:00Z"
    }


@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "message": "PaperGraph Intelligence Orchestrator",
        "version": "0.1.0",
        "docs": "/docs"
    }


@app.get("/test")
async def test_llm():
    """Test LLM generation"""
    if not llm:
        raise HTTPException(status_code=503, detail="LLM service not available")

    try:
        response = llm.generate("Say hello in 5 words")
        return {
            "prompt": "Say hello in 5 words",
            "response": response,
            "provider": os.getenv("LLM_PROVIDER", "ollama")
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"LLM generation failed: {str(e)}")


@app.get("/stats")
async def get_stats():
    """Get system statistics (stub)"""
    return {
        "papers": 0,
        "entities": 0,
        "relationships": 0,
        "entity_breakdown": {
            "Organism": 0,
            "Disease": 0,
            "Molecule": 0
        },
        "relationship_breakdown": {
            "PRODUCES": 0,
            "AFFECTS": 0,
            "INHIBITS": 0,
            "CORRELATES_WITH": 0
        },
        "average_confidence": 0.0,
        "pre_seeded_papers": 0,
        "user_added_papers": 0,
        "last_updated": "2024-01-01T00:00:00Z"
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)