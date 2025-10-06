# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PaperGraph Intelligence is a multi-agent literature analysis system that builds a knowledge graph to identify research gaps and guide scientific discovery. The system uses intelligent agents to fetch papers from PubMed, extract entities/relationships, and synthesize answers based on accumulated knowledge.

**Current Focus:** Ollama-powered local MVP. Free to use, no cloud dependencies required.

## System Architecture

### Core Components
- **Ollama**: Local LLM service running on host (llama3.1:8b primary, mistral:7b backup) - completely free
- **Neo4j**: Knowledge graph database storing papers, entities, and relationships
- **Redis**: Message queue and state management between agents
- **Orchestrator**: FastAPI service coordinating agent workflows
- **Streamlit UI**: User interface for queries and visualization

**Important:** Ollama runs on the host machine at `localhost:11434`. Docker services connect via `host.docker.internal:11434`. AWS Bedrock support exists as stub for future enterprise use, but Ollama is the primary focus.

### Multi-Agent System
```
User Question → Query Strategist → Fetcher Agent → KG Builder → Orchestrator → UI
```

1. **Query Strategist**: Analyzes knowledge graph gaps and generates targeted PubMed queries (stub)
2. **Fetcher Agent**: Retrieves papers from PubMed using Biopython E-utilities API
3. **KG Builder Agent**: Extracts entities and relationships using LLM and spaCy (stub)
4. **Orchestrator**: Coordinates workflow and synthesizes final answers (partial)

## Project Structure

```
papergraph/
├── docker-compose.yml          # Service orchestration
├── .env.example                 # Environment template
├── orchestrator/               # FastAPI coordination service
│   ├── pyproject.toml          # UV dependencies
│   ├── Dockerfile              # Container definition
│   └── main.py                 # FastAPI application
├── agents/
│   ├── shared/                 # Common utilities
│   │   ├── __init__.py
│   │   ├── llm_service.py      # LLM abstraction layer
│   │   └── neo4j_client.py     # Knowledge graph client
│   ├── fetcher/                # PubMed paper retrieval
│   │   ├── pyproject.toml
│   │   ├── Dockerfile
│   │   └── agent.py
│   ├── kg_builder/             # Entity and relationship extraction
│   │   ├── pyproject.toml
│   │   ├── Dockerfile
│   │   └── agent.py
│   └── query_strategist/       # Knowledge gap analysis
│       ├── pyproject.toml
│       ├── Dockerfile
│       └── agent.py
├── ui/                         # Streamlit interface
│   ├── pyproject.toml
│   ├── Dockerfile
│   └── app.py
├── seeds/                      # Pre-built knowledge graphs
└── scripts/
    ├── setup.sh                # Automated setup script
    ├── test_system.py          # System verification
    └── prepare_seeds.py        # Seed data generation
```

## Development Commands

### Quick Setup
```bash
# Automated setup (recommended for first time)
./scripts/setup.sh

# Manual setup
cp .env.example .env                    # Configure environment (add PUBMED_EMAIL)
docker-compose up -d                    # Start all services
docker-compose exec ollama ollama pull llama3.1:8b   # Pull primary LLM model
docker-compose exec ollama ollama pull mistral:7b    # Pull backup model
python scripts/init_schema.py           # Initialize Neo4j schema (optional)
```

### Service Management
```bash
docker-compose up -d                    # Start all services
docker-compose up --build [service]     # Rebuild and start specific service
docker-compose down                     # Stop all services
docker-compose down -v                  # Stop and remove volumes (full reset)
docker-compose restart [service]        # Restart specific service
docker-compose logs -f [service]        # View service logs (follow mode)
docker-compose ps                       # Check service status
docker-compose exec [service] bash      # Access service shell
```

### Local Development with UV
```bash
# Work on orchestrator
cd orchestrator
uv sync                                 # Install dependencies (creates .venv)
uv run uvicorn main:app --reload        # Run with hot reload on port 8000
uv run pytest                           # Run tests
uv run pytest -v                        # Run tests with verbose output
uv run ruff check .                     # Lint code
uv run black .                          # Format code

# Work on agents (similar pattern for fetcher, kg_builder, query_strategist)
cd agents/fetcher
uv sync                                 # Install dependencies
uv run python agent.py                  # Run agent directly (standalone mode)

# Note: Each service has its own pyproject.toml and isolated dependencies
# The shared/ module is copied into agent containers at build time
```

### Testing and Verification
```bash
# System health check
python scripts/test_system.py

# Service health endpoints
curl http://localhost:8000/health        # Orchestrator
curl http://localhost:7474               # Neo4j
curl http://localhost:11434/api/tags     # Ollama
curl http://localhost:8501/_stcore/health # Streamlit

# Check Redis
docker-compose exec redis redis-cli ping
```

### Service URLs
- Streamlit UI: http://localhost:8501
- Orchestrator API: http://localhost:8000
- Neo4j Browser: http://localhost:7474 (neo4j/papergraph123)
- Ollama API: http://localhost:11434

## Knowledge Graph Schema

### Nodes
- `Paper {pmid, title, abstract, year, pre_seeded, analyzed}`
- `Organism {name, ncbi_id}`
- `Disease {name, mesh_id}`
- `Molecule {name, chebi_id}`

### Relationships
- `-[:PRODUCES {confidence, paper_pmid}]->`
- `-[:AFFECTS {confidence, paper_pmid}]->`
- `-[:INHIBITS {confidence, paper_pmid}]->`
- `-[:CORRELATES_WITH {effect_size, p_value, paper_pmid}]->`

## Environment Variables

Required in `.env`:
```bash
PUBMED_EMAIL=your_email@example.com     # Required by NCBI
LLM_PROVIDER=ollama                     # "ollama" or "bedrock"
LLM_MODEL=llama3.1:8b                   # Model name
NEO4J_PASSWORD=papergraph123            # Neo4j auth
```

## Testing Approach

### Unit Testing
- Test individual agent functions
- Mock external services (PubMed, LLM)
- Test knowledge graph operations

### Integration Testing
- End-to-end query flow
- Agent communication via Redis
- Data persistence in Neo4j

### Manual Testing
```bash
# Use automated system test
python scripts/test_system.py

# Test individual services
curl -X POST http://localhost:8000/query \
  -H "Content-Type: application/json" \
  -d '{"question": "How do gut bacteria influence obesity?"}'

# Test Redis queues
docker-compose exec redis redis-cli
> LPUSH test_queue '{"message": "test"}'
> LPOP test_queue

# Verify data in Neo4j browser (http://localhost:7474)
# Username: neo4j, Password: papergraph123
MATCH (n) RETURN count(n)  # Count all nodes
```

## Key Development Patterns

### LLM Service
The `agents/shared/llm_service.py` provides LLM abstraction with Ollama as the primary provider. Bedrock exists as a future stub. Always use `get_llm_provider()` which defaults to Ollama for the MVP.

**Usage Pattern:**
```python
from shared.llm_service import get_llm_provider, generate_text

# Factory pattern (recommended)
llm = get_llm_provider()  # Reads LLM_PROVIDER env var
response = llm.generate(prompt, system="You are a helpful assistant")

# Quick helper
response = generate_text(prompt, system="...")
```

### Redis Communication
Agents communicate via Redis queues using JSON messages. Use descriptive queue names like `fetch_queue`, `analysis_queue`.

**Message Pattern:**
```python
import redis
import json

r = redis.from_url(os.getenv("REDIS_URL"))
r.lpush("queue_name", json.dumps({"key": "value"}))  # Producer
message = json.loads(r.brpop("queue_name")[1])       # Consumer (blocking)
```

### Neo4j Operations
Use the `KnowledgeGraph` client class (`agents/shared/neo4j_client.py`) for all database operations. It provides helper methods for creating entities, relationships, and querying coverage.

**Usage Pattern:**
```python
from shared.neo4j_client import KnowledgeGraph

kg = KnowledgeGraph()  # Reads NEO4J_* env vars
kg.create_paper(pmid, title, abstract, year)
kg.create_entity("E. coli", "Organism", ncbi_id="562")
kg.create_relationship(subj_name, subj_type, "PRODUCES", obj_name, obj_type,
                       confidence=0.9, paper_pmid=pmid)
```

### Pre-seeded Knowledge
The system supports 6 pre-built topics (Gut Microbiome, Protein Engineering, etc.) with 20-30 papers each for faster initial queries.

## Development Workflows

### Adding New Entity Types
1. Update Neo4j schema in `agents/shared/neo4j_client.py`
2. Modify entity extraction prompts in `agents/kg_builder/agent.py`
3. Update UI visualization logic in `ui/app.py`
4. Test with `python scripts/test_system.py`

### Adding New Relationship Types
1. Update relationship extraction in `agents/kg_builder/agent.py`
2. Test extraction with sample papers
3. Verify in Neo4j browser at http://localhost:7474

### Debugging Issues
```bash
# Check service logs
docker-compose logs -f orchestrator
docker-compose logs -f fetcher
docker-compose logs -f kg_builder

# Check all services at once
docker-compose logs --tail=50

# Test individual components
cd orchestrator && uv run pytest -v
python scripts/test_system.py

# Test connections
curl http://localhost:8000/health          # Orchestrator
curl http://localhost:11434/api/tags       # Ollama (check models)
docker-compose exec redis redis-cli ping   # Redis
docker-compose exec neo4j cypher-shell -u neo4j -p papergraph123 "RETURN 1"  # Neo4j

# Reset environment (complete cleanup)
docker-compose down -v
docker-compose up --build

# Rebuild single service after code changes
docker-compose up --build orchestrator
```

### Swapping LLM Providers
Change `LLM_PROVIDER` environment variable to switch between Ollama and Bedrock without code changes. The `agents/shared/llm_service.py` provides the abstraction layer.

## Performance Targets

- Query processing: <2 minutes end-to-end
- Paper fetching: <30 seconds for 20 papers
- Entity extraction: <10 seconds per paper
- Neo4j queries: <5 seconds for coverage analysis

## Architecture Decisions

- **Docker Compose**: Simple multi-container orchestration for development
- **Python 3.11**: Team expertise and rich scientific library ecosystem
- **py2neo**: Pythonic Neo4j client with good abstraction
- **Biopython**: Standard library for PubMed integration
- **spaCy**: Fast NLP with pre-trained models for entity recognition
- **FastAPI**: Modern async API framework with auto-documentation
- **UV**: Fast Python package manager for all services
- **Multi-container**: Each service has its own pyproject.toml and Dockerfile
- **Shared Module**: `agents/shared/` contains common code (llm_service.py, neo4j_client.py) copied into agent containers at build time
- Use the .ideas directory to read or write documentations related to implementation and roadmap