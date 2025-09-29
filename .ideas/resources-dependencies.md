# PaperGraph Intelligence - Resources & Dependencies

**Version:** 1.0  
**Last Updated:** September 2024  
**Package Manager:** UV + Docker

---

## Table of Contents

1. [Python Package Manager - UV](#python-package-manager---uv)
2. [Core Python Dependencies](#core-python-dependencies)
3. [LLM & AI Libraries](#llm--ai-libraries)
4. [Database Clients](#database-clients)
5. [Docker Images & Services](#docker-images--services)
6. [Development Tools](#development-tools)
7. [UV Project Setup](#uv-project-setup)
8. [Quick Reference](#quick-reference)

---

## Python Package Manager - UV

### What is UV?

UV is a modern, extremely fast Python package and project manager written in Rust. It's 10-100x faster than pip and handles virtual environments automatically.

**Official Resources:**
- 📚 **Documentation**: https://docs.astral.sh/uv/
- 🐙 **GitHub**: https://github.com/astral-sh/uv
- 📦 **PyPI**: https://pypi.org/project/uv/
- 🎓 **Getting Started**: https://docs.astral.sh/uv/getting-started/

### Installation

```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows (PowerShell)
irm https://astral.sh/uv/install.ps1 | iex

# With pip (if needed)
pip install uv

# Verify installation
uv --version
```

### Key Features

- ⚡ **10-100x faster** than pip
- 🔒 **Automatic virtual environment** management
- 📦 **Lock files** for reproducible builds (uv.lock)
- 🔄 **Drop-in pip replacement** (supports pyproject.toml)
- 🐳 **Works great with Docker**

### Basic Commands

```bash
# Create new project
uv init papergraph

# Add dependency
uv add fastapi

# Add dev dependency
uv add --dev pytest

# Install all dependencies
uv sync

# Run command in venv
uv run python script.py

# Update dependencies
uv lock --upgrade
```

---

## Core Python Dependencies

### Web Framework

#### FastAPI
**Purpose**: Orchestrator API, async endpoints  
**Version**: 0.104.1+  
**Resources**:
- 📚 Docs: https://fastapi.tiangolo.com/
- 🐙 GitHub: https://github.com/tiangolo/fastapi
- 🎓 Tutorial: https://fastapi.tiangolo.com/tutorial/
- 📖 Advanced: https://fastapi.tiangolo.com/advanced/

```bash
uv add "fastapi[all]>=0.104.1"
```

**Key Features for This Project**:
- Async/await support (crucial for agent coordination)
- Auto-generated OpenAPI docs
- Pydantic validation
- WebSocket support (future)

#### Uvicorn
**Purpose**: ASGI server for FastAPI  
**Version**: 0.24.0+  
**Resources**:
- 📚 Docs: https://www.uvicorn.org/
- 🐙 GitHub: https://github.com/encode/uvicorn

```bash
uv add "uvicorn[standard]>=0.24.0"
```

### UI Framework

#### Streamlit
**Purpose**: User interface  
**Version**: 1.28.0+  
**Resources**:
- 📚 Docs: https://docs.streamlit.io/
- 🐙 GitHub: https://github.com/streamlit/streamlit
- 🎓 Gallery: https://streamlit.io/gallery
- 🎨 Components: https://streamlit.io/components

```bash
uv add streamlit>=1.28.0
```

**Key Streamlit Resources**:
- Session State: https://docs.streamlit.io/library/api-reference/session-state
- Caching: https://docs.streamlit.io/library/advanced-features/caching
- Custom Components: https://docs.streamlit.io/library/components

---

## LLM & AI Libraries

### Ollama (Local LLM)

**Purpose**: Local LLM inference  
**Resources**:
- 🌐 Website: https://ollama.com/
- 📚 Docs: https://github.com/ollama/ollama/tree/main/docs
- 🐙 GitHub: https://github.com/ollama/ollama
- 🐋 Docker: https://hub.docker.com/r/ollama/ollama
- 📦 Models: https://ollama.com/library
- 🔌 API: https://github.com/ollama/ollama/blob/main/docs/api.md

**Python Client**:
```bash
uv add ollama
```

**Docker Usage**:
```yaml
services:
  ollama:
    image: ollama/ollama:latest
    ports:
      - "11434:11434"
    volumes:
      - ollama_data:/root/.ollama
```

**Model Management**:
```bash
# In Docker
docker exec -it ollama ollama pull llama3.1:8b
docker exec -it ollama ollama pull mistral:7b
docker exec -it ollama ollama list

# Direct
ollama pull llama3.1:8b
ollama list
ollama run llama3.1:8b
```

**API Examples**: https://github.com/ollama/ollama/blob/main/docs/api.md

### spaCy (NLP)

**Purpose**: Entity extraction  
**Version**: 3.7.2+  
**Resources**:
- 📚 Docs: https://spacy.io/
- 🐙 GitHub: https://github.com/explosion/spaCy
- 📖 Models: https://spacy.io/models
- 🎓 Course: https://course.spacy.io/

```bash
uv add spacy>=3.7.2

# Download models
uv run python -m spacy download en_core_web_sm
uv run python -m spacy download en_core_web_md  # Better accuracy
```

**Model Comparison**: https://spacy.io/models/en#en_core_web_sm

### Biopython

**Purpose**: PubMed API access  
**Version**: 1.81+  
**Resources**:
- 📚 Docs: https://biopython.org/
- 📖 Tutorial: https://biopython.org/DIST/docs/tutorial/Tutorial.html
- 🐙 GitHub: https://github.com/biopython/biopython
- 🔬 Entrez: https://biopython.org/docs/1.81/api/Bio.Entrez.html

```bash
uv add biopython>=1.81
```

**Entrez API Guide**: https://www.ncbi.nlm.nih.gov/books/NBK25501/  
**E-utilities**: https://www.ncbi.nlm.nih.gov/books/NBK25500/

---

## Database Clients

### Neo4j (Knowledge Graph)

**Purpose**: Graph database for knowledge storage  
**Version**: 5.15+ (Community)  
**Resources**:
- 🌐 Website: https://neo4j.com/
- 📚 Docs: https://neo4j.com/docs/
- 🐋 Docker: https://hub.docker.com/_/neo4j
- 🎓 GraphAcademy: https://graphacademy.neo4j.com/
- 📖 Cypher: https://neo4j.com/docs/cypher-manual/current/

**Python Driver (py2neo)**:
```bash
uv add py2neo>=2021.2.3
```

**Resources**:
- py2neo Docs: https://py2neo.org/
- py2neo GitHub: https://github.com/py2neo-org/py2neo

**Docker Usage**:
```yaml
services:
  neo4j:
    image: neo4j:5.15-community
    ports:
      - "7474:7474"  # Browser
      - "7687:7687"  # Bolt
    environment:
      - NEO4J_AUTH=neo4j/papergraph123
      - NEO4J_PLUGINS=["apoc"]
      - NEO4J_server_memory_heap_initial__size=512m
      - NEO4J_server_memory_heap_max__size=2G
    volumes:
      - neo4j_data:/data
```

**Important Neo4j Docker Resources**:
- Docker Guide: https://neo4j.com/docs/operations-manual/current/docker/
- Configuration: https://neo4j.com/docs/operations-manual/current/configuration/
- APOC Plugin: https://neo4j.com/docs/apoc/current/

**Cypher Query Language**:
- Cheat Sheet: https://neo4j.com/docs/cypher-cheat-sheet/
- Best Practices: https://neo4j.com/developer/cypher-style-guide/

### Redis (State Management)

**Purpose**: Message queue, caching  
**Version**: 7.x  
**Resources**:
- 🌐 Website: https://redis.io/
- 📚 Docs: https://redis.io/docs/
- 🐋 Docker: https://hub.docker.com/_/redis
- 🎓 University: https://university.redis.com/
- 📖 Commands: https://redis.io/commands/

**Python Client**:
```bash
uv add redis>=5.0.1
```

**Resources**:
- redis-py Docs: https://redis-py.readthedocs.io/
- redis-py GitHub: https://github.com/redis/redis-py

**Docker Usage**:
```yaml
services:
  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"
    volumes:
      - redis_data:/data
    command: redis-server --appendonly yes --maxmemory 256mb --maxmemory-policy allkeys-lru
```

**Important Redis Docker Resources**:
- Docker Guide: https://redis.io/docs/getting-started/install-stack/docker/
- Configuration: https://redis.io/docs/management/config/
- Persistence: https://redis.io/docs/management/persistence/

**Redis Patterns for This Project**:
- Pub/Sub: https://redis.io/docs/manual/pubsub/
- Lists (Queues): https://redis.io/docs/data-types/lists/
- Key Expiration: https://redis.io/commands/expire/

---

## Docker Images & Services

### Official Docker Documentation

**Core Resources**:
- 📚 Docs: https://docs.docker.com/
- 🐋 Hub: https://hub.docker.com/
- 📖 Compose: https://docs.docker.com/compose/
- 🎓 Get Started: https://docs.docker.com/get-started/

### Docker Compose

**Purpose**: Multi-container orchestration  
**Version**: 3.8+  
**Resources**:
- Compose File: https://docs.docker.com/compose/compose-file/
- Networking: https://docs.docker.com/compose/networking/
- Volumes: https://docs.docker.com/storage/volumes/
- Best Practices: https://docs.docker.com/develop/dev-best-practices/

**Key Resources for This Project**:
- Health Checks: https://docs.docker.com/engine/reference/builder/#healthcheck
- Depends On: https://docs.docker.com/compose/startup-order/
- Environment Variables: https://docs.docker.com/compose/environment-variables/

### Project-Specific Images

#### Python Base
```dockerfile
FROM python:3.11-slim
```
- Image: https://hub.docker.com/_/python
- Tags: https://hub.docker.com/_/python/tags

#### Ollama
```yaml
ollama:
  image: ollama/ollama:latest
```
- Image: https://hub.docker.com/r/ollama/ollama
- Docs: https://github.com/ollama/ollama/blob/main/docs/docker.md

#### Neo4j
```yaml
neo4j:
  image: neo4j:5.15-community
```
- Image: https://hub.docker.com/_/neo4j
- Docs: https://neo4j.com/docs/operations-manual/current/docker/

#### Redis
```yaml
redis:
  image: redis:7-alpine
```
- Image: https://hub.docker.com/_/redis
- Docs: https://redis.io/docs/getting-started/install-stack/docker/

---

## Development Tools

### Data Visualization

#### Plotly
**Purpose**: Graph visualization in UI  
**Version**: 5.17.0+  
**Resources**:
- 📚 Docs: https://plotly.com/python/
- 🐙 GitHub: https://github.com/plotly/plotly.py
- 📊 Examples: https://plotly.com/python/basic-charts/

```bash
uv add plotly>=5.17.0
```

### HTTP Client

#### Requests
**Purpose**: HTTP requests, API calls  
**Version**: 2.31.0+  
**Resources**:
- 📚 Docs: https://requests.readthedocs.io/
- 🐙 GitHub: https://github.com/psf/requests

```bash
uv add requests>=2.31.0
```

### Validation

#### Pydantic
**Purpose**: Data validation (included with FastAPI)  
**Version**: 2.5.0+  
**Resources**:
- 📚 Docs: https://docs.pydantic.dev/
- 🐙 GitHub: https://github.com/pydantic/pydantic

```bash
uv add pydantic>=2.5.0
```

---

## UV Project Setup

### Converting from requirements.txt to UV

#### Step 1: Install UV

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

#### Step 2: Create pyproject.toml

```bash
# Initialize UV project
uv init papergraph
cd papergraph

# Or convert existing project
uv init
```

#### Step 3: Add Dependencies

**From requirements.txt**:
```bash
# Add each dependency
cat requirements.txt | xargs -I {} uv add {}

# Or in bulk (if formatted correctly)
uv add $(cat requirements.txt)
```

**Or manually**:
```bash
# Core dependencies
uv add "fastapi[all]>=0.104.1"
uv add "uvicorn[standard]>=0.24.0"
uv add streamlit>=1.28.0
uv add redis>=5.0.1
uv add py2neo>=2021.2.3
uv add biopython>=1.81
uv add spacy>=3.7.2
uv add plotly>=5.17.0
uv add requests>=2.31.0
uv add pydantic>=2.5.0

# Dev dependencies
uv add --dev pytest>=7.4.0
uv add --dev black>=23.0.0
uv add --dev ruff>=0.1.0
```

### Project Structure with UV

```
papergraph/
├── pyproject.toml          # Project config + dependencies
├── uv.lock                 # Lock file (like package-lock.json)
├── .python-version         # Python version
├── orchestrator/
│   ├── main.py
│   └── ...
├── agents/
│   ├── shared/
│   ├── fetcher/
│   ├── kg_builder/
│   └── query_strategist/
└── ui/
    └── app.py
```

### pyproject.toml Example

```toml
[project]
name = "papergraph"
version = "1.0.0"
description = "Multi-agent literature analysis with knowledge graph"
readme = "README.md"
requires-python = ">=3.11"
dependencies = [
    "fastapi[all]>=0.104.1",
    "uvicorn[standard]>=0.24.0",
    "streamlit>=1.28.0",
    "redis>=5.0.1",
    "py2neo>=2021.2.3",
    "biopython>=1.81",
    "spacy>=3.7.2",
    "plotly>=5.17.0",
    "requests>=2.31.0",
    "pydantic>=2.5.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.4.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
]

[tool.uv]
dev-dependencies = [
    "pytest>=7.4.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
```

### Docker Integration with UV

#### Dockerfile with UV

```dockerfile
FROM python:3.11-slim

WORKDIR /app

# Install UV
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Copy dependency files
COPY pyproject.toml uv.lock ./

# Install dependencies
RUN uv sync --frozen --no-cache

# Copy application code
COPY . .

# Run with UV
CMD ["uv", "run", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
```

#### Alternative: Pre-built venv

```dockerfile
FROM python:3.11-slim

WORKDIR /app

# Install UV
RUN pip install uv

# Copy files
COPY pyproject.toml uv.lock ./
COPY . .

# Create venv and install
RUN uv venv
RUN uv sync --frozen

# Activate venv and run
ENV PATH="/app/.venv/bin:$PATH"
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
```

### UV Commands for Development

```bash
# Install dependencies
uv sync

# Add new dependency
uv add <package>

# Add dev dependency
uv add --dev <package>

# Update dependencies
uv lock --upgrade

# Run Python scripts
uv run python script.py

# Run commands in venv
uv run pytest
uv run black .
uv run uvicorn main:app --reload

# Show installed packages
uv pip list

# Export to requirements.txt (if needed)
uv pip compile pyproject.toml -o requirements.txt
```

---

## Quick Reference

### Essential Commands

```bash
# UV Commands
uv sync                    # Install all dependencies
uv add <package>           # Add dependency
uv add --dev <package>     # Add dev dependency
uv run <command>           # Run in venv
uv lock --upgrade          # Update lock file

# Docker Commands
docker-compose up -d       # Start all services
docker-compose ps          # Check status
docker-compose logs -f     # View logs
docker-compose down        # Stop services
docker-compose restart     # Restart services

# Neo4j Commands
docker exec -it neo4j cypher-shell -u neo4j -p password
MATCH (n) RETURN count(n); # Count nodes

# Redis Commands
docker exec -it redis redis-cli
PING                       # Test connection
KEYS *                     # List all keys
GET key                    # Get value

# Ollama Commands
docker exec -it ollama ollama list        # List models
docker exec -it ollama ollama pull <model> # Download model
```

### Port Reference

| Service | Port | URL |
|---------|------|-----|
| Orchestrator API | 8000 | http://localhost:8000 |
| Streamlit UI | 8501 | http://localhost:8501 |
| Neo4j Browser | 7474 | http://localhost:7474 |
| Neo4j Bolt | 7687 | bolt://localhost:7687 |
| Redis | 6379 | redis://localhost:6379 |
| Ollama API | 11434 | http://localhost:11434 |

### Environment Variables

```bash
# Database
NEO4J_URI=bolt://neo4j:7687
NEO4J_USER=neo4j
NEO4J_PASSWORD=papergraph123
REDIS_URL=redis://redis:6379

# LLM
OLLAMA_URL=http://ollama:11434
LLM_PROVIDER=ollama
LLM_MODEL=llama3.1:8b

# PubMed
PUBMED_EMAIL=your_email@example.com
PUBMED_API_KEY=optional_api_key

# Application
MAX_PAPERS_PER_QUERY=20
```

### Debugging Tips

**Neo4j Not Starting:**
```bash
# Check logs
docker-compose logs neo4j

# Common fix: Clear data
docker-compose down -v
docker-compose up -d neo4j
```

**Ollama Model Issues:**
```bash
# Verify models
docker exec -it ollama ollama list

# Re-download
docker exec -it ollama ollama pull llama3.1:8b
```

**Redis Connection:**
```bash
# Test connection
docker exec -it redis redis-cli ping

# Check memory
docker exec -it redis redis-cli info memory
```

**UV Issues:**
```bash
# Clear cache
uv cache clean

# Reinstall
rm -rf .venv uv.lock
uv sync
```

---

## Learning Resources

### Tutorials & Courses

**Neo4j**:
- Neo4j Fundamentals: https://graphacademy.neo4j.com/courses/neo4j-fundamentals/
- Cypher Fundamentals: https://graphacademy.neo4j.com/courses/cypher-fundamentals/
- Building Neo4j Applications with Python: https://graphacademy.neo4j.com/courses/app-python/

**FastAPI**:
- Full Tutorial: https://fastapi.tiangolo.com/tutorial/
- Advanced User Guide: https://fastapi.tiangolo.com/advanced/
- Real Python FastAPI: https://realpython.com/fastapi-python-web-apis/

**Docker**:
- Docker Get Started: https://docs.docker.com/get-started/
- Docker Compose Tutorial: https://docs.docker.com/compose/gettingstarted/
- Best Practices: https://docs.docker.com/develop/dev-best-practices/

**Streamlit**:
- 30 Days of Streamlit: https://30days.streamlit.app/
- Advanced Features: https://docs.streamlit.io/library/advanced-features

### Community & Support

- **Neo4j Community**: https://community.neo4j.com/
- **FastAPI Discord**: https://discord.gg/VQjSZaeJmf
- **Streamlit Forum**: https://discuss.streamlit.io/
- **Docker Forums**: https://forums.docker.com/
- **UV Discussions**: https://github.com/astral-sh/uv/discussions

---

## Version Compatibility Matrix

| Component | Minimum Version | Tested Version | Python Version |
|-----------|----------------|----------------|----------------|
| Python | 3.11 | 3.11 | - |
| UV | 0.1.0 | latest | - |
| FastAPI | 0.104.1 | 0.104.1 | 3.11+ |
| Streamlit | 1.28.0 | 1.28.0 | 3.8-3.12 |
| Neo4j | 5.15 | 5.15 | - |
| Redis | 7.0 | 7.2 | - |
| Ollama | - | latest | - |
| spaCy | 3.7.2 | 3.7.2 | 3.7-3.12 |
| Biopython | 1.81 | 1.81 | 3.7+ |

---

## Migration Checklist: pip → UV

- [ ] Install UV: `curl -LsSf https://astral.sh/uv/install.sh | sh`
- [ ] Create pyproject.toml: `uv init`
- [ ] Add dependencies: `uv add <packages>`
- [ ] Update Dockerfiles to use UV
- [ ] Test build: `uv sync`
- [ ] Update CI/CD to use UV
- [ ] Update documentation
- [ ] Train team on UV commands
- [ ] Verify lock file in version control
- [ ] Remove old requirements.txt (optional)

---

## Additional Resources

### Official Links

- **PubMed E-utilities**: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **Ollama Models**: https://ollama.com/library
- **Neo4j APOC**: https://neo4j.com/docs/apoc/current/
- **Redis Commands**: https://redis.io/commands/
- **FastAPI Best Practices**: https://github.com/zhanymkanov/fastapi-best-practices
- **Docker Multi-stage Builds**: https://docs.docker.com/build/building/multi-stage/

### Monitoring & Debugging

- **Neo4j Query Profiling**: https://neo4j.com/docs/cypher-manual/current/query-tuning/
- **Redis Monitoring**: https://redis.io/docs/management/optimization/
- **FastAPI Debugging**: https://fastapi.tiangolo.com/tutorial/debugging/
- **Docker Logs**: https://docs.docker.com/config/containers/logging/

---

**Last Updated**: September 27, 2024  
**Maintainer**: Development Team  
**UV Version**: latest  
**Docker Compose Version**: 3.8
