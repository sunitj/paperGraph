# PR #1: Project Scaffolding

**Sprint:** Week 1, Days 1-2  
**Goal:** Set up Docker environment and basic project structure

## What to Build

- Docker Compose with all services
- Basic directory structure
- Requirements.txt for each component
- README with setup instructions

## Deliverable

```bash
docker-compose up -d
# All services start successfully
```

## Directory Structure

```
papergraph/
├── docker-compose.yml
├── README.md
├── .env.example
├── orchestrator/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── main.py
├── agents/
│   ├── shared/
│   │   ├── __init__.py
│   │   ├── llm_service.py
│   │   └── neo4j_client.py
│   ├── fetcher/
│   │   ├── Dockerfile
│   │   ├── requirements.txt
│   │   └── agent.py
│   ├── kg_builder/
│   │   ├── Dockerfile
│   │   ├── requirements.txt
│   │   └── agent.py
│   └── query_strategist/
│       ├── Dockerfile
│       ├── requirements.txt
│       └── agent.py
├── ui/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── app.py
├── seeds/
│   └── .gitkeep
└── scripts/
    └── prepare_seeds.py
```

## docker-compose.yml

```yaml
version: '3.8'

services:
  ollama:
    image: ollama/ollama:latest
    ports:
      - "11434:11434"
    volumes:
      - ollama_data:/root/.ollama
    command: serve

  neo4j:
    image: neo4j:5.15-community
    ports:
      - "7474:7474"
      - "7687:7687"
    environment:
      - NEO4J_AUTH=neo4j/papergraph123
      - NEO4J_PLUGINS=["apoc"]
    volumes:
      - neo4j_data:/data

  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"
    volumes:
      - redis_data:/data

  orchestrator:
    build: ./orchestrator
    ports:
      - "8000:8000"
    environment:
      - NEO4J_URI=bolt://neo4j:7687
      - NEO4J_USER=neo4j
      - NEO4J_PASSWORD=papergraph123
      - REDIS_URL=redis://redis:6379
      - OLLAMA_URL=http://ollama:11434
    depends_on:
      - neo4j
      - redis
      - ollama

  fetcher:
    build:
      context: .
      dockerfile: agents/fetcher/Dockerfile
    environment:
      - REDIS_URL=redis://redis:6379
      - NEO4J_URI=bolt://neo4j:7687
      - NEO4J_USER=neo4j
      - NEO4J_PASSWORD=papergraph123
      - PUBMED_EMAIL=${PUBMED_EMAIL}
    depends_on:
      - redis
      - neo4j

  kg_builder:
    build:
      context: .
      dockerfile: agents/kg_builder/Dockerfile
    environment:
      - REDIS_URL=redis://redis:6379
      - NEO4J_URI=bolt://neo4j:7687
      - NEO4J_USER=neo4j
      - NEO4J_PASSWORD=papergraph123
      - OLLAMA_URL=http://ollama:11434
    depends_on:
      - redis
      - neo4j
      - ollama

  query_strategist:
    build:
      context: .
      dockerfile: agents/query_strategist/Dockerfile
    environment:
      - REDIS_URL=redis://redis:6379
      - NEO4J_URI=bolt://neo4j:7687
      - NEO4J_USER=neo4j
      - NEO4J_PASSWORD=papergraph123
      - OLLAMA_URL=http://ollama:11434
    depends_on:
      - redis
      - neo4j
      - ollama

  ui:
    build: ./ui
    ports:
      - "8501:8501"
    environment:
      - ORCHESTRATOR_URL=http://orchestrator:8000
    depends_on:
      - orchestrator

volumes:
  ollama_data:
  neo4j_data:
  redis_data:
```

## Basic Requirements Files

### orchestrator/requirements.txt
```
fastapi==0.104.1
uvicorn==0.24.0
redis==5.0.1
py2neo==2021.2.3
pydantic==2.5.0
```

### agents/fetcher/requirements.txt
```
biopython==1.81
redis==5.0.1
py2neo==2021.2.3
```

### agents/kg_builder/requirements.txt
```
spacy==3.7.2
redis==5.0.1
py2neo==2021.2.3
requests==2.31.0
```

### agents/query_strategist/requirements.txt
```
redis==5.0.1
py2neo==2021.2.3
requests==2.31.0
```

### ui/requirements.txt
```
streamlit==1.28.0
requests==2.31.0
py2neo==2021.2.3
plotly==5.17.0
```

## Dockerfiles

### orchestrator/Dockerfile
```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY orchestrator/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY orchestrator/ .
COPY agents/shared/ /app/shared/

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
```

### agents/fetcher/Dockerfile
```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY agents/fetcher/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY agents/fetcher/ .
COPY agents/shared/ /app/shared/

CMD ["python", "agent.py"]
```

### agents/kg_builder/Dockerfile
```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY agents/kg_builder/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN python -m spacy download en_core_web_sm

COPY agents/kg_builder/ .
COPY agents/shared/ /app/shared/

CMD ["python", "agent.py"]
```

### agents/query_strategist/Dockerfile
```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY agents/query_strategist/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY agents/query_strategist/ .
COPY agents/shared/ /app/shared/

CMD ["python", "agent.py"]
```

### ui/Dockerfile
```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY ui/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY ui/ .

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
```

## README.md

```markdown
# PaperGraph Intelligence

Multi-agent literature analysis system with knowledge graph learning.

## Quick Start

1. Clone repository
2. Copy `.env.example` to `.env` and add your email:
   ```
   PUBMED_EMAIL=your_email@example.com
   ```

3. Start services:
   ```bash
   docker-compose up -d
   ```

4. Pull LLM model:
   ```bash
   docker exec -it papergraph-ollama-1 ollama pull llama3.1:8b
   ```

5. Access UI:
   - Streamlit: http://localhost:8501
   - Neo4j Browser: http://localhost:7474

## Services

- **Ollama**: Local LLM service (port 11434)
- **Neo4j**: Knowledge graph (ports 7474, 7687)
- **Redis**: State management (port 6379)
- **Orchestrator**: API gateway (port 8000)
- **Agents**: Fetcher, KG Builder, Query Strategist
- **UI**: Streamlit interface (port 8501)
```

## Testing

```bash
# Verify all services are running
docker-compose ps

# Check logs
docker-compose logs -f

# Test Neo4j connection
curl http://localhost:7474

# Test Ollama
curl http://localhost:11434/api/tags
```

## Acceptance Criteria

- [ ] All services start without errors
- [ ] Neo4j browser accessible
- [ ] Ollama responds to health check
- [ ] Redis accepts connections
- [ ] README documentation complete
