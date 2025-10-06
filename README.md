# PaperGraph Intelligence

Multi-agent literature analysis system with knowledge graph learning. Uses AI agents to discover research gaps and guide scientific discovery through intelligent paper analysis.

## 🚀 Quick Start

> **Current Status**: MVP implementation with working infrastructure and stub agent implementations. Ollama LLM integration and Neo4j knowledge graph are functional. Agent logic uses placeholder implementations for rapid development.

### Prerequisites
- Docker and Docker Compose
- [UV](https://docs.astral.sh/uv/) for local Python development (optional)

### Setup

1. **Clone and configure**
   ```bash
   git clone <repository-url>
   cd paperGraph
   cp .env.example .env
   ```

2. **Configure environment**
   Edit `.env` and add your email (required for PubMed):
   ```bash
   PUBMED_EMAIL=your_email@example.com
   ```

3. **Install Ollama (if not already installed)**
   ```bash
   # The system uses host Ollama installation, not containerized
   # Download from https://ollama.ai or use:
   curl -fsSL https://ollama.com/install.sh | sh

   # Pull required models
   ollama pull llama3.1:8b
   ollama pull mistral:7b
   ```

4. **Start services**
   ```bash
   docker-compose up -d
   ```

5. **Access interfaces**
   - **Streamlit UI**: http://localhost:8501
   - **API Documentation**: http://localhost:8000/docs
   - **Neo4j Browser**: http://localhost:7474 (neo4j/papergraph123)

## 📋 System Architecture

### Services
- **Ollama**: Local LLM service on host (port 11434) - services connect via `host.docker.internal`
- **Neo4j**: Knowledge graph database (ports 7474, 7687)
- **Redis**: Message queue and state management (port 6379)
- **Orchestrator**: FastAPI coordination service (port 8000)
- **Agents**: Specialized processing containers
  - Fetcher: PubMed paper retrieval
  - KG Builder: Entity and relationship extraction
  - Query Strategist: Knowledge gap analysis
- **UI**: Streamlit interface (port 8501)

**Note:** Ollama runs on the host machine, not in a container. Services access it via `host.docker.internal:11434`.

### Data Flow
```
User Question → Query Strategist → Fetcher Agent → KG Builder → Orchestrator → UI
```

## 🔬 PubMed API Setup

The fetcher agent uses NCBI's E-utilities API to retrieve papers. **You must provide an email address** (required by NCBI policy).

### Rate Limits
- **Without API key:** 3 requests per second
- **With API key:** 10 requests per second (recommended for heavy use)

### Get an API Key (Optional)
1. Create NCBI account: https://www.ncbi.nlm.nih.gov/account/
2. Generate API key in account settings
3. Add to `.env`: `PUBMED_API_KEY=your_key_here`

### Configure
```bash
# .env
PUBMED_EMAIL=your_email@example.com  # Required
PUBMED_API_KEY=your_api_key          # Optional but recommended
```

---

## 🛠️ Development

### Local Development with UV

Install dependencies for a specific service:
```bash
cd orchestrator
uv sync
uv run uvicorn main:app --reload
```

### Using Docker for Development

```bash
# Build and start specific service
docker-compose up --build orchestrator

# View logs
docker-compose logs -f fetcher

# Execute commands in container
docker-compose exec orchestrator bash
```

### Testing

```bash
# Check all services are running
docker-compose ps

# Test service health
curl http://localhost:8000/health
curl http://localhost:7474
curl http://localhost:11434/api/tags

# Test fetcher agent with real PubMed query
python scripts/test_fetcher.py

# View logs
docker-compose logs -f
```

## 📊 Usage

### 1. Load Pre-seeded Knowledge (Optional)
Available topics:
- Gut Microbiome & Obesity
- Protein Engineering
- Cancer Immunotherapy
- CRISPR Gene Editing
- Neurodegeneration
- Antibiotic Resistance

### 2. Ask Questions
Examples:
- "How do gut bacteria influence obesity?"
- "What are the mechanisms of CRISPR gene editing?"
- "How does machine learning improve protein design?"

### 3. Explore Results
- **Synthesized Answer**: Evidence-based response with citations
- **Knowledge Gaps**: Identified areas needing more research
- **Papers**: Source papers with abstracts and links
- **Graph**: Network visualization of entities and relationships

## 🔧 Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `PUBMED_EMAIL` | Your email (required by NCBI) | - |
| `LLM_PROVIDER` | "ollama" or "bedrock" | ollama |
| `LLM_MODEL` | Model name | llama3.1:8b |
| `MAX_PAPERS_PER_QUERY` | Papers per search | 20 |

### Performance Tuning

- **LLM Models**: llama3.1:8b (faster) vs mistral:7b (lighter)
- **Memory**: Minimum 8GB RAM recommended
- **Storage**: ~10GB for models and data

## 🏗️ Project Structure

```
papergraph/
├── docker-compose.yml          # Service orchestration
├── orchestrator/               # FastAPI coordination service
│   ├── pyproject.toml         # UV dependencies
│   ├── Dockerfile             # UV-optimized container
│   └── main.py                # FastAPI application
├── agents/
│   ├── shared/                # Common utilities
│   │   ├── llm_service.py     # LLM abstraction layer
│   │   └── neo4j_client.py    # Knowledge graph client
│   ├── fetcher/               # PubMed paper retrieval
│   ├── kg_builder/            # Entity extraction
│   └── query_strategist/      # Gap analysis
├── ui/                        # Streamlit interface
├── seeds/                     # Pre-built knowledge graphs
└── scripts/                   # Setup utilities
```

## 🚨 Troubleshooting

### Common Issues

1. **Services not running**
   ```bash
   # Check if services are up
   docker-compose ps

   # If not running, start them
   docker-compose up -d
   ```

2. **Ollama not responding**
   ```bash
   docker-compose restart ollama
   docker-compose exec ollama ollama list
   ```

3. **Neo4j authentication failed**
   ```bash
   # Verify .env has NEO4J_PASSWORD=papergraph123
   # If changed, reset Neo4j
   docker-compose down -v  # Remove volumes
   docker-compose up -d neo4j
   ```

4. **Services won't start**
   ```bash
   docker-compose down
   docker system prune -f
   docker-compose up --build
   ```

5. **Port conflicts**
   ```bash
   # Check what's using the ports
   lsof -i :8000  # Orchestrator
   lsof -i :8501  # Streamlit
   lsof -i :7474  # Neo4j
   lsof -i :11434 # Ollama

   # Either stop conflicting services or edit docker-compose.yml port mappings
   ```

6. **Missing UV lock files**
   ```bash
   # Lock files are generated on first sync
   cd orchestrator
   uv sync  # Creates uv.lock
   ```

### Logs and Debugging

```bash
# View all logs
docker-compose logs

# Follow specific service
docker-compose logs -f orchestrator

# Check container status
docker-compose ps
```

## 📚 API Documentation

Interactive API documentation available at:
- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

Key endpoints:
- `GET /health` - System health check
- `GET /stats` - Knowledge graph statistics
- `POST /query` - Process research questions

## 📊 Implementation Status

### ✅ Completed (MVP)
- **Infrastructure**: Docker Compose orchestration with all services
- **Ollama Integration**: Local LLM service with llama3.1:8b and mistral:7b
- **Neo4j Setup**: Knowledge graph database with schema
- **Redis**: Message queue for agent communication
- **Orchestrator**: FastAPI service with health checks and stats endpoints
- **UI Framework**: Streamlit interface with system monitoring
- **Shared Modules**: LLM abstraction layer and Neo4j client
- **Fetcher Agent**: Real PubMed integration with Biopython E-utilities API, rate limiting, error handling

### 🚧 Stub Implementations (Ready for Enhancement)
- **KG Builder Agent**: Placeholder entity extraction (needs spaCy + LLM prompts)
- **Query Strategist**: Placeholder gap analysis (needs query generation logic)
- **End-to-End Flow**: Basic orchestration (needs full workflow implementation)

### 🔜 Future Enhancements
- Real PubMed API integration with Biopython
- LLM-powered entity and relationship extraction
- Knowledge graph coverage analysis
- Intelligent query generation based on gaps
- Pre-seeded knowledge graphs for common domains
- Advanced visualization and citation tracking

## 🤝 Contributing

See implementation status above for development priorities.

## 📄 License

[License information to be added]

## 🆘 Support

For issues and questions:
- Check the troubleshooting section above
- Review service logs with `docker-compose logs`
- Ensure all environment variables are set correctly