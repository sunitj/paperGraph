# PaperGraph Intelligence

Multi-agent literature analysis system with knowledge graph learning. Uses AI agents to discover research gaps and guide scientific discovery through intelligent paper analysis.

## 🚀 Quick Start

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

3. **Start services**
   ```bash
   docker-compose up -d
   ```

4. **Pull LLM models** (required for first setup)
   ```bash
   # Wait for Ollama to start (~30 seconds)
   docker exec papergraph-ollama-1 ollama pull llama3.1:8b
   docker exec papergraph-ollama-1 ollama pull mistral:7b
   ```

5. **Access interfaces**
   - **Streamlit UI**: http://localhost:8501
   - **API Documentation**: http://localhost:8000/docs
   - **Neo4j Browser**: http://localhost:7474 (neo4j/papergraph123)

## 📋 System Architecture

### Services
- **Ollama**: Local LLM service (port 11434)
- **Neo4j**: Knowledge graph database (ports 7474, 7687)
- **Redis**: Message queue and state management (port 6379)
- **Orchestrator**: FastAPI coordination service (port 8000)
- **Agents**: Specialized processing containers
  - Fetcher: PubMed paper retrieval
  - KG Builder: Entity and relationship extraction
  - Query Strategist: Knowledge gap analysis
- **UI**: Streamlit interface (port 8501)

### Data Flow
```
User Question → Query Strategist → Fetcher Agent → KG Builder → Orchestrator → UI
```

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

1. **Ollama not responding**
   ```bash
   docker-compose restart ollama
   docker exec papergraph-ollama-1 ollama list
   ```

2. **Neo4j authentication failed**
   ```bash
   docker-compose down -v  # Remove volumes
   docker-compose up -d neo4j
   ```

3. **Services won't start**
   ```bash
   docker-compose down
   docker system prune -f
   docker-compose up --build
   ```

4. **Port conflicts**
   Edit `docker-compose.yml` to change port mappings

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

## 🤝 Contributing

This is the MVP implementation (PR1: Project Scaffolding). Future PRs will add:
- Real LLM integration (PR2)
- Neo4j schema and operations (PR3)
- PubMed integration (PR4)
- Entity extraction (PR5)
- Gap analysis (PR6)
- Complete orchestration (PR7)
- Full UI implementation (PR8)

## 📄 License

[License information to be added]

## 🆘 Support

For issues and questions:
- Check the troubleshooting section above
- Review service logs with `docker-compose logs`
- Ensure all environment variables are set correctly