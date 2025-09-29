# PaperGraph Intelligence - API Endpoint Specifications

**Version:** 1.0.0  
**Base URL:** `http://localhost:8000`  
**API Type:** REST  
**Authentication:** None (MVP)  
**Content-Type:** application/json

---

## Table of Contents

1. [Health & Status](#health--status)
2. [Knowledge Base Management](#knowledge-base-management)
3. [Query & Search](#query--search)
4. [System Information](#system-information)
5. [Error Codes](#error-codes)
6. [Rate Limits](#rate-limits)

---

## Health & Status

### GET /health

Health check endpoint to verify service status.

**Description**: Returns the operational status of all system components including LLM service, Redis, and Neo4j.

#### Request

```http
GET /health HTTP/1.1
Host: localhost:8000
```

#### Response

**Status Code**: `200 OK`

```json
{
  "status": "healthy",
  "llm_provider": "ollama",
  "llm_status": "healthy",
  "redis": "connected",
  "neo4j": "connected",
  "timestamp": "2024-09-27T10:30:00Z"
}
```

#### Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | Overall system status: "healthy", "degraded", "unhealthy" |
| `llm_provider` | string | Active LLM provider: "ollama" or "bedrock" |
| `llm_status` | string | LLM service status: "healthy", "unhealthy", "unknown" |
| `redis` | string | Redis connection status: "connected", "disconnected" |
| `neo4j` | string | Neo4j connection status: "connected", "disconnected" |
| `timestamp` | string | ISO 8601 timestamp of health check |

#### Example

```bash
curl http://localhost:8000/health
```

```json
{
  "status": "healthy",
  "llm_provider": "ollama",
  "llm_status": "healthy",
  "redis": "connected",
  "neo4j": "connected"
}
```

#### Error Responses

**Status Code**: `503 Service Unavailable`

```json
{
  "status": "unhealthy",
  "error": "Neo4j connection failed",
  "details": "Connection refused"
}
```

---

## Knowledge Base Management

### GET /seeds

List all available seed knowledge graphs.

**Description**: Returns metadata about available pre-built knowledge graphs that can be loaded.

#### Request

```http
GET /seeds HTTP/1.1
Host: localhost:8000
```

#### Response

**Status Code**: `200 OK`

```json
{
  "seeds": [
    {
      "name": "Gut Microbiome Obesity",
      "description": "Explores how gut bacteria influence weight, metabolism, and obesity through various mechanisms including SCFA production and inflammation.",
      "papers": 30,
      "filename": "gut_microbiome_obesity.json"
    },
    {
      "name": "Protein Engineering",
      "description": "Rational design and directed evolution of proteins for therapeutics, enzymes, and biosensors including machine learning approaches.",
      "papers": 28,
      "filename": "protein_engineering.json"
    }
  ],
  "total": 6
}
```

#### Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `seeds` | array | Array of seed objects |
| `seeds[].name` | string | Display name of the seed topic |
| `seeds[].description` | string | Brief description of the research area |
| `seeds[].papers` | integer | Number of papers in the seed |
| `seeds[].filename` | string | Seed file identifier |
| `total` | integer | Total number of available seeds |

#### Example

```bash
curl http://localhost:8000/seeds | jq
```

---

### POST /load_seed

Load a pre-built seed knowledge graph into the system.

**Description**: Loads papers and relationships from a seed file into Neo4j. This populates the knowledge graph with existing research for faster querying.

#### Request

```http
POST /load_seed HTTP/1.1
Host: localhost:8000
Content-Type: application/json

{
  "seed": "Protein Engineering"
}
```

#### Request Body

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `seed` | string | Yes | Name of the seed to load (must match name from /seeds) |

#### Response

**Status Code**: `200 OK`

```json
{
  "status": "success",
  "seed": "Protein Engineering",
  "papers_count": 28,
  "relationships_count": 145,
  "entities_created": 67,
  "load_time_seconds": 23.4,
  "message": "Seed loaded successfully"
}
```

#### Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | Operation status: "success" or "error" |
| `seed` | string | Name of the loaded seed |
| `papers_count` | integer | Number of papers loaded |
| `relationships_count` | integer | Number of relationships created |
| `entities_created` | integer | Number of entities (organisms, diseases, molecules) created |
| `load_time_seconds` | float | Time taken to load seed |
| `message` | string | Human-readable status message |

#### Example

```bash
curl -X POST http://localhost:8000/load_seed \
  -H "Content-Type: application/json" \
  -d '{"seed": "Protein Engineering"}'
```

#### Error Responses

**Status Code**: `404 Not Found` - Seed doesn't exist

```json
{
  "detail": "Seed 'Invalid Name' not found",
  "available_seeds": ["Gut Microbiome Obesity", "Protein Engineering", ...]
}
```

**Status Code**: `500 Internal Server Error` - Loading failed

```json
{
  "detail": "Failed to load seed",
  "error": "Neo4j connection timeout"
}
```

---

## Query & Search

### POST /query

Process a natural language research question.

**Description**: Main endpoint for querying the system. Analyzes the knowledge graph, identifies gaps, fetches relevant papers, extracts entities, and synthesizes an answer.

#### Request

```http
POST /query HTTP/1.1
Host: localhost:8000
Content-Type: application/json

{
  "question": "How do gut bacteria influence obesity?",
  "max_papers": 20
}
```

#### Request Body

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `question` | string | Yes | Natural language research question (10-500 chars) |
| `max_papers` | integer | No | Maximum papers to fetch (default: 20, range: 5-50) |

#### Response

**Status Code**: `200 OK`

```json
{
  "question_id": "a3f8d2c1",
  "answer": "Gut bacteria influence obesity through multiple mechanisms. Firmicutes bacteria produce short-chain fatty acids (SCFAs) [1,2] which affect energy metabolism and inflammation. Studies show that the Firmicutes to Bacteroidetes ratio is elevated in obese individuals [3]. Akkermansia muciniphila has been identified as protective against obesity by strengthening the gut barrier [4].\n\nRecent research indicates that SCFAs can modulate insulin resistance and adipocyte function [5]. However, the exact mechanisms linking specific bacterial species to metabolic outcomes remain partially understood, particularly for understudied species like Akkermansia.",
  "papers": [
    {
      "pmid": "35789123",
      "title": "Firmicutes produce SCFAs that modulate metabolic pathways",
      "year": 2023,
      "abstract": "Recent studies demonstrate...",
      "url": "https://pubmed.ncbi.nlm.nih.gov/35789123"
    }
  ],
  "gaps": [
    {
      "type": "weak_evidence",
      "subject": "Akkermansia",
      "predicate": "AFFECTS",
      "object": "obesity",
      "confidence": 0.35,
      "priority": "high",
      "description": "Limited evidence (confidence: 0.35) for Akkermansia's effect on obesity"
    },
    {
      "type": "understudied",
      "subject": "diet",
      "predicate": "MODULATES",
      "object": "microbiome composition",
      "paper_count": 2,
      "priority": "medium",
      "description": "Only 2 papers cover diet's effect on microbiome"
    },
    {
      "type": "missing_entity",
      "entity": "Bifidobacterium",
      "priority": "high",
      "description": "No information found about Bifidobacterium in current knowledge base"
    }
  ],
  "entities": [
    "gut bacteria",
    "Firmicutes",
    "Akkermansia",
    "obesity",
    "short-chain fatty acids"
  ],
  "subgraph": {
    "nodes": [
      {"id": "Firmicutes", "type": "Organism", "label": "Firmicutes"},
      {"id": "SCFA", "type": "Molecule", "label": "short-chain fatty acids"},
      {"id": "obesity", "type": "Disease", "label": "obesity"}
    ],
    "edges": [
      {
        "source": "Firmicutes",
        "target": "SCFA",
        "type": "PRODUCES",
        "confidence": 0.85,
        "papers": 15
      },
      {
        "source": "SCFA",
        "target": "obesity",
        "type": "AFFECTS",
        "confidence": 0.78,
        "papers": 8
      }
    ]
  },
  "metadata": {
    "processing_time_seconds": 87.3,
    "papers_fetched": 18,
    "papers_analyzed": 18,
    "new_entities_created": 12,
    "new_relationships_created": 24,
    "queries_generated": 3
  }
}
```

#### Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `question_id` | string | Unique identifier for this query |
| `answer` | string | Synthesized answer with citations [1], [2] |
| `papers` | array | Array of paper objects analyzed |
| `papers[].pmid` | string | PubMed ID |
| `papers[].title` | string | Paper title |
| `papers[].year` | integer | Publication year |
| `papers[].abstract` | string | Paper abstract (truncated) |
| `papers[].url` | string | PubMed URL |
| `gaps` | array | Identified knowledge gaps |
| `gaps[].type` | string | Gap type: "weak_evidence", "understudied", "missing_entity" |
| `gaps[].priority` | string | Priority: "high", "medium", "low" |
| `gaps[].description` | string | Human-readable gap description |
| `entities` | array | Entities extracted from question |
| `subgraph` | object | Knowledge graph subset relevant to question |
| `subgraph.nodes` | array | Entity nodes |
| `subgraph.edges` | array | Relationships between entities |
| `metadata` | object | Processing metadata and statistics |

#### Example

```bash
curl -X POST http://localhost:8000/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "How does CRISPR gene editing work?",
    "max_papers": 15
  }' | jq
```

#### Error Responses

**Status Code**: `408 Request Timeout` - Processing timeout

```json
{
  "detail": "Query processing timeout (>120 seconds)",
  "question_id": "a3f8d2c1",
  "partial_results": {
    "analysis_complete": true,
    "papers_fetched": 0
  }
}
```

**Status Code**: `400 Bad Request` - Invalid request

```json
{
  "detail": "Validation error",
  "errors": [
    {
      "field": "question",
      "message": "Question must be between 10 and 500 characters"
    }
  ]
}
```

**Status Code**: `429 Too Many Requests` - Rate limit exceeded

```json
{
  "detail": "Rate limit exceeded",
  "retry_after": 60,
  "message": "Maximum 10 queries per minute"
}
```

---

## System Information

### GET /stats

Get knowledge graph statistics.

**Description**: Returns current statistics about the knowledge graph including paper count, entity count, and relationship count.

#### Request

```http
GET /stats HTTP/1.1
Host: localhost:8000
```

#### Response

**Status Code**: `200 OK`

```json
{
  "papers": 248,
  "entities": 892,
  "relationships": 1547,
  "entity_breakdown": {
    "Organism": 342,
    "Disease": 178,
    "Molecule": 372
  },
  "relationship_breakdown": {
    "PRODUCES": 456,
    "AFFECTS": 523,
    "INHIBITS": 234,
    "CORRELATES_WITH": 334
  },
  "average_confidence": 0.73,
  "pre_seeded_papers": 180,
  "user_added_papers": 68,
  "last_updated": "2024-09-27T15:45:00Z"
}
```

#### Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `papers` | integer | Total papers in knowledge graph |
| `entities` | integer | Total entities (organisms, diseases, molecules) |
| `relationships` | integer | Total relationships between entities |
| `entity_breakdown` | object | Count by entity type |
| `relationship_breakdown` | object | Count by relationship type |
| `average_confidence` | float | Average confidence score (0.0-1.0) |
| `pre_seeded_papers` | integer | Papers from loaded seeds |
| `user_added_papers` | integer | Papers added via queries |
| `last_updated` | string | ISO 8601 timestamp of last update |

#### Example

```bash
curl http://localhost:8000/stats | jq
```

---

### GET /query/{question_id}

Get details of a previous query.

**Description**: Retrieve full details of a previously processed query including analysis and results.

#### Request

```http
GET /query/a3f8d2c1 HTTP/1.1
Host: localhost:8000
```

#### Path Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `question_id` | string | Question ID from original query response |

#### Response

**Status Code**: `200 OK`

```json
{
  "question_id": "a3f8d2c1",
  "question": "How do gut bacteria influence obesity?",
  "timestamp": "2024-09-27T14:30:00Z",
  "status": "complete",
  "answer": "Gut bacteria influence obesity through...",
  "papers": [...],
  "gaps": [...],
  "entities": [...],
  "subgraph": {...},
  "metadata": {...}
}
```

#### Error Responses

**Status Code**: `404 Not Found`

```json
{
  "detail": "Query 'a3f8d2c1' not found or expired",
  "message": "Query results expire after 1 hour"
}
```

---

### GET /entities/{entity_name}

Get information about a specific entity.

**Description**: Retrieve all relationships and papers mentioning a specific entity.

#### Request

```http
GET /entities/Firmicutes HTTP/1.1
Host: localhost:8000
```

#### Path Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `entity_name` | string | Name of the entity (URL encoded) |

#### Response

**Status Code**: `200 OK`

```json
{
  "name": "Firmicutes",
  "type": "Organism",
  "properties": {
    "ncbi_id": "1239"
  },
  "relationships": [
    {
      "predicate": "PRODUCES",
      "target": "short-chain fatty acids",
      "target_type": "Molecule",
      "confidence": 0.85,
      "paper_count": 15
    },
    {
      "predicate": "AFFECTS",
      "target": "obesity",
      "target_type": "Disease",
      "confidence": 0.72,
      "paper_count": 8
    }
  ],
  "papers": [
    {
      "pmid": "35789123",
      "title": "Firmicutes and metabolic health",
      "year": 2023
    }
  ],
  "statistics": {
    "total_relationships": 12,
    "total_papers": 23,
    "average_confidence": 0.78
  }
}
```

#### Error Responses

**Status Code**: `404 Not Found`

```json
{
  "detail": "Entity 'InvalidName' not found in knowledge graph"
}
```

---

## Error Codes

### Standard Error Response Format

All error responses follow this structure:

```json
{
  "detail": "Human-readable error message",
  "error_code": "ERROR_CODE",
  "timestamp": "2024-09-27T10:30:00Z",
  "path": "/query",
  "request_id": "abc123"
}
```

### Error Code Reference

| HTTP Status | Error Code | Description |
|-------------|-----------|-------------|
| 400 | INVALID_REQUEST | Request validation failed |
| 400 | INVALID_QUESTION | Question format invalid |
| 400 | INVALID_SEED_NAME | Seed name not found |
| 404 | NOT_FOUND | Resource not found |
| 404 | SEED_NOT_FOUND | Specified seed doesn't exist |
| 404 | QUERY_NOT_FOUND | Query ID not found or expired |
| 404 | ENTITY_NOT_FOUND | Entity not in knowledge graph |
| 408 | TIMEOUT | Operation timeout |
| 408 | QUERY_TIMEOUT | Query processing exceeded time limit |
| 429 | RATE_LIMIT_EXCEEDED | Too many requests |
| 500 | INTERNAL_ERROR | Server error |
| 500 | LLM_ERROR | LLM service error |
| 500 | NEO4J_ERROR | Database error |
| 500 | REDIS_ERROR | Cache/queue error |
| 503 | SERVICE_UNAVAILABLE | Service temporarily unavailable |
| 503 | LLM_UNAVAILABLE | LLM service not responding |

---

## Rate Limits

### Current Limits (MVP)

| Endpoint | Limit | Window |
|----------|-------|--------|
| `/query` | 10 requests | 1 minute |
| `/load_seed` | 5 requests | 5 minutes |
| `/stats` | 60 requests | 1 minute |
| `/health` | Unlimited | - |
| `/seeds` | 60 requests | 1 minute |

### Rate Limit Headers

Responses include rate limit information:

```http
HTTP/1.1 200 OK
X-RateLimit-Limit: 10
X-RateLimit-Remaining: 7
X-RateLimit-Reset: 1695825600
```

### Rate Limit Exceeded Response

**Status Code**: `429 Too Many Requests`

```json
{
  "detail": "Rate limit exceeded",
  "error_code": "RATE_LIMIT_EXCEEDED",
  "limit": 10,
  "window": "1 minute",
  "retry_after": 45,
  "message": "Please wait 45 seconds before retrying"
}
```

---

## Request/Response Examples

### Complete Query Flow

```bash
# 1. Check system health
curl http://localhost:8000/health

# 2. List available seeds
curl http://localhost:8000/seeds

# 3. Load a seed
curl -X POST http://localhost:8000/load_seed \
  -H "Content-Type: application/json" \
  -d '{"seed": "Protein Engineering"}'

# 4. Check statistics
curl http://localhost:8000/stats

# 5. Ask a question
curl -X POST http://localhost:8000/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "How is machine learning used in protein design?",
    "max_papers": 20
  }' > response.json

# 6. Extract question_id
QUESTION_ID=$(jq -r '.question_id' response.json)

# 7. Get query details later
curl http://localhost:8000/query/$QUESTION_ID

# 8. Explore specific entity
curl http://localhost:8000/entities/directed%20evolution
```

### Python Client Example

```python
import requests

class PaperGraphClient:
    def __init__(self, base_url="http://localhost:8000"):
        self.base_url = base_url
    
    def health(self):
        return requests.get(f"{self.base_url}/health").json()
    
    def list_seeds(self):
        return requests.get(f"{self.base_url}/seeds").json()
    
    def load_seed(self, seed_name):
        return requests.post(
            f"{self.base_url}/load_seed",
            json={"seed": seed_name}
        ).json()
    
    def query(self, question, max_papers=20):
        return requests.post(
            f"{self.base_url}/query",
            json={"question": question, "max_papers": max_papers}
        ).json()
    
    def get_stats(self):
        return requests.get(f"{self.base_url}/stats").json()
    
    def get_entity(self, entity_name):
        return requests.get(
            f"{self.base_url}/entities/{entity_name}"
        ).json()

# Usage
client = PaperGraphClient()

# Check health
print(client.health())

# Load seed
client.load_seed("Protein Engineering")

# Query
result = client.query(
    "How does directed evolution improve enzymes?",
    max_papers=15
)

print(f"Answer: {result['answer']}")
print(f"Papers: {len(result['papers'])}")
print(f"Gaps: {len(result['gaps'])}")
```

### JavaScript Client Example

```javascript
class PaperGraphClient {
  constructor(baseUrl = 'http://localhost:8000') {
    this.baseUrl = baseUrl;
  }

  async health() {
    const response = await fetch(`${this.baseUrl}/health`);
    return response.json();
  }

  async listSeeds() {
    const response = await fetch(`${this.baseUrl}/seeds`);
    return response.json();
  }

  async loadSeed(seedName) {
    const response = await fetch(`${this.baseUrl}/load_seed`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ seed: seedName })
    });
    return response.json();
  }

  async query(question, maxPapers = 20) {
    const response = await fetch(`${this.baseUrl}/query`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ question, max_papers: maxPapers })
    });
    return response.json();
  }

  async getStats() {
    const response = await fetch(`${this.baseUrl}/stats`);
    return response.json();
  }
}

// Usage
const client = new PaperGraphClient();

(async () => {
  // Check health
  const health = await client.health();
  console.log('System status:', health.status);

  // Load seed
  await client.loadSeed('Protein Engineering');

  // Query
  const result = await client.query(
    'How is machine learning used in protein design?'
  );
  
  console.log('Answer:', result.answer);
  console.log('Papers:', result.papers.length);
  console.log('Gaps:', result.gaps.length);
})();
```

---

## API Versioning

### Current Version: v1

Base URL includes version: `http://localhost:8000/v1` (optional in MVP)

### Version Header

```http
API-Version: 1.0.0
```

### Deprecation Notice

Deprecated endpoints include warning header:

```http
X-API-Deprecated: true
X-API-Deprecation-Date: 2024-12-31
X-API-Sunset-Date: 2025-03-31
```

---

## WebSocket Support (Future)

Real-time query updates via WebSocket (planned for Phase 2):

```javascript
const ws = new WebSocket('ws://localhost:8000/ws/query');

ws.send(JSON.stringify({
  action: 'query',
  question: 'How do gut bacteria affect obesity?'
}));

ws.onmessage = (event) => {
  const update = JSON.parse(event.data);
  console.log('Status:', update.status);
  // Updates: "analyzing", "fetching", "extracting", "synthesizing", "complete"
};
```

---

## OpenAPI/Swagger Specification

Full OpenAPI 3.0 specification available at:
- Interactive docs: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc
- JSON spec: http://localhost:8000/openapi.json

---

## Support & Feedback

**Issues**: GitHub Issues  
**Email**: support@papergraph.example.com  
**Documentation**: https://docs.papergraph.example.com  
**Status Page**: https://status.papergraph.example.com

---

**Last Updated**: September 27, 2024  
**API Version**: 1.0.0  
**Document Version**: 1.0
