# PaperGraph Intelligence - PR Roadmap Status

**Last Updated:** October 4, 2024
**Current Sprint:** Week 2

---

## PR Completion Status

### ✅ PR #1: Project Scaffolding (COMPLETE)
**Status:** Complete
**Completed:** September 29, 2024
**Commit:** 3decadf, 45c3ea8

**Delivered:**
- Docker Compose multi-service architecture
- All service directories with Dockerfiles and pyproject.toml
- UV-based dependency management
- README with setup instructions
- .env.example configuration template

**Acceptance Criteria:** All AC-1.x met
- Services start successfully
- Docker Compose orchestration working
- Documentation complete

**Notes:**
- Ollama runs on host instead of containerized (practical decision for port conflict)
- Used UV instead of requirements.txt for modern Python dependency management

---

### ✅ PR #2: Ollama Integration & LLM Service (COMPLETE)
**Status:** Complete
**Completed:** September 29, 2024
**Commit:** 3decadf

**Delivered:**
- LLM abstraction layer (`agents/shared/llm_service.py`)
- OllamaProvider with health checks and text generation
- BedrockProvider stub for future AWS support
- Factory pattern for provider selection
- Integration with llama3.1:8b and mistral:7b models

**Acceptance Criteria:** All AC-2.x met
- LLM generates text responses
- Health check implemented
- Provider switching via environment variable

---

### ✅ PR #3: Neo4j Setup & Basic Schema (COMPLETE)
**Status:** Complete
**Completed:** September 29, 2024
**Commit:** 3decadf

**Delivered:**
- KnowledgeGraph client class (`agents/shared/neo4j_client.py`)
- Schema with unique constraints (Paper.pmid, Organism.name, Disease.name, Molecule.name)
- Indexes for performance optimization
- CRUD operations (create_paper, create_entity, create_relationship)
- Query methods (get_entity_coverage, get_subgraph, mark_paper_analyzed)
- init_schema() for database initialization

**Acceptance Criteria:** All AC-3.x met
- Schema constraints created
- All CRUD operations working
- Query methods functional

---

### ✅ PR #4: Fetcher Agent - PubMed Integration (COMPLETE)
**Status:** Complete
**Completed:** October 4, 2024
**Commit:** [Pending]

**Delivered:**
- Real PubMed integration with Biopython E-utilities API
- Entrez.esearch() for searching PubMed by query
- Entrez.efetch() for retrieving paper details
- MEDLINE/XML parsing for PMID, title, abstract, year
- Rate limiting (3 req/sec without API key, 10 req/sec with key)
- Error handling for PubMed API failures
- Test script (`scripts/test_fetcher.py`)
- Documentation updates (README, CLAUDE.md)

**Acceptance Criteria:** All AC-4.x met
- ✅ AC-4.1: Agent starts and connects
- ✅ AC-4.2: Fetches real papers from PubMed
- ✅ AC-4.3: Test script created
- ✅ AC-4.7: Uses Biopython Entrez
- ✅ AC-4.8: Respects rate limits
- ✅ AC-4.9: Parses all required fields

**Testing:**
- Manual testing via Docker Compose
- Integration test script provided
- Rate limiting validated

**Notes:**
- Replaced stub implementation returning fake papers (STUB000001, etc.)
- Now fetches real papers from NCBI PubMed database
- Handles structured abstracts and date parsing
- Ready for PR #5 (KG Builder can now process real papers)

---

## 🚧 In Progress

### ⏳ PR #5: KG Builder Agent - Entity Extraction
**Status:** Stub Implementation
**Target:** Week 2, Days 4-5

**Current State:**
- Agent structure complete
- Stub entity extraction
- Needs: spaCy NER integration
- Needs: LLM-based relationship extraction

**Dependencies:**
- ✅ PR #4 complete (real papers available)

**Next Steps:**
1. Implement spaCy named entity recognition
2. Design LLM prompts for relationship extraction
3. Parse entities into Organism, Disease, Molecule types
4. Extract relationships (PRODUCES, AFFECTS, INHIBITS, CORRELATES_WITH)
5. Store in Neo4j with confidence scores

---

### ⏳ PR #6: Query Strategist Agent - Gap Analysis
**Status:** Stub Implementation
**Target:** Week 2, Days 6-7

**Current State:**
- Agent structure complete
- Stub gap analysis
- Needs: Knowledge graph coverage analysis
- Needs: Query generation logic

**Dependencies:**
- ⏳ PR #5 in progress (needs real knowledge graph)

---

### ⏳ PR #7: Orchestrator - Agent Coordination
**Status:** Partial Implementation
**Target:** Week 3, Days 1-2

**Current State:**
- ✅ FastAPI application structure
- ✅ Health endpoint
- ✅ Stats endpoint
- ❌ `/query` endpoint (main workflow)
- ❌ `/load_seed` endpoint
- ❌ Agent coordination logic

**Dependencies:**
- ⏳ PR #5 and PR #6 in progress

---

### ⏳ PR #8: Streamlit UI - User Interface
**Status:** Stub Implementation
**Target:** Week 3, Days 3-4

**Current State:**
- ✅ Streamlit app structure
- ❌ Query input interface
- ❌ Results display
- ❌ Graph visualization

**Dependencies:**
- ⏳ PR #7 in progress (needs `/query` endpoint)

---

### ⏳ PR #9: KG Seeding System
**Status:** Not Started
**Target:** Week 3, Days 5-7

**Current State:**
- Empty `seeds/` directory
- Needs: 6 pre-seeded topic JSON files
- Needs: `scripts/prepare_seeds.py`

**Dependencies:**
- None (independent feature, can start anytime)

---

## 📋 Next Actions

### Immediate (Next PR)
**PR #5: KG Builder Agent - Entity Extraction**

**Estimated Duration:** 3-4 days

**Implementation Steps:**
1. Add spaCy integration to `agents/kg_builder/agent.py`
2. Download and configure `en_core_web_sm` model
3. Design LLM prompts for relationship extraction
4. Implement entity classification (Organism, Disease, Molecule)
5. Extract relationships with confidence scores
6. Store in Neo4j with paper citations
7. Create test script
8. Update documentation

**Why Next:**
- Real papers now available from PR #4
- Core value proposition of the system
- Enables knowledge graph building
- Unblocks PR #6 (Query Strategist needs real graph)

---

### Subsequent PRs (In Order)

**PR #6: Query Strategist** (Week 2, Days 6-7)
- Gap analysis with real knowledge graph
- Query generation logic
- **Depends on:** PR #5

**PR #7: Orchestrator** (Week 3, Days 1-2)
- `/query` endpoint implementation
- Agent coordination via Redis
- Answer synthesis
- **Depends on:** PR #5, PR #6

**PR #8: Streamlit UI** (Week 3, Days 3-4)
- Query interface
- Results visualization
- Graph rendering
- **Depends on:** PR #7

**PR #9: Seeding System** (Week 3, Days 5-7)
- Can be done in parallel
- 6 pre-seeded topics
- Independent of other PRs

---

## 🎯 Sprint Goals

### Week 2 (Current)
- [x] PR #4: Real PubMed Integration
- [ ] PR #5: Entity Extraction with spaCy + LLM
- [ ] PR #6: Gap Analysis and Query Generation

### Week 3
- [ ] PR #7: Orchestrator Coordination
- [ ] PR #8: Streamlit UI
- [ ] PR #9: Seeding System

### Week 4 (Integration & Polish)
- [ ] End-to-end testing
- [ ] Performance optimization
- [ ] Documentation review
- [ ] Deployment preparation

---

## 📝 Notes

### Key Decisions Made
1. **Ollama on Host**: Decided to run Ollama on host machine instead of containerized to avoid port conflicts and simplify setup
2. **UV for Dependencies**: Chose UV over pip/requirements.txt for faster, more reliable dependency management
3. **Real PubMed First**: Prioritized real data fetching before entity extraction to enable better testing

### Lessons Learned
1. Docker build contexts need to be at project root for shared module copying
2. UV sync without lock files works better for rapid iteration
3. Host Ollama via `host.docker.internal` works well for development

### Technical Debt
- None significant at this stage
- All stub implementations are clearly marked and tracked
- Documentation kept up to date with actual implementation

---

## 🚀 MVP Readiness

**Current Completion:** ~40%

**Functional PRs:** 4/9 (44%)
**Stub PRs:** 5/9 (56%)

**Path to MVP:**
- PR #5 (Entity Extraction): +20% → 60% complete
- PR #6 (Query Strategist): +15% → 75% complete
- PR #7 (Orchestrator): +15% → 90% complete
- PR #8 (UI): +5% → 95% complete
- PR #9 (Seeding): +5% → 100% complete

**Estimated Completion:** 2-3 weeks from current date
