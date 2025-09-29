# PaperGraph Intelligence - Product Requirements Document

**Version:** 1.0  
**Date:** September 2024  
**Status:** MVP Development  
**Owner:** Sunit Jain

---

## Executive Summary

PaperGraph Intelligence is a multi-agent literature analysis system that learns from a knowledge graph to identify research gaps and guide scientific discovery. Unlike traditional literature search tools that treat each query independently, PaperGraph builds persistent knowledge that improves with each paper analyzed, enabling intelligent gap-driven research.

### Vision
Enable researchers to discover novel connections in scientific literature through AI agents that collaboratively build and query a living knowledge graph.

### Target Users
- Computational biologists and bioinformatics researchers
- PhD students and postdocs conducting literature reviews
- R&D teams in biotech/pharma companies
- Research scientists exploring new domains

---

## Problem Statement

### Current State
Researchers face significant challenges in literature analysis:
1. **Information Overload**: 3+ million papers published annually in PubMed alone
2. **Isolated Queries**: Each search starts from scratch with no learning
3. **Missing Connections**: Novel relationships buried across disparate papers
4. **Knowledge Gaps**: No systematic way to identify understudied areas
5. **Manual Synthesis**: Hours spent connecting findings across papers

### Impact
- Missed research opportunities and novel connections
- Duplicated effort across research teams
- Slow hypothesis generation
- Difficulty identifying research gaps for grant proposals

---

## Goals and Objectives

### Primary Goals
1. **Enable Gap-Driven Discovery**: Automatically identify and fill knowledge gaps in research literature
2. **Build Persistent Knowledge**: Create a learning system that improves with each paper analyzed
3. **Accelerate Research**: Reduce literature review time from days to minutes
4. **Surface Novel Connections**: Find relationships across papers that manual review would miss

### Success Metrics (3 months post-launch)
- **Utility**: 50+ researchers use the system weekly
- **Learning**: Knowledge graph grows to 1000+ papers, 5000+ relationships
- **Discovery**: Users report 5+ novel connections found
- **Efficiency**: 70% reduction in literature review time (self-reported)
- **Engagement**: 60% weekly return rate

### MVP Success Criteria (Launch)
- ✅ System analyzes papers and builds knowledge graph
- ✅ Gap analysis identifies weak evidence and missing links
- ✅ 6 pre-seeded research topics available
- ✅ End-to-end query flow works in <2 minutes
- ✅ Graph visualization shows knowledge network
- ✅ System learns from each paper processed

---

## User Personas

### Persona 1: Dr. Sarah Chen - Computational Biologist
**Background**: 3 years postdoc, researching gut microbiome  
**Goals**: Find novel bacterial-disease connections, publish high-impact papers  
**Pain Points**: Overwhelmed by literature, missing cross-domain connections  
**Use Case**: "I need to know what's NOT known about Akkermansia's role in obesity"

### Persona 2: Alex Kumar - PhD Student
**Background**: 2nd year PhD in protein engineering  
**Goals**: Complete literature review for dissertation, identify research gap  
**Pain Points**: Manual review takes weeks, hard to identify gaps for thesis  
**Use Case**: "I need to find understudied areas in ML-guided protein design"

### Persona 3: Dr. James Liu - Pharma R&D Scientist
**Background**: 10 years industry, leading cancer immunotherapy research  
**Goals**: Validate hypotheses, find competitive intelligence, identify targets  
**Pain Points**: Need comprehensive view quickly, track emerging mechanisms  
**Use Case**: "Show me all known PD-1 resistance mechanisms and what's missing"

---

## User Stories

### Epic 1: Knowledge Base Setup
- **US-1.1**: As a researcher, I want to load a pre-seeded topic so I can start with existing knowledge
- **US-1.2**: As a user, I want to see system statistics so I know the knowledge base size
- **US-1.3**: As a researcher, I want to explore the knowledge graph so I can understand relationships

### Epic 2: Intelligent Querying
- **US-2.1**: As a researcher, I want to ask natural language questions so I don't need to learn query syntax
- **US-2.2**: As a user, I want the system to identify knowledge gaps so I can find research opportunities
- **US-2.3**: As a scientist, I want targeted papers fetched so I get relevant recent research
- **US-2.4**: As a researcher, I want to see which entities are understudied so I can prioritize investigations

### Epic 3: Answer Synthesis
- **US-3.1**: As a user, I want synthesized answers with citations so I can trust the information
- **US-3.2**: As a researcher, I want to see confidence scores so I know evidence strength
- **US-3.3**: As a scientist, I want to access original papers so I can verify claims
- **US-3.4**: As a user, I want gap identification explained so I understand what's missing

### Epic 4: Knowledge Building
- **US-4.1**: As a researcher, I want papers automatically added to the graph so knowledge persists
- **US-4.2**: As a user, I want entity extraction automated so I don't do manual work
- **US-4.3**: As a scientist, I want relationship confidence tracked so I know evidence quality
- **US-4.4**: As a researcher, I want the system to learn from my queries so it gets smarter

---

## Features and Requirements

### F1: Pre-Seeded Knowledge Graphs
**Priority**: P0 (Must Have)

**Requirements**:
- 6 research topics pre-loaded (Gut Microbiome, Cancer Immunotherapy, Antibiotic Resistance, CRISPR, Neurodegeneration, Protein Engineering)
- 20-30 papers per topic with extracted entities and relationships
- One-click loading via UI
- Topics cover diverse biological domains

**Acceptance**:
- Seed files exist for all 6 topics
- Load time <30 seconds per seed
- Knowledge graph visible in Neo4j after loading
- UI shows paper count and description

### F2: Natural Language Query Interface
**Priority**: P0 (Must Have)

**Requirements**:
- Text input for questions
- No special syntax required
- Support for biological domain questions
- Real-time feedback during processing

**Acceptance**:
- Users can type free-form questions
- System extracts entities from questions
- Processing status visible to user
- Works with questions about mechanisms, relationships, evidence

### F3: Knowledge Gap Analysis
**Priority**: P0 (Must Have)

**Requirements**:
- Identify weak evidence (confidence <0.5)
- Find understudied relationships (papers <3)
- Detect missing entities
- Prioritize gaps by severity

**Acceptance**:
- System queries existing knowledge graph
- Gaps categorized by type (weak/missing/understudied)
- Gaps displayed with priority levels
- Gap analysis completes in <10 seconds

### F4: Intelligent Paper Fetching
**Priority**: P0 (Must Have)

**Requirements**:
- Generate targeted PubMed queries based on gaps
- Fetch 10-20 relevant papers
- Filter by recency and relevance
- Handle PubMed rate limits

**Acceptance**:
- Queries are specific to identified gaps
- Papers retrieved within 30 seconds
- Abstracts captured for analysis
- No duplicate papers fetched

### F5: Entity and Relationship Extraction
**Priority**: P0 (Must Have)

**Requirements**:
- Extract organisms, diseases, molecules
- Identify relationships (produces, affects, inhibits)
- Assign confidence scores
- Create knowledge graph nodes/edges

**Acceptance**:
- Extracts 5-10 entities per paper
- Identifies 3-5 relationships per paper
- Entities correctly classified by type
- Relationships have confidence scores

### F6: Answer Synthesis
**Priority**: P0 (Must Have)

**Requirements**:
- Generate coherent answers using LLM
- Cite specific papers
- Incorporate knowledge graph context
- Note gaps and uncertainties

**Acceptance**:
- Answers are 2-3 paragraphs
- Papers cited with [number] format
- Answers incorporate both new papers and existing KG
- Gaps mentioned in response

### F7: Graph Visualization
**Priority**: P1 (Should Have)

**Requirements**:
- Display knowledge network
- Show entity relationships
- Interactive exploration
- Link to Neo4j browser

**Acceptance**:
- Graph renders in UI
- Nodes and edges visible
- Neo4j query examples provided
- Link to full browser works

### F8: Multi-Agent Coordination
**Priority**: P0 (Must Have)

**Requirements**:
- Query Strategist analyzes gaps
- Fetcher retrieves papers
- KG Builder extracts entities
- Orchestrator coordinates flow

**Acceptance**:
- Agents communicate via Redis
- Each agent processes independently
- State shared via Redis/Neo4j
- Full pipeline completes in <2 minutes

---

## Technical Architecture

### System Components

**Core Services**:
1. **Ollama** - Local LLM service (llama3.1:8b)
2. **Neo4j** - Knowledge graph database
3. **Redis** - State management and message queue
4. **Orchestrator** - FastAPI coordination service
5. **Agents** - Specialized processing containers
6. **Streamlit UI** - User interface

**Agent Architecture**:
```
User Question
    ↓
Query Strategist → Analyzes KG gaps → Generates queries
    ↓
Fetcher Agent → Retrieves papers → Stores in Neo4j
    ↓
KG Builder Agent → Extracts entities → Updates graph
    ↓
Orchestrator → Synthesizes answer → Returns to UI
```

**Data Flow**:
1. User submits question via UI
2. Orchestrator sends to Query Strategist queue
3. Strategist analyzes KG, identifies gaps, generates queries
4. Queries sent to Fetcher queue
5. Fetcher retrieves papers, stores in Neo4j, marks for analysis
6. KG Builder processes papers, extracts entities/relationships
7. Orchestrator collects results, synthesizes answer
8. UI displays answer, papers, gaps, graph

### Knowledge Graph Schema

**Nodes**:
- `Paper {pmid, title, abstract, year, analyzed}`
- `Organism {name, ncbi_id}`
- `Disease {name, mesh_id}`
- `Molecule {name, chebi_id}`

**Relationships**:
- `-[:PRODUCES {confidence, paper_pmid}]->`
- `-[:AFFECTS {confidence, paper_pmid}]->`
- `-[:INHIBITS {confidence, paper_pmid}]->`
- `-[:CORRELATES_WITH {effect_size, p_value, paper_pmid}]->`

### Technology Stack

| Component | Technology | Rationale |
|-----------|-----------|-----------|
| **LLM** | Ollama (llama3.1:8b) | Free, local, fast; easy Bedrock swap |
| **Knowledge Graph** | Neo4j Community | Best graph database, free tier |
| **Message Queue** | Redis | Simple, fast, pub/sub support |
| **API** | FastAPI | Modern, async, auto-docs |
| **UI** | Streamlit | Fastest prototyping |
| **Orchestration** | Docker Compose | Simple multi-container |
| **Language** | Python 3.11 | Team expertise, library ecosystem |

---

## User Experience Flow

### Primary Flow: Ask a Question

1. **Load Knowledge Base** (optional)
   - User selects pre-seeded topic from sidebar
   - Clicks "Load Seed"
   - System loads 20-30 papers into graph
   - Stats update showing paper count

2. **Ask Question**
   - User enters natural language question
   - Example: "How do gut bacteria influence obesity?"
   - Clicks "Search" button

3. **Processing** (60-90 seconds)
   - Progress indicator shows: "Analyzing knowledge graph..."
   - System identifies gaps in background
   - Fetches targeted papers
   - Extracts entities and relationships
   - Updates knowledge graph

4. **Results Display**
   - **Answer Section**: 2-3 paragraph synthesis with citations
   - **Gaps Section**: Color-coded gaps (weak evidence, missing entities)
   - **Papers Section**: Expandable cards with abstracts and PubMed links
   - **Graph Section**: Network visualization of entities and relationships

5. **Follow-up** (optional)
   - User asks related question
   - System uses expanded knowledge graph
   - Faster response due to cached knowledge

---

## Out of Scope (Post-MVP)

### Phase 2 Features
- PDF parsing and full-text analysis
- Real-time paper monitoring and alerts
- Multi-user support and shared workspaces
- Advanced graph algorithms (PageRank, community detection)
- Custom entity types and ontology integration
- Citation network analysis
- Export functionality (reports, presentations)

### Phase 3 Features
- Browser extension for in-context analysis
- API for programmatic access
- Collaboration features (annotations, sharing)
- Advanced visualizations (temporal analysis, topic evolution)
- Integration with reference managers (Zotero, Mendeley)
- Fine-tuned biological NER models

### Explicitly NOT Building
- PDF storage or full-text repository
- User authentication (single-user MVP)
- Payment/subscription system
- Mobile applications
- Real-time collaboration
- Data export to external systems

---

## Technical Constraints

### MVP Limitations
1. **Abstracts Only**: No PDF parsing (time constraint)
2. **Single User**: No auth or multi-tenancy
3. **Local Deployment**: Docker Compose only (no cloud scaling)
4. **Simple Entities**: 3 types only (Organism, Disease, Molecule)
5. **Basic NER**: spaCy + heuristics (no BioBERT fine-tuning)
6. **Synchronous Flow**: No real-time updates (polling-based)

### Performance Targets
- Query processing: <2 minutes end-to-end
- Paper fetching: <30 seconds for 20 papers
- Entity extraction: <10 seconds per paper
- Graph query: <5 seconds for coverage analysis
- Seed loading: <30 seconds for 30 papers
- UI response: <1 second for interactions

### Scalability Considerations
- Redis queue can handle 100+ concurrent queries
- Neo4j can store 100K+ papers (community edition)
- Ollama processes 1-2 papers/second
- Docker Compose suitable for single-user deployment

---

## Success Metrics and KPIs

### Launch Metrics (Week 1)
- System uptime: >95%
- Query success rate: >90%
- Average processing time: <2 minutes
- Knowledge graph size: 150+ papers (6 seeds)

### Adoption Metrics (Month 1)
- Active users: 20+
- Queries per week: 100+
- Seeds loaded: All 6 used at least once
- Return users: 50%+

### Quality Metrics (Ongoing)
- Entity extraction accuracy: >70% (manual validation sample)
- Relationship relevance: >60% (user feedback)
- Answer quality: 4/5 average rating
- Gap identification: 3+ gaps per query

### Learning Metrics (Month 3)
- Knowledge graph growth: 1000+ papers
- Relationship count: 5000+
- Average confidence increase: +10% for re-queried entities
- Novel connections found: 5+ reported by users

---

## Risk Assessment

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| LLM hallucination in entity extraction | High | Medium | Confidence scores, manual validation sample |
| PubMed API rate limits | Medium | Medium | Respect limits, implement backoff, use API key |
| Neo4j performance with large graphs | Low | High | Use indexes, limit query depth, optimize schema |
| Docker resource constraints | Medium | Low | Document min requirements, optimize containers |
| Entity extraction quality | High | Medium | Use established keywords, provide examples to LLM |

### Product Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Users need more than 6 seed topics | High | Low | Document how to create custom seeds |
| Processing time too slow | Medium | High | Set expectations, show progress, optimize |
| Answers not comprehensive enough | Medium | Medium | Improve prompts, fetch more papers, tune synthesis |
| Users need PDF analysis | High | Medium | Clearly scope as post-MVP, document workaround |
| Gap identification not useful | Low | High | Validate with users early, iterate on logic |

---

## Development Timeline

### Week 1: Foundation (Sept 27 - Oct 3)
- **Days 1-2**: Project scaffolding, Docker setup
- **Days 3-4**: Ollama integration, LLM service
- **Days 5-7**: Neo4j setup, schema, basic operations

### Week 2: Core Agents (Oct 4 - Oct 10)
- **Days 1-3**: Fetcher Agent (PubMed integration)
- **Days 4-5**: KG Builder Agent (entity extraction)
- **Days 6-7**: Query Strategist Agent (gap analysis)

### Week 3: Integration & Launch (Oct 11 - Oct 17)
- **Days 1-2**: Orchestrator (coordination, synthesis)
- **Days 3-4**: Streamlit UI (user interface)
- **Days 5-7**: Seed generation, testing, documentation

### Post-Launch: Iteration
- **Week 4**: Bug fixes, performance tuning
- **Week 5-6**: User feedback, feature refinement
- **Week 7-8**: Plan Phase 2 based on usage

---

## Launch Readiness Checklist

### Technical Readiness
- [ ] All services start with `docker-compose up`
- [ ] All 6 seeds load successfully
- [ ] End-to-end query flow works
- [ ] Error handling for common failures
- [ ] Logging implemented for debugging
- [ ] Performance meets targets (<2 min queries)

### Documentation Readiness
- [ ] README with setup instructions
- [ ] User guide with example queries
- [ ] Architecture documentation
- [ ] Troubleshooting guide
- [ ] API documentation (Swagger)

### Quality Readiness
- [ ] All PR acceptance criteria met
- [ ] Integration tests passing
- [ ] Sample queries validated
- [ ] Knowledge graph accuracy spot-checked
- [ ] UI/UX reviewed

### Launch Readiness
- [ ] Demo video created
- [ ] Blog post written
- [ ] GitHub repo public (if applicable)
- [ ] Feedback mechanism in place
- [ ] Support plan defined

---

## Support and Maintenance

### MVP Support Model
- **User Support**: GitHub issues or email
- **Bug Priority**: Critical (24h), High (3d), Medium (1w)
- **Updates**: Weekly during first month
- **Monitoring**: Manual log review daily

### Key Metrics to Monitor
- System uptime and errors
- Query processing time
- Agent failure rates
- Knowledge graph growth
- User feedback and issues

---

## Appendix

### A. Example Queries by Seed

**Gut Microbiome Obesity**:
- "How do gut bacteria influence obesity?"
- "What is the role of Akkermansia?"
- "How do SCFAs affect metabolism?"

**Protein Engineering**:
- "How is machine learning used in protein design?"
- "What is directed evolution?"
- "How to engineer thermostable enzymes?"

**Cancer Immunotherapy**:
- "How do PD-1 inhibitors work?"
- "What causes immunotherapy resistance?"
- "CAR-T therapy mechanisms?"

### B. Technical Dependencies

**Python Packages**:
- fastapi==0.104.1
- streamlit==1.28.0
- biopython==1.81
- spacy==3.7.2
- redis==5.0.1
- py2neo==2021.2.3
- plotly==5.17.0

**External Services**:
- PubMed E-utilities API
- Ollama API
- Neo4j Bolt protocol

### C. Glossary

- **Knowledge Graph (KG)**: Network database of entities and relationships
- **Entity**: Biological concept (organism, disease, molecule)
- **Relationship**: Connection between entities (produces, affects, inhibits)
- **Gap**: Weak evidence or missing knowledge in the graph
- **Seed**: Pre-loaded set of papers on a research topic
- **Agent**: Autonomous processing component with specific responsibility
- **Confidence**: Numeric score (0-1) indicating evidence strength

---

**Document Status**: APPROVED FOR DEVELOPMENT  
**Next Review**: Post-MVP (Week 4)  
**Feedback**: GitHub Issues or sunit@example.com
