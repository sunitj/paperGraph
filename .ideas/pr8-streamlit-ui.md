# PR #8: Streamlit UI - User Interface

**Sprint:** Week 3, Days 3-4  
**Goal:** Build interactive user interface

## What to Build

1. Question input interface
2. Seed selection sidebar
3. Graph visualization
4. Results display

## Deliverable

```bash
# Visit http://localhost:8501
# Select seed → Ask question → View results + graph
```

## ui/app.py

```python
import streamlit as st
import requests
import json
import plotly.graph_objects as go
from datetime import datetime

# Page config
st.set_page_config(
    page_title="PaperGraph Intelligence",
    page_icon="🔬",
    layout="wide"
)

# API URL
API_URL = "http://orchestrator:8000"

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        margin-bottom: 2rem;
    }
    .gap-badge {
        background-color: #ff4b4b;
        color: white;
        padding: 0.25rem 0.5rem;
        border-radius: 0.25rem;
        font-size: 0.8rem;
    }
</style>
""", unsafe_allow_html=True)


def load_seeds():
    """Load available seeds"""
    try:
        response = requests.get(f"{API_URL}/seeds", timeout=5)
        if response.status_code == 200:
            return response.json()['seeds']
    except:
        pass
    return []


def load_seed(seed_name):
    """Load a seed into the knowledge graph"""
    try:
        response = requests.post(
            f"{API_URL}/load_seed",
            json={"seed": seed_name},
            timeout=30
        )
        return response.json()
    except Exception as e:
        st.error(f"Error loading seed: {e}")
        return None


def query_system(question):
    """Send query to orchestrator"""
    with st.spinner("🔍 Analyzing knowledge graph and searching literature..."):
        try:
            response = requests.post(
                f"{API_URL}/query",
                json={"question": question},
                timeout=120
            )
            
            if response.status_code == 200:
                return response.json()
            else:
                st.error(f"Error: {response.status_code}")
                return None
        
        except Exception as e:
            st.error(f"Query error: {e}")
            return None


def get_stats():
    """Get system statistics"""
    try:
        response = requests.get(f"{API_URL}/stats", timeout=5)
        if response.status_code == 200:
            return response.json()
    except:
        pass
    return None


def create_network_graph(subgraph_data):
    """Create network visualization"""
    
    # For MVP, create a simple visualization
    # In production, properly parse Neo4j paths
    
    fig = go.Figure()
    
    # Add placeholder nodes
    fig.add_trace(go.Scatter(
        x=[0, 1, 2, 1],
        y=[0, 1, 0, -1],
        mode='markers+text',
        text=['Entity 1', 'Entity 2', 'Entity 3', 'Entity 4'],
        textposition="top center",
        marker=dict(size=20, color='lightblue'),
        hoverinfo='text'
    ))
    
    # Add placeholder edges
    edges_x = [0, 1, None, 1, 2, None, 1, 1, None]
    edges_y = [0, 1, None, 1, 0, None, 1, -1, None]
    
    fig.add_trace(go.Scatter(
        x=edges_x,
        y=edges_y,
        mode='lines',
        line=dict(width=1, color='#888'),
        hoverinfo='none'
    ))
    
    fig.update_layout(
        showlegend=False,
        hovermode='closest',
        margin=dict(b=0, l=0, r=0, t=0),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=400
    )
    
    return fig


# Main UI
st.markdown('<h1 class="main-header">🔬 PaperGraph Intelligence</h1>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Multi-agent literature analysis with knowledge graph learning</p>', unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.header("🌱 Knowledge Base")
    
    # Load seeds
    seeds = load_seeds()
    
    if seeds:
        seed_options = ["None (Start Fresh)"] + [s['name'] for s in seeds]
        selected_seed = st.selectbox(
            "Select seed topic:",
            seed_options,
            help="Pre-loaded research topics with 20-30 papers"
        )
        
        # Show seed info
        if selected_seed != "None (Start Fresh)":
            seed_info = next(s for s in seeds if s['name'] == selected_seed)
            st.info(f"📄 {seed_info['papers']} papers\n\n{seed_info.get('description', '')}")
            
            if st.button("🚀 Load Seed", use_container_width=True):
                with st.spinner("Loading knowledge graph..."):
                    result = load_seed(selected_seed)
                    if result:
                        st.success(f"✅ Loaded {result['papers_count']} papers!")
                        st.rerun()
    
    # System stats
    st.divider()
    st.subheader("📊 System Stats")
    
    stats = get_stats()
    if stats:
        col1, col2 = st.columns(2)
        col1.metric("Papers", stats['papers'])
        col2.metric("Entities", stats['entities'])
        st.metric("Relationships", stats['relationships'])
    
    # Info
    st.divider()
    st.markdown("""
    ### How it works
    1. Select a seed topic or start fresh
    2. Ask research questions
    3. System identifies knowledge gaps
    4. Agents fetch and analyze papers
    5. Get answers with context
    """)

# Main content
tab1, tab2 = st.tabs(["💬 Ask Questions", "📈 Knowledge Graph"])

with tab1:
    # Question input
    question = st.text_input(
        "Ask a research question:",
        placeholder="How do gut bacteria influence obesity?",
        help="Ask about mechanisms, relationships, or evidence in scientific literature"
    )
    
    col1, col2, col3 = st.columns([2, 1, 1])
    with col1:
        search_button = st.button("🔍 Search", type="primary", use_container_width=True)
    with col2:
        max_papers = st.number_input("Max papers", 5, 50, 20, 5)
    
    if search_button and question:
        # Query the system
        result = query_system(question)
        
        if result:
            # Display answer
            st.markdown("### 💡 Answer")
            st.markdown(result['answer'])
            
            # Display gaps
            if result.get('gaps'):
                st.markdown("### 🔍 Knowledge Gaps Identified")
                for gap in result['gaps'][:3]:
                    gap_type = gap.get('type', 'unknown')
                    
                    if gap_type == 'weak_evidence':
                        st.warning(f"**Weak Evidence**: {gap.get('subject')} → {gap.get('object')} (confidence: {gap.get('confidence', 0):.2f})")
                    elif gap_type == 'understudied':
                        st.info(f"**Understudied**: {gap.get('subject')} → {gap.get('object')} ({gap.get('paper_count', 0)} papers)")
                    elif gap_type == 'missing_entity':
                        st.error(f"**Missing**: No information found about {gap.get('entity')}")
            
            # Display papers
            st.markdown("### 📚 Papers Analyzed")
            
            for i, paper in enumerate(result['papers'][:5], 1):
                with st.expander(f"{i}. {paper['title']} ({paper.get('year', 'N/A')})"):
                    st.markdown(f"**PMID**: [{paper['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']})")
                    if 'abstract' in paper:
                        st.markdown(paper['abstract'][:300] + "...")
            
            # Display entities
            if result.get('entities'):
                st.markdown("### 🏷️ Entities Identified")
                st.markdown(", ".join([f"`{e}`" for e in result['entities']]))

with tab2:
    st.markdown("### 📊 Knowledge Network")
    
    stats = get_stats()
    if stats and stats['relationships'] > 0:
        # Show network graph
        fig = create_network_graph({})
        st.plotly_chart(fig, use_container_width=True)
        
        # Query examples
        st.markdown("### 🔎 Explore in Neo4j Browser")
        st.code("""
# View all relationships
MATCH (n)-[r]->(m) RETURN n, r, m LIMIT 50

# Find weak evidence (gaps)
MATCH (a)-[r]->(b)
WHERE r.confidence < 0.5
RETURN a.name, type(r), b.name, r.confidence
ORDER BY r.confidence

# Find most connected entities
MATCH (n)-[r]-()
RETURN n.name, count(r) as connections
ORDER BY connections DESC
LIMIT 10
        """, language="cypher")
        
        st.info("💡 Open Neo4j Browser at http://localhost:7474 (neo4j/papergraph123)")
    else:
        st.info("📭 No knowledge graph data yet. Load a seed or ask a question to get started.")

# Footer
st.divider()
st.markdown("""
<div style='text-align: center; color: #666; padding: 1rem;'>
    <p>PaperGraph Intelligence | Multi-Agent Literature Analysis</p>
    <p style='font-size: 0.8rem;'>Powered by Ollama LLMs + Neo4j Knowledge Graph</p>
</div>
""", unsafe_allow_html=True)
```

## ui/requirements.txt

```txt
streamlit==1.28.0
requests==2.31.0
plotly==5.17.0
```

## Testing

```bash
# Start UI
docker-compose up -d ui

# Open browser
open http://localhost:8501

# Or from CLI
streamlit run ui/app.py
```

## Features

### 1. Seed Selection
- Load pre-built knowledge graphs
- Shows paper count and description
- One-click loading

### 2. Question Interface
- Natural language input
- Adjustable paper count
- Real-time processing feedback

### 3. Results Display
- Synthesized answer
- Identified gaps (color-coded)
- Paper references with links
- Entity extraction

### 4. Graph Visualization
- Network graph (placeholder in MVP)
- Neo4j query examples
- Direct link to Neo4j browser

### 5. System Stats
- Live paper count
- Entity count
- Relationship count

## Acceptance Criteria

- [ ] UI loads successfully
- [ ] Can select and load seeds
- [ ] Question input works
- [ ] Results display correctly
- [ ] Papers show with PubMed links
- [ ] Gaps are highlighted
- [ ] System stats update
- [ ] Graph visualization renders
- [ ] Responsive design
- [ ] Error handling works
