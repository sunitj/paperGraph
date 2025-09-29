#!/usr/bin/env python3
"""
PaperGraph Intelligence UI
Streamlit user interface (stub implementation)
"""

import streamlit as st
import requests
import os
import json


# Configure Streamlit page
st.set_page_config(
    page_title="PaperGraph Intelligence",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Configuration
ORCHESTRATOR_URL = os.getenv("ORCHESTRATOR_URL", "http://localhost:8000")


def check_orchestrator_health():
    """Check if orchestrator is healthy"""
    try:
        response = requests.get(f"{ORCHESTRATOR_URL}/health", timeout=5)
        return response.status_code == 200, response.json()
    except:
        return False, None


def main():
    """Main UI application"""

    # Title
    st.title("🔬 PaperGraph Intelligence")
    st.markdown("Multi-agent literature analysis with knowledge graph learning")

    # Sidebar
    with st.sidebar:
        st.header("System Status")

        # Health check
        healthy, health_data = check_orchestrator_health()

        if healthy:
            st.success("✅ System Online")
            if health_data:
                st.json(health_data)
        else:
            st.error("❌ System Offline")
            st.write("Orchestrator not reachable")

    # Main content tabs
    tab1, tab2, tab3 = st.tabs(["Query", "Knowledge Graph", "Settings"])

    with tab1:
        st.header("Ask a Research Question")

        # Question input
        question = st.text_area(
            "Enter your research question:",
            placeholder="e.g., How do gut bacteria influence obesity?",
            height=100
        )

        col1, col2 = st.columns([1, 3])

        with col1:
            max_papers = st.number_input("Max Papers", min_value=5, max_value=50, value=20)

        with col2:
            if st.button("🔍 Search", type="primary", disabled=not question.strip()):
                if question.strip():
                    with st.spinner("Processing your question..."):
                        # Stub implementation
                        st.info("Query processing is not yet implemented (stub UI)")

                        # Show fake results
                        st.success("Query completed!")

                        st.subheader("Answer")
                        st.write("""
                        This is a stub answer. In the full implementation, this would contain
                        a synthesized response based on the knowledge graph and newly fetched papers.
                        """)

                        st.subheader("Papers Found")
                        st.write("5 papers would be listed here")

                        st.subheader("Knowledge Gaps")
                        st.write("Identified gaps would be shown here")

    with tab2:
        st.header("Knowledge Graph Visualization")
        st.info("Graph visualization will be implemented in future PRs")

        # Stub statistics
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Papers", "0")

        with col2:
            st.metric("Entities", "0")

        with col3:
            st.metric("Relationships", "0")

        with col4:
            st.metric("Avg Confidence", "0.0")

    with tab3:
        st.header("System Settings")

        st.subheader("Pre-seeded Topics")
        st.info("Seed loading will be implemented in future PRs")

        topics = [
            "Gut Microbiome & Obesity",
            "Protein Engineering",
            "Cancer Immunotherapy",
            "CRISPR Gene Editing",
            "Neurodegeneration",
            "Antibiotic Resistance"
        ]

        for topic in topics:
            col1, col2 = st.columns([3, 1])
            with col1:
                st.write(f"📚 {topic}")
            with col2:
                st.button("Load", key=f"load_{topic}", disabled=True)


if __name__ == "__main__":
    main()