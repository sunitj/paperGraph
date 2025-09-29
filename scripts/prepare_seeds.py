#!/usr/bin/env python3
"""
Prepare seed data for knowledge graph
(Stub implementation for PR1)
"""

import json
import os


def generate_stub_seed(topic_name: str, description: str, paper_count: int = 25) -> dict:
    """Generate a stub seed file for a research topic"""

    seed_data = {
        "topic": topic_name,
        "description": description,
        "version": "0.1.0",
        "created": "2024-01-01T00:00:00Z",
        "papers": [],
        "entities": [],
        "relationships": []
    }

    # Generate fake papers
    for i in range(1, paper_count + 1):
        paper = {
            "pmid": f"SEED{topic_name.replace(' ', '').upper()}{i:03d}",
            "title": f"Study {i}: {topic_name} research paper",
            "abstract": f"This is a seed abstract for paper {i} about {topic_name}. "
                       f"It contains relevant information about {description.lower()}.",
            "year": 2023,
            "pre_seeded": True
        }
        seed_data["papers"].append(paper)

    # Generate fake entities
    entities = [
        {"name": "entity_1", "type": "Organism"},
        {"name": "entity_2", "type": "Disease"},
        {"name": "entity_3", "type": "Molecule"},
    ]
    seed_data["entities"] = entities

    # Generate fake relationships
    relationships = [
        {
            "subject": "entity_1",
            "subject_type": "Organism",
            "predicate": "AFFECTS",
            "object": "entity_2",
            "object_type": "Disease",
            "confidence": 0.8,
            "paper_pmid": seed_data["papers"][0]["pmid"]
        }
    ]
    seed_data["relationships"] = relationships

    return seed_data


def create_all_seeds():
    """Create all stub seed files"""

    seeds = [
        ("Gut Microbiome Obesity", "How gut bacteria influence weight and metabolism", 30),
        ("Protein Engineering", "Rational design and directed evolution of proteins", 28),
        ("Cancer Immunotherapy", "Immune system approaches to cancer treatment", 32),
        ("CRISPR Gene Editing", "CRISPR-Cas9 and related gene editing technologies", 26),
        ("Neurodegeneration", "Alzheimer's, Parkinson's and related diseases", 29),
        ("Antibiotic Resistance", "Bacterial resistance mechanisms and solutions", 27)
    ]

    # Create seeds directory if it doesn't exist
    seeds_dir = os.path.join(os.path.dirname(__file__), "..", "seeds")
    os.makedirs(seeds_dir, exist_ok=True)

    for topic, description, count in seeds:
        print(f"Creating seed for: {topic}")

        seed_data = generate_stub_seed(topic, description, count)

        # Save to JSON file
        filename = topic.lower().replace(" ", "_").replace("'", "") + ".json"
        filepath = os.path.join(seeds_dir, filename)

        with open(filepath, 'w') as f:
            json.dump(seed_data, f, indent=2)

        print(f"  Saved: {filename} ({count} papers)")

    print(f"\n✅ Created {len(seeds)} seed files in {seeds_dir}")


if __name__ == "__main__":
    print("PaperGraph Seed Preparation (Stub)")
    print("=" * 40)
    create_all_seeds()