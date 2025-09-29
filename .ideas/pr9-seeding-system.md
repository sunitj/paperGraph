# PR #9: KG Seeding System

**Sprint:** Week 3, Days 5-7  
**Goal:** Create pre-seeded knowledge graphs for 6 research topics

## What to Build

1. Seed data preparation script
2. 6 seed JSON files with 20-30 papers each
3. Seed loading integration
4. Documentation

## Topics

1. Gut Microbiome & Obesity
2. Cancer Immunotherapy
3. Antibiotic Resistance Mechanisms
4. CRISPR Gene Editing
5. Neurodegenerative Diseases
6. **Protein Engineering** ✨ (NEW)

## Deliverable

```bash
# Create seeds
python scripts/prepare_seeds.py

# Load via API
curl -X POST http://localhost:8000/load_seed -d '{"seed":"Protein Engineering"}'
```

## scripts/prepare_seeds.py

```python
#!/usr/bin/env python3
"""
Prepare seed knowledge graphs from PubMed
"""

import json
import os
from Bio import Entrez
from pathlib import Path

Entrez.email = os.getenv("PUBMED_EMAIL", "papergraph@example.com")


def fetch_papers(query: str, num_papers: int = 30) -> list:
    """Fetch papers from PubMed"""
    print(f"  Fetching papers for: {query}")
    
    try:
        # Search
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=num_papers,
            sort="relevance"
        )
        record = Entrez.read(handle)
        handle.close()
        
        pmids = record.get("IdList", [])
        print(f"    Found {len(pmids)} papers")
        
        if not pmids:
            return []
        
        # Fetch details
        handle = Entrez.efetch(
            db="pubmed",
            id=pmids,
            rettype="abstract",
            retmode="xml"
        )
        articles = Entrez.read(handle)
        handle.close()
        
        papers = []
        for article in articles.get('PubmedArticle', []):
            try:
                medline = article.get('MedlineCitation', {})
                pmid = str(medline.get('PMID', ''))
                
                article_data = medline.get('Article', {})
                title = article_data.get('ArticleTitle', '')
                
                abstract_data = article_data.get('Abstract', {})
                abstract_text = abstract_data.get('AbstractText', [])
                
                if isinstance(abstract_text, list):
                    abstract = ' '.join([str(t) for t in abstract_text])
                else:
                    abstract = str(abstract_text)
                
                pub_date = article_data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
                year = pub_date.get('Year', '2024')
                
                if '-' in str(year):
                    year = str(year).split('-')[0]
                
                try:
                    year = int(year)
                except:
                    year = 2024
                
                if title and abstract and len(abstract) > 100:
                    papers.append({
                        'pmid': pmid,
                        'title': title,
                        'abstract': abstract,
                        'year': year
                    })
            
            except Exception as e:
                print(f"      Error parsing article: {e}")
                continue
        
        print(f"    Processed {len(papers)} papers")
        return papers
    
    except Exception as e:
        print(f"    Error: {e}")
        return []


def create_seed(topic: str, description: str, query: str, num_papers: int = 30):
    """Create a seed file"""
    print(f"\n📚 Creating seed: {topic}")
    
    papers = fetch_papers(query, num_papers)
    
    if not papers:
        print(f"  ⚠️  No papers found for {topic}")
        return
    
    seed = {
        'topic': topic,
        'description': description,
        'query': query,
        'papers': papers,
        'created_at': '2024-09-26'
    }
    
    # Save to file
    filename = topic.lower().replace(" ", "_").replace("&", "and")
    filepath = Path('seeds') / f"{filename}.json"
    
    os.makedirs('seeds', exist_ok=True)
    
    with open(filepath, 'w') as f:
        json.dump(seed, f, indent=2)
    
    print(f"  ✅ Created: {filepath} ({len(papers)} papers)")


# Seed definitions
SEEDS = [
    {
        'topic': 'Gut Microbiome Obesity',
        'description': 'Explores how gut bacteria influence weight, metabolism, and obesity through various mechanisms including SCFA production and inflammation.',
        'query': 'gut microbiome obesity mechanisms metabolic',
        'num_papers': 30
    },
    {
        'topic': 'Cancer Immunotherapy',
        'description': 'Research on immune checkpoint inhibitors, CAR-T cells, and other immunotherapy approaches for cancer treatment.',
        'query': 'cancer immunotherapy PD-1 CTLA-4 checkpoint inhibitor',
        'num_papers': 30
    },
    {
        'topic': 'Antibiotic Resistance',
        'description': 'Mechanisms of bacterial antibiotic resistance including efflux pumps, beta-lactamases, and horizontal gene transfer.',
        'query': 'antibiotic resistance mechanisms bacteria efflux pump',
        'num_papers': 30
    },
    {
        'topic': 'CRISPR Gene Editing',
        'description': 'CRISPR-Cas9 and other gene editing technologies, including applications, off-target effects, and base editing.',
        'query': 'CRISPR Cas9 gene editing therapeutic applications',
        'num_papers': 30
    },
    {
        'topic': 'Neurodegenerative Diseases',
        'description': 'Mechanisms of Alzheimer\'s, Parkinson\'s, and other neurodegenerative diseases including protein aggregation and neuroinflammation.',
        'query': 'neurodegeneration Alzheimer Parkinson amyloid tau protein',
        'num_papers': 30
    },
    {
        'topic': 'Protein Engineering',
        'description': 'Rational design and directed evolution of proteins for therapeutics, enzymes, and biosensors including machine learning approaches.',
        'query': 'protein engineering design directed evolution therapeutic enzyme',
        'num_papers': 30
    }
]


def main():
    """Create all seeds"""
    print("="*60)
    print("PaperGraph Seed Generation")
    print("="*60)
    print(f"\nCreating {len(SEEDS)} seed knowledge graphs...")
    
    for seed_config in SEEDS:
        create_seed(**seed_config)
    
    print("\n" + "="*60)
    print("✅ All seeds created!")
    print("\nSeeds available:")
    for seed_config in SEEDS:
        print(f"  - {seed_config['topic']}")
    print("\nTo load a seed:")
    print("  curl -X POST http://localhost:8000/load_seed \\")
    print("    -H 'Content-Type: application/json' \\")
    print("    -d '{\"seed\":\"Protein Engineering\"}'")
    print("="*60)


if __name__ == "__main__":
    main()
```

## Seed File Format

```json
{
  "topic": "Protein Engineering",
  "description": "Rational design and directed evolution of proteins...",
  "query": "protein engineering design directed evolution therapeutic enzyme",
  "created_at": "2024-09-26",
  "papers": [
    {
      "pmid": "38123456",
      "title": "Machine Learning-Guided Protein Engineering for Enhanced Stability",
      "abstract": "Recent advances in machine learning have revolutionized protein engineering...",
      "year": 2024
    }
  ]
}
```

## Seed Descriptions

### 1. Gut Microbiome Obesity
**Focus**: Firmicutes/Bacteroidetes ratio, SCFA production, inflammation  
**Key entities**: Organisms (bacteria), Molecules (SCFAs), Disease (obesity)  
**Typical queries**: "How do gut bacteria affect weight?", "What is the role of SCFAs?"

### 2. Cancer Immunotherapy
**Focus**: PD-1, CTLA-4, CAR-T cells, tumor microenvironment  
**Key entities**: Molecules (checkpoint inhibitors), Disease (cancer), Mechanisms  
**Typical queries**: "How do checkpoint inhibitors work?", "What is CAR-T therapy?"

### 3. Antibiotic Resistance
**Focus**: Efflux pumps, beta-lactamases, horizontal gene transfer  
**Key entities**: Organisms (bacteria), Genes, Mechanisms  
**Typical queries**: "How does E. coli resist antibiotics?", "What are efflux pumps?"

### 4. CRISPR Gene Editing
**Focus**: Cas9, guide RNA, off-target effects, base editing  
**Key entities**: Proteins (Cas9), Molecules (RNA), Applications  
**Typical queries**: "How does CRISPR work?", "What are off-target effects?"

### 5. Neurodegenerative Diseases
**Focus**: Amyloid-beta, tau protein, neuroinflammation, mitochondria  
**Key entities**: Proteins, Disease (Alzheimer's, Parkinson's), Pathways  
**Typical queries**: "What causes Alzheimer's?", "Role of tau protein?"

### 6. Protein Engineering ✨
**Focus**: Rational design, directed evolution, ML-guided design, enzyme engineering  
**Key entities**: Proteins, Methods (design, evolution), Applications (therapeutic, industrial)  
**Typical queries**: "How to engineer stable proteins?", "What is directed evolution?"

## Usage

```bash
# 1. Generate all seeds
python scripts/prepare_seeds.py

# 2. View generated seeds
ls -lh seeds/
cat seeds/protein_engineering.json | jq '.topic, .papers | length'

# 3. Load in UI
# Go to http://localhost:8501
# Select "Protein Engineering" from sidebar
# Click "Load Seed"

# 4. Or load via API
curl -X POST http://localhost:8000/load_seed \
  -H "Content-Type: application/json" \
  -d '{"seed":"Protein Engineering"}'

# 5. Query the seed
curl -X POST http://localhost:8000/query \
  -H "Content-Type: application/json" \
  -d '{"question":"How is machine learning used in protein engineering?"}'
```

## Example Questions per Seed

### Gut Microbiome Obesity
- "How do gut bacteria influence obesity?"
- "What is the role of Akkermansia in metabolic health?"
- "How do SCFAs affect inflammation?"

### Cancer Immunotherapy
- "How do PD-1 inhibitors work?"
- "What is the difference between PD-1 and CTLA-4?"
- "How does CAR-T therapy target tumors?"

### Antibiotic Resistance
- "How does E. coli develop resistance?"
- "What are the main resistance mechanisms?"
- "How does horizontal gene transfer spread resistance?"

### CRISPR Gene Editing
- "How accurate is CRISPR gene editing?"
- "What are the therapeutic applications of CRISPR?"
- "How to minimize off-target effects?"

### Neurodegenerative Diseases
- "What causes amyloid plaque formation?"
- "How does tau protein contribute to Alzheimer's?"
- "What is the role of neuroinflammation?"

### Protein Engineering
- "How to design thermostable enzymes?"
- "What is the difference between rational design and directed evolution?"
- "How is machine learning used in protein engineering?"
- "What are best practices for antibody engineering?"

## Testing

```bash
# Test seed creation
python scripts/prepare_seeds.py

# Verify files
ls seeds/*.json

# Test loading
curl http://localhost:8000/seeds | jq

# Load and query
curl -X POST http://localhost:8000/load_seed -H "Content-Type: application/json" -d '{"seed":"Protein Engineering"}'

curl -X POST http://localhost:8000/query -H "Content-Type: application/json" -d '{"question":"How does directed evolution improve protein function?"}'
```

## Acceptance Criteria

- [ ] All 6 seeds created successfully
- [ ] Each seed has 20-30 papers
- [ ] Protein Engineering seed included
- [ ] Seeds load via API
- [ ] Seeds visible in UI dropdown
- [ ] Papers have abstracts
- [ ] JSON format is valid
- [ ] Documentation complete
- [ ] Example queries provided
