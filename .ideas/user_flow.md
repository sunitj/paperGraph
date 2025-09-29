1. User selects seed topic (e.g., "Gut Microbiome & Obesity") 
   → Loads pre-built KG with 30 papers
   
2. User asks: "How do gut bacteria influence obesity?"
   → Query Strategist checks KG for gaps
   → Fetches targeted papers to fill gaps
   → Agents extract entities and relationships
   → Updates KG with new knowledge
   
3. User asks follow-up: "What about Akkermansia?"
   → System knows Akkermansia is underexplored (from KG)
   → Searches specifically for Akkermansia mechanisms
   → Shows novel pathway in graph visualization
   
4. System suggests: "Gap found: diet→microbiome link missing"
   → User can explore suggested gaps