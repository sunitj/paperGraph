#!/bin/bash

echo "Setting up Ollama..."

# Wait for Ollama to be ready
echo "Waiting for Ollama service..."
until curl -s http://localhost:11434/api/tags > /dev/null; do
    sleep 2
done

echo "✅ Ollama is ready"

# Pull models
echo "Pulling llama3.1:8b (this may take a few minutes)..."
docker exec papergraph-ollama-1 ollama pull llama3.1:8b

echo "Pulling mistral:7b as backup..."
docker exec papergraph-ollama-1 ollama pull mistral:7b

echo "✅ Models ready!"
echo ""
echo "Available models:"
docker exec papergraph-ollama-1 ollama list