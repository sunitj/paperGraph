#!/bin/bash

# PaperGraph Intelligence Setup Script
# Automates initial setup and verification

set -e

echo "🔬 PaperGraph Intelligence Setup"
echo "================================="
echo

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper functions
log_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

log_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

log_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

log_error() {
    echo -e "${RED}❌ $1${NC}"
}

# Check prerequisites
check_prerequisites() {
    log_info "Checking prerequisites..."

    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed. Please install Docker first."
        exit 1
    fi

    if ! command -v docker-compose &> /dev/null; then
        log_error "Docker Compose is not installed. Please install Docker Compose first."
        exit 1
    fi

    log_success "Prerequisites check passed"
}

# Setup environment file
setup_env() {
    log_info "Setting up environment configuration..."

    if [ ! -f .env ]; then
        if [ -f .env.example ]; then
            cp .env.example .env
            log_success "Created .env from .env.example"
            log_warning "Please edit .env and add your PUBMED_EMAIL"
        else
            log_error ".env.example not found"
            exit 1
        fi
    else
        log_success ".env already exists"
    fi
}

# Start services
start_services() {
    log_info "Starting Docker services..."

    docker-compose down --remove-orphans
    docker-compose up -d

    log_success "Services started"
}

# Wait for services
wait_for_services() {
    log_info "Waiting for services to be ready..."

    # Wait for Ollama
    log_info "Waiting for Ollama..."
    for i in {1..60}; do
        if curl -f http://localhost:11434/api/tags &> /dev/null; then
            log_success "Ollama is ready"
            break
        fi
        if [ $i -eq 60 ]; then
            log_error "Ollama failed to start within 60 seconds"
            exit 1
        fi
        sleep 1
    done

    # Wait for Neo4j
    log_info "Waiting for Neo4j..."
    for i in {1..60}; do
        if curl -f http://localhost:7474 &> /dev/null; then
            log_success "Neo4j is ready"
            break
        fi
        if [ $i -eq 60 ]; then
            log_error "Neo4j failed to start within 60 seconds"
            exit 1
        fi
        sleep 1
    done

    # Wait for Redis
    log_info "Waiting for Redis..."
    for i in {1..30}; do
        if docker-compose exec -T redis redis-cli ping &> /dev/null; then
            log_success "Redis is ready"
            break
        fi
        if [ $i -eq 30 ]; then
            log_error "Redis failed to start within 30 seconds"
            exit 1
        fi
        sleep 1
    done

    # Wait for Orchestrator
    log_info "Waiting for Orchestrator..."
    for i in {1..60}; do
        if curl -f http://localhost:8000/health &> /dev/null; then
            log_success "Orchestrator is ready"
            break
        fi
        if [ $i -eq 60 ]; then
            log_error "Orchestrator failed to start within 60 seconds"
            exit 1
        fi
        sleep 1
    done
}

# Pull LLM models
pull_models() {
    log_info "Pulling LLM models (this may take several minutes)..."

    log_info "Pulling llama3.1:8b..."
    docker-compose exec ollama ollama pull llama3.1:8b

    log_info "Pulling mistral:7b as backup..."
    docker-compose exec ollama ollama pull mistral:7b

    log_success "LLM models ready"
}

# Generate seed data
generate_seeds() {
    log_info "Generating seed data..."

    if command -v python3 &> /dev/null; then
        python3 scripts/prepare_seeds.py
        log_success "Seed data generated"
    else
        log_warning "Python3 not found, skipping seed generation"
    fi
}

# Verify setup
verify_setup() {
    log_info "Verifying setup..."

    # Check service health
    services=("ollama:11434/api/tags" "neo4j:7474" "redis" "orchestrator:8000/health" "ui:8501")

    for service in "${services[@]}"; do
        if [[ $service == *":"* ]]; then
            name=$(echo $service | cut -d: -f1)
            endpoint=$(echo $service | cut -d: -f2-)

            if [[ $name == "redis" ]]; then
                if docker-compose exec -T redis redis-cli ping &> /dev/null; then
                    log_success "$name is healthy"
                else
                    log_error "$name is not responding"
                fi
            else
                if curl -f "http://localhost:$endpoint" &> /dev/null; then
                    log_success "$name is healthy"
                else
                    log_error "$name is not responding"
                fi
            fi
        fi
    done
}

# Show access information
show_access_info() {
    echo
    log_success "Setup complete! 🎉"
    echo
    echo "📋 Access Information:"
    echo "┌─────────────────────────────────────────────┐"
    echo "│ Service              │ URL                  │"
    echo "├─────────────────────────────────────────────┤"
    echo "│ Streamlit UI         │ http://localhost:8501│"
    echo "│ API Documentation    │ http://localhost:8000│"
    echo "│ Neo4j Browser        │ http://localhost:7474│"
    echo "│ Ollama API           │ http://localhost:11434│"
    echo "└─────────────────────────────────────────────┘"
    echo
    echo "🔐 Neo4j Credentials:"
    echo "   Username: neo4j"
    echo "   Password: papergraph123"
    echo
    echo "📚 Next Steps:"
    echo "   1. Open the Streamlit UI"
    echo "   2. Try asking a research question"
    echo "   3. Explore the knowledge graph in Neo4j"
    echo
}

# Main execution
main() {
    check_prerequisites
    setup_env
    start_services
    wait_for_services
    pull_models
    generate_seeds
    verify_setup
    show_access_info
}

# Run main function
main "$@"