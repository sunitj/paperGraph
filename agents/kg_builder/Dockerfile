FROM python:3.11-slim

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

# Set work directory
WORKDIR /app

# Copy project files
COPY agents/kg_builder/pyproject.toml agents/kg_builder/uv.lock* ./
COPY agents/shared/ /app/shared/

# Install dependencies
RUN uv sync --frozen --no-cache

# Download spaCy model
RUN uv run python -m spacy download en_core_web_sm

# Copy application code
COPY agents/kg_builder/ .

# Run with uv
CMD ["uv", "run", "python", "agent.py"]