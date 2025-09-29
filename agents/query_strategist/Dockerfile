FROM python:3.11-slim

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

# Set work directory
WORKDIR /app

# Copy project files
COPY agents/query_strategist/pyproject.toml agents/query_strategist/uv.lock* ./
COPY agents/shared/ /app/shared/

# Install dependencies
RUN uv sync --frozen --no-cache

# Copy application code
COPY agents/query_strategist/ .

# Run with uv
CMD ["uv", "run", "python", "agent.py"]