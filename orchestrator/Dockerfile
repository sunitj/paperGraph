FROM python:3.11-slim

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

# Set work directory
WORKDIR /app

# Copy project files
COPY orchestrator/pyproject.toml orchestrator/uv.lock* ./

# Install dependencies
RUN uv sync --no-cache

# Copy application code
COPY orchestrator/ .

# Copy shared modules
COPY agents/shared/ /app/shared/

# Expose port
EXPOSE 8000

# Run with uv
CMD ["uv", "run", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]