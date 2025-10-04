FROM python:3.11-slim

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

# Set work directory
WORKDIR /app

# Copy project files
COPY ui/pyproject.toml ui/uv.lock* ./

# Install dependencies
RUN uv sync --no-cache

# Copy application code
COPY ui/ .

# Expose port
EXPOSE 8501

# Run with uv
CMD ["uv", "run", "streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]