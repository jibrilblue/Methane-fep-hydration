FROM python:3.10-slim

# Exact pins required by the platform
RUN pip install --no-cache-dir pytest==8.4.1 pytest-json-ctrf==0.3.5

# The verifier re-runs MBAR on the submitted reduced potentials
RUN pip install --no-cache-dir numpy==1.26.4 pymbar==4.0.3

COPY ./tests/ /tests/
RUN mkdir -p /app/output
RUN chmod 600 /tests/ground_truth.json
