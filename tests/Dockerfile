FROM python:3.10-slim

RUN pip install --no-cache-dir pytest==8.4.1 pytest-json-ctrf==0.3.5 numpy==1.26.4 scipy==1.13.1 pymbar==4.0.3

COPY ./tests/ /tests/
RUN mkdir -p /app/output
RUN chmod 600 /tests/ground_truth.json
