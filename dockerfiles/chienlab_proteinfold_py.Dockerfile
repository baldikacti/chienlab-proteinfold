FROM python:3.13-slim

WORKDIR /app

RUN pip install --no-cache-dir pandas requests

COPY bin/tsv2json.py .

RUN chmod +x tsv2json.py

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

CMD ["bash"]