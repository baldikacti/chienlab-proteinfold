FROM python:3.13-slim

WORKDIR /app

RUN pip install --no-cache-dir pandas requests

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

CMD ["bash"]