FROM ubuntu:24.04

RUN DEBIAN_FRONTEND=noninteractive \
    apt-get update --quiet \
    && apt-get install --yes --quiet python3 python3-pip python3-venv python3-dev \
    && apt-get install --yes --quiet git wget gcc g++ make

RUN python3 -m venv /boltz_venv
ENV PATH="/boltz_venv/bin:$PATH"

RUN pip3 install --no-cache-dir --upgrade pip
RUN pip3 install --no-cache-dir boltz[cuda] -U

CMD [ "bash" ]