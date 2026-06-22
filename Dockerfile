# syntax=docker/dockerfile:1

ARG NEXTCLADE_VERSION=3.21.2
FROM nextstrain/nextclade:${NEXTCLADE_VERSION}-debian

ARG NEXTCLADE_VERSION=3.21.2
ARG NEXTCLADE_DATASET_NAME=nextstrain/sars-cov-2/wuhan-hu-1/orfs
ARG NEXTCLADE_DATASET_TAG=2026-04-21--09-39-50Z

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    NEXTCLADE_EXE=nextclade \
    NEXTCLADE_DATASET=/opt/nextclade/datasets/sars-cov-2 \
    NEXTCLADE_DATASET_NAME=${NEXTCLADE_DATASET_NAME} \
    NEXTCLADE_DATASET_TAG=${NEXTCLADE_DATASET_TAG} \
    SARS2_PIPELINE_CONTAINER_IMAGE=sars2-bioinformatics:nextclade-${NEXTCLADE_VERSION}

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        python3 \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt

RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir -r /app/requirements.txt

ENV PATH=/opt/venv/bin:${PATH}

RUN mkdir -p "${NEXTCLADE_DATASET}" \
    && nextclade dataset get \
        --name "${NEXTCLADE_DATASET_NAME}" \
        --tag "${NEXTCLADE_DATASET_TAG}" \
        --output-dir "${NEXTCLADE_DATASET}" \
    && nextclade --version \
    && test -f "${NEXTCLADE_DATASET}/pathogen.json"

COPY . /app

CMD ["bash", "-c", "/opt/venv/bin/python main.py sars2 && /opt/venv/bin/python scripts/generate_report.py"]
