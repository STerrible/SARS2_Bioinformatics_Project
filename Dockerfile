# syntax=docker/dockerfile:1

ARG NEXTCLADE_VERSION=3.21.2
FROM nextstrain/nextclade:${NEXTCLADE_VERSION}-debian

ARG NEXTCLADE_VERSION=3.21.2
ARG NEXTCLADE_DATASET_NAME=nextstrain/sars-cov-2/wuhan-hu-1/orfs
ARG NEXTCLADE_DATASET_TAG=2026-04-21--09-39-50Z
ARG BIOCONDA_PREFIX=/opt/bioconda

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    NEXTCLADE_EXE=nextclade \
    NEXTCLADE_DATASET=/opt/nextclade/datasets/sars-cov-2 \
    NEXTCLADE_DATASET_NAME=${NEXTCLADE_DATASET_NAME} \
    NEXTCLADE_DATASET_TAG=${NEXTCLADE_DATASET_TAG} \
    MAMBA_ROOT_PREFIX=/opt/micromamba \
    SARS2_PIPELINE_CONTAINER_IMAGE=sars2-bioinformatics:full-${NEXTCLADE_VERSION}

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        curl \
        python3 \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt

RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir -r /app/requirements.txt

RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba \
    && micromamba create -y -p "${BIOCONDA_PREFIX}" -c conda-forge -c bioconda \
        snippy=4.6.0 \
        fasttree=2.1.11 \
    && micromamba clean --all --yes

ENV PATH=${BIOCONDA_PREFIX}/bin:/opt/venv/bin:${PATH}

RUN printf 'export PATH=%s/bin:/opt/venv/bin:$PATH\n' "${BIOCONDA_PREFIX}" > /etc/profile.d/project-path.sh

RUN mkdir -p "${NEXTCLADE_DATASET}" \
    && nextclade dataset get \
        --name "${NEXTCLADE_DATASET_NAME}" \
        --tag "${NEXTCLADE_DATASET_TAG}" \
        --output-dir "${NEXTCLADE_DATASET}" \
    && nextclade --version \
    && test -f "${NEXTCLADE_DATASET}/pathogen.json"

RUN snippy --check \
    && FastTree -help >/tmp/fasttree-help.txt 2>&1 || test -s /tmp/fasttree-help.txt

COPY . /app

CMD ["bash", "-c", "/opt/venv/bin/python main.py sars2 && /opt/venv/bin/python scripts/generate_report.py"]
