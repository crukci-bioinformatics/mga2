FROM rust:1.51-slim-buster

LABEL authors="Matt Eldridge" \
      version="2.0.2" \
      description="Multi-Genome Analysis (MGA) contaminant screen for genomic sequence data"

RUN apt-get update && \
    apt-get install -y git curl procps && \
    apt-get clean

ARG CONDA_VERSION=py39_4.11.0
ARG CONDA_MD5=4e2f31e0b2598634c80daa12e4981647

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -o miniconda3.sh && \
    echo "${CONDA_MD5}  miniconda3.sh" > miniconda3.md5 && \
    mkdir -p /opt && \
    sh miniconda3.sh -b -p /opt/conda && \
    rm miniconda3.sh miniconda3.md5

COPY conda.yml .

RUN /opt/conda/bin/conda env create -f conda.yml && /opt/conda/bin/conda clean -a

ENV PATH /opt/mga2/bin:/opt/conda/envs/mga2/bin:$PATH

# ensures that a new build is run if commits have been made since the last build
ADD "https://api.github.com/repos/crukci-bioinformatics/mga2/commits?per_page=1" /opt/.latest_mga_commit

ARG VERSION=2.0.2

RUN git clone --branch ${VERSION} --depth 1 https://github.com/crukci-bioinformatics/mga2 /opt/mga2 && \
    cargo install --path /opt/mga2 --root /opt/mga2

