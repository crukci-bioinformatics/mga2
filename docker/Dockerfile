FROM rust:1.51-slim-buster

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk" \
           version "2.0.0" \
           description "Multi-Genome Analysis (MGA) contaminant screen for genomic sequence data"

RUN apt-get update && \
    apt-get install -y git curl procps && \
    apt-get clean

ARG CONDA_VERSION=py38_4.9.2
ARG CONDA_MD5=122c8c9beb51e124ab32a0fa6426c656

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda3.sh && \
    echo "${CONDA_MD5}  miniconda3.sh" > miniconda3.md5 && \
    mkdir -p /opt && \
    sh miniconda3.sh -b -p /opt/conda && \
    rm miniconda3.sh miniconda3.md5

COPY conda.yml .

RUN /opt/conda/bin/conda env create -f conda.yml && /opt/conda/bin/conda clean -a

ENV PATH /opt/conda/envs/mga2/bin:$PATH

# ensures that a new build is run if commits have been made since the last build
ADD "https://api.github.com/repos/crukci-bioinformatics/mga2/commits?per_page=1" /opt/mga2/.latest_commit

ARG VERSION=master
RUN git clone --branch ${VERSION} --depth 1 https://github.com/crukci-bioinformatics/mga2 /opt/mga2/build && \
    cargo install --path /opt/mga2/build

