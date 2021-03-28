FROM rust:1.51-slim-buster

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk" \
           version "2.0.0" \
           description "Multi-Genome Analysis (MGA) contaminant screen for genomic sequence data"

RUN apt update && apt install -y git curl procps

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /opt/miniconda3.sh && \
    sh /opt/miniconda3.sh -b -p /opt/conda && \
    rm /opt/miniconda3.sh

ENV PATH ${PATH}:/opt/conda/bin

COPY conda.yml .

RUN conda env create -f conda.yml && conda clean -a

ENV PATH /opt/conda/envs/mga2/bin:$PATH

ADD "https://api.github.com/repos/crukci-bioinformatics/mga2/commits?per_page=1" /opt/mga2/.latest_commit

RUN git clone https://github.com/crukci-bioinformatics/mga2 /opt/mga2/build && \
    cargo install --path /opt/mga2/build

