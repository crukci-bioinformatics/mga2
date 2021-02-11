FROM rust:1.49-slim-buster

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk" \
           version "1.0" \
           description "Multi-Genome Analysis (MGA) contaminant screen for genomic sequence data"

RUN apt update && apt install -y \
    git \
    curl

RUN git clone https://github.com/crukci-bioinformatics/mga2 /opt/mga2 && \
    cargo install --path /opt/mga2

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /opt/miniconda3.sh && \
    sh /opt/miniconda3.sh -b -p /opt/conda && \
    rm /opt/miniconda3.sh

ENV PATH /opt/conda/bin:$PATH

COPY conda.yml .

RUN conda env update -n root -f conda.yml && conda clean -a

