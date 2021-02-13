FROM rust:1.49-slim-buster

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk" \
           version "2.0.0" \
           description "Multi-Genome Analysis (MGA) contaminant screen for genomic sequence data"

RUN apt update && apt install -y \
    git \
    curl \
    gdebi-core \
    libxml2-dev \
    libssl-dev

RUN git clone https://github.com/crukci-bioinformatics/mga2 /opt/mga2 && \
    cargo install --path /opt/mga2

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /opt/miniconda3.sh && \
    sh /opt/miniconda3.sh -b -p /opt/conda && \
    rm /opt/miniconda3.sh

ENV PATH /opt/conda/bin:${PATH}

COPY conda.yml .

RUN conda env update -n root -f conda.yml && conda clean -a

ENV R_VERSION=4.0.3

RUN curl -O https://cdn.rstudio.com/r/debian-10/pkgs/r-${R_VERSION}_1_amd64.deb && \
    gdebi -n r-${R_VERSION}_1_amd64.deb && \
    rm r-${R_VERSION}_1_amd64.deb

ENV PATH /opt/R/${R_VERSION}/bin:${PATH}

RUN R -e 'install.packages("tidyverse", repos = "https://cloud.r-project.org")'

RUN R -e 'install.packages("optparse", repos = "https://cloud.r-project.org")'

