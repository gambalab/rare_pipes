# $ git clone https://github.com/gambalab/honey_pipes
# $ cd honey_pipes
# $ sudo docker build -f ./Dockerfile -t gambalab/rare .

# stage 1: build GLnexus
FROM ubuntu:18.04 AS builder
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive
ARG build_type=Release

# dependencies
RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     curl wget ca-certificates git-core less netbase \
     g++ cmake autoconf make file valgrind \
     libjemalloc-dev libzip-dev libsnappy-dev libbz2-dev zlib1g-dev liblzma-dev libzstd-dev \
     python3-pyvcf bcftools pv

# Copy in the local source tree / build context
ADD ./GLnexus /GLnexus
WORKDIR /GLnexus

# compile GLnexus
RUN cmake -DCMAKE_BUILD_TYPE=$build_type . && make -j4

# Stage 2: miniconda envs and tools
FROM continuumio/miniconda3:latest AS base
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        file \
        default-jre \
        wget \
        libjemalloc2 \
        bcftools \
        tabix \
        pv \
        littler && \
        apt-get autoremove && \
        apt-get clean && \
    rm -rf /var/lib/apt/lists


RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

RUN conda create -y -n bio \
                    bioconda::bcftools=1.20 \
                    bioconda::tabix=0.2.6 \
                    bioconda::bedtools=2.31.1 \
                    bioconda::vcf2tsvpy=0.6.1 \
                    bioconda::rtg-tools=3.12.1 \
                    && conda clean -a
# Install Maverick
RUN conda create -y -n maverick python=3.7 \
        && conda init \
        && . ~/.bashrc \
        && conda activate maverick \
        && conda install -y pip \
        && pip install pandas matplotlib scikit-learn scipy biopython tensorflow==2.7 transformers tf-models-official==2.7 \
        && pip install numpy --upgrade \
        && conda deactivate \
        && conda clean -a

WORKDIR /opt/maverick
COPY ./Maverick_resources.tar.gz /opt/
RUN tar -xvzf /opt/Maverick_resources.tar.gz -C /opt/maverick/
RUN rm /opt/Maverick_resources.tar.gz
COPY ./Maverick /opt/maverick/Maverick
RUN rm /opt/maverick/Maverick/InferenceScripts/*.sh
COPY ./scripts/*.* /opt/maverick/Maverick/InferenceScripts
RUN chmod +x /opt/maverick/Maverick/InferenceScripts/*.sh
RUN chmod +x /opt/maverick/Maverick/InferenceScripts/postprocess_results.R

# Install snpEff
COPY ./snpEff /opt/snpEff
RUN chmod +x /opt/snpEff/scripts/filter_variants.sh

# Install Picard
WORKDIR /opt/picard
RUN wget https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar
RUN chmod +x /opt/picard/picard.jar

# Copy Rdata
COPY ./Rdata /opt/rare/Rdata

# install CADA
COPY ./CADA /opt/CADA
WORKDIR /opt/CADA
RUN conda create -y -n cada python=3.7 \
        && conda init \
        && . ~/.bashrc \
        && conda activate cada \
        && pip install -e . \
        && conda deactivate \
        && conda clean -a
COPY ./scripts/prioritizing.py /opt/CADA/src/CADA/
RUN chmod +x /opt/maverick/Maverick/InferenceScripts/add_CADA.R
RUN chmod +x /opt/maverick/Maverick/InferenceScripts/clean_HPO.R

# Setup GLnexus
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2
COPY --from=builder /GLnexus/glnexus_cli /usr/local/bin/
ADD https://github.com/mlin/spVCF/releases/download/v1.0.0/spvcf /usr/local/bin/
RUN chmod +x /usr/local/bin/spvcf

ENV PATH="${PATH}":/opt/maverick/Maverick/InferenceScripts/:/opt/snpEff/scripts/:opt/picard