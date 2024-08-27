# $ git clone https://github.com/gambalab/honey_pipes
# $ cd honey_pipes
# $ sudo docker build -f ./Dockerfile -t gambalab/rare .

# Stage miniconda envs
FROM continuumio/miniconda3:latest AS base

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        file \
        default-jre \
        wget \
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

ENV PATH="${PATH}":/opt/maverick/Maverick/InferenceScripts/:/opt/snpEff/scripts/:opt/picard