# hash:sha256:319c6123ec71ed086b4a6ce6f20f0b833843eae9b12bfbfba95f87bd788e3d96
FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bowtie2=2.3.4.1-1 \
        build-essential=12.4ubuntu1 \
        default-jdk=2:1.11-68ubuntu1~18.04.1 \
        default-jre=2:1.11-68ubuntu1~18.04.1 \
        libgit2-dev=0.26.0+dfsg.1-1.1ubuntu0.2 \
        libssl-dev=1.1.1-1ubuntu2.1~18.04.6 \
        pkg-config=0.29.1-0ubuntu2 \
        r-base=3.4.4-1ubuntu1 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y --channel=bioconda \
        biopython==1.72 \
        fastqc==0.11.7 \
        pandas==0.23.4 \
        plotly==3.2.1 \
        scipy==1.1.0 \
        trim-galore==0.6.4 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    cutadapt==2.4
