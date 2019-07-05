FROM continuumio/anaconda3
MAINTAINER olga.botvinnik@czbiohub.org

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="e.g. https://github.com/czbiohub/nf-kmer-similarity"


WORKDIR /home

USER root

# Add user "main" because that's what is expected by this image
RUN useradd -ms /bin/bash main


ENV PACKAGES zlib1g git g++ make ca-certificates gcc zlib1g-dev libc6-dev procps

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN conda install --yes Cython bz2file pytest numpy matplotlib scipy sphinx alabaster

RUN cd /home && \
    git clone https://github.com/dib-lab/khmer.git -b master && \
    cd khmer && \
    python3 setup.py install

# Check that khmer was installed properly
RUN trim-low-abund.py --help
RUN trim-low-abund.py --version

# RUN conda install --channel bioconda --yes sourmash

# Required for multiprocessing of 10x bam file
# RUN pip install pathos bamnostic

ENV SOURMASH_VERSION 'olgabot/dayhoff'
RUN cd /home && \
    git clone --branch $SOURMASH_VERSION https://github.com/czbiohub/sourmash.git && \
    cd sourmash && \
    python3 setup.py install

RUN which -a python3
RUN python3 --version
RUN sourmash info
COPY docker/sysctl.conf /etc/sysctl.conf
