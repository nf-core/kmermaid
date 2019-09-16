FROM continuumio/anaconda3
MAINTAINER olga.botvinnik@czbiohub.org

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="e.g. https://github.com/czbiohub/nf-kmer-similarity"
ARG BUILD_DATE

WORKDIR /home

USER root

ENV PATH="/opt/conda/bin:${PATH}"

# Add user "main" because that's what is expected by this image
RUN useradd -ms /bin/bash main


ENV PACKAGES zlib1g git g++ make ca-certificates gcc zlib1g-dev libc6-dev procps libbz2-dev libcurl4-openssl-dev libssl-dev

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN conda install --yes Cython bz2file pytest numpy matplotlib scipy sphinx alabaster

RUN cd /home && \
    git clone https://github.com/dib-lab/khmer.git -b master && \
    cd khmer && \
    python3 setup.py install && \
    trim-low-abund.py --version && \
    trim-low-abund.py --help

RUN which -a python3
RUN python3 --version
# Required for multiprocessing of 10x bam file
RUN pip install pathos pysam

# ENV SOURMASH_VERSION master
RUN cd /home && \
    git clone https://github.com/dib-lab/sourmash.git && \
    cd sourmash && \
    python3 setup.py install

RUN which -a python3
RUN python3 --version
RUN sourmash info
COPY docker/sysctl.conf /etc/sysctl.conf
