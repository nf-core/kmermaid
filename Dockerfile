FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/kmermaid pipeline"

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
ARG BUILD_DATE
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="e.g. https://github.com/nf-core/kmermaid"

WORKDIR /home

ENV PACKAGES zlib1g git g++ make ca-certificates gcc zlib1g-dev libc6-dev procps libbz2-dev libcurl4-openssl-dev libssl-dev

# Sourmash requires last 3 of the above apt packages for installation of htslib in pysam

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-kmermaid-0.1dev/bin:$PATH

RUN which -a pip
RUN which -a python
RUN cd /home && \
    git clone https://github.com/czbiohub/sourmash.git && \
    cd sourmash && \
    python setup.py install

RUN which -a sourmash

RUN sourmash info
COPY docker/sysctl.conf /etc/sysctl.conf