FROM nfcore/base
LABEL description="Docker image containing all requirements for czbiohub/nf-kmer-similarity pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-crisprvar-1.0dev/bin:$PATH

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="e.g. https://github.com/czbiohub/nf-kmer-similarity"


WORKDIR /home


ENV SOURMASH_VERSION 'olgabot/dayhoff'
RUN cd /home && \
    git clone --branch $SOURMASH_VERSION https://github.com/czbiohub/sourmash.git && \
    cd sourmash && \
    python3 setup.py install

RUN sourmash info
COPY docker/sysctl.conf /etc/sysctl.conf
