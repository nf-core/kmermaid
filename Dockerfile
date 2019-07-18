FROM nfcore/base
LABEL description="Docker image containing all requirements for czbiohub/nf-kmer-similarity pipeline"

RUN conda update -n base -c defaults conda
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/czbiohub-nf-kmer-similarity-0.1dev/bin:$PATH

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="e.g. https://github.com/czbiohub/nf-kmer-similarity"


RUN which -a sourmash

RUN which -a python3
RUN python3 --version
RUN sourmash info
COPY docker/sysctl.conf /etc/sysctl.conf
