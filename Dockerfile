FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/kmermaid pipeline"

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
ARG BUILD_DATE
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="e.g. https://github.com/nf-core/kmermaid"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-kmermaid-0.1dev/bin:$PATH

RUN sourmash info
COPY docker/sysctl.conf /etc/sysctl.conf
# TODO: Change to be format below
# LABEL authors="Olga Botvinnik" \
#       description="Docker image containing all requirements for nf-core/kmermaid pipeline"

# COPY environment.yml /
# RUN conda env create -f /environment.yml && conda clean -a
# ENV PATH /opt/conda/envs/nf-core-kmermaid-1.1dev/bin:$PATH
