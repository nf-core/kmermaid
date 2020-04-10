FROM nfcore/base:1.7
LABEL description="Docker image containing all requirements for nf-core/kmermaid pipeline"

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
ARG BUILD_DATE
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="e.g. https://github.com/nf-core/kmermaid"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-kmermaid-1.0.0dev/bin:$PATH

RUN sourmash info
RUN bam2fasta info
RUN khtools index --help
RUN khtools translate --help
COPY docker/sysctl.conf /etc/sysctl.conf
