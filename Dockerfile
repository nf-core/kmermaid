<<<<<<< HEAD
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
=======
FROM nfcore/base:dev
LABEL authors="Olga Botvinnik" \
      description="Docker image containing all software requirements for the nf-core/kmermaid pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-kmermaid-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-kmermaid-1.0dev > nf-core-kmermaid-1.0dev.yml
>>>>>>> upstream/TEMPLATE
