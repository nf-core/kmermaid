FROM nfcore/base
MAINTAINER olga.botvinnik@czbiohub.org

# Suggested tags from https://microbadger.com/labels
ARG VCS_REF
LABEL org.label-schema.vcs-ref=$VCS_REF \
org.label-schema.vcs-url="https://github.com/czbiohub/nf-kmer-similarity"


WORKDIR /home

USER root

# Add user "main" because that's what is expected by this image
# RUN useradd -ms /bin/bash main


# ENV PACKAGES zlib1g git g++ make ca-certificates gcc zlib1g-dev libc6-dev procps

### don't modify things below here for version updates etc.

WORKDIR /home

# RUN apt-get update && \
#     apt-get install -y --no-install-recommends ${PACKAGES} && \
#     apt-get clean

# Set always yes
RUN conda config --set always_yes yes --set changeps1 no

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN conda update --yes -n base -c defaults conda && conda init $(basename $SHELL) && exec $SHELL
# Use conda to install khmer and sourmash scientific python dependencies
# RUN conda install --yes Cython bz2file pytest numpy matplotlib scipy sphinx alabaster khmer

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/czbiohub-nf-kmer-similarity-0.1/bin:$PATH

# Check that khmer was installed properly
RUN trim-low-abund.py --help
RUN trim-low-abund.py --version

# Required for multiprocessing of 10x bam file
# RUN pip install pathos bamnostic

# ENV SOURMASH_VERSION master
RUN cd /home && \
    git clone https://github.com/dib-lab/sourmash.git && \
    cd sourmash && \
    pip install .

RUN which -a python3
RUN python3 --version
RUN sourmash info
COPY docker/sysctl.conf /etc/sysctl.conf

# Copy utility scripts to docker image
COPY bin/* /usr/local/bin/
