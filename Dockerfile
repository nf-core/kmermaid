FROM nfcore/base
LABEL authors="Olga Botvinnik" \
      description="Docker image containing all requirements for nf-core/kmermaid pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-kmermaid-1.1dev/bin:$PATH
