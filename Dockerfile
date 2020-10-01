FROM nfcore/base:1.10.2
LABEL authors="Olga Botvinnik" \
      description="Docker image containing all software requirements for the nf-core/kmermaid pipeline"

COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-kmermaid-1.0.0dev > nf-core-kmermaid-1.0.0dev.yml
ENV PATH /opt/conda/envs/nf-core-kmermaid-1.0.0dev/bin:$PATH
RUN sourmash info
RUN bam2fasta info
RUN sencha index --help
RUN sencha translate --help
