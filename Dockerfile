FROM nfcore/base:1.12.1
LABEL authors="Olga Botvinnik" \
      description="Docker image containing all software requirements for the nf-core/kmermaid pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda install -c conda-forge mamba
RUN mamba env create -f /environment.yml && mamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-kmermaid-0.1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-kmermaid-0.1.0dev > nf-core-kmermaid-0.1.0dev.yml

# Install super fast rust code to remove nuisance hashes (e.g. ribosomal) from signatures
RUN git clone https://github.com/luizirber/2021-01-27-olga-remove-protein/ 
RUN cd 2021-01-27-olga-remove-protein  && cargo build --release 
# Add "subtract" command to path
ENV PATH $HOME/2021-01-27-olga-remove-protein/target/release:$PATH

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
