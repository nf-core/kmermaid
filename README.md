
# nf-core/kmermaid

**Compare DNA/RNA/protein sequences on k-mer content**.

[![Build Status](https://travis-ci.com/nf-core/kmermaid.svg?branch=master)](https://travis-ci.com/nf-core/kmermaid)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/kmermaid.svg)](https://hub.docker.com/r/nfcore/kmermaid)

## Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.


## Documentation
The nf-core/kmermaid pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Usage

### With a samples.csv file:

```
nextflow run nf-core/kmermaid --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv
```

### With R1, R2 read pairs:

```
nextflow run nf-core/kmermaid --outdir s3://olgabot-maca/nf-kmer-similarity/ \
  --read_pairs 's3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq/*{R1,R2}*.fastq.gz,s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/*{1,2}.fastq.gz'
```

### With SRA ids:

```
nextflow run nf-core/kmermaid --outdir s3://olgabot-maca/nf-kmer-similarity/ --sra SRP016501
```

### With fasta files:

```
nextflow run nf-core/kmermaid --outdir s3://olgabot-maca/nf-kmer-similarity/ \
  --fastas '*.fasta'
```

### With bam file

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ \
  --bam 'possorted_genome_bam.bam'
```

## Credits
nf-core/kmermaid was originally written by Olga Botvinnik with contributions from Pranathi Vemuri.

