# nf-kmer-similarity

This is a [Nextflow](nextflow.io) workflow for running k-mer similarity

[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/czbiohub/nf-kmer-similarity.svg)](https://cloud.docker.com/u/czbiohub/repository/docker/czbiohub/nf-kmer-similarity)

## Usage

### With a samples.csv file:

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv
```

### With multitple s3 directories:

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --directories s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq,s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/
```

### With SRA ids:

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --sra SRP016501
```
