# nf-kmer-similarity

This is a [Nextflow](nextflow.io) workflow for running k-mer similarity

[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/czbiohub/nf-kmer-similarity.svg)](https://cloud.docker.com/u/czbiohub/repository/docker/czbiohub/nf-kmer-similarity)

## Usage

### With a samples.csv file:

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv
```

### With R1, R2 read pairs:

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ \
  --read_pairs 's3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq/*{R1,R2}*.fastq.gz,s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/*{1,2}.fastq.gz'
```

### With SRA ids:

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --sra SRP016501
```

### With fasta files:

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ \
  --fastas '*.fasta'
```

### With bam file

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ \
  --bam 'possorted_genome_bam.bam' --one_signature_per_record
```