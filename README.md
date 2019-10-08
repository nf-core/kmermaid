# nf-kmer-similarity

This is a [Nextflow](nextflow.io) workflow for running k-mer similarity.

[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/czbiohub/nf-kmer-similarity.svg)](https://cloud.docker.com/u/czbiohub/repository/docker/czbiohub/nf-kmer-similarity)

## Usage

By default, this pipeline creates a [MinHash](https://en.wikipedia.org/wiki/MinHash) sketch of sequencing reads using [sourmash](https://sourmash.readthedocs.io), then compares them all using a [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) . Here are the default parameters:

- log2 sketch sizes of 10, 12, 14, 16  (as if `--log2_sketch_sizes 10,12,14,16` was specified on the command line), so 2^10, 2^12, 2^14, 2^16 = 1024, 4096, 16 384, 65 536 hashed k-mers for each sample
- Compute both DNA and protein signatures (as if `--molecules dna,protein` was specified on the command line). The protein k-mers are obtained by doing [six-frame translation](https://en.wikipedia.org/wiki/Reading_frame#/media/File:Open_reading_frame.jpg) on the DNA k-mers
- K-mer sizes of 21, 27, 33, 51 (as if `--ksizes 21,27,33,51` was specified on the command line).
  - If using the `--splitKmer` option, keep in mind that the k-mer size in this case is the two halves of the split k-mer, which you can visualize as `[---ksize---]N[---ksize---]`. So the default k-mer sizes for `--splitKmer` is 9 and 15, for a total sequence unit size of `2*15+1 = 31` and `2*9+1 = 19` which is as if you specified on the command line `--splitKmer --ksize 9,15`. Additionally k-mer sizes with `--splitKmer` must be divisible by 3 (yes, this is inconvenient)

### With a samples.csv file:

This is where you'd have a csv file with a `sample_id,read1,read2` header containing the sample id and paths to each of your R1 and R2 read files.

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

### With Split Kmer Analysis [SKA](https://github.com/simonrharris/SKA):

Note: the meaning of `ksize` is different with split k-mers, so now the value specified by `--ksize` is just under half of the total sampled sequence size, where the middle base can be any base (`N`) `[---ksize---]N[---ksize---]`. Note that `--splitKmer` can only work with DNA sequence and does not work with `protein` specified in `--molecules`.

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv --splitKmer --subsample 1000
```

### With Split Kmer Analysis [SKA](https://github.com/simonrharris/SKA) and fastq subsampling with [seqtk](https://github.com/lh3/seqtk):

The `subsample` command is often necessary because the `ska` tool uses ALL the reads rather than a MinHash subsampling of them. If your input files are rather big, then the `ska` sketching command (`ska fastq`) runs out of memory, or it takes so long that it's untenable. The `--subsample` command specifies the number of reads to be used.

```
nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv --splitKmer --subsample 1000
```
