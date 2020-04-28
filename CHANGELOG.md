# nf-core/kmermaid: Changelog

## v1.0.0dev

### Documentation updates

### Pipeline enhancements & Fixes

* Add option to use Dayhoff encoding for sourmash

### Dependency Updates

* Add samtools, screed, tqdm to dependencies

## v1.0 - 6 March 2019

* Add option to use Dayhoff encoding for sourmash.
* Add `bam2fasta` process to kmermaid pipeline and flags involved.
* Add `extract_coding` and `peptide_bloom_filter` process and flags involved.
* Add `track_abundance` feature to keep track of hashed kmer frequency.
* Add social preview image
* Add `fastp` process for trimming reads
* Add option to use compressed `.tgz` file containing output from 10X Genomics' `cellranger count` outputs, including `possorted_genome_bam.bam` and `barcodes.tsv` files
* Add samtools_fastq_unaligned and samtools_fastq_aligned process for converting bam to per cell
barcode fastq
* Remove `one_signature_per_record` flag and add bam2fasta count_umis_percell and make_fastqs_percell instead of bam2fasta sharding method
* Update renaming of `khtools` commands to `sencha`
* Make sure `samtools_fastq_aligned` outputs ALL aligned reads, regardless of mapping quality or primary alignment status

### Dependency updates

* Add `ska` and `seqtk` to container dependencies
* Add `fastp` to container requirements
* Add `fastqc` to environment.yml
* Add [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo to environment.yml (now renamed to [czbiohub/sencha](https://github.com/czbiohub/sencha/))
* Update Dockerfile with sourmash compute bam input dependencies
* Add `track_abundance` feature to keep track of hashed kmer frequency.
* Add [czbiohub/sencha](https://github.com/czbiohub/kh-tools/) repo to environment.yml
* Add [`czbiohub/bam2fasta`](https://github.com/czbiohub/bam2fasta/) repo to environment.yml
* Update sourmash to version 3.2.2

## v0.1.0 - 6 March 2019

Initial release of nf-core/kmermaid, created with the [nf-core](http://nf-co.re/) template.
