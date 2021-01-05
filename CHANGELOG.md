# nf-core/kmermaid: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - [date]

Initial release of nf-core/kmermaid, created with the [nf-core](https://nf-co.re/) template.

### `Added`

* Add option to use Dayhoff encoding for sourmash.
* Add `bam2fasta` process to kmermaid pipeline and flags involved.
* Add `extract_coding` and `peptide_bloom_filter` process and flags involved.
* Add `track_abundance` feature to keep track of hashed kmer frequency.
* Add social preview image
* Add `fastp` process for trimming reads
* Add option to use compressed `.tgz` file containing output from 10X Genomics' `cellranger count` outputs, including `possorted_genome_bam.bam` and `barcodes.tsv` files
* Add samtools_fastq_unaligned and samtools_fastq_aligned process for converting bam to per cell
barcode fastq
* Add version printing for sencha, bam2fasta, and sourmash in Dockerfile, update versions in environment.yml
* For processes translate, sourmash compute  add cpus=1 as they are only serial ([#107](https://github.com/nf-core/kmermaid/pull/107))

### `Fixed`

* Remove `one_signature_per_record` flag and add bam2fasta count_umis_percell and make_fastqs_percell instead of bam2fasta sharding method
* Update renaming of `khtools` commands to `sencha`
* Make sure `samtools_fastq_aligned` outputs ALL aligned reads, regardless of mapping quality or primary alignment status
* Add `--protein_fastas` option for translated protein input
* Rename splitkmer to `split_kmer` and add `--skip_compare option` to skip `sourmash_compare_sketches` process
* Increase CPUs in `high_memory_long` profile from 1 to 10
* add `--skip_compute option` to skip `sourmash_compute_sketch_*`
* add option to write non-coding nucleotide sequences fasta files while doing sencha translate
* Used `.combine()` instead of `each` to do cartesian product of all possible molecules, ksizes, and sketch values
* Use ripgrep instead of bam2fasta to make per-cell fastq, which will hopefully make resuming long-running pipelines on bams much faster
* Fix the use of `skip_multiqc` flag condition with if and not when
* Updated sencha=1.0.3 to fix the bug in memory errors possibly with the numpy array on unique filenames ([PR #96 on sencha](https://github.com/czbiohub/leaftea/pull/96))
* Do `sourmash compute` on all input ksizes, and all peptide molecule types, at once to save disk reading/writing efforts

### `Dependencies`

### `Deprecated`
