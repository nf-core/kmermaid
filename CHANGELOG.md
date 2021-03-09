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
* Add `sourmash sig merge` for aligned/unaligned signatures from bam files, and add `--skip_sig_merge` option to turn it off
* Add `--protein_fastas` option for creating sketches of already-translated protein sequences
* Add `--skip_compare option` to skip `sourmash_compare_sketches` process
* Add merging of aligned/unaligned parts of single-cell data ([#117](https://github.com/nf-core/kmermaid/pull/117))
* Add renamed package dependency orpheum (used to be known as sencha)

### `Fixed`

#### Resources

* Increase CPUs in `high_memory_long` profile from 1 to 10

#### Naming

* Rename splitkmer to `split_kmer`

#### Per-cell fastqs and bams

* Remove `one_signature_per_record` flag and add bam2fasta count_umis_percell and make_fastqs_percell instead of bam2fasta sharding method
* Use ripgrep instead of bam2fasta to make per-cell fastq, which will hopefully make resuming long-running pipelines on bams much faster
* Make sure `samtools_fastq_aligned` outputs ALL aligned reads, regardless of mapping quality or primary alignment status

#### Sourmash

* add `--skip_compute option` to skip `sourmash_compute_sketch_*`
* Used `.combine()` instead of `each` to do cartesian product of all possible molecules, ksizes, and sketch values
* Do `sourmash compute` on all input ksizes, and all peptide molecule types, at once to save disk reading/writing efforts

#### Translate

* Updated sencha=1.0.3 to fix the bug in memory errors possibly with the numpy array on unique filenames ([PR #96 on orpheum](https://github.com/czbiohub/orpheum/pull/96))
* Add option to write non-coding nucleotide sequences fasta files while doing sencha translate
* Don't save translate csvs and jsons by default, add separate `--save_translate_json` and `--save_translate_csv`
* Updated `sencha translate` default parameters to be `--ksize 8 --jaccard-threshold 0.05` because those were the most successful
* Update renaming of `khtools` commands to `sencha`

#### MultiQC

* Fix the use of `skip_multiqc` flag condition with if and not when

### `Dependencies`

### `Deprecated`

* Removed ability to specify multiple `--scaled` or `--num-hashes` values to enable merging of signatures
