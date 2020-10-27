
# nf-core/kmermaid: Output

This document describes the output produced by the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Sourmash sketch](#sourmash-sketch) - Compute a k-mer sketch of each sample
* [Sourmash compare](#sourmash-compare) - Compare all samples on k-mer sketches
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## Sourmash Sketch

[Sourmash](https://sourmash.readthedocs.io/en/latest/) is a tool to compute MinHash sketches on nucleotide (DNA/RNA) and protein sequences. It allows for fast comparisons of sequences based on their nucleotide content.

**Output directory: `results/sourmash/sketches`**

For each sample and provided `molecules`, `ksizes` and `sketch_num_hashes_log2`, a file is created:

* `sample_molecule-${molecule}__ksize-${ksize}__${sketch_value}__track_abundance-${track_abundance}.sig`

## Sourmash Compare

**Output directory: `results/compare_sketches`**

For each provided `molecules`, `ksizes` and `sketch_num_hashes_log2`, a file is created containing a symmetric matrix of the similarity between all samples, written as a comma-separated variable file:

* `similarities_molecule-${molecule}__ksize-${ksize}__${sketch_value}__track_abundance-${track_abundance}.csv`
