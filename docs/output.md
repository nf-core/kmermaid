# nf-core/nf-kmer-similarity: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.


## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Sourmash](#sourmash) - MinHash to subset the reads before comparing samples
  * [Sourmash sketch](#sourmash-sketch) - Compute a k-mer sketch of each sample
  * [Sourmash compare](#sourmash-compare) - Compare all samples on k-mer sketches
* [Split K-mer Analysis (SKA)](#split-k-mer-analysis-ska)
  * [SKA sketch](#ska-sketch) - Compute a k-mer sketch of each sample
  * [SKA compare](#ska-compare) - Compare all samples on k-mer sketches
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline


## Sourmash

[Sourmash](https://sourmash.readthedocs.io/en/latest/) is a tool to compute MinHash sketches on nucleotide (DNA/RNA) and protein sequences. It allows for fast comparisons of sequences based on their nucleotide content.

### Sourmash sketch

**Output directory: `results/sourmash/sketches`**

For each sample and provided `molecule`, `ksize` and `log2_sketch_size`, a file is created:

* `sample_molecule-$molecule_ksize-$ksize_log2sketchsize-$log2_sketch_size.sig`

For example:

```
SRR4050379_molecule-dayhoff_ksize-3_log2sketchsize-2.sig
SRR4050379_molecule-dayhoff_ksize-3_log2sketchsize-4.sig
SRR4050379_molecule-dayhoff_ksize-9_log2sketchsize-2.sig
SRR4050379_molecule-dayhoff_ksize-9_log2sketchsize-4.sig
SRR4050379_molecule-dna_ksize-3_log2sketchsize-2.sig
SRR4050379_molecule-dna_ksize-3_log2sketchsize-4.sig
SRR4050379_molecule-dna_ksize-9_log2sketchsize-2.sig
SRR4050379_molecule-dna_ksize-9_log2sketchsize-4.sig
SRR4050379_molecule-protein_ksize-3_log2sketchsize-2.sig
SRR4050379_molecule-protein_ksize-3_log2sketchsize-4.sig
SRR4050379_molecule-protein_ksize-9_log2sketchsize-2.sig
SRR4050379_molecule-protein_ksize-9_log2sketchsize-4.sig
```

### Sourmash compare

**Output directory: `results/sourmash/compare`**

For each provided `molecule`, `ksize` and `log2_sketch_size`, a file is created containing a symmetric matrix of the similarity between all samples, written as a comma-separated variable file:

* `molecule-$molecule_ksize-$ksize_log2sketchsize-$log2_sketch_size.csv`

For example,

```
similarities_molecule-dna_ksize-3_log2sketchsize-2.csv
similarities_molecule-dna_ksize-3_log2sketchsize-4.csv
similarities_molecule-dna_ksize-9_log2sketchsize-2.csv
similarities_molecule-dna_ksize-9_log2sketchsize-4.csv
similarities_molecule-protein_ksize-3_log2sketchsize-2.csv
similarities_molecule-protein_ksize-3_log2sketchsize-4.csv
similarities_molecule-protein_ksize-9_log2sketchsize-2.csv
similarities_molecule-protein_ksize-9_log2sketchsize-4.csv
```

## Split K-mer Analysis (SKA)

[Split K-mer analysis (SKA)](https://github.com/simonrharris/SKA) is a program to take ALL the reads from a sample and find split k-mers.

### SKA sketch

**Output directory: `results/ska/sketches`**



### SKA compare

**Output directory: `results/ska/compare`**



## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info
