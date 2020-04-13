
# nf-core/kmermaid: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Sourmash sketch](#sourmash-sketch) - Compute a k-mer sketch of each sample
* [Sourmash compare](#sourmash-compare) - Compare all samples on k-mer sketches
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## Sourmash Sketch

[Sourmash](https://sourmash.readthedocs.io/en/latest/) is a tool to compute MinHash sketches on nucleotide (DNA/RNA) and protein sequences. It allows for fast comparisons of sequences based on their nucleotide content.

**Output directory: `results/sourmash/sketches`**

For each sample and provided `molecule`, `ksize` and `log2_sketch_size`, a file is created:

* `sample_molecule-$molecule_ksize-$ksize_log2sketchsize-$log2_sketch_size.sig`

For example:

```bash
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

## Sourmash Compare

**Output directory: `results/sourmash`**

For each provided `molecule`, `ksize` and `log2_sketch_size`, a file is created containing a symmetric matrix of the similarity between all samples, written as a comma-separated variable file:

* `molecule-$molecule_ksize-$ksize_log2sketchsize-$log2_sketch_size.csv`

For example,

```bash
similarities_molecule-dna_ksize-3_log2sketchsize-2.csv
similarities_molecule-dna_ksize-3_log2sketchsize-4.csv
similarities_molecule-dna_ksize-9_log2sketchsize-2.csv
similarities_molecule-dna_ksize-9_log2sketchsize-4.csv
similarities_molecule-protein_ksize-3_log2sketchsize-2.csv
similarities_molecule-protein_ksize-3_log2sketchsize-4.csv
similarities_molecule-protein_ksize-9_log2sketchsize-2.csv
similarities_molecule-protein_ksize-9_log2sketchsize-4.csv
```

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
