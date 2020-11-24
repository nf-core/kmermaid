# nf-core/kmermaid: Output

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/kmermaid/output](https://nf-co.re/kmermaid/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [MultiQC](#multiqc) - Aggregate report describing results from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
* [Sourmash Sketch](#sourmash-sketch) - Compute a k-mer sketch of each sample
* [Sourmash Compare](#sourmash-compare) - Compare all samples on k-mer sketches

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `fastp` directory.

**Output files:**

* `fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
* `fastqc/zips/`
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

## Sourmash Sketch

[Sourmash](https://sourmash.readthedocs.io/en/latest/) is a tool to compute MinHash sketches on nucleotide (DNA/RNA) and protein sequences. It allows for fast comparisons of sequences based on their nucleotide content.

**Output directory: `results/sourmash/sketches`**

For each sample and provided `molecules`, `ksizes` and `sketch_num_hashes_log2`, a file is created:

* `sample_molecule-${molecule}__ksize-${ksize}__${sketch_value}__track_abundance-${track_abundance}.sig`

For example:

```bash
SRR4050379_molecule-dayhoff_ksize-3_sketch_num_hashes_log2-2.sig
SRR4050379_molecule-dayhoff_ksize-3_sketch_num_hashes_log2-4.sig
SRR4050379_molecule-dayhoff_ksize-9_sketch_num_hashes_log2-2.sig
SRR4050379_molecule-dayhoff_ksize-9_sketch_num_hashes_log2-4.sig
SRR4050379_molecule-dna_ksize-3_sketch_num_hashes_log2-2.sig
SRR4050379_molecule-dna_ksize-3_sketch_num_hashes_log2-4.sig
SRR4050379_molecule-dna_ksize-9_sketch_num_hashes_log2-2.sig
SRR4050379_molecule-dna_ksize-9_sketch_num_hashes_log2-4.sig
SRR4050379_molecule-protein_ksize-3_sketch_num_hashes_log2-2.sig
SRR4050379_molecule-protein_ksize-3_sketch_num_hashes_log2-4.sig
SRR4050379_molecule-protein_ksize-9_sketch_num_hashes_log2-2.sig
SRR4050379_molecule-protein_ksize-9_sketch_num_hashes_log2-4.sig
```

## Sourmash Compare

**Output directory: `results/compare_sketches`**

For each provided `molecules`, `ksizes` and `sketch_num_hashes_log2`, a file is created containing a symmetric matrix of the similarity between all samples, written as a comma-separated variable file:

* `similarities_molecule-${molecule}__ksize-${ksize}__${sketch_value}__track_abundance-${track_abundance}.csv`
For example,

```bash
similarities_molecule-dna_ksize-9_sketch_num_hashes_log2-4.csv
similarities_molecule-protein_ksize-3_sketch_num_hashes_log2-2.csv
similarities_molecule-protein_ksize-3_sketch_num_hashes_log2-4.csv
similarities_molecule-protein_ksize-9_sketch_num_hashes_log2-2.csv
similarities_molecule-protein_ksize-9_sketch_num_hashes_log2-4.csv
```

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.
