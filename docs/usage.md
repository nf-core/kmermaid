# nf-core/kmermaid: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
        * [`docker`](#docker)
        * [`awsbatch`](#awsbatch)
        * [`standard`](#standard)
        * [`none`](#none)
    * [Read inputs](#read-inputs)
        * [`--read_pairs`](#--read_pairs)
        * [`--read_singles`](#--read_singles)
        * [`--csv_pairs`](#--csv_pairs)
        * [`--csv_singles`](#--csv_singles)
        * [`--fastas`](#--fastas)
        * [`--sra`](#--sra)
    * [Sketch parameters](#sketch-parameters)
        * [`--molecule`](#--molecule)
        * [`--ksize`](#--ksize)
        * [`--log2_sketch_size`](#--log2_sketch_size)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
    * [`-awsbatch`](#-awsbatch)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_emails`](#--plaintext_emails)
    * [`--sampleLevel`](#--sampleLevel)


## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run nf-core/nf-kmer-similarity --reads '*_R{1,2}.fastq.gz' -profile standard,docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/nf-kmer-similarity
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/nf-kmer-similarity releases page](https://github.com/nf-core/nf-kmer-similarity/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main Arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile standard,docker` - the order of arguments is important!

* `standard`
    * The default profile, used if `-profile` is not specified at all.
    * Runs locally and expects all software to be installed and available on the `PATH`.
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`nfcore/nf-kmer-similarity`](http://hub.docker.com/r/nfcore/nf-kmer-similarity/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `awsbatch`
    * A generic configuration profile to be used with AWS Batch.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

## Read inputs

This pipeline can take a large variety of input data, ranging from single-end or paired-end FASTQ files (`fastq.gz` totally cool, too), [SRA](https://www.ncbi.nlm.nih.gov/sra) ids or fasta files.

### `--read_pairs`
Use this to specify the location of your input *paired-end* FastQ files. Multiple paths can be separated by seimcolons (`;`). For example:

```bash
--read_pairs 'path/to/data/sample_*_{1,2}.fastq.gz;more/data/sample_*_{1,2}.fastq.gz'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--read_singles`

Use this to specify the location of your input *single-end* FastQ files.  Multiple paths can be separated by seimcolons (`;`). For example:

```bash
--read_singles 'path/to/data/sample_*.fastq;more/data/sample*.fastq.gz'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character

If left unspecified, no reads are used.

### `--csv_pairs`
Use this to specify the location of a csv containing the columns `sample_id`, `read1`, `read2` to your input *paired-end* FastQ files. For example:

```bash
--csv_pairs samples.csv
```

Please note the following requirements:

1. The CSV contains a header, typically `sample_id`, `read1`, `read2`, but the names aren't used.
1. The first column contains a sample ID
1. The second column contains the full path to the R1 read
1. The third column contains the full path to the R2 read

If left unspecified, no reads are used.

### `--csv_singles`
Use this to specify the location of a csv containing the columns `sample_id`, `read1` to your input *single-end* FastQ files. For example:

```bash
--csv_singles samples.csv
```

Please note the following requirements:

1. The CSV contains a header, typically `sample_id`, `read1` but the names aren't used.
1. The first column contains a sample ID
1. The second column contains the full path to the R1 read

If left unspecified, no reads are used.

### `--fastas`

Use this to specify the location of fasta sequence files. Multiple inputs can be separated by seimcolons (`;`). For example:

```bash
--fastas 'path/to/data/elephant.fasta;more/data/*.fasta'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path *may* have at one or more `*` wildcard character

If left unspecified, no samples are used.

### `--sra`

Use this to specify the location of fasta sequence files. Multiple inputs can be separated by seimcolons (`;`). For example:

```bash
--sra 'SRR4050379;SRR4050380;SRP016501'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. Any of the `SRR`, `SRP`, or `PRJNA` ids can be used

If left unspecified, no samples are used.

## Sketch parameters

[K-mer](https://en.wikipedia.org/wiki/K-mer) [MinHash](https://en.wikipedia.org/wiki/MinHash) sketches are defined by three parameters:

1. The [molecule](#--molecule) used to create k-mer words from each sample
1. The [ksize](#--ksize) used to extract k-mer words from each sample
1. The number of k-mer words specified by the [log2 sketch size](#--log2_sketch_size)

### `--molecule`

The molecule can be either `dna`, `protein`, or `dayhoff`, and if all of them are desired, then they can be separated by columns.

* `dna` indicates to use the raw nucleotide sequence from each input file to create k-mers
* `protein` indicates to translate each DNA k-mer into protein using [6-frame translation](https://en.wikipedia.org/wiki/Reading_frame#/media/File:Open_reading_frame.jpg), and hash the translated peptide fragment k-mers
* `dayhoff`, developed by [Margaret Oakley Dayhoff](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff) is an extension of `protein`, where in addition to translating each amino acid k-mer in six frames, each amino acid is remapped to a degenerate amino acid encoding. This degenerate encoding doesn't penalize large evolutionary distances as a k-mer changing. For example, an small residue change of an Alanine (`A`) to a Glycine (`G`) doesn't hash to the same value and thus do not match, but a Dayhoff-encoded amino acid would allow for these small changes.

| Amino acid    | Property              | Dayhoff |
|---------------|-----------------------|---------|
| C             | Sulfur polymerization | a       |
| A, G, P, S, T | Small                 | b       |
| D, E, N, Q    | Acid and amide        | c       |
| H, K, R       | Basic                 | d       |
| I, L, M, V    | Hydrophobic           | e       |
| F, W, Y       | Aromatic              | f       |

**Example parameters**

* Default:
  * `--molecule dna,protein,dayhoff`
* DNA only:
  * `--molecule dna`


### `--ksize`

The fundamental unit of the sketch is a [hashed](https://en.wikipedia.org/wiki/Hash_function) [k-mer](https://en.wikipedia.org/wiki/K-mer). The `--ksize` parameter determines how long of a DNA word to use to create the k-mers. When the  molecule is `protein` or `dayhoff`, then the `ksize/3` is used to create each k-mer.

*NB: if either `protein` or `dayhoff` is specified, the k-mer size must be divisible by 3*

**Example parameters**

* Default:
  * `--ksize 21,27,33,51`
* k-mer size of 21 only:
  * `--ksize 21`

### `--log2_sketch_size`

The log2 sketch size specifies the number of k-mers to use for the sketch. We use the log2 of the sketch size instead of the raw number of k-mers to be compatible for comparison with [`dashing`](https://github.com/dnbaker/dashing) that uses HyperLogLog instead of MinHash.

**Example parameters**

* Default:
  * `--log2_sketch_size 10,12,14,16`
* Only a log2 sketch size of 8 (2^8 = 256):
  * `--log2_sketch_size 8`

### `--fasta`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.


**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, you can specify a config file using `-c` that contains the following:

```nextflow

```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`
Set to disable colourful command line output and live life in monochrome.
