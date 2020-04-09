# nf-core/kmermaid: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
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
        * [`--bam`](#--bam)
        * [`--barcodes_file`](#--barcodes_file)
        * [`--rename_10x_barcodes`](#--rename_10x_barcodes)
        * [`--save_fastas`](#--save_fastas)
    * [K-merization/Sketching program options](#k-merization-sketching-program-options)
        * [`--splitKmer`](#--splitKmer)

    * [Sketch parameters](#sketch-parameters)
        * [`--molecule`](#--molecule)
        * [`--ksize`](#--ksize)
        * [`--log2_sketch_size`](#--log2_sketch_size)
        * [`--track_abundance`](#--track_abundance)
    * [Bam optional parameters](#bam-optional-parameters)
        * [`--min_umi_per_barcode`](#--min_umi_per_barcode)
        * [`--write_barcode_meta_csv`](#--write_barcode_meta_csv)
        * [`--line_count`](#--line_count)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/kmermaid --reads '*_R{1,2}.fastq.gz' -profile docker
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
nextflow pull nf-core/kmermaid
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/kmermaid releases page](https://github.com/nf-core/kmermaid/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/kmermaid`](http://hub.docker.com/r/nfcore/kmermaid/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/kmermaid`](http://hub.docker.com/r/nfcore/kmermaid/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

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

## K-merization/Sketching program Options

By default, the k-merization and sketch creation program is [sourmash](https://sourmash.readthedocs.io).

### `--splitKmer`

If `--splitKmer` is specified, then the [Split K-mer Analysis (SKA)](https://github.com/simonrharris/SKA) program ([publication](https://www.biorxiv.org/content/10.1101/453142v1)) is used to obtain k-mers from the data. This allows for a SNP to be present in the middle of a k-mer which can be advantageous for metagenomic analyses or with single-cell ATAC-seq data.

#### What does `--ksize` mean when `--splitKmer` is set

The meaning of `ksize` is different with split k-mers, so now the value specified by `--ksize` is just under half of the total sampled sequence size, where the middle base can be any base (`N`) `[---ksize---]N[---ksize---]`. When `--splitKmer` is set, then the default k-mer sizes are 9 and 15, for a total sequence unit size of `2*15+1 = 31` and `2*9+1 = 19` which is as if you specified on the command line `--splitKmer --ksize 9,15`. Additionally k-mer sizes with `--splitKmer` must be divisible by 3 (yes, this is inconvenient) and between 3 and 60 (inclusive). So the "total" `2*k+1` sizes can be:

* k = 3 --> 2*3 + 1 = 7 total length
* k = 6 --> 2*6 + 1 = 13 total length
* k = 9 --> 2*9 + 1 = 18 total length
* k = 12 --> 2*12 + 1 = 25 total length
* k = 15 --> 2*15 + 1 = 31 total length
* ...
* k = 60 --> 2*60 + 1 = 121 total length

#### `--subsample` reads when `--splitKmer` is set

The `subsample` command is often necessary because the `ska` tool uses ALL the reads rather than a MinHash subsampling of them. If your input files are rather big, then the `ska` sketching command (`ska fastq`) runs out of memory, or it takes so long that it's untenable. The `--subsample` command specifies the number of reads to be used. When e.g. `--subsample 1000` is set, then 1000 reads (or read pairs) are randomly subsampled from the data using [seqtk](https://github.com/lh3/seqtk).


#### Which `--molecules` are valid when `--splitKmer` is set

Currently, `--splitKmer` only works with DNA sequence and not protein sequence, and thus will fail if `protein` or `dayhoff` is specified in `--molecules`.

### `--bam`
For bam/10x files, Use this to specify the location of the bam file. For example:

```bash
--bam /path/to/data/10x-example/possorted_genome_bam
```
### `--barcodes_file`
For bam/10x files, Use this to specify the location of tsv (tab separated file) containing cell barcodes. For example:

```bash
--barcodes_file /path/to/data/10x-example/barcodes.tsv
```

If left unspecified, barcodes are derived from bam are used.

### `--rename_10x_barcodes`
For bam/10x files, Use this to specify the location of your tsv (tab separated file) containing map of cell barcodes and their corresponding new names(e.g row in the tsv file: AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1).
For example:

```bash
--rename_10x_barcodes /path/to/data/10x-example/barcodes_renamer.tsv
```
If left unspecified, barcodes in bam as given in barcodes_file are not renamed.

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


### `--track_abundance`

* Tracking abundance - add this parameter if we want to keep track of the number of times a hashed kmer appears.
  * `--track_abundance`



### `--log2_sketch_size`

The log2 sketch size specifies the number of k-mers to use for the sketch. We use the log2 of the sketch size instead of the raw number of k-mers to be compatible for comparison with [`dashing`](https://github.com/dnbaker/dashing) that uses HyperLogLog instead of MinHash.

**Example parameters**

* Default:
  * `--log2_sketch_size 10,12,14,16`
* Only a log2 sketch size of 8 (2^8 = 256):
  * `--log2_sketch_size 8`


### `--save_fastas`

1. The [save_fastas ](#--save_fastas ) used to save the sequences of each unique barcode in the bam file. It is a path relative to outdir to save unique barcodes to files namely {CELL_BARCODE}.fasta. These fastas are computed once for one permutation of ksize, molecule, and log2_sketch_size, further used to compute the signatures and compare signature matrix for all permutations of ksizes, molecules, and log2_sketch_size. This is done to save the time on saving the computational time and storage in obtaining unique barcodes, sharding the bam file.

**Example parameters**

* Default: Save fastas in a directory called fastsas inside outdir:
  * `--save_fastas "fastas"`


## Bam optional parameters

### `--write_barcode_meta_csv`
This creates a CSV containing the number of reads and number of UMIs per barcode, written in a path relative to `${params.outdir}/barcode_metadata`. This csv file is empty with just header when the min_umi_per_barcode is zero i.e reads and UMIs per barcode are calculated only when the barcodes are filtered based on [min_umi_per_barcode](#--min_umi_per_barcode)

**Example parameters**

* Default: barcode metadata is not saved
* Save fastas in a file cinside outdir/barcode/metadata:
  * `--write_barcode_meta_csv "barcodes_counts.csv"`


### `--min_umi_per_barcode`
The parameter `--min_umi_per_barcode` ensures that a barcode is only considered a valid barcode read and its sketch is written if number of unique molecular identifiers (UMIs, aka molecular barcodes) are greater than the value specified.

**Example parameters**

* Default: min_umi_per_barcode is 0
* Set minimum UMI per cellular barcode as 10:
  * `--min_umi_per_barcode 10`


### `--line_count`
The parameter `--line_count` specifies the number of alignments/lines in each bam shard.
**Example parameters**

* Default: line_count is 350
* Save fastas in a directory called fastas inside outdir:
  * `--line_count 400`


## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

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

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
