# nf-core/kmermaid: Changelog

## v1.0.0dev

* Add option to use Dayhoff encoding for sourmash
* Add bam file process to kmermaid pipeline and optional params for the same
* Add social preview image

### Dependency updates

* Add `ska` and `seqtk` to container dependencies
* Add `fastp` to container requirements
* Add `fastqc` to environment.yml
* Add [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo to environment.yml
* Update Dockerfile with sourmash compute bam input dependencies
* Add "fastp" to container requirements
* Update Dockerfile with sourmash compute bam input dependencies
* Add `track_abundance` feature to keep track of hashed kmer frequency.
* Add [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo to environment.yml
* Add [`czbiohub/bam2fasta`](https://github.com/czbiohub/bam2fasta/) repo to environment.yml

## v0.1.0 - 6 March 2019

Initial release of nf-core/kmermaid, created with the [nf-core](http://nf-co.re/) template.

## v1.0dev - 6 March 2019

Initial release of nf-core/nf-kmer-similarity, created with the [nf-core](http://nf-co.re/) template.
