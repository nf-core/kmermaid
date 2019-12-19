# nf-core/kmermaid: Changelog

## v1.0.0dev

* Add option to use Dayhoff encoding for sourmash
* Add bam file process to kmermaid pipeline and optional params
  for the same
* Add `track_abundance` feature to keep track of hashed kmer frequency.

### Dependency updates

* Add fastqc to environment.yml
* Add [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo to environment.yml
* Update Dockerfile with sourmash compute bam input dependencies
* Add [`bam2fasta`](https://pypi.org/project/bam2fasta/) to environment.yml
* Use Olga's branch of sourmash
* Add `ska` and `seqtk` to container dependencies
* Add `fastp` to container requirements
* Add `ska` and `seqtk` to container dependencies

## v0.1.0 - 6 March 2019

Initial release of nf-core/kmermaid, created with the [nf-core](http://nf-co.re/) template.


## v1.0dev - 6 March 2019
Initial release of nf-core/nf-kmer-similarity, created with the [nf-core](http://nf-co.re/) template.
