# nf-core/nf-kmer-similarity: Changelog

## v1.1dev

* Add option to use Dayhoff encoding for sourmash
* Add bam file process to kmermaid pipeline and optional params
  for the same

### Dependency updates
* Add fastqc to environment.yml
* Add [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo to environment.yml
* Update Dockerfile with sourmash compute bam input dependencies
* Add "fastp" to container requirements
* Update Dockerfile with sourmash compute bam input dependencies
* Update Dockerfile with [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo

## v1.0 - 6 March 2019

Initial release of nf-core/nf-kmer-similarity, created with the [nf-core](http://nf-co.re/) template.
