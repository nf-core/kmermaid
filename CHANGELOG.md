# nf-core/nf-kmer-similarity: Changelog

## v1.1dev

* Add option to use Dayhoff encoding for sourmash
* Add bam2fasta process to kmermaid pipeline and optional params
  for the same
* Update dockerfile with bam2fasta dependencies
* Add "fastp" to container requirements

### Dependency updates
* Add fastqc to environment.yml
* Add [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo to environment.yml
* Update Dockerfile with sourmash compute bam input dependencies
* Add "fastp" to container requirements
* Update Dockerfile with sourmash compute bam input dependencies
>>>>>>> 4a04a895b5a68abe28303eb45987e9c955decee0

## v1.0 - 6 March 2019

Initial release of nf-core/nf-kmer-similarity, created with the [nf-core](http://nf-co.re/) template.
