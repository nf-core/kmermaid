# nf-core/kmermaid: Changelog

## v1.1dev

* Add option to use Dayhoff encoding for sourmash
* Add bam file process to kmermaid pipeline and optional params
  for the same

### Dependency updates

* Add fastqc to environment.yml
* Add [czbiohub/khtools](https://github.com/czbiohub/kh-tools/) repo to environment.yml
* Update Dockerfile with sourmash compute bam input dependencies
* Add [`bam2fasta`](https://pypi.org/project/bam2fasta/) to environment.yml
* Use Olga's branch of sourmash
* Add `ska` and `seqtk` to container dependencies
* Add `fastp` to container requirements
* Add `track_abundance` feature to to keep track of hashed kmer frequency.

## v1.0dev - 6 March 2019
Initial release of nf-core/nf-kmer-similarity, created with the [nf-core](http://nf-co.re/) template.
