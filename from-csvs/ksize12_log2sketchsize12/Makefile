run:
	nextflow run ../main.nf --samples ../samples.csv --log2_sketch_size 12 --ksize 12 --molecule 'protein' --outdir s3://olgabot-maca/nf-kmer-similarity/human_mouse_zebrafish/ \
		-work-dir s3://olgabot-maca/nextflow-workdir-test/ \
              -bucket-dir s3://olgabot-maca/nextflow-bucket-dir-test/ \
              -with-trace -with-timeline -with-dag -with-report -latest -resume
