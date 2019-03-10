run_aws:
	nextflow run main.nf \
		-work-dir s3://olgabot-maca/nextflow-workdir-test/ \
		-bucket-dir s3://olgabot-maca/nextflow-bucket-dir-test/ \
		-with-trace -with-timeline -with-dag -with-report -latest -resume

run:
	nextflow run main.nf

run_ndnd:
	sudo nextflow run main.nf -work-dir ${HOME}/pure-scratch/nextflow/

run_aws_csv:
	pushd from-csvs
	nextflow run main.nf \
		-work-dir s3://olgabot-maca/nextflow-workdir-test/ \
		-bucket-dir s3://olgabot-maca/nextflow-bucket-dir-test/ \
		-with-trace -with-timeline -with-dag -with-report -latest -resume
	popd

run_ndnd_local:
	sudo nextflow run main.nf -work-dir ${HOME}/pure-scratch/nextflow/  -process.executor='local'
