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
	sudo nextflow run main.nf -work-dir ${HOME}/pure-scratch/nextflow/ \
		-process.executor='local'


test_sra:
		nextflow run main.nf --sra "SRP016501" -profile aws

test_samplescsv:
		nextflow run main.nf --ksizes 21 --log2_sketch_sizes 10 \
			--molecules dna \
			--samples testing/samples.csv -profile local

test_sra_samplescsv:
	nextflow run main.nf \
		-profile aws \
		--sra SRP016501 \
		--samples testing/samples.csv


test: test_sra test_samplescsv test_sra_samplescsv
