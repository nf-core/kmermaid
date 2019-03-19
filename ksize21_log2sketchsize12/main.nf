
params.samples = "samples.csv"
params.log2_sketch_size = 12
params.ksize = 21
params.molecule = 'protein'
params.outdir = "s3://olgabot-maca/nf-kmer-similarity/human_mouse_zebrafish/"

sketch_id = "molecule-${params.molecule}_ksize-${params.ksize}_log2sketchsize-${params.log2_sketch_size}"


Channel
	.fromPath(params.samples)
	.splitCsv(header:true)
	.map{ row -> tuple(row.sample_id, file(row.read1), file(row.read2))}
	.set{ samples_ch }


//AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

process sourmash_compute_sketch {
	tag "${sample_id}_${sketch_id}"
	publishDir "${params.outdir}/sketches/${sketch_id}", mode: 'copy'
	container 'czbiohub/kmer-hashing'

	// If job fails, try again with more memory
	memory { 2.GB * task.attempt }
	errorStrategy 'retry'

	input:
	set sample_id, file(read1), file(read2) from samples_ch

	output:
	file "${sample_id}.sig" into sourmash_sketches

	script:
	"""
	sourmash compute \
		--num-hashes \$((2**$params.log2_sketch_size)) \
		--ksizes $params.ksize \
		--$params.molecule \
		--output ${sample_id}.sig \
		--merge '$sample_id' $read1 $read2
	"""
}


process sourmash_compare_sketches {
	tag "fromcsvs_${sketch_id}"

	container 'czbiohub/kmer-hashing'
	publishDir "${params.outdir}/", mode: 'copy'
	memory { 64.GB * task.attempt }
	errorStrategy 'retry'

	input:
	file ("sketches/${sketch_id}/*") from sourmash_sketches.collect()

	output:
	file "similarities_${sketch_id}.csv"

	script:
	"""
	sourmash compare \
        --ksize $params.ksize \
        --$params.molecule \
        --csv similarities_${sketch_id}.csv \
        --traverse-directory .
	"""

}
