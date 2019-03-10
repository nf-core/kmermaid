
params.samples = "samples.csv"
params.log2_sketch_size = 10
params.ksize = 15
params.molecule = 'protein'
params.outdir = "s3://olgabot-maca/nf-kmer-similarity/human_mouse_zebrafish/"

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
	tag "${sample_id}_molecule-${params.molecule}_ksize-${params.ksize}_log2sketchsize-${params.log2_sketch_size}"
	publishDir "${params.outdir}/sketches/${sketch_id}", mode: 'copy'
	container 'czbiohub/kmer-hashing'
	memory '2 GB'

	input:
	set val(name), file(reads) from samples_ch

	output:
	file "${name}.sig" into sourmash_sketches

	script:
	"""
	sourmash compute \
		--num-hashes \$((2**$log2_sketch_size)) \
		--ksizes $ksize \
		--$molecule \
		--output ${name}.sig \
		--merge '$name' $reads
	"""
}


process sourmash_compare_sketches {
	container 'czbiohub/kmer-hashing'
	publishDir "${params.outdir}/", mode: 'copy'
	memory '32 GB'

	input:
	file ("sketches/${sketch_id}/*") from sourmash_sketches.collect()

	output:
	file "similarities_${sketch_id}.csv"

	script:
	"""
	sourmash compare \
        --ksize $ksize \
        --$molecule \
        --csv similarities_ksize=${ksize}_molecule=${molecule}.csv \
        --traverse-directory .
	"""

}
