
	params.samples = "samples.csv"
// params.log2_sketch_size = 12
// params.ksize = 15
// params.molecule = 'protein'
params.outdir = "s3://olgabot-maca/nf-kmer-similarity/human_mouse_zebrafish/"

// sketch_id = "molecule-${params.molecule}_ksize-${params.ksize}_log2sketchsize-${params.log2_sketch_size}"

// if (params.molecule == "protein") {
// 	other_molecule = "dna"
// } else {
// 	other_molecule = "protein"
// }


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


ksizes = Channel.from([15, 21, 27, 33, 51])
molecules = Channel.from(['protein', 'dna'])
log2_sketch_sizes = Channel.from([10, 12, 14, 16])

parameters = molecules
	.combine(ksizes)
	.combine(log2_sketch_sizes)


process sourmash_compute_sketch {
	tag "${sample_id}_molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"
	publishDir "${params.outdir}/sketches/molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}", mode: 'copy'
	container 'czbiohub/kmer-hashing'

	// If job fails, try again with more memory
	memory { 2.GB * task.attempt }
	errorStrategy 'retry'

	input:
	each ksizes from ksizes
	each molecule from molecules
	each log2_sketch_sizes from log2_sketch_sizes
	set sample_id, file(read1), file(read2) from samples_ch

	output:
	file "${sample_id}.sig" into sourmash_sketches

	script:
	"""
	sourmash compute \
		--num-hashes \$((2**$log2_sketch_size)) \
		--ksizes $ksize \
		--$molecule \
		--output ${sample_id}.sig \
		--merge '$sample_id' $read1 $read2
	"""
}


process sourmash_compare_sketches {
	tag "fromcsvs_molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"

	container 'czbiohub/kmer-hashing'
	publishDir "${params.outdir}/", mode: 'copy'
	memory { 1024.GB * task.attempt }
	// memory { sourmash_sketches.size() < 100 ? 8.GB :
	// 	sourmash_sketches.size() * 100.MB * task.attempt}
	errorStrategy 'retry'

	input:
	each ksizes from ksizes
	each molecule from molecules
	each log2_sketch_sizes from log2_sketch_sizes
	file ("sketches/molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}/*") from sourmash_sketches.collect()

	output:
	file "similarities_molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}.csv"

	script:
	"""
	sourmash compare \
        --ksize $ksize \
        --$molecule \
        --csv similarities_molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}.csv \
        --traverse-directory .
	"""

}
