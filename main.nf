Channel
    .fromFilePairs('s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/SRR*_{1,2}.fastq.gz')
    .set{ samples_ch }


ksize = 15
log2_sketch_size = 10
molecule = 'protein'

//AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

process sourmash_compute_sketch {
	tag "$name"

	container 'czbiohub/kmer-hashing'

	input:
	set val(name), file(reads) from samples_ch

	output:
	file "${name}.sig" into sourmash_sketches

	script:
	"""
	sourmash compute \
		--num-hashes \$((2**$log2_sketch_size)) \
		--ksizes $ksize \
		--output ${name}.sig \
		--merge '$name' $reads
	"""
}

process sourmash_compare_sketches {
	container 'czbiohub/kmer-hashing'

	input:
	file sketches from sourmash_sketches.collect()

	output:
	file "similarities_ksize=${ksize}_molecule=${molecule}.csv"

	script:
	"""
	sourmash compare \
        --ksize $ksize \
        --$molecule \
        --csv similarities_ksize=${ksize}_molecule=${molecule}.csv \
        $sketches
	"""

}