
println "asdlfkijasdlkf"

def paths = file('directories.txt').readLines().findAll { it.size()>0 }

println paths


// // --- Try both glob patterns ---
// // This is the Illumina style R1/R2 specification
// Channel
// 	.from(paths)
// 	.map{ it + "**_{R1,R2}_*.fastq.gz"}
// 	.set{ illumina_style }

// // This is the SRA style R1/R2 specifiation
// Channel
// 	.from(paths)
// 	.map{ it + "**_{1,2}.fastq.gz"}
// 	.set{ sra_style }


// paths_concatenated = illumina_style.concat(sra_style)


// Channel
// 	.from(paths_concatenated)
// 	.map{ fromFilePairs(it) }
// 	.println()
// 	// .subscribe onNext: { println it }, onComplete: { println 'Done' }

// 	// .set{ samples_ch }


// // Channel
// // 	.fromFilePairs(  "${paths}**_{1,2}.fastq.gz" )
// //     .set{ samples_ch }

// // Channel
// // 	.fromFilePairs(  "${paths}**_{R1,R2}_*.fastq.gz" )
// // 	.combine(sra_style)
// //     .set{ samples_ch }


// println "samples_ch" samples_ch

Channel
	.fromPath("s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq/")
	.println()


Channel
	.fromFilePairs("s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq/**_{1,2}.fastq.gz")
	.set{ human_samples }
println human_samples

Channel
	.fromFilePairs("s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/**_{1,2}.fastq.gz")
	.set{ zebrafish_samples }

println zebrafish_samples

samples_ch = human_samples.concat(zebrafish_samples)

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
		--$molecule \
		--output ${name}.sig \
		--merge '$name' $reads
	"""
}

process sourmash_compare_sketches {
	container 'czbiohub/kmer-hashing'

	memory '32 GB'

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
