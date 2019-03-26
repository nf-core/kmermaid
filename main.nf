
def helpMessage() {
    log.info """
    ==============================================================
      _                            _       _ _         _ _
     | |_ ___ _____ ___ ___    ___|_|_____|_| |___ ___|_| |_ _ _
     | '_|___|     | -_|  _|  |_ -| |     | | | .'|  _| |  _| | |
     |_,_|   |_|_|_|___|_|    |___|_|_|_|_|_|_|__,|_| |_|_| |_  |
                                                            |___|
    ==============================================================

    Usage:

    The typical command for running the pipeline is as follows.

    With a samples.csv file:

      nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv

    With multitple s3 directories:

      nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --directories s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq,s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/

    With SRA ids:

      nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --sra SRP016501


    Mandatory Arguments:
      --outdir                      Local or S3 directory to output the comparison matrix to

    Sample Arguments -- One or more of:
      --samples                     CSV file with columns id, read1, read2 for each sample
      --directories                 Local or s3 directories containing *R{1,2}*.fastq.gz
                                    files, separated by commas
      --sra                         SRR, ERR, SRP IDs representing a project. Only compatible with
                                    Nextflow 19.03-edge or greater


    Options:
      --ksizes                      Which nucleotide k-mer sizes to use. Multiple are
                                    separated by commas. Default is '21,27,33,51'
      --molecules                   Which molecule to compare on. Default is both DNA
                                    and protein, i.e. 'dna,protein'
      --log2-sketch-sizes           Which log2 sketch sizes to use. Multiple are separated
                                    by commas. Default is '10,12,14,16'
    """.stripIndent()
}



// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Get all possible input samples
reads_ch = Channel.create()

// Provided SRA ids
if (params.sra){
  Channel
      .fromSRA( params.sra )
      .set{ sra_ch }
  reads_ch = reads_ch.concat(sra_ch)
}
// Provided a samples.csv file
if (params.samples){
  Channel
  	.fromPath(params.samples)
  	.splitCsv(header:true)
  	.map{ row -> tuple(row.sample_id, file(row.read1), file(row.read2))}
  	.set{ samples_ch }
  reads_ch = reads_ch.concat(samples_ch)
}
// Provided s3 or local directories
if (params.directories){
  Channel
    .fromFilePairs(params.directories.splitCsv())
    .set{ directories_ch }
  reads_ch = reads_ch.concat(directories_ch)
}


// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}


params.ksizes = [21, 27, 33, 51]
params.molecules =  ['dna', 'protein']
params.log2_sketch_sizes = [10, 12, 14, 16]


process sourmash_compute_sketch {
	tag "${sample_id}_molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"
	publishDir "${params.outdir}/sketches/molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}", mode: 'copy'
	container 'czbiohub/kmer-hashing'

	// If job fails, try again with more memory
	memory { 2.GB * task.attempt }
	errorStrategy 'retry'

	input:
	each ksize from params.ksizes
	each molecule from params.molecules
	each log2_sketch_size from params.log2_sketch_sizes
	set sample_id, file(read1), file(read2) from reads_ch

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
	each ksize from ksizes
	each molecule from molecules
	each log2_sketch_size from log2_sketch_sizes
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
