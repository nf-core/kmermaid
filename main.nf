
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

    With a samples.csv file containing the columns sample_id,read1,read2:

      nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv

    With one or more s3 directories:

      nextflow run czbiohub/nf-kmer-similarity --outdir s3://olgabot-maca/nf-kmer-similarity/ --directories s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq,s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/

    With SRA ids (requires nextflow v19.03-edge or greater):

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
      --log2_sketch_sizes           Which log2 sketch sizes to use. Multiple are separated
                                    by commas. Default is '10,12,14,16'
      --one_signature_per_record    Make a k-mer signature for each record in the FASTQ/FASTA files.
                                    Useful for comparing e.g. assembled transcriptomes or metagenomes.
                                    (Not typically used for raw sequencing data as this would create
                                    a k-mer signature for each read!)
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

//  // Samples from SRA
//  sra_ch = Channel.create()
//  // R1, R2 pairs from a samples.csv file
//  samples_ch = Channel.create()
//  // Extract R1, R2 pairs from a directory
//  directories_ch = Channel.create()
//
//  // Provided SRA ids
//  if (params.sra){
//    sra_ch = Channel
//        .fromSRA( params.sra?.toString()?.tokenize(',') )
//  }
//  // Provided a samples.csv file
//  if (params.samples){
   // samples_ch = Channel
   // 	.fromPath(params.samples)
   // 	.splitCsv(header:true)
   // 	.map{ row -> tuple(row.sample_id, tuple(row.read1, row.read2))}
//  }
//  // Provided s3 or local directories
//  if (params.directories){
//    directories_ch = Channel
//      .from(params.directories?.toString()?.tokenize(','))
//      .map(it + "*{R1,R2}*.fastq.gz")
//      .fromFilePairs()
//  }
//
// sra_ch.concat(samples_ch, directories_ch)
//   .set{ reads_ch }

Channel
  .fromPath(params.samples)
  .splitCsv(header:true)
  .map{ row -> tuple(row.sample_id, tuple(row.read1, row.read2))}
  .set{ reads_ch }


// Channel
//   .fromSRA( params.sra?.toString()?.tokenize(',') )
//   .set{ reads_ch }

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}


params.ksizes = '21,27,33,51'
params.molecules =  'dna,protein'
params.log2_sketch_sizes = '10,12,14,16'

// Parse the parameters
ksizes = params.ksizes?.toString().tokenize(',')
molecules = params.molecules?.toString().tokenize(',')
log2_sketch_sizes = params.log2_sketch_sizes?.toString().tokenize(',')


process sourmash_compute_sketch {
	tag "${sample_id}_molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"
	publishDir "${params.outdir}/sketches/molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}", mode: 'copy'
	container 'czbiohub/nf-kmer-similarity'

	// If job fails, try again with more memory
	memory { 2.GB * task.attempt }
	errorStrategy 'retry'
  maxRetries 5

	input:
	each ksize from ksizes
	each molecule from molecules
	each log2_sketch_size from log2_sketch_sizes
	set sample_id, file(reads) from reads_ch

	output:
	file "${sample_id}.sig" into sourmash_sketches

	script:
  if ( params.one_signature_per_record ){
    """
    sourmash compute \
      --num-hashes \$((2**$log2_sketch_size)) \
      --ksizes $ksize \
      --$molecule \
      --output ${sample_id}.sig \
      $read1 $read2
    """
  } else {
    """
    sourmash compute \
      --num-hashes \$((2**$log2_sketch_size)) \
      --ksizes $ksize \
      --$molecule \
      --output ${sample_id}.sig \
      --merge '$sample_id' $reads
    """
  }

}


process sourmash_compare_sketches {
	tag "molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"

	container 'czbiohub/nf-kmer-similarity'
	publishDir "${params.outdir}/", mode: 'copy'
	memory { 8.GB * task.attempt }
	errorStrategy 'retry'
  maxRetries 5

	input:
	each ksize from ksizes
	each molecule from molecules
	each log2_sketch_size from params.log2_sketch_sizes
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
