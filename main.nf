
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

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv


    With read pairs in one or more semicolon-separated s3 directories:

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ \
        --read_pairs s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq/*{R1,R2}*.fastq.gz;s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/*{R1,R2}*.fastq.gz


    With plain ole fastas in one or more semicolon-separated s3 directories:

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/choanoflagellates_richter2018/ \
        --fastas /home/olga/data/figshare/choanoflagellates_richter2018/1_choanoflagellate_transcriptomes/*.fasta


    With SRA ids (requires nextflow v19.03-edge or greater):

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ --sra SRP016501


    Mandatory Arguments:
      --outdir                      Local or S3 directory to output the comparison matrix to

    Sample Arguments -- One or more of:
      --read_pairs                  Local or s3 directories containing *R{1,2}*.fastq.gz
                                    files, separated by commas
      --read_singles                Local or s3 directories of single-end read files, separated by commas
      --csv_pairs                   CSV file with columns id, read1, read2 for each sample
      --csv_singles                 CSV file with columns id, read1, read2 for each sample
      --fastas
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

// Samples from SRA
sra_ch = Channel.empty()

// R1, R2 pairs from a samples.csv file
samples_ch = Channel.empty()

// Single-enede reads from a samples.csv file
csv_singles_ch = Channel.empty()

// Extract R1, R2 pairs from a directory
read_pairs_ch = Channel.empty()

// Extract single-ended from a directory
read_singles_ch = Channel.empty()

// vanilla fastas
fastas_ch = Channel.empty()


 // Provided SRA ids
 if (params.sra){
   sra_ch = Channel
       .fromSRA( params.sra?.toString()?.tokenize(';') )
       .ifEmpty { exit 1, "params.sra ${params.sra} was not found - no input files supplied" }
 }
 // Provided a samples.csv file of read pairs
 if (params.csv_pairs){
   samples_ch = Channel
    .fromPath(params.csv_pairs)
    .splitCsv(header:true)
    .map{ row -> tuple(row[0], tuple(file(row[1]), file(row[2])))}
    .ifEmpty { exit 1, "params.csv_pairs was empty - no input files supplied" }
}

 // Provided a samples.csv file of single-ended reads
 if (params.csv_singles){
   csv_singles_ch = Channel
    .fromPath(params.csv_singles)
    .splitCsv(header:true)
    .map{ row -> tuple(row[0], tuple(file(row[1])))}
    .ifEmpty { exit 1, "params.csv_singles was empty - no input files supplied" }
}

 // Provided fastq gz read pairs
 if (params.read_pairs){
   read_pairs_ch = Channel
     .fromFilePairs(params.read_pairs?.toString()?.tokenize(';'))
     .ifEmpty { exit 1, "params.read_pairs was empty - no input files supplied" }
 }
 // Provided fastq gz read pairs
 if (params.read_singles){
   read_singles_ch = Channel
     .fromFilePairs(params.read_singles?.toString()?.tokenize(';'), size: 1)
     .ifEmpty { exit 1, "params.read_singles was empty - no input files supplied" }
}
 // Provided vanilla fastas
 if (params.fastas){
   fastas_ch = Channel
     .fromPath(params.fastas?.toString()?.tokenize(';'))
     .map{ f -> tuple(f.baseName, tuple(file(f))) }
     .ifEmpty { exit 1, "params.fastas was empty - no input files supplied" }
 }

// Parameters for testing
 if(params.read_paths_singles){
     read_paths_single_end_ch = Channel
         .from(params.read_paths_singles)
         .map { row -> [ row[0], [file(row[1][0])]] }
         .ifEmpty { exit 1, "params.read_paths_single_end was empty - no input files supplied" }
 }
if (params.read_paths) {
     read_paths_ch = Channel
         .from(params.read_paths)
         .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
         .ifEmpty { exit 1, "params.read_paths was empty - no input files supplied" }
 }


 sra_ch.concat(samples_ch, csv_singles_ch, read_pairs_ch, read_singles_ch, fastas_ch, read_paths_single_end_ch, read_paths_ch)
  .set{ reads_ch }


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

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


// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.read_pairs)     summary['Read Pairs']                 = params.read_pairs
if(params.read_singles)     summary['Single-end reads']         = params.read_singles
if(params.csv_pairs) summary['Paired-end samples.csv']            = params.csv_pairs
if(params.csv_singles) summary['Single-end samples.csv']    = params.csv_singles
if(params.sra)       summary['SRA']                             = params.sra
if(params.fasta)     summary["FASTAs"]                          = params.fasta
if(params.read_paths) summary['Read paths (paired-end)']            = params.read_paths
if(params.read_paths_singles) summary['Read paths (single-end)']    = params.read_paths_singles
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[0m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-kmer-similarity-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/kmer-similarity Workflow Summary'
    section_href: 'https://github.com/nf-core/kmer-similarity'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


process sourmash_compute_sketch {
	tag "${sample_id}_${sketch_id}"
	publishDir "${params.outdir}/sketches", mode: 'copy'
	container 'czbiohub/nf-kmer-similarity'

	// If job fails, try again with more memory
	// memory { 8.GB * task.attempt }
	errorStrategy 'retry'
  maxRetries 3

	input:
	each ksize from ksizes
	each molecule from molecules
	each log2_sketch_size from log2_sketch_sizes
	set sample_id, file(reads) from reads_ch

	output:
  set val(sketch_id), val(molecule), val(ksize), val(log2_sketch_size), file("${sample_id}_${sketch_id}.sig") into sourmash_sketches

	script:
  sketch_id = "molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"
  molecule = molecule
  not_dna = molecule == 'dna' ? '' : '--no-dna'
  ksize = ksize
  if ( params.one_signature_per_record ){
    """
    sourmash compute \\
      --num-hashes \$((2**$log2_sketch_size)) \\
      --ksizes $ksize \\
      --$molecule \\
      $not_dna \\
      --output ${sample_id}_${sketch_id}.sig \\
      $reads
    """
  } else {
    """
    sourmash compute \\
      --num-hashes \$((2**$log2_sketch_size)) \\
      --ksizes $ksize \\
      --$molecule \\
      $not_dna \\
      --output ${sample_id}_${sketch_id}.sig \\
      --merge '$sample_id' $reads
    """
  }

}

// sourmash_sketches.println()
// sourmash_sketches.groupTuple(by: [0,3]).println()

process sourmash_compare_sketches {
	tag "${sketch_id}"

	container 'czbiohub/nf-kmer-similarity'
	publishDir "${params.outdir}/", mode: 'copy'
	errorStrategy 'retry'
  maxRetries 3

	input:
  set val(sketch_id), val(molecule), val(ksize), val(log2_sketch_size), file ("sketches/*.sig") \
    from sourmash_sketches.groupTuple(by: [0, 3])

	output:
	file "similarities_${sketch_id}.csv"

	script:
	"""
	sourmash compare \\
        --ksize ${ksize[0]} \\
        --${molecule[0]} \\
        --csv similarities_${sketch_id}.csv \\
        --traverse-directory .
	"""

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/kmer-similarity v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
