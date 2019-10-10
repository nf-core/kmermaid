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


    With BAM file:

      nextflow run main.nf \
      --outdir ./results \
      --bam possorted_genome_bam.bam

    Mandatory Arguments:
      --outdir                      Local or S3 directory to output the comparison matrix to

    Sample Arguments -- One or more of:
      --read_pairs                  Local or s3 directories containing *R{1,2}*.fastq.gz
                                    files, separated by commas
      --read_singles                Local or s3 directories of single-end read files, separated by commas
      --csv_pairs                   CSV file with columns id, read1, read2 for each sample
      --csv_singles                 CSV file with columns id, read1, read2 for each sample
      --fastas                      Path to FASTA sequence files. Can be semi-colon-separated
      --bam                         Path to 10x BAM file
      --bai                         Path to 10x BAI file
      --save_fastas                 For bam files, Path relative to outdir to save unique barcodes to {CELL_BARCODE}.fasta
      --sra                         SRR, ERR, SRP IDs representing a project. Only compatible with
                                    Nextflow 19.03-edge or greater


    Options:
      --ksizes                      Which nucleotide k-mer sizes to use. Multiple are
                                    separated by commas. Default is '21,27,33,51'
      --molecules                   Which molecule to compare on. Default is both DNA
                                    and protein, i.e. 'dna,protein,dayhoff'
      --log2_sketch_sizes           Which log2 sketch sizes to use. Multiple are separated
                                    by commas. Default is '10,12,14,16'
      --one_signature_per_record    Make a k-mer signature for each record in the FASTQ/FASTA files.
                                    Useful for comparing e.g. assembled transcriptomes or metagenomes.
                                    (Not typically used for raw sequencing data as this would create
                                    a k-mer signature for each read!)
      --write_barcode_meta_csv      For bam files, Csv file name relative to outdir/barcode_metadata to write number of reads and number of umis per barcode.
                                    This csv file is empty with just header when the min_umi_per_barcode is zero i.e
                                    Reads and umis per barcode are calculated only when the barcodes are filtered
                                    based on min_umi_per_barcode
      --min_umi_per_barcode           A barcode is only considered a valid barcode read
                                    and its signature is written if number of umis are greater than min_umi_per_barcode
      --barcodes_file               For bam files, Optional absolute path to a .tsv barcodes file if the input is unfiltered 10x bam file
      --rename_10x_barcodes         For bam files, Optional absolute path to a .tsv Tab-separated file mapping 10x barcode name
                                    to new name, e.g. with channel or cell annotation label
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

// read_paths is only used for testing
read_paths_ch = Channel.empty()

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

// Parameters for testing
if (params.read_paths) {
     read_paths_ch = Channel
         .from(params.read_paths)
         .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
         .ifEmpty { exit 1, "params.read_paths (${params.read_paths}) was empty - no input files supplied" }
 } else {
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
      .ifEmpty { exit 1, "params.csv_pairs (${params.csv_pairs}) was empty - no input files supplied" }
  }

   // Provided a samples.csv file of single-ended reads
   if (params.csv_singles){
     csv_singles_ch = Channel
      .fromPath(params.csv_singles)
      .splitCsv(header:true)
      .map{ row -> tuple(row[0], tuple(file(row[1])))}
      .ifEmpty { exit 1, "params.csv_singles (${params.csv_singles}) was empty - no input files supplied" }
  }

   // Provided fastq gz read pairs
   if (params.read_pairs){
     read_pairs_ch = Channel
       .fromFilePairs(params.read_pairs?.toString()?.tokenize(';'))
       .ifEmpty { exit 1, "params.read_pairs (${params.read_pairs}) was empty - no input files supplied" }
   }
   // Provided fastq gz read singles
   if (params.read_singles){
     read_singles_ch = Channel
       .fromFilePairs(params.read_singles?.toString()?.tokenize(';'), size: 1)
       .ifEmpty { exit 1, "params.read_singles (${params.read_singles}) was empty - no input files supplied" }
  }
   // Provided vanilla fastas
   if (params.fastas){
     fastas_ch = Channel
       .fromPath(params.fastas?.toString()?.tokenize(';'))
       .map{ f -> tuple(f.baseName, tuple(file(f))) }
       .ifEmpty { exit 1, "params.fastas (${params.fastas}) was empty - no input files supplied" }
   }

  if (params.bam) {
  Channel.fromPath(params.bam, checkIfExists: true)
        .ifEmpty { exit 1, "Bam file not found: ${params.bam}" }
        .map{ f -> tuple(f.baseName, tuple(file(f))) }
        .set{bam_ch}
  Channel.fromPath(params.bai, checkIfExists: true)
        .ifEmpty { exit 1, "Bai file not found in ${params.bai}" }
        .set{bai_ch}
  }

  // If barcodes is as expected, check if it exists and set channel
  if (params.barcodes_file) {
     Channel.fromPath(params.barcodes_file, checkIfExists: true)
        .ifEmpty { exit 1, "Barcodes file not found: ${params.barcodes_file}" }
        .set{barcodes_ch}
  }
  else {
    Channel.from(false)
        .set{barcodes_ch}
  }

  // If renamer barcode file is as expected, check if it exists and set channel
  if (params.rename_10x_barcodes) {
     Channel.fromPath(params.rename_10x_barcodes, checkIfExists: true)
        .ifEmpty { exit 1, "Barcodes file not found: ${params.rename_10x_barcodes}" }
        .set{rename_10x_barcodes_ch}
  }
  else {
    Channel.from(false)
        .set{rename_10x_barcodes_ch}
  }
}

if (!params.bam) { 
sra_ch.concat(samples_ch, csv_singles_ch, read_pairs_ch,
 read_singles_ch, fastas_ch, read_paths_ch)
 .ifEmpty{ exit 1, "No reads provided! Check read input files"}
 .set{ reads_ch }
}


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


// Parse the parameters
ksizes = params.ksizes?.toString().tokenize(',')
ksize = params.ksizes[0]
molecules = params.molecules?.toString().tokenize(',')
molecule = molecules[0]
log2_sketch_sizes = params.log2_sketch_sizes?.toString().tokenize(',')
log2_sketch_size = params.log2_sketch_sizes[0]


// For bam files, set a folder name to save the optional barcode metadata csv
if (!params.write_barcode_meta_csv) {
  barcode_metadata_folder = ""
}
else {
  barcode_metadata_folder = "barcode_metadata"
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// Input reads
if(params.read_pairs)   summary['Read Pairs']                 = params.read_pairs
if(params.read_singles) summary['Single-end reads']         = params.read_singles
if(params.csv_pairs)    summary['Paired-end samples.csv']            = params.csv_pairs
if(params.csv_singles)  summary['Single-end samples.csv']    = params.csv_singles
if(params.sra)          summary['SRA']                             = params.sra
if(params.fastas)       summary["FASTAs"]                          = params.fastas
if(params.bam)          summary["BAM"]                             = params.bam
if(params.bam)          summary["BAI"]                             = params.bai
if(params.barcodes_file)          summary["Barcodes"]              = params.barcodes_file
if(params.rename_10x_barcodes)    summary["Renamer barcodes"]      = params.rename_10x_barcodes
if(params.read_paths)   summary['Read paths (paired-end)']         = params.read_paths
// Sketch parameters
summary['K-mer sizes']            = params.ksizes
summary['Molecule']               = params.molecules
summary['Log2 Sketch Sizes']      = params.log2_sketch_sizes
summary['One Sig per Record']         = params.one_signature_per_record
// 10x parameters
if(params.min_umi_per_barcode) summary['Count valid reads'] = params.min_umi_per_barcode
if(params.save_fastas) summary['Saved Fastas '] = params.save_fastas
if(params.write_barcode_meta_csv) summary['Barcode umi read metadata'] = params.write_barcode_meta_csv
// Resource information
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


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.txt"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    sourmash info &> v_sourmash.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

if (params.bam) {
  process sourmash_compute_sketch_bam {
    tag "${sample_id}_${sketch_id}"
    publishDir "${params.outdir}/sketches", pattern: '*.sig', mode: 'copy'
    publishDir "${params.outdir}/${params.save_fastas}", pattern: '*.fasta', saveAs: { filename -> "${params.outdir}/${params.save_fastas}/${filename.replace("|", "-")}"}
    publishDir "${params.outdir}/${barcode_metadata_folder}", pattern: '*.csv', mode: 'copy'


    // If job fails, try again with more memory
    // memory { 8.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

    input:
    ksize
    molecule
    log2_sketch_size
    set sample_id, file(bam) from bam_ch
    file(barcodes_file) from barcodes_ch
    file(bai) from bai_ch
    file(rename_10x_barcodes) from rename_10x_barcodes_ch

    output:
    set val(sample_id), file("*.fasta") into reads_ch
    // https://github.com/nextflow-io/patterns/blob/master/docs/optional-output.adoc
    file("${params.write_barcode_meta_csv}") optional true

    script:
    sketch_id = "molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"
    molecule = molecule
    not_dna = molecule != 'dna' ? '--no-dna' : ''
    ksize = ksize

    min_umi_per_barcode = params.min_umi_per_barcode ? "--count-valid-reads ${params.min_umi_per_barcode}" : ''
    metadata = params.write_barcode_meta_csv ? "--write-barcode-meta-csv ${params.write_barcode_meta_csv}": ''
    save_fastas = "--save-fastas ${params.save_fastas}"

    def barcodes_file = params.barcodes_file ? "--barcodes-file ${barcodes_file.baseName}.tsv": ''
    def rename_10x_barcodes = params.rename_10x_barcodes ? "--rename-10x-barcodes ${rename_10x_barcodes.baseName}.tsv": ''
    """
      sourmash compute \\
        --ksize $ksize \\
        --$molecule \\
        $not_dna \\
        --num-hashes \$((2**$log2_sketch_size)) \\
        --processes ${task.cpus} \\
        $min_umi_per_barcode \\
        $rename_10x_barcodes \\
        $barcodes_file \\
        $save_fastas \\
        $metadata \\
        --output ${sample_id}_${sketch_id}.sig \\
        --input-is-10x $bam
      find . -type f -name "*.fasta" | while read src; do if [[ \$src == *"|"* ]]; then mv "\$src" \$(echo "\$src" | tr "|" "_"); fi done
    """
  }
}

Channel.from(reads_ch)
process sourmash_compute_sketch_fastx {
  tag "${sample_id}_${sketch_id}"
  publishDir "${params.outdir}/sketches", mode: 'copy'

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

  if ( params.one_signature_per_record ) {
    """
    sourmash compute \\
      --num-hashes \$((2**$log2_sketch_size)) \\
      --ksizes $ksize \\
      --$molecule \\
      $not_dna \\
      --output ${sample_id}_${sketch_id}.sig \\
      $reads
    """
  }
  else {
    """
    sourmash compute \\
      --num-hashes \$((2**$log2_sketch_size)) \\
      --ksizes $ksize \\
      --$molecule \\
      $not_dna \\
      --output ${sample_id}_${sketch_id}.sig \\
      --merge '$sample_id' \\
      $reads
    """
  }
}


process sourmash_compare_sketches {
  tag "${sketch_id}"

  container "$workflow.container"
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
        --processes ${task.cpus} \\
        --csv similarities_${sketch_id}.csv \\
        --traverse-directory .
  """

}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/kmermaid] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/kmermaid] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/kmermaid] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/kmermaid] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = file( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = file( "${output_d}/pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = file( "${output_d}/pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/kmermaid]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/kmermaid]${c_red} Pipeline completed with errors${c_reset}"
    }

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
