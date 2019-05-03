
def helpMessage() {
    log.info """
    =========================================
     czbiohub/nf-kmer-similarity v${workflow.manifest.version}
    =========================================

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
      --samples                     CSV file with columns id, read1, read2 for each sample
      --fastas
      --read_pairs                 Local or s3 directories containing *R{1,2}*.fastq.gz
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
      --minlength                   Minimum length of reads after trimming, default 100

    Other Options:
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}



// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.fastas = false
params.sra = false
params.samples = false
params.read_pairs = false
params.one_signature_per_record = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.minlength = 100
params.outdir = "nextflow-output"

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * SET UP CONFIGURATION VARIABLES
 */

 // Samples from SRA
 sra_ch = Channel.empty()
 // R1, R2 pairs from a samples.csv file
 samples_ch = Channel.empty()
 // Extract R1, R2 pairs from a directory
 read_pairs_ch = Channel.empty()
 // vanilla fastas
 fastas_ch = Channel.empty()

 // Provided SRA ids
 if (params.sra){
   sra_ch = Channel
       .fromSRA( params.sra?.toString()?.tokenize(';') )
 }
 // Provided a samples.csv file
 if (params.samples){
   samples_ch = Channel
    .fromPath(params.samples)
    .splitCsv(header:true)
    .map{ row -> tuple(row.sample_id, tuple(file(row.read1), file(row.read2)))}
 }
 // Provided fastq gz read pairs
 if (params.read_pairs){
   read_pairs_ch = Channel
     .fromFilePairs(params.read_pairs?.toString()?.tokenize(';'))
     .println()
 }
 // Provided vanilla fastas
 if (params.fastas){
   fastas_ch = Channel
     .fromPath(params.fastas?.toString()?.tokenize(';'))
     .map{ f -> tuple(f.baseName, tuple(file(f))) }
 }

 sra_ch.concat(samples_ch, read_pairs_ch, fastas_ch)
  .into{ read_files_fastqc; read_files_trimming }


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
log.info """==================================================================
      _                            _       _ _         _ _
     | |_ ___ _____ ___ ___    ___|_|_____|_| |___ ___|_| |_ _ _
     | '_|___|     | -_|  _|  |_ -| |     | | | .'|  _| |  _| | |
     |_,_|   |_|_|_|___|_|    |___|_|_|_|_|_|_|__,|_| |_|_| |_  |
                                                            |___|

czbiohub/nf-kmer-similarity v${workflow.manifest.version}"
=================================================================="""
def summary = [:]
summary['Pipeline Name']       = 'czbiohub/nf-kmer-similarity'
summary['Pipeline Version']    = workflow.manifest.version
summary['Run Name']            = custom_runName ?: workflow.runName
summary['Read Pairs']          = params.read_pairs
summary['Input fastas']        = params.fastas
summary['SRA/ENA IDs']         = params.sra
summary['Molecules']           = params.molecules
summary['K-mer sizes']         = params.ksizes
summary['Log2 sketch sizes']   = params.log2_sketch_sizes
summary["Make one signature per record?"] = params.one_signature_per_record

summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-large-assembly-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'czbiohub/nf-kmer-simlarity Workflow Summary'
    section_href: 'https://github.com/czbiohub/nf-kmer-simlarity'
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
    container 'czbiohub/nf-kmer-similarity'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    conda activate multiqc && multiqc --version > v_multiqc.txt
    conda activate fastqc && fastqc --version > v_fastqc.txt
    conda deactivate
    sourmash info > v_sourmash.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    conda activate fastqc
    fastqc -q $reads
    """
}



/*
 * STEP 2 - trim reads - Fastp
 */
process fastp {
    tag "$name"
    publishDir "${params.outdir}/fastp", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_trimming

    output:
    file "*_fastp.{zip,html}" into fastp_results
    val name into sample_ids
    file "${name}_R1_fastp_trimmed.fastq.gz" into read1_trimmed
    file "${name}_R2_fastp_trimmed.fastq.gz" into read2_trimmed

    script:
    read1 = reads[0]
    read2 = reads[1]
    """
    conda activate fastp
    fastp --in1 $read1 --in2 $read2 \
      --length_required ${params.minlength} \
      --thread ${task.cpus} \
      --overrepresentation_analysis \
      --out1 ${name}_R1_fastp_trimmed.fastq.gz \
      --out2 ${name}_R2_fastp_trimmed.fastq.gz \
      -h ${name}_fastp.html \
      -j ${name}_fastp.json
    """
}



process sourmash_compute_sketch {
	tag "${sample_id}_${sketch_id}"
	publishDir "${params.outdir}/sketches", mode: 'copy'

	// If job fails, try again with more memory
	// memory { 8.GB * task.attempt }
	errorStrategy 'retry'
  maxRetries 3

	input:
	each ksize from ksizes
	each molecule from molecules
	set sample_id from sample_ids.collect()
  set read1 from read1_trimmed.collect()
  set read2 from read2_trimmed.collect()

	output:
  set val(sketch_id), val(molecule), val(ksize), val(log2_sketch_size), file("${sample_id}_${sketch_id}.sig") into sourmash_sketches

	script:
  sketch_id = "molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}"
  molecule = molecule
  ksize = ksize
  if ( params.one_signature_per_record ){
    """
    sourmash compute \
      --num-hashes \$((2**$log2_sketch_size)) \
      --ksizes $ksize \
      --$molecule \
      --output ${sample_id}_${sketch_id}.sig \
      $read1 $read2
    """
  } else {
    """
    sourmash compute \
      --num-hashes \$((2**$log2_sketch_size)) \
      --ksizes $ksize \
      --$molecule \
      --output ${sample_id}_${sketch_id}.sig \
      --merge '$sample_id' $reads
    """
  }

}

// sourmash_sketches.println()
// sourmash_sketches.groupTuple(by: [0,3]).println()

process sourmash_compare_sketches {
	tag "${sketch_id}"

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
	sourmash compare \
        --ksize ${ksize[0]} \
        --${molecule[0]} \
        --csv similarities_${sketch_id}.csv \
        --traverse-directory .
	"""

}


/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    conda activate r
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/test] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/test] FAILED: $workflow.runName"
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
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/test] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/test] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/test] Pipeline Complete"

}
