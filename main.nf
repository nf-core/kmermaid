#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/kmermaid
========================================================================================
 nf-core/kmermaid Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/kmermaid
----------------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info nfcoreHeader()
    log.info """
    =========================================
     nf-core/kmermaid v${workflow.manifest.version}
    =========================================

    Usage:

    The typical command for running the pipeline is as follows.

    With a samples.csv file containing the columns sample_id,read1,read2:

      nextflow run nf-core/kmermaid \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv


    With read pairs in one or more semicolon-separated s3 directories:

      nextflow run nf-core/kmermaid \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ \
        --read_pairs s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq/*{R1,R2}*.fastq.gz;s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/*{R1,R2}*.fastq.gz


    With plain ole fastas in one or more semicolon-separated s3 directories:

      nextflow run nf-core/kmermaid \
        --outdir s3://olgabot-maca/nf-kmer-similarity/choanoflagellates_richter2018/ \
        --fastas /home/olga/data/figshare/choanoflagellates_richter2018/1_choanoflagellate_transcriptomes/*.fasta


    With SRA ids (requires nextflow v19.03-edge or greater):

      nextflow run nf-core/kmermaid \
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
      --save_fastas                 For bam files, Path relative to outdir to save unique barcodes to {CELL_BARCODE}.fasta
      --save_intermediate_files     save temporary fastas and chunks of bam files
                                    in the absolute path given by this flag
                                    By default, they are saved in temp directory.
                                    An important note is This might cause
                                    not enough space on the device left depending on the size of your bam file and harddisk space allocated for tmp folder on your machine, so its better to specify a directory.
                                    These files are deleted automatically at the end of the program.
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
      --track_abundance             Track abundance of each hashed k-mer, could be useful for cancer RNA-seq or ATAC-seq analyses
      --skip_trimming               If provided, skip fastp trimming of reads

    Split K-mer options:
      --splitKmer                   If provided, use SKA to compute split k-mer sketches instead of
                                    sourmash to compute k-mer sketches
      --subsample                   Integer value to subsample reads from input fastq files

    Bam file options:
      --write_barcode_meta_csv      For bam files, Csv file name relative to outdir/barcode_metadata to write number of reads and number of umis per barcode.
                                    This csv file is empty with just header when the tenx_min_umi_per_cell is zero i.e
                                    Reads and umis per barcode are calculated only when the barcodes are filtered
                                    based on tenx_min_umi_per_cell
      --tenx_min_umi_per_cell         A barcode is only considered a valid barcode read
                                    and its signature is written if number of umis are greater than tenx_min_umi_per_cell
      --shard_size                  Number of alignment to contain in each sharded bam file
      --barcodes_file               For bam files, Optional absolute path to a .tsv barcodes file if the input is unfiltered 10x bam file
      --rename_10x_barcodes         For bam files, Optional absolute path to a .tsv Tab-separated file mapping 10x barcode name
                                    to new name, e.g. with channel or cell annotation label

    Extract protein-coding options:
      --peptide_fasta               Path to a well-curated fasta file of protein sequences. Used to filter for coding reads
      --bloomfilter_tablesize       Maximum table size for bloom filter creation

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

// read_paths is only used for testing
read_paths_ch = Channel.empty()

// Samples from SRA
sra_ch = Channel.empty()

// R1, R2 pairs from a samples.csv file
csv_pairs_ch = Channel.empty()

// Single-enede reads from a samples.csv file
csv_singles_ch = Channel.empty()

// Extract R1, R2 pairs from a directory
read_pairs_ch = Channel.empty()

// Extract single-ended from a directory
read_singles_ch = Channel.empty()

// vanilla fastas
fastas_ch = Channel.empty()

// 10X Genomics .tgz file containing possorted_genome_bam file
tenx_tgz_ch = Channel.empty()


// Parameters for testing
if (params.read_paths) {
     read_paths_ch = Channel
        .from(params.read_paths)
        .map { row -> if (row[1].size() == 2) [ row[0], [file(row[1][0]), file(row[1][1])]]
              else [row[0], [file(row[1][0])]]}
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
     csv_pairs_ch = Channel
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

   // Provided fastq gz paired-end reads
   if (params.read_pairs){
     read_pairs_ch = Channel
       .fromFilePairs(params.read_pairs?.toString()?.tokenize(';'), size: 2)
       .ifEmpty { exit 1, "params.read_pairs (${params.read_pairs}) was empty - no input files supplied" }
   }
   // Provided fastq gz single-end reads
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
       .dump( tag: 'bam' )
       .set{ tenx_bam_for_unaligned_fastq_ch}
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

  if (params.tenx_tgz) {
    Channel.fromPath(params.tenx_tgz, checkIfExists: true)
       .dump(tag: 'tenx_tgz_before_mri_filter')
       .filter{ ~/.+[^mri]\.tgz/ }
       .ifEmpty { exit 1, "10X .tgz file not found: ${params.tenx_tgz}" }
       .dump(tag: 'tenx_tgz_after_mri_filter')
       .set{ tenx_tgz_ch }
  }
}

if (params.peptide_fasta) {
Channel.fromPath(params.peptide_fasta, checkIfExists: true)
     .ifEmpty { exit 1, "Peptide fasta file not found: ${params.peptide_fasta}" }
     .set{ ch_peptide_fasta }
}

if (params.subsample) {
  if (params.bam){
     exit 1, "Cannot provide both a bam file with --bam and specify --subsample"
  } else {
    if (params.skip_trimming){
      sra_ch.concat(csv_pairs_ch, csv_singles_ch, read_pairs_ch,
        read_singles_ch, fastas_ch, read_paths_ch)
        .ifEmpty{ exit 1, "No reads provided! Check read input files"}
        .set{ subsample_reads_ch }
    } else {
      sra_ch.concat(
          csv_pairs_ch, csv_singles_ch, read_pairs_ch,
          read_singles_ch, read_paths_ch)
        .ifEmpty{ exit 1, "No reads provided! Check read input files"}
        .set{ ch_read_files_trimming }
    }
  }
} else {
  if (!(params.tenx_tgz || params.bam)) {
    if(params.skip_trimming){
      sra_ch.concat(
          csv_pairs_ch, csv_singles_ch, read_pairs_ch,
          read_singles_ch, fastas_ch, read_paths_ch)
       .ifEmpty{ exit 1, "No reads provided! Check read input files"}
       .set{ reads_ch }
    } else {
      if (params.fastas) {
        // With fasta files - combine everything that can be trimmed
        sra_ch.concat(
            csv_pairs_ch, csv_singles_ch, read_pairs_ch,
            read_singles_ch, read_paths_ch)
          .set{ ch_read_files_trimming }
      } else {
        // No fasta files - combine everything and error out
        sra_ch.concat(
            csv_pairs_ch, csv_singles_ch, read_pairs_ch,
            read_singles_ch, read_paths_ch)
          .ifEmpty{ exit 1, "No reads provided! Check read input files"}
          .set{ ch_read_files_trimming }
      }
    }
  } else {
//   Do nothing - can't combine the fastq files and bam files (yet)
      sra_ch.concat(
          csv_pairs_ch, csv_singles_ch, read_pairs_ch,
          read_singles_ch, read_paths_ch)
        .set{ ch_non_bam_reads }
    }
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

if (workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}


if (params.splitKmer){
    params.ksizes = '15,9'
    params.molecules = 'dna'
} else {
    params.ksizes = '21,27,33,51'
}


// Parse the parameters

ksizes = params.ksizes?.toString().tokenize(',')
molecules = params.molecules?.toString().tokenize(',')
peptide_molecules = molecules.findAll { it != "dna" }
log2_sketch_sizes = params.log2_sketch_sizes?.toString().tokenize(',')

int bloomfilter_tablesize = Math.round(Float.valueOf(params.bloomfilter_tablesize))

peptide_ksize = params.translate_peptide_ksize
peptide_molecule = params.translate_peptide_molecule
jaccard_threshold = params.translate_jaccard_threshold
track_abundance = params.track_abundance

if (params.bam){
  // Extract the fasta just once using sourmash
  single_ksize = ksizes[0]
  single_molecule = molecules[0]
  single_log2_sketch_size = log2_sketch_sizes[0]
}

// Tenx parameters
tenx_tags = params.tenx_tags
tenx_cell_barcode_pattern = params.tenx_cell_barcode_pattern
tenx_molecular_barcode_pattern = params.tenx_molecular_barcode_pattern
tenx_min_umi_per_cell = params.tenx_min_umi_per_cell

if (params.splitKmer && 'protein' in molecules){
  exit 1, "Cannot specify 'protein' in `--molecules` if --splitKmer is set"
}


// For bam files, set a folder name to save the optional barcode metadata csv
if (!params.write_barcode_meta_csv) {
  barcode_metadata_folder = ""
}
else {
  barcode_metadata_folder = "barcode_metadata"
}

// For bam files, one_signature_per_record is true
if (params.bam) {
  one_signature_per_record = true
}
else {
  one_signature_per_record = params.one_signature_per_record
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
if(params.barcodes_file)          summary["Barcodes"]              = params.barcodes_file
if(params.rename_10x_barcodes)    summary["Renamer barcodes"]      = params.rename_10x_barcodes
if(params.read_paths)   summary['Read paths (paired-end)']         = params.read_paths
// Sketch parameters
summary['Skip trimming?'] = params.skip_trimming
summary['K-mer sizes']            = params.ksizes
summary['Molecule']               = params.molecules
summary['Log2 Sketch Sizes']      = params.log2_sketch_sizes
summary['One Sig per Record']     = params.one_signature_per_record
summary['Track Abundance']        = params.track_abundance
// 10x parameters
if(params.tenx_tgz) summary["10x .tgz"] = params.tenx_tgz
if(params.tenx_tgz) summary["10x SAM tags"] = params.tenx_tags
if(params.tenx_tgz) summary["10x Cell pattern"] = params.tenx_cell_barcode_pattern
if(params.tenx_tgz) summary["10x UMI pattern"] = params.tenx_molecular_barcode_pattern
if(params.tenx_tgz) summary['Min UMI/cell'] = params.tenx_min_umi_per_cell
// Extract coding parameters
if(params.peptide_fasta) summary["Peptide fasta"] = params.peptide_fasta
if(params.peptide_fasta) summary['Peptide ksize'] = params.translate_peptide_ksize
if(params.peptide_fasta) summary['Peptide molecule'] = params.translate_peptide_molecule
if(params.peptide_fasta) summary['Bloom filter table size'] = params.bloomfilter_tablesize
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
    id: 'nf-core-kmermaid-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/kmermaid Workflow Summary'
    section_href: 'https://github.com/nf-core/kmermaid'
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
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    sourmash info &> v_sourmash.txt
    bam2fasta info &> v_bam2fasta.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

if (params.peptide_fasta){
  process make_protein_index {
    tag "${peptides}__${bloom_id}"
    label "low_memory"

    publishDir "${params.outdir}/protein_index", mode: 'copy'

    input:
    file(peptides) from ch_peptide_fasta
    peptide_ksize
    peptide_molecule

    output:
    set val(bloom_id), val(peptide_molecule), file("${peptides.simpleName}__${bloom_id}.bloomfilter") into ch_khtools_bloom_filter

    script:
    bloom_id = "molecule-${peptide_molecule}_ksize-${peptide_ksize}"
    """
    khtools index \\
      --tablesize ${bloomfilter_tablesize} \\
      --molecule ${peptide_molecule} \\
      --peptide-ksize ${peptide_ksize} \\
      --save-as ${peptides.simpleName}__${bloom_id}.bloomfilter \\
      ${peptides}
    """
  }
}

if (params.bam) {
  process bam2fasta {
    tag "$sample_id"
    label "high_memory"
    publishDir "${params.outdir}/${params.save_fastas}/", pattern: '*.fastq.gz', mode: 'copy'
    publishDir "${params.outdir}/${barcode_metadata_folder}", pattern: '*.csv', mode: 'copy'


    // If job fails, try again with more memory
    // memory { 8.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 1

    input:
    file(barcodes_file) from barcodes_ch
    file(bam) from tenx_bam_for_unaligned_fastq_ch
    file(rename_10x_barcodes) from rename_10x_barcodes_ch

    output:
    set val(sample_id), file("${params.save_fastas}/*.fastq.gz") into per_cell_fastqs_ch
    // https://github.com/nextflow-io/patterns/blob/master/docs/optional-output.adoc
    file("${params.write_barcode_meta_csv}") optional true

    script:
    sample_id = "${bam.simpleName}"
    tenx_min_umi_per_cell = params.tenx_min_umi_per_cell ? "--min-umi-per-barcode ${params.tenx_min_umi_per_cell}" : ''
    shard_size = params.shard_size ? "--shard-size ${params.shard_size}" : ''
    metadata = params.write_barcode_meta_csv ? "--write-barcode-meta-csv ${params.write_barcode_meta_csv}": ''
    save_fastas = "--save-fastas ${params.save_fastas}"
    save_intermediate_files = "--save-intermediate-files ${params.save_intermediate_files}"
    processes = "--processes ${params.max_cpus}"
    output_format = "--output-format fastq.gz"
    channel_id = "--channel-id ${sample_id}"

    def barcodes_file = params.barcodes_file ? "--barcodes-file ${barcodes_file.baseName}.tsv": ''
    def rename_10x_barcodes = params.rename_10x_barcodes ? "--rename-10x-barcodes ${rename_10x_barcodes.baseName}.tsv": ''
    """
      bam2fasta percell \\
        $processes \\
        $tenx_min_umi_per_cell \\
        $shard_size \\
        $rename_10x_barcodes \\
        $barcodes_file \\
        $save_fastas \\
        $output_format \\
        $channel_id \\
        $save_intermediate_files \\
        $metadata \\
        --filename ${bam}
        find ${params.save_fastas}/ -type f -name "*.fastq.gz" | while read src; do if [[ \$src == *"|"* ]]; then mv "\$src" \$(echo "\$src" | tr "|" "_"); fi done
    """
  }

  if (params.skip_trimming) {
    reads_ch = ch_non_bam_reads.concat(per_cell_fastqs_ch)
  } else {
    ch_read_files_trimming = ch_non_bam_reads.concat(per_cell_fastqs_ch)
  }
}

if (params.tenx_tgz) {
  process tenx_tgz_extract_bam {
    tag "$sample_id"
    publishDir "${params.outdir}/10x-bams", mode: 'copy'

    input:
    file(tenx_tgz) from tenx_tgz_ch

    output:
    set val(sample_id), file(bam) into tenx_bam_for_unaligned_fastq_ch, tenx_bam_for_aligned_fastq_ch
    file(bai)
    set val(sample_id), file(barcodes) into tenx_bam_barcodes_ch

    script:
    sample_id = "${tenx_tgz.simpleName}"
    bam = "${sample_id}__possorted_genome_bam.bam"
    bai = "${sample_id}__possorted_genome_bam.bam.bai"
    barcodes = "${sample_id}__barcodes.tsv"
    """
    tar xzvf ${tenx_tgz} \\
      ${sample_id}/outs/possorted_genome_bam.bam.bai \\
      ${sample_id}/outs/possorted_genome_bam.bam \\
      ${sample_id}/outs/filtered_gene_bc_matrices
    # Rename the files so there aren't conflicting duplicate filenames for the future
    mv ${sample_id}/outs/possorted_genome_bam.bam ${bam}
    mv ${sample_id}/outs/possorted_genome_bam.bam.bai ${bai}
    mv ${sample_id}/outs/filtered_gene_bc_matrices/*/barcodes.tsv ${barcodes}
    """
  }
}

if (params.tenx_tgz) {
  process samtools_fastq_aligned {
    tag "${channel_id}"
    publishDir "${params.outdir}/10x-fastqs/per-channel/aligned", mode: 'copy'
    label "mid_cpu"

    input:
    set val(channel_id), file(bam) from tenx_bam_for_unaligned_fastq_ch

    output:
    set val(channel_id), val("aligned"), file(reads) into tenx_reads_aligned_counting_ch, tenx_reads_aligned_concatenation_ch

    script:
    reads = "${channel_id}__aligned.fastq.gz"
    """
    samtools view -ub -F 256 -q 255 ${bam} \\
        | samtools fastq --threads ${task.cpus} -T ${tenx_tags} - \\
        | gzip -c - \\
          > ${reads}
    """
  }

  process samtools_fastq_unaligned {
    tag "${channel_id}"
    publishDir "${params.outdir}/10x-fastqs/per-channel/unaligned", mode: 'copy'
    label "mid_cpu"

    input:
    set val(channel_id), file(bam) from tenx_bam_for_aligned_fastq_ch

    output:
    set val(channel_id), val("unaligned"), file(reads) into tenx_reads_unaligned_unfiltered_ch

    script:
    reads = "${channel_id}__unaligned.fastq.gz"
    """
    samtools view -f4 ${bam} \\
      | grep -E '${tenx_cell_barcode_pattern}' \\
      | samtools fastq --threads ${task.cpus} -T ${tenx_tags} - \\
      | gzip -c - \\
        > ${reads} \\
      || touch ${reads}
    """
    // The '||' means that if anything in the previous step fails, do the next thing
    // It's bash magic from: https://stackoverflow.com/a/3822649/1628971
  }
  // Remove empty files
  // it[0] = channel_id
  // it[1] = "unaligned"
  // it[2] = read file
  // gzipped files are 20 bytes when "empty" due to the header
  tenx_reads_unaligned_unfiltered_ch.filter{ it -> it[2].size() > 20 }
    .set{ tenx_reads_unaligned_ch }

  // Put fastqs from aligned and unaligned reads into a single channel
  tenx_reads_aligned_concatenation_ch.mix(tenx_reads_unaligned_ch)
    .dump(tag: "tenx_reads_ch")
    .set{ tenx_reads_ch }

  if ((params.tenx_min_umi_per_cell > 0) || !params.barcodes_file) {
    process count_umis_per_cell {
      tag "${is_aligned_channel_id}"
      label 'low_memory_long'

      publishDir "${params.outdir}/10x-fastqs/umis-per-cell/", mode: 'copy'

      input:
      set val(channel_id), val(is_aligned), file(reads) from tenx_reads_aligned_counting_ch

      output:
      file(umis_per_cell)
      set val(channel_id), file(good_barcodes) into good_barcodes_unfiltered_ch

      script:
      is_aligned_channel_id = "${channel_id}__${is_aligned}"
      umis_per_cell = "${is_aligned_channel_id}__n_umi_per_cell.csv"
      good_barcodes = "${is_aligned_channel_id}__barcodes.tsv"
      """
      bam2fasta count_umis_percell \\
          --filename ${reads} \\
          --min-umi-per-barcode ${tenx_min_umi_per_cell} \\
          --cell-barcode-pattern '${tenx_cell_barcode_pattern}' \\
          --molecular-barcode-pattern '${tenx_molecular_barcode_pattern}' \\
          --write-barcode-meta-csv ${umis_per_cell} \\
          --barcodes-significant-umis-file ${good_barcodes}
      """
    }
    // Make sure good barcodes file is nonempty so next step doesn't start
    // it[0] = channel id
    // it[1] = good_barcodes file
    good_barcodes_unfiltered_ch.filter{ it -> it[1].size() > 0 }
      .ifEmpty{ exit 1, "No cell barcodes found with at least ${tenx_min_umi_per_cell} molecular barcodes (UMIs) per cell"}
      .set{ good_barcodes_ch }

  } else {
    // Use barcodes extracted from the tenx .tgz file
    good_barcodes_ch = tenx_bam_barcodes_ch
  }

  tenx_reads_ch.combine(good_barcodes_ch, by: 0)
    .dump(tag: 'tenx_reads_ch__combine__good_barcodes_ch')
    .set{ tenx_reads_with_good_barcodes_ch }

  process extract_per_cell_fastqs {
    tag "${is_aligned_channel_id}"
    label "high_memory_long"
    publishDir "${params.outdir}/10x-fastqs/per-cell/${channel_id}/", mode: 'copy', pattern: '*.fastq', saveAs: { filename -> "${filename.replace("|", "-")}"}

    input:
    // Example input:
    // ['mouse_lung', 'aligned', mouse_lung__aligned.fastq.gz, mouse_lung__aligned__barcodes.tsv]
    set val(channel_id), val(is_aligned), file(reads), file(barcodes) from tenx_reads_with_good_barcodes_ch

    output:
    file('*.fastq.gz') into per_channel_cell_reads_ch

    script:
    is_aligned_channel_id = "${channel_id}__${is_aligned}"
    """
    bam2fasta make_fastqs_percell \\
        --filename ${reads} \\
        --barcodes-significant-umis-fil ${barcodes} \\
        --cell-barcode-pattern '${tenx_cell_barcode_pattern}' \\
        --channel-id ${is_aligned_channel_id} \\
        --output-format 'fastq.gz'

    # Decoy file just in case there are no reads found,
    # to prevent this process from erroring out
    touch empty.fastq.gz
    """
  }
  // Make per-cell fastqs into a flat channel that matches the read channels of yore
  per_channel_cell_reads_ch
    .dump(tag: 'per_channel_cell_reads_ch')
    .flatten()
    .filter{ it -> it.size() > 0 }   // each item is just a single file, no need to do it[1]
    .map{ it -> tuple(it.simpleName, file(it)) }
    .dump(tag: 'per_cell_fastqs_ch')
    .set{ per_cell_fastqs_ch }

  if (params.skip_trimming) {
    reads_ch = ch_non_bam_reads.concat(per_cell_fastqs_ch)
  } else {
    ch_read_files_trimming = ch_non_bam_reads.concat(per_cell_fastqs_ch)
  }
}



if (!params.skip_trimming){
  process fastp {
      label 'process_low'
      tag "$name"
      publishDir "${params.outdir}/fastp", mode: 'copy',
        saveAs: {filename ->
                    if (filename.indexOf(".fastq.gz") == -1) "logs/$filename"
                    else if (reads[1] == null) "single_end/$filename"
                    else if (reads[1] != null) "paired_end/$filename"
                    else null
                }

      input:
      set val(name), file(reads) from ch_read_files_trimming

      output:
      set val(name), file("*trimmed.fastq.gz") into ch_reads_trimmed
      file "*fastp.json" into ch_fastp_results
      file "*fastp.html" into ch_fastp_html

      script:
      // One set of reads --> single end
      if (reads[1] == null) {
          """
          fastp \\
              --in1 ${reads} \\
              --out1 ${name}_R1_trimmed.fastq.gz \\
              --json ${name}_fastp.json \\
              --html ${name}_fastp.html
          """
      } else if (reads[1] != null ){
        // More than one set of reads --> paired end
          """
          fastp \\
              --in1 ${reads[0]} \\
              --in2 ${reads[1]} \\
              --out1 ${name}_R1_trimmed.fastq.gz \\
              --out2 ${name}_R2_trimmed.fastq.gz \\
              --json ${name}_fastp.json \\
              --html ${name}_fastp.html
          """
      } else {
        """
        echo name ${name}
        echo reads: ${reads}
        echo "Number of reads is not equal to 1 or 2 --> don't know how to trim non-paired-end and non-single-end reads"
        """
      }
  }
  // Concatenate trimmed fastq files with fastas
  if (params.subsample){
    // Concatenate trimmed reads with fastas for subsequent subsampling
    subsample_reads_ch = ch_reads_trimmed.concat(fastas_ch)
  } else {
    // Concatenate trimmed reads with fastas for signature generation
    reads_ch = ch_reads_trimmed.concat(fastas_ch)
  }
}

if (params.subsample) {
    process subsample_input {
	tag "${id}_subsample"
	publishDir "${params.outdir}/seqtk/", mode: 'copy'

	input:
	set id, file(reads) from subsample_reads_ch

	output:

	set val(id), file("*_${params.subsample}.fastq.gz") into reads_ch

	script:
	read1 = reads[0]
	read2 = reads[1]
	read1_prefix = read1.name.minus(".fastq.gz") // TODO: change to RE to match fasta as well?
	read2_prefix = read2.name.minus(".fastq.gz")

  """
  seqtk sample -s100 ${read1} ${params.subsample} > ${read1_prefix}_${params.subsample}.fastq.gz
  seqtk sample -s100 ${read2} ${params.subsample} > ${read2_prefix}_${params.subsample}.fastq.gz
  """
  }
}



if (params.peptide_fasta){
  process translate {
    tag "${sample_id}"
    label "low_memory_long"
    publishDir "${params.outdir}/translate/", mode: 'copy'

    input:
    set bloom_id, molecule, file(bloom_filter) from ch_khtools_bloom_filter.collect()
    set sample_id, file(reads) from reads_ch

    output:
    // TODO also extract nucleotide sequence of coding reads and do sourmash compute using only DNA on that?
    set val(sample_id), file("${sample_id}__coding_reads_peptides.fasta") into ch_coding_peptides
    set val(sample_id), file("${sample_id}__coding_reads_nucleotides.fasta") into ch_coding_nucleotides
    set val(sample_id), file("${sample_id}__coding_scores.csv") into ch_coding_scores_csv
    set val(sample_id), file("${sample_id}__coding_summary.json") into ch_coding_scores_json

    script:
    """
    khtools translate \\
      --molecule ${molecule} \\
      --coding-nucleotide-fasta ${sample_id}__coding_reads_nucleotides.fasta \\
      --csv ${sample_id}__coding_scores.csv \\
      --json-summary ${sample_id}__coding_summary.json \\
      --jaccard-threshold ${jaccard_threshold} \\
      --peptides-are-bloom-filter \\
      ${bloom_filter} \\
      ${reads} > ${sample_id}__coding_reads_peptides.fasta
    """
  }
  // Remove empty files
  // it[0] = sample id
  // it[1] = sequence fasta file
  ch_coding_nucleotides_nonempty = ch_coding_nucleotides.filter{ it[1].size() > 0 }
  ch_coding_peptides_nonempty = ch_coding_peptides.filter{ it[1].size() > 0 }

} else {
  // Send reads directly into coding/noncoding
  reads_ch
    .set{ ch_coding_nucleotides_nonempty }
}

if (params.splitKmer){
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     CREATE SKA SKETCH                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  process ska_compute_sketch {
      tag "${sketch_id}"
      publishDir "${params.outdir}/ska/sketches/", mode: 'copy'
      errorStrategy 'retry'
      maxRetries 3


    input:
    each ksize from ksizes
    set id, file(reads) from reads_ch

    output:
    set val(ksize), file("${sketch_id}.skf") into ska_sketches

    script:
    sketch_id = "${id}_ksize_${ksize}"

      """
      ska fastq \\
        -k $ksize \\
        -o ${sketch_id} \\
        ${reads}
      """

    }
} else {
process sourmash_compute_sketch_fastx_nucleotide {
  tag "${sample_id}_${sketch_id}"
  label "mid_memory"
  publishDir "${params.outdir}/sketches_nucleotide/${sketch_id}", mode: 'copy'

  input:
  each ksize from ksizes
  each log2_sketch_size from log2_sketch_sizes
  set sample_id, file(reads) from ch_coding_nucleotides_nonempty

  output:
  set val(sketch_id), val("dna"), val(ksize), val(log2_sketch_size), file("${sample_id}_${sketch_id}.sig") into sourmash_sketches_nucleotide

  script:
  // Don't calculate DNA signature if this is protein, to minimize disk,
  // memory and IO requirements in the future
  ksize = ksize
  sketch_id = "molecule-dna_ksize-${ksize}_log2sketchsize-${log2_sketch_size}_trackabundance-${params.track_abundance}"
  track_abundance_flag = track_abundance ? '--track-abundance' : ''

  if ( params.one_signature_per_record ){
    """
    sourmash compute \\
      --num-hashes \$((2**$log2_sketch_size)) \\
      --ksizes $ksize \\
      --dna \\
      $track_abundance_flag \\
      --output ${sample_id}_${sketch_id}.sig \\
      $reads
    """
  }
  else {
    """
    sourmash compute \\
      --num-hashes \$((2**$log2_sketch_size)) \\
      --ksizes $ksize \\
      --dna \\
      $track_abundance_flag \\
      --output ${sample_id}_${sketch_id}.sig \\
      --merge '$sample_id' \\
      $reads
    """
    }
  }
}


if (params.peptide_fasta){
  process sourmash_compute_sketch_fastx_peptide {
    tag "${sample_id}_${sketch_id}"
    label "low_memory"
    publishDir "${params.outdir}/sketches_peptide/${sketch_id}", mode: 'copy'

    input:
    each ksize from ksizes
    each molecule from peptide_molecules
    each log2_sketch_size from log2_sketch_sizes
    set sample_id, file(reads) from ch_coding_peptides_nonempty

    output:
    set val(sketch_id), val(molecule), val(ksize), val(log2_sketch_size), file("${sample_id}_${sketch_id}.sig") into sourmash_sketches_peptide

    script:
    sketch_id = "molecule-${molecule}_ksize-${ksize}_log2sketchsize-${log2_sketch_size}_trackabundance-${params.track_abundance}"
    molecule = molecule
    ksize = ksize
    track_abundance_flag = track_abundance ? '--track-abundance' : ''

    if ( params.one_signature_per_record ) {
      """
      sourmash compute \\
        --num-hashes \$((2**$log2_sketch_size)) \\
        --ksizes $ksize \\
        --input-is-protein \\
        --$molecule \\
        --no-dna \\
        $track_abundance_flag \\
        --output ${sample_id}_${sketch_id}.sig \\
        $reads
      """
    }
    else {
      """
      sourmash compute \\
        --num-hashes \$((2**$log2_sketch_size)) \\
        --ksizes $ksize \\
        --input-is-protein \\
        --$molecule \\
        --no-dna \\
        $track_abundance_flag \\
        --output ${sample_id}_${sketch_id}.sig \\
        --merge '$sample_id' \\
        $reads
      """
    }
  }
} else {
  sourmash_sketches_peptide = Channel.empty()
}


if (params.splitKmer){
     process ska_compare_sketches {
    tag "${sketch_id}"
    publishDir "${params.outdir}/ska/compare/", mode: 'copy'

    input:
    set val(ksize), file (sketches) from ska_sketches.groupTuple()

    output:
     // uploaded distances, clusters, and graph connecting (dot) file
    file "ksize_${ksize}*"

    script:
    """
    ska distance -o ksize_${ksize} -s 25 -i 0.95 ${sketches}
    """

    }
  } else {
  process sourmash_compare_sketches {
    // Combine peptide and nucleotide sketches
    sourmash_sketches = sourmash_sketches_peptide.concat(sourmash_sketches_nucleotide)
    tag "${sketch_id}"
    publishDir "${params.outdir}/sourmash/compare", mode: 'copy'

    input:
    set val(sketch_id), val(molecule), val(ksize), val(log2_sketch_size), file ("sourmash/sketches/*.sig") \
      from sourmash_sketches.groupTuple(by: [0, 3])

    output:
    file "similarities_${sketch_id}.csv"

    script:
    processes = "--processes ${task.cpus}"
    """
    sourmash compare \\
          --ksize ${ksize[0]} \\
          --${molecule[0]} \\
          --csv similarities_${sketch_id}.csv \\
          --traverse-directory .
    """

  }
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
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision //
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
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( "${output_d}/pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( "${output_d}/pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
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
    ${c_purple}  nf-core/kmermaid v${workflow.manifest.version}${c_reset}
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
