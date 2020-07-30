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
      --protein_fastas              Path to protein fasta inputs
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
      --track_abundance             Track abundance of each hashed k-mer, could be useful for cancer RNA-seq or ATAC-seq analyses
      --skip_trimming               If provided, skip fastp trimming of reads
      --skip_compare                If provided, skip comparison of hashes using sourmash compare
      --skip_compute                If provided, skip computing of signatures using sourmash compute

     Sketch size options:
      --sketch_num_hashes                 Number of hashes to use for making the sketches.
                                    Mutually exclusive with --sketch_num_hashes_log2
      --sketch_num_hashes_log2            Which log2 sketch sizes to use. Multiple are separated by commas.
                                    Default is '10,12,14,16'. Mutually exclusive with --sketch_num_hashes
      --sketch_scaled               Observe every 1/N hashes per sample, rather than a "flat rate" of N hashes
                                    per sample. This way, the number of hashes scales by the sequencing depth.
                                    Mutually exclusive with --sketch_scaled_log2
      --sketch_scaled_log2          Same as --sketch_scaled, but instead of specifying the true number of hashes,
                                    specify the power to take 2 to. Mutually exlusive with --sketch_scaled


    Split K-mer options:
      --split_kmer                   If provided, use SKA to compute split k-mer sketches instead of
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

    Translate RNA-seq reads into protein-coding sequences options:
      --reference_proteome_fasta    Path to a well-curated fasta file of protein sequences. Used to filter for coding reads
      --translate_peptide_ksize     K-mer size to use for translating RNA into protein.
                                    Default: 9, which is good for 'protein'. If using dayhoff, suggest 15
      --translate_peptide_molecule  Which molecular encoding to use for translating. Default: "protein"
                                    If your reference proteome is quite different from your species of interest,
                                    suggest using "dayhoff" encoding
      --translate_jaccard_threshold Minimum fraction of overlapping translated k-mers from the read to match to the reference. Default: 0.95
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

// Boolean for if an nucleotide input exists anywhere
have_nucleotide_fasta_input = params.fastas || params.fasta_paths
have_nucleotide_input = params.read_paths || params.sra || params.csv_pairs || params.csv_singles || params.read_pairs || params.read_singles || have_nucleotide_fasta_input || params.bam || params.tenx_tgz

if (!params.split_kmer){
  have_sketch_num_hashes = params.sketch_num_hashes || params.sketch_num_hashes_log2 || params.sketch_scaled || params.sketch_scaled_log2
  if (!have_sketch_num_hashes) {
    exit 1, "Must provide one of --sketch_num_hashes, --sketch_num_hashes_log2, --sketch_scaled, --sketch_scaled_log2 for Sourmash!"
  }
}

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
       .dump ( tag: 'fastas_ch' )
       .ifEmpty { exit 1, "params.fastas (${params.fastas}) was empty - no input files supplied" }
   } else if (params.fasta_paths) {
     fastas_ch = Channel
       .from(params.fasta_paths)
       .map { row -> if (row[1].size() == 2) [ row[0], [file(row[1][0]), file(row[1][1])]]
             else [row[0], [file(row[1][0])]]}
       .dump ( tag: 'fastas_ch' )
       .ifEmpty { exit 1, "params.fasta_paths (${params.fastas}) was empty - no input files supplied" }
   }

  if (params.bam) {
  Channel.fromPath(params.bam, checkIfExists: true)
        .map{ f -> tuple(f.baseName, tuple(file(f))) }
       .ifEmpty { exit 1, "Bam file not found: ${params.bam}" }
       .dump( tag: 'bam' )
       .into{ tenx_bam_for_unaligned_fastq_ch; tenx_bam_for_aligned_fastq_ch}
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
    have_nucleotide_input = true
    Channel.fromPath(params.tenx_tgz, checkIfExists: true)
       .dump(tag: 'tenx_tgz_before_mri_filter')
       .filter{ ~/.+[^mri]\.tgz/ }
       .ifEmpty { exit 1, "10X .tgz file not found: ${params.tenx_tgz}" }
       .dump(tag: 'tenx_tgz_after_mri_filter')
       .set{ tenx_tgz_ch }
  }
}

////////////////////////////////////////////////////
/* --          Parse protein fastas            -- */
////////////////////////////////////////////////////
if (params.protein_fastas){
  Channel.fromPath(params.protein_fastas?.toString()?.tokenize(';'))
      .map{ f -> tuple(f.baseName, tuple(file(f))) }
      .ifEmpty { exit 1, "params.protein_fastas was empty - no input files supplied" }
      .set { ch_protein_fastas }
} else if (params.protein_fasta_paths){
  Channel
    .from(params.protein_fasta_paths)
    .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true)] ] }
    .ifEmpty { exit 1, "params.protein_fasta_paths was empty - no input files supplied" }
    .dump(tag: "protein_fasta_paths")
    .set { ch_protein_fastas }
} else {
  ch_protein_fastas = Channel.empty()
}

if (params.reference_proteome_fasta) {
Channel.fromPath(params.reference_proteome_fasta, checkIfExists: true)
     .ifEmpty { exit 1, "Reference proteome file not found: ${params.reference_proteome_fasta}" }
     .set{ ch_reference_proteome_fasta }
}

////////////////////////////////////////////////////
/* --    Concatenate all nucleotide inputs     -- */
////////////////////////////////////////////////////
// Add _unchecked suffix because have not yet checked if these files are not empty
if (params.subsample) {
  if (params.bam){
     exit 1, "Cannot provide both a bam file with --bam and specify --subsample"
  } else {
    if (params.skip_trimming){
      sra_ch.concat(csv_pairs_ch, csv_singles_ch, read_pairs_ch,
        read_singles_ch, fastas_ch, read_paths_ch)
        .set{ subsample_reads_ch_unchecked }
    } else {
      sra_ch.concat(
          csv_pairs_ch, csv_singles_ch, read_pairs_ch,
          read_singles_ch, read_paths_ch)
        // .ifEmpty{ exit 1, "No reads provided! Check read input files"}
        .set{ ch_read_files_trimming_unchecked }
    }
  }
} else {
  if (!(params.tenx_tgz || params.bam)) {
    if(params.skip_trimming){
      sra_ch.concat(
          csv_pairs_ch, csv_singles_ch, read_pairs_ch,
          read_singles_ch, fastas_ch, read_paths_ch)
       .set{ reads_ch_unchecked }
    } else {
      if (have_nucleotide_fasta_input) {
        // With fasta files - combine everything that can be trimmed
        sra_ch.concat(
            csv_pairs_ch, csv_singles_ch, read_pairs_ch,
            read_singles_ch, read_paths_ch)
          .dump ( tag: 'ch_read_files_trimming_unchecked__with_fastas' )
          .into { ch_read_files_trimming_to_trim; ch_read_files_trimming_to_check_size }
      } else {
        // No fasta files - combine everything and error out
        sra_ch.concat(
            csv_pairs_ch, csv_singles_ch, read_pairs_ch,
            read_singles_ch, read_paths_ch)
          .dump ( tag: 'ch_read_files_trimming_unchecked__no_fastas' )
          .set{ ch_read_files_trimming_unchecked }
      }
    }
  } else {
      sra_ch.concat(
          csv_pairs_ch, csv_singles_ch, read_pairs_ch,
          read_singles_ch, read_paths_ch)
        .dump ( tag: 'ch_non_bam_reads_unchecked__concatenated' )
        .set{ ch_non_bam_reads_unchecked }
    }
}

protein_input = params.protein_fastas || params.protein_fasta_paths
if (!protein_input) {
  if (params.subsample && params.skip_trimming ) {
    subsample_reads_ch_unchecked
      .ifEmpty{  exit 1, "No reads provided! Check read input files" }
      .set { subsample_ch_reads_for_ribosomal_removal }
  }
  if (params.skip_trimming && !(params.bam || params.tenx_tgz)) {
    reads_ch_unchecked
      .ifEmpty{ exit 1, "No reads provided! Check read input files" }
      .set { ch_reads_for_ribosomal_removal }
    ch_read_files_trimming_to_check_size = Channel.empty()
  } else if (params.bam || params.tenx_tgz) {
    ch_non_bam_reads_unchecked
      // No need to check if empty since there is bam input
      .set { ch_non_bam_reads }
  } else if (!have_nucleotide_fasta_input) {
      // if no fastas, then definitely trimming the remaining reads
      ch_read_files_trimming_unchecked
        .ifEmpty{ exit 1, "No reads provided! Check read input files" }
        .into { ch_read_files_trimming_to_trim; ch_read_files_trimming_to_check_size }
  }
} else {
  // Since there exists protein input, don't check if these are empty
  if (params.subsample) {
    subsample_reads_ch_unchecked
      .set { subsample_ch_reads_for_ribosomal_removal }
  }
  if (params.skip_trimming) {
    reads_ch_unchecked
      .set { ch_reads_for_ribosomal_removal }
    ch_read_files_trimming_to_check_size = Channel.empty()
  } else if (!have_nucleotide_fasta_input) {
    ch_read_files_trimming_unchecked
      .into { ch_read_files_trimming_to_trim; ch_read_files_trimming_to_check_size }
  }
  if (params.bam) {
    ch_non_bam_reads_unchecked
      .set { ch_non_bam_reads }
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


if (params.split_kmer){
    params.ksizes = '15,9'
    params.molecules = 'dna'
} else {
    params.ksizes = '21,27,33,51'
}

// Get rRNA databases
// Default is set to bundled DB list in `assets/rrna-db-defaults.txt`

rRNA_database = file(params.rrna_database_manifest)
if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}
Channel
    .from( rRNA_database.readLines() )
    .map { row -> file(row) }
    .set { sortmerna_fasta }


// Parse the parameters

ksizes = params.ksizes?.toString().tokenize(',')
molecules = params.molecules?.toString().tokenize(',')
peptide_molecules = molecules.findAll { it != "dna" }


def make_sketch_id (molecule, ksize, sketch_value, track_abundance, sketch_style) {
  if (sketch_style == 'size') {
    sketch_value = "num_hashes-${sketch_value}"
  } else {
    sketch_value = "scaled-${sketch_value}"
  }

  sketch_id = "molecule-${molecule}__ksize-${ksize}__${sketch_value}__track_abundance-${track_abundance}"
  return sketch_id
}

// Create the --num-hashes or --scaled flag for sourmash
def make_sketch_value_flag(sketch_style, sketch_value) {
  if (sketch_style == "size") {
    number_flag = "--num-hashes ${sketch_value}"
  } else if (sketch_style == "scaled" ) {
    number_flag = "--scaled ${sketch_value}"
  } else {
    exit 1, "${sketch_style} is not a valid sketch counting style! Only 'scaled' and 'size' are valid"
  }
  return number_flag
}

int bloomfilter_tablesize = Math.round(Float.valueOf(params.bloomfilter_tablesize))

peptide_ksize = params.translate_peptide_ksize
peptide_molecule = params.translate_peptide_molecule
jaccard_threshold = params.translate_jaccard_threshold
track_abundance = params.track_abundance

// Parse sketch value and style parameters
sketch_num_hashes = params.sketch_num_hashes
sketch_num_hashes_log2 = params.sketch_num_hashes_log2
sketch_scaled = params.sketch_scaled
sketch_scaled_log2 = params.sketch_scaled_log2
have_sketch_value = params.sketch_num_hashes || params.sketch_num_hashes_log2 || params.sketch_scaled || params.sketch_scaled_log2

if (!have_sketch_value && !params.split_kmer) {
  exit 1, "None of --sketch_num_hashes, --sketch_num_hashes_log2, --sketch_scaled, --sketch_scaled_log2 was provided! Provide one (1) and only one to specify the style and amount of hashes per sourmash sketch"
}


// Tenx parameters
tenx_tags = params.tenx_tags
tenx_cell_barcode_pattern = params.tenx_cell_barcode_pattern
tenx_molecular_barcode_pattern = params.tenx_molecular_barcode_pattern
tenx_min_umi_per_cell = params.tenx_min_umi_per_cell

if (params.split_kmer && 'protein' in molecules){
  exit 1, "Cannot specify 'protein' in `--molecules` if --split_kmer is set"
}


// For bam files, set a folder name to save the optional barcode metadata csv
if (!params.write_barcode_meta_csv) {
  barcode_metadata_folder = ""
}
else {
  barcode_metadata_folder = "barcode_metadata"
}


// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

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
if(params.protein_fastas) summary["Protein FASTAs"]                = params.protein_fastas
if(params.bam)          summary["BAM"]                             = params.bam
if(params.barcodes_file)          summary["Barcodes"]              = params.barcodes_file
if(params.rename_10x_barcodes)    summary["Renamer barcodes"]      = params.rename_10x_barcodes
if(params.read_paths)   summary['Read paths (paired-end)']         = params.read_paths
// Sketch parameters
summary['Skip trimming?'] = params.skip_trimming
summary['Skip compare?'] = params.skip_compare
summary['Skip compute?'] = params.skip_compute
summary['Skip multiqc?'] = params.skip_multiqc
summary['K-mer sizes']            = params.ksizes
summary['Molecule']               = params.molecules
summary['Track Abundance']        = params.track_abundance
// -- Sketch size parameters --
if (params.sketch_num_hashes) summary['Sketch Sizes']                  = params.sketch_num_hashes
if (params.sketch_num_hashes_log2) summary['Sketch Sizes (log2)']      = params.sketch_num_hashes_log2
if (params.sketch_scaled) summary['Sketch scaled']               = params.sketch_scaled
if (params.sketch_scaled_log2) summary['Sketch scaled (log2)']   = params.sketch_scaled_log2
// 10x parameters
if(params.tenx_tgz) summary["10x .tgz"] = params.tenx_tgz
if(params.tenx_tgz) summary["10x SAM tags"] = params.tenx_tags
if(params.tenx_tgz) summary["10x Cell pattern"] = params.tenx_cell_barcode_pattern
if(params.tenx_tgz) summary["10x UMI pattern"] = params.tenx_molecular_barcode_pattern
if(params.tenx_tgz) summary['Min UMI/cell'] = params.tenx_min_umi_per_cell
// Extract coding parameters
if(params.reference_proteome_fasta) summary["Peptide fasta"] = params.reference_proteome_fasta
if(params.reference_proteome_fasta) summary['Peptide ksize'] = params.translate_peptide_ksize
if(params.reference_proteome_fasta) summary['Peptide molecule'] = params.translate_peptide_molecule
if(params.reference_proteome_fasta) summary['Bloom filter table size'] = params.bloomfilter_tablesize
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

if ( !params.split_kmer && have_sketch_value ) {
  // Only use this for sourmash sketches, not split k-mer sketches
  /*
   * Validate sketch sizes
   */
  process validate_sketch_values {
      publishDir "${params.outdir}/pipeline_info", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf(".txt") > 0) filename
          else null
      }
      input:
      val sketch_num_hashes
      val sketch_num_hashes_log2
      val sketch_scaled
      val sketch_scaled_log2

      output:
      file sketch_values into ch_sketch_values
      file sketch_style into ch_sketch_style

      script:
      sketch_style = "sketch_style.txt"
      sketch_values = 'sketch_values.txt'
      """
      validate_sketch_values.py \\
        --sketch_num_hashes ${sketch_num_hashes} \\
        --sketch_num_hashes_log2 ${sketch_num_hashes_log2} \\
        --sketch_scaled ${sketch_scaled} \\
        --sketch_scaled_log2 ${sketch_scaled_log2} \\
        --output ${sketch_values} \\
        --sketch_style ${sketch_style}
      """
  }

  // Parse sketch style into value
  ch_sketch_style
    .splitText()
    .dump ( tag: 'ch_sketch_style' )
    .map { it -> it.replaceAll('\\n', '' ) }
    // .first()
    .dump ( tag: 'sketch_style_parsed' )
    .into { sketch_style_for_nucleotides; sketch_style_for_proteins }
  // sketch_style = sketch_styles[0]
  // println "sketch_style_parsed: ${sketch_style_parsed}"
  // println "sketch_style: ${sketch_style}"

  // Parse file into values
  ch_sketch_values
    .splitText()
    .map { it -> it.replaceAll('\\n', '')}
    .dump ( tag : 'sketch_values_parsed' )
    .collect()
    .dump ( tag : 'sketch_values_parsed__collected' )
    .set { sketch_values }

}


Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-kmermaid-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/kmermaid Workflow Summary'
    section_href: 'https://github.com/nf-core/kmermaid'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }



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
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    bam2fasta info &> v_bam2fasta.txt
    fastp --version &> v_fastp.txt
    samtools --version &> v_samtools.txt
    ska version &> v_ska.txt
    sortmerna --version &> v_sortmerna.txt
    sourmash -v &> v_sourmash.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

if (params.reference_proteome_fasta){
  process make_protein_index {
    tag "${peptides}__${bloom_id}"
    label "low_memory"

    publishDir "${params.outdir}/protein_index", mode: 'copy'

    input:
    file(peptides) from ch_reference_proteome_fasta
    peptide_ksize
    peptide_molecule

    output:
    set val(bloom_id), val(peptide_molecule), file("${peptides.simpleName}__${bloom_id}.bloomfilter") into ch_sencha_bloom_filter

    script:
    bloom_id = "molecule-${peptide_molecule}_ksize-${peptide_ksize}"
    """
    sencha index \\
      --tablesize ${bloomfilter_tablesize} \\
      --molecule ${peptide_molecule} \\
      --peptide-ksize ${peptide_ksize} \\
      --save-as ${peptides.simpleName}__${bloom_id}.bloomfilter \\
      ${peptides}
    """
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

if (params.tenx_tgz || params.bam) {
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
    samtools view -ub -F 4 ${bam} \\
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
    set val(channel_id), val("unaligned"), file(reads) into tenx_reads_unaligned_ch

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

  // Put fastqs from aligned and unaligned reads into a single channel
  tenx_reads_aligned_concatenation_ch.mix(tenx_reads_unaligned_ch)
    .dump(tag: "tenx_ch_reads_for_ribosomal_removal")
    .set{ tenx_ch_reads_for_ribosomal_removal }

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

  } else if (params.barcodes) {
    good_barcodes_ch = barcodes_ch
  }
  else {
    // Use barcodes extracted from the tenx .tgz file
    good_barcodes_ch = tenx_bam_barcodes_ch
  }

  tenx_ch_reads_for_ribosomal_removal.combine(good_barcodes_ch, by: 0)
    .dump(tag: 'tenx_ch_reads_for_ribosomal_removal__combine__good_barcodes_ch')
    .set{ tenx_reads_with_good_barcodes_ch }

  process extract_per_cell_fastqs {
    tag "${is_aligned_channel_id}"
    label "high_memory_long"
    publishDir "${params.outdir}/10x-fastqs/per-cell/${channel_id}/", mode: 'copy', pattern: '*.fastq.gz', saveAs: { filename -> "${filename.replace("|", "-")}"}

    input:
    // Example input:
    // ['mouse_lung', 'aligned', mouse_lung__aligned.fastq.gz, mouse_lung__aligned__barcodes.tsv]
    set val(channel_id), val(is_aligned), file(reads), file(barcodes) from tenx_reads_with_good_barcodes_ch

    output:
    file('*.fastq.gz') into per_channel_cell_ch_reads_for_ribosomal_removal

    script:
    is_aligned_channel_id = "${channel_id}__${is_aligned}"
    def rename_10x_barcodes = params.rename_10x_barcodes ? "--rename-10x-barcodes ${rename_10x_barcodes.baseName}.tsv": ''
    processes = "--processes ${task.cpus}"
    """
    bam2fasta make_fastqs_percell \\
        $processes \\
        --filename ${reads} \\
        --barcodes-significant-umis-fil ${barcodes} \\
        $rename_10x_barcodes \\
        --cell-barcode-pattern '${tenx_cell_barcode_pattern}' \\
        --channel-id ${is_aligned_channel_id}__ \\
        --output-format 'fastq.gz'
    # Decoy file just in case there are no reads found,
    # to prevent this process from erroring out
    touch empty.fastq.gz
    """
  }
  // Make per-cell fastqs into a flat channel that matches the read channels of yore
  // Filtering out fastq.gz files less than 200 bytes (arbitary number)
  // ~200 bytes is about the size of a file with a single read or less
  // We can't use .size() > 0 because it's fastq.gz is gzipped content
  per_channel_cell_ch_reads_for_ribosomal_removal
    .dump(tag: 'per_channel_cell_ch_reads_for_ribosomal_removal')
    .flatten()
    .filter{ it -> it.size() > 200 }   // each item is just a single file, no need to do it[1]
    .map{ it -> tuple(it.simpleName, file(it)) }
    .dump(tag: 'per_cell_fastqs_ch')
    .set{ per_cell_fastqs_ch }

  if (params.skip_trimming) {
    ch_reads_for_ribosomal_removal = ch_non_bam_reads.concat(per_cell_fastqs_ch)
  } else {
    ch_non_bam_reads
      .mix ( per_cell_fastqs_ch )
      .dump ( tag: 'ch_non_bam_reads__per_cell_fastqs_ch' )
      .into{ ch_read_files_trimming_to_trim; ch_read_files_trimming_to_check_size }
  }
}


if ( have_nucleotide_input ) {
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
        set val(name), file(reads) from ch_read_files_trimming_to_trim

        output:
        set val(name), file("*trimmed.fastq.gz") into ch_reads_all_trimmed
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

  // Filtering out fastq.gz files less than 200 bytes (arbitary number)
  // ~200 bytes is about the size of a file with a single read or less
  // We can't use .size() > 0 because it's fastq.gz is gzipped content
  ch_reads_all_trimmed
    .dump ( tag: 'ch_reads_all_trimmed' )
    .branch {
      // Paired is a tuple of two reads
      paired: it[1].size() == 2
      single: true
    }
    .set { ch_reads_trimmed_branched }

  ch_reads_trimmed_branched.paired
    .filter{ it -> it[1][0].size() > 200 }
    .dump ( tag: 'ch_reads_trimmed_paired' )
    .set{ ch_reads_trimmed_paired }

  ch_reads_trimmed_branched.single
    .filter{ it -> it[1].size() > 200 }
    .dump ( tag: 'ch_reads_trimmed_single' )
    .set{ ch_reads_trimmed_single }

  ch_reads_trimmed_single
    .mix ( ch_reads_trimmed_paired )
    .set { ch_reads_trimmed }

  // Concatenate trimmed fastq files with fastas
  if (params.subsample){
    // Concatenate trimmed reads with fastas for subsequent subsampling
    ch_reads_trimmed
      .concat( fastas_ch )
      .dump ( tag: 'trimmed_reads__concat_fastas' )
      .set { subsample_ch_reads_for_ribosomal_removal }
  } else {
    // Concatenate trimmed reads with fastas for signature generation
    ch_reads_for_ribosomal_removal = ch_reads_trimmed.concat(fastas_ch)
  }
} else {
  ch_fastp_results = Channel.from(false)
}

if (params.subsample) {
  process subsample_input {
  	tag "${id}_subsample"
  	publishDir "${params.outdir}/seqtk/", mode: 'copy'

  	input:
  	set val(id), file(reads) from subsample_ch_reads_for_ribosomal_removal

  	output:
  	set val(id), file("*_${params.subsample}.fastq.gz") into ch_reads_for_ribosomal_removal

  	script:
  	read1 = reads[0]
  	read2 = reads[1]
  	read1_prefix = read1.simpleName
  	read2_prefix = read2.simpleName

    """
    seqtk sample -s100 ${read1} ${params.subsample} > ${read1_prefix}_${params.subsample}.fastq.gz
    seqtk sample -s100 ${read2} ${params.subsample} > ${read2_prefix}_${params.subsample}.fastq.gz
    """
    }
  }

/*
 * STEP 2+ - SortMeRNA - remove rRNA sequences on request
 */
if (!params.remove_ribo_rna) {
    ch_reads_for_ribosomal_removal
        .set { reads_ch }
    sortmerna_logs = Channel.empty()
} else {
    process sortmerna_index {
        label 'low_memory'
        tag "${fasta.baseName}"

        input:
        file(fasta) from sortmerna_fasta

        output:
        val("${fasta.baseName}") into sortmerna_db_name
        file("$fasta") into sortmerna_db_fasta
        file("${fasta.baseName}*") into sortmerna_db

        script:
        """
        indexdb_rna --ref $fasta,${fasta.baseName} -m 3072 -v
        """
    }

    process sortmerna {
        label 'low_memory'
        tag "$name"
        publishDir "${params.outdir}/SortMeRNA", mode: "${params.publish_dir_mode}",
            saveAs: {filename ->
                if (filename.indexOf("_rRNA_report.txt") > 0) "logs/$filename"
                else if (params.save_non_rrna_reads) "reads/$filename"
                else null
            }

        input:
        set val(name), file(reads) from ch_reads_for_ribosomal_removal
        val(db_name) from sortmerna_db_name.collect()
        file(db_fasta) from sortmerna_db_fasta.collect()
        file(db) from sortmerna_db.collect()

        output:
        set val(name), file("*.fq.gz") into reads_ch
        file "*_rRNA_report.txt" into sortmerna_logs


        script:
        //concatenate reference files: ${db_fasta},${db_name}:${db_fasta},${db_name}:...
        def Refs = ''
        for (i=0; i<db_fasta.size(); i++) { Refs+= ":${db_fasta[i]},${db_name[i]}" }
        Refs = Refs.substring(1)

        // One set of reads --> single end
        if (reads[1] == null) {
            """
            gzip -d --force < ${reads} > all-reads.fastq
            sortmerna --ref ${Refs} \
                --reads all-reads.fastq \
                --num_alignments 1 \
                -a ${task.cpus} \
                --fastx \
                --aligned rRNA-reads \
                --other non-rRNA-reads \
                --log -v
            gzip --force < non-rRNA-reads.fastq > ${name}.fq.gz
            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
        } else {
            """
            gzip -d --force < ${reads[0]} > reads-fw.fq
            gzip -d --force < ${reads[1]} > reads-rv.fq
            merge-paired-reads.sh reads-fw.fq reads-rv.fq all-reads.fastq
            sortmerna --ref ${Refs} \
                --reads all-reads.fastq \
                --num_alignments 1 \
                -a ${task.cpus} \
                --fastx --paired_in \
                --aligned rRNA-reads \
                --other non-rRNA-reads \
                --log -v
            unmerge-paired-reads.sh non-rRNA-reads.fastq non-rRNA-reads-fw.fq non-rRNA-reads-rv.fq
            gzip < non-rRNA-reads-fw.fq > ${name}-fw.fq.gz
            gzip < non-rRNA-reads-rv.fq > ${name}-rv.fq.gz
            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
        }
    }
  }



  if (params.reference_proteome_fasta){
    process translate {
      tag "${sample_id}"
      label "low_memory_long"
      publishDir "${params.outdir}/translate/", mode: 'copy'

      input:
      set bloom_id, molecule, file(bloom_filter) from ch_sencha_bloom_filter.collect()
      set sample_id, file(reads) from reads_ch

      output:
      // TODO also extract nucleotide sequence of coding reads and do sourmash compute using only DNA on that?
      set val(sample_id), file("${sample_id}__coding_reads_peptides.fasta") into ch_translated_protein_seqs
      set val(sample_id), file("${sample_id}__coding_reads_nucleotides.fasta") into ch_translatable_nucleotide_seqs
      set val(sample_id), file("${sample_id}__coding_scores.csv") into ch_coding_scores_csv
      set val(sample_id), file("${sample_id}__coding_summary.json") into ch_coding_scores_json

    script:
    """
    sencha translate \\
      --molecule ${molecule} \\
      --coding-nucleotide-fasta ${sample_id}__coding_reads_nucleotides.fasta \\
      --csv ${sample_id}__coding_scores.csv \\
      --json-summary ${sample_id}__coding_summary.json \\
      --jaccard-threshold ${jaccard_threshold} \\
      --peptide-ksize ${peptide_ksize} \\
      --peptides-are-bloom-filter \\
      ${bloom_filter} \\
      ${reads} > ${sample_id}__coding_reads_peptides.fasta
    """
    }

    // Remove empty files
    // it[0] = sample id
    // it[1] = sequence fasta file
    ch_translatable_nucleotide_seqs_nonempty = ch_translatable_nucleotide_seqs.filter{ it[1].size() > 0 }
    ch_translated_protein_seqs.filter{ it[1].size() > 0 }
      .mix ( ch_protein_fastas )
      .set { ch_translated_protein_seqs_nonempty }

  } else {
    // Send reads directly into coding/noncoding
    reads_ch
      .set{ ch_translatable_nucleotide_seqs_nonempty }
  }

  if (params.split_kmer){
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
  } else if (!params.skip_compute) {
    process sourmash_compute_sketch_fastx_nucleotide {
      tag "${sig_id}"
      label "low_memory"
      publishDir "${params.outdir}/sourmash/sketches_nucleotide", mode: "${params.publish_dir_mode}",
          saveAs: {filename ->
              if (filename.indexOf(".csv") > 0) "description/$filename"
              else if (filename.indexOf(".sig") > 0) "sigs/$filename"
              else null
          }

      input:
      val sketch_style from sketch_style_for_nucleotides.collect()
      each ksize from ksizes
      each sketch_value from sketch_values
      val track_abundance
      set sample_id, file(reads) from ch_translatable_nucleotide_seqs_nonempty

      output:
      file(csv) into ch_sourmash_sig_describe_nucleotides
      set val(sketch_id), val("dna"), val(ksize), val(sketch_value), file(sig) into sourmash_sketches_all_nucleotide

      script:
      sketch_style = sketch_style[0]
      // Don't calculate DNA signature if this is protein, to minimize disk,
      // memory and IO requirements in the future
      ksize = ksize
      sketch_id = make_sketch_id("dna", ksize, sketch_value, track_abundance, sketch_style)
      sketch_value_flag = make_sketch_value_flag(sketch_style, sketch_value)
      track_abundance_flag = track_abundance ? '--track-abundance' : ''
      processes = "--processes ${task.cpus}"
      sig_id = "${sample_id}__${sketch_id}"
      sig = "${sig_id}.sig"
      csv = "${sig_id}.csv"
      """
        sourmash compute \\
          ${sketch_value_flag} \\
          --ksizes $ksize \\
          --dna \\
          $processes \\
          $track_abundance_flag \\
          --output ${sig} \\
          $reads
        sourmash sig describe --csv ${csv} ${sig}
      """
    }
    sourmash_sketches_nucleotide = sourmash_sketches_all_nucleotide.filter{ it[4].size() > 0 }
  }
} else {
  sourmash_sketches_nucleotide = Channel.empty()
  ch_protein_fastas
    .set { ch_translated_protein_seqs_nonempty }
}



if ((!have_nucleotide_input) || params.skip_trimming || have_nucleotide_fasta_input) {
  // Only protein input or skip trimming, or fastas which can't be trimmed.
  ch_fastp_results = Channel.from(false)
}

if (!have_nucleotide_input) {
  // Only protein input, can't do sortMeRNA
  sortmerna_logs = Channel.empty()
}


if (!params.skip_compute && (protein_input || params.reference_proteome_fasta)){
  process sourmash_compute_sketch_fastx_peptide {
    tag "${sig_id}"
    label "low_memory"
    publishDir "${params.outdir}/sourmash/sketches_peptide", mode: "${params.publish_dir_mode}",
        saveAs: {filename ->
            if (filename.indexOf(".csv") > 0) "description/$filename"
            else if (filename.indexOf(".sig") > 0) "sigs/$filename"
            else null
        }

    input:
    val sketch_style from sketch_style_for_proteins.collect()
    each ksize from ksizes
    each molecule from peptide_molecules
    each sketch_value from sketch_values
    val track_abundance
    set sample_id, file(reads) from ch_translated_protein_seqs_nonempty

    output:
    file(csv) into ch_sourmash_sig_describe_peptides
    set val(sketch_id), val(molecule), val(ksize), val(sketch_value), file(sig) into sourmash_sketches_all_peptide

    script:
    sketch_style = sketch_style[0]
    molecule = molecule
    ksize = ksize
    sketch_id = make_sketch_id(molecule, ksize, sketch_value, track_abundance, sketch_style)
    sketch_value_flag = make_sketch_value_flag(sketch_style, sketch_value)
    track_abundance_flag = track_abundance ? '--track-abundance' : ''
    processes = "--processes ${task.cpus}"
    sig_id = "${sample_id}__${sketch_id}"
    sig = "${sig_id}.sig"
    csv = "${sig_id}.csv"
    """
      sourmash compute \\
        ${sketch_value_flag} \\
        --ksizes $ksize \\
        --input-is-protein \\
        --$molecule \\
        --name ${sample_id} \\
        --no-dna \\
        $processes \\
        $track_abundance_flag \\
        --output ${sig} \\
        $reads
      sourmash sig describe --csv ${csv} ${sig}
    """
    }
  sourmash_sketches_peptide = sourmash_sketches_all_peptide.filter{ it[4].size() > 0 }
} else {
  sourmash_sketches_peptide = Channel.empty()
}

if (params.split_kmer){
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
  }
// If skip_compute is true, skip compare must be specified as true as well
if (!params.split_kmer && !params.skip_compare && !params.skip_compute) {
  process sourmash_compare_sketches {
    // Combine peptide and nucleotide sketches
    sourmash_sketches = sourmash_sketches_peptide.concat(sourmash_sketches_nucleotide)
    tag "${sketch_id}"
    publishDir "${params.outdir}/sourmash/compare", mode: 'copy'

    input:
    set val(sketch_id), val(molecule), val(ksize), val(sketch_value), file ("sourmash/sketches/*.sig") \
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
          $processes \\
          --traverse-directory .
    """

  }
}


/*
 * STEP 16 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: "${params.publish_dir_mode}"

    when:
    !params.skip_multiqc

    input:
    file multiqc_config from ch_multiqc_config
    file ('fastp/*') from ch_fastp_results.collect().ifEmpty([])
    file ('sortmerna/*') from sortmerna_logs.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $multiqc_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file \\
        -m custom_content \\
        -m fastp \\
        -m sortmerna
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
