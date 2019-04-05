params.samples = "testing/samples.csv"
params.sra = "SRP016501"
params.directories = "s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/*{R1,R2}*.fastq.gz"

// Samples from SRA
sra_ch = Channel.empty()
// R1, R2 pairs from a samples.csv file
samples_ch = Channel.empty()
// Extract R1, R2 pairs from a directory
directories_ch = Channel.empty()

// Provided SRA ids
if (params.sra){
  sra_ch = Channel
      .fromSRA( params.sra?.toString()?.tokenize(',') )
}
// Provided a samples.csv file
if (params.samples){
  samples_ch = Channel
   .fromPath(params.samples)
   .splitCsv(header:true)
   .map{ row -> tuple(row.sample_id, tuple(row.read1, row.read2))}
}
// Provided s3 or local directories
if (params.directories){
  directories_ch = Channel
    .fromFilePairs(params.directories?.toString()?.tokenize(','))
}

sra_ch.concat(samples_ch, directories_ch)
 .set{ reads_ch }

println reads_ch
