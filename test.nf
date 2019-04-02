

params.sra = "SRP016501"
params.samples = "testing/samples.csv"

// Samples from SRA
sra_ch = Channel.create()
// R1, R2 pairs from a samples.csv file
samples_ch = Channel.create()
// Extract R1, R2 pairs from a directory
directories_ch = Channel.create()

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
    .from(params.directories?.toString()?.tokenize(','))
    .map(it + "*{R1,R2}*.fastq.gz")
    .fromFilePairs()
}

// sra_ch.concat(samples_ch, directories_ch).subscribe { println it }

sra_ch.concat(samples_ch, directories_ch)
  .set { reads_ch }
println reads_ch
