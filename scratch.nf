housekeeping_protein_fasta = false
housekeeping_rna_fasta = true

ch_refseq_moltype_to_fasta = Channel.from(["protein", housekeeping_protein_fasta], ["rna", housekeeping_rna_fasta])
ch_refseq_moltype_to_fasta
    // filter if the second item, the fasta is false
    .filter{ !it[1] }
    // Take only the first item, the molecule type
    .map{ it[0] }
    .println()
