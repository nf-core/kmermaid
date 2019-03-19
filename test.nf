
ksizes = Channel.from([15, 21, 27, 33, 51])
molecules = Channel.from(['protein', 'dna'])
log2_sketch_sizes = Channel.from([10, 12, 14, 16]

parameters = molecules
	.combine(ksizes)
	.combine(log2_sketch_sizes)
	.println()
