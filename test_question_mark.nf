ksizes = params.ksizes.splitCsv() ?: [21, 27, 33, 51]

Channel
  .from(ksizes)
  .println()
