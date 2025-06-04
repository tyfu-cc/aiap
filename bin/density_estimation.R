#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]

if (file.info(infile)$size == 0) {
  file.create(outfile)
} else {
  if (length(expanded_counts) < 2) {
    file.create(outfile)
  } else {
    dens <- density(expanded_counts, bw = "SJ")
    write.table(
      data.frame(dens$x, dens$y),
      file = outfile,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
}
