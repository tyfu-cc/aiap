#!/usr/bin/env nextflow

process S2_2_COMPLEXITY {
  publishDir "${params.outdir}-${meta.id}/2_2_complexity", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(aligned_reads)
  )

  output:
  path("step2_3_yield*.result")

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  preseq lc_extrap -o "step2_3_yield${suffix}.result" -B ${aligned_reads}
  """
}
