#!/usr/bin/env nextflow

process PRESEQ {
  publishDir "${params.outdir}/${meta.id}/preseq_results", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(aligned_reads)
  )

  output:
  path("*yield.result"), emit: preseq_yield

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  preseq lc_extrap -o "${meta.id}${suffix}_yield.result" -B ${aligned_reads}
  """
}
