#!/usr/bin/env nextflow

process S4_5_VISUAL {
  publishDir "${params.outdir}-${meta.id}/plots", mode: "copy"

  input:
  val(meta)

  output:
  path("*.png"), optional: true 

  script:
  """
  plots.R
  """
}
