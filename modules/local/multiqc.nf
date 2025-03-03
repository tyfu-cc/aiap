#!/usr/bin/env nextflow

process MULTIQC {
  publishDir "${params.outdir}/multiqc", mode: "copy"
  
  input:
  path(multiqc_config)
  path("*")

  output:
  path("multiqc_report.html"), emit: report
  path("multiqc_report_data"), emit: data

  script:
  """
  multiqc . -n multiqc_report.html --config ${multiqc_config}
  """
}
