#!/usr/bin/env nextflow

process S4_6_SUMMARIZE {
  // publishDir "${params.outdir}-${meta.id}", mode: "copy"
  publishDir "${params.outdir}-${meta.id}", mode: "copy"

  input:
  val(meta)
  path(trimlog)
  path(methylqa_results)
  path(files)
  path(s1_1_json)
  path(s3_1_json)
  path(s4_1_json)
  path(s4_2_json)
  path(s4_4_json)

  output:
  path("*.json"), optional: true
  path("*.html")
  
  script:
  """
  jq -s "add" ${s1_1_json} ${s3_1_json} ${s4_1_json} ${s4_2_json} ${s4_4_json} > data_collection.json
  multiqc .
  """
}
