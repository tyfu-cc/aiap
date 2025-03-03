#!/usr/bin/env nextflow

process CUTADAPT_TRIM {
  publishDir "${params.outdir}/${meta.id}/trimmed_results", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(raw_reads)
  )
  tuple(
    val(adapter_1),
    val(adapter_2)
  )

  output:
  tuple(
    val(meta), 
    path("*trimmed*.fastq"), 
  )
  path("${meta.id}*.trimlog"), emit: trimlog

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  if (meta.paired)
    """
    cutadapt -a ${adapter_1} -A ${adapter_2} \
      --quality-cutoff=15,10 --minimum-length=36 \
      -o ${meta.id}${suffix}_trimmed_1.fastq -p ${meta.id}${suffix}_trimmed_2.fastq \
      ${raw_reads[0]} ${raw_reads[1]} \
      > ${meta.id}${suffix}.trimlog
    """
  else
    """
    cutadapt -a ${adapter_1} \
      --quality-cutoff=15,10 --minimum-length=36 \
      -o ${meta.id}${suffix}_trimmed_1.fastq \
      ${raw_reads} \
      > ${meta.id}${suffix}.trimlog
    """ 
}
