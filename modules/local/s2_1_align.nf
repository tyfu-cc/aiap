#!/usr/bin/env nextflow

process S2_1_ALIGN {
  publishDir "${params.outdir}-${meta.id}/2_1_align", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(trimmed)
  )
  path(ref_files)

  output:
  tuple(
    val(meta), 
    path("step2_1_aligned*.bam")
  )

  script:
  def ref_fa = ref_files.find{ it.name.endsWith(".fa") } ?: ref_files.find{ it.name.endsWith(".sa") }.name[0..-4]
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  def bwa_cmd = meta.paired ?
    "bwa mem -t ${params.threads} -v 1 ${ref_fa} ${trimmed[0]} ${trimmed[1]}" : 
    "bwa mem -t ${params.threads} -v 1 ${ref_fa} ${trimmed}"
  """
  ${bwa_cmd} \
    | samtools view -bS -@ ${params.threads} - \
    | samtools sort - -o step2_1_aligned${suffix}.bam -t temp_aln -@ ${params.threads}
  """ 
}
