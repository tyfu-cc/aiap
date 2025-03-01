#!/usr/bin/env nextflow

process S1_1_TRIM {
  // container "${projectDir}/images/test-s.sif"

  publishDir "${params.outdir}-${meta.id}/1_1_trim", mode: "copy"

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
    path("step1_1_trimmed_*.fastq"), 
  )
  path("step1_1*.trimlog"), emit: trimlog
  path("data_collection_s1_1*.json"), emit: json

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  if (meta.paired)
    """
    cutadapt -a ${adapter_1} -A ${adapter_2} \
      --quality-cutoff=15,10 --minimum-length=36 \
      -o step1_1_trimmed_1${suffix}.fastq -p step1_1_trimmed_2${suffix}.fastq \
      ${raw_reads[0]} ${raw_reads[1]} \
      > step1_1${suffix}.trimlog

    total_reads=\$(awk '/Total read pairs processed:/ {gsub(",", ""); print \$NF}' step1_1${suffix}.trimlog)
    written_reads=\$(awk '/Pairs written/ {gsub(",", ""); print \$(NF-1)}' step1_1${suffix}.trimlog)

    jq --argjson total \$total_reads --argjson written_reads \$written_reads \
      '.total = \$total | .written_reads = \$written_reads' -n > data_collection_s1_1${suffix}.json
    """
  else
    """
    cutadapt -a ${adapter_1} \
      --quality-cutoff=15,10 --minimum-length=36 \
      -o step1_1_trimmed_1${suffix}.fastq \
      ${raw_reads} \
      > step1_1${suffix}.trimlog

    raw_reads=\$(awk '/Total reads processed:/ {gsub(",", ""); print \$NF}' step1_1${suffix}.trimlog
    written_reads=\$(awk '/Reads written/ {gsub(",", ""); print \$(NF-1)}' step1_1${suffix}.trimlog

    jq --arg total \$total_reads --arg written_reads \$written_reads \
      '.total = \$total | .written_reads = \$written_reads' -n > data_collection_s1_1${suffix}.json
    """ 
}
