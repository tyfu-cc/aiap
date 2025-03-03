#!/usr/bin/env nextflow

process MACS_PEAK_CALLING {
  publishDir "${params.outdir}/${meta.id}/peak_calling_results", mode: "copy"
  
  input:
  tuple(
    val(meta),
    path(input_bed)
  )
  path(black_list)

  output:
  tuple(
    val(meta), 
    path("*_bl_removed_case.open.bed"),
    path("*_peaks.narrowPeak"),
    path("*_peaks.xls"),
    path("*_summits.bed"),
  )
  path("*_peak_len_dist.tsv"), emit: peak_len_dist
  path("*_peak_len_dist_mqc.tsv"), emit: peak_len_dist_mqc

  script:
  def case_bed = input_bed.find { it.name.contains("case") }
  def ctrl_bed = input_bed.find { it.name.contains("ctrl") }
  if (meta.cc)
    """
    awk '{if (\$2 > \$3) sub(\$2, 0); print}' OFS="\\t" ${case_bed} \
      | bedtools intersect -iobuf 200M -a - -b ${black_list} -v \
      > ${meta.id}_bl_removed_case.open.bed
    awk '{if (\$2 > \$3) sub(\$2, 0); print}' OFS="\\t" ${ctrl_bed} \
      | bedtools intersect -iobuf 200M -a - -b ${black_list} -v \
      > ${meta.id}_bl_removed_ctrl.open.bed

    # Peak calling using MACS3
    # q-value cutoff: 0.01
    macs3 callpeak -t ${meta.id}_bl_removed_case.open.bed \
      -c ${meta.id}_bl_removed_ctrl.open.bed \
      -g ${params.macs_genome} \
      -q 0.01 \
      -n ${meta.id} \
      --keep-dup 1000 --nomodel --shift 0 --extsize 150

    # Peak length distribution
    # The length of each peak is calculated
    # Counts of unique peak lengths are stored
    awk '{print \$3 - \$2 + 1}' ${meta.id}_peaks.narrowPeak \
      | sort -n \
      | uniq -c \
      | awk 'BEGIN {OFS="\t"} {print \$2, \$1}' \
      > ${meta.id}_peak_len_dist.tsv

    density_estimation.R ${meta.id}_peak_len_dist.tsv ${meta.id}_peak_len_dist_mqc.tsv
    """
  else
    """
    awk '{if (\$2 > \$3) sub(\$2, 0); print}' OFS="\\t" ${input_bed} \
      | bedtools intersect -iobuf 200M -a - -b ${black_list} -v \
      > ${meta.id}_bl_removed_case.open.bed

    # Peak calling using MACS3
    # q-value cutoff: 0.01
    macs3 callpeak -t ${meta.id}_bl_removed_case.open.bed \
      -g ${params.macs_genome} \
      -q 0.01 \
      -n ${meta.id} \
      --keep-dup 1000 --nomodel --shift 0 --extsize 150

    # Peak length distribution
    # The length of each peak is calculated
    # Counts of unique peak lengths are stored
    awk '{print \$3 - \$2 + 1}' ${meta.id}_peaks.narrowPeak \
      | sort -n \
      | uniq -c \
      | awk 'BEGIN {OFS="\t"} {print \$2, \$1}' \
      > ${meta.id}_peak_len_dist.tsv

    density_estimation.R ${meta.id}_peak_len_dist.tsv ${meta.id}_peak_len_dist_mqc.tsv
    """
}
