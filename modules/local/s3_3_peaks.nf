#!/usr/bin/env nextflow

process S3_3_PEAKS {
  publishDir "${params.outdir}-${meta.id}/3_3_peaks", mode: "copy"
  
  input:
  tuple(
    val(meta),
    path(input_bed)
  )
  path(black_list)

  output:
  tuple(
    val(meta), 
    path("step3_3_rmbl_case.open.bed"),
    path("step3_3_peakcall_peaks.narrowPeak"),
    path("step3_3_peakcall_peaks.xls"),
    path("step3_3_peakcall_summits.bed"),
    path("step3_3_peak_length_distri.result")
  )
  path("*.png"), emit: plot, optional: true

  script:
  def case_bed = input_bed.find { it.name.contains("case") }
  def ctrl_bed = input_bed.find { it.name.contains("ctrl") }
  if (meta.cc)
    """
    awk '{if (\$2 > \$3) sub(\$2, 0); print}' OFS="\\t" ${case_bed} \
      | bedtools intersect -iobuf 200M -a - -b ${black_list} -v \
      > step3_3_rmbl_case.open.bed
    awk '{if (\$2 > \$3) sub(\$2, 0); print}' OFS="\\t" ${ctrl_bed} \
      | bedtools intersect -iobuf 200M -a - -b ${black_list} -v \
      > step3_3_rmbl_ctrl.open.bed

    # Peak calling using MACS3
    # q-value cutoff: 0.01
    macs3 callpeak -t step3_3_rmbl_case.open.bed \
      -c step3_3_rmbl_ctrl.open.bed \
      -g ${params.macs_genome} \
      -q 0.01 \
      -n step3_3_peakcall \
      --keep-dup 1000 --nomodel --shift 0 --extsize 150

    # Peak length distribution
    # The length of each peak is calculated
    # Counts of unique peak lengths are stored
    awk '{print \$3 - \$2 + 1}' step3_3_peakcall_peaks.narrowPeak \
      | sort -n \
      | uniq -c \
      | awk '{print \$2, \$1}' \
      > step3_3_peak_length_distri.result

    # plots.R s3_3
    """
  else
    """
    awk '{if (\$2 > \$3) sub(\$2, 0); print}' OFS="\\t" ${input_bed} \
      | bedtools intersect -iobuf 200M -a - -b ${black_list} -v \
      > step3_3_rmbl_case.open.bed

    # Peak calling using MACS3
    # q-value cutoff: 0.01
    macs3 callpeak -t step3_3_rmbl_case.open.bed \
      -g ${params.macs_genome} \
      -q 0.01 \
      -n step3_3_peakcall \
      --keep-dup 1000 --nomodel --shift 0 --extsize 150

    # Peak length distribution
    # The length of each peak is calculated
    # Counts of unique peak lengths are stored
    awk '{print \$3 - \$2 + 1}' step3_3_peakcall_peaks.narrowPeak \
      | sort -n \
      | uniq -c \
      | awk '{print \$2, \$1}' \
      > step3_3_peak_length_distri.result

    # plots.R s3_3 step3_3_peak_length_distri.result
    """
}


