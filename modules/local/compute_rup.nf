#!/usr/bin/env nextflow

process COMPUTE_RUP {
  publishDir "${params.outdir}/${meta.id}/rup_results", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(input_bed),
    path(input_peak)
  )
  path(nochrM_chrom_size)

  output:
  tuple(
    val(meta),
    path("*insertion_site.bedGraph"),
    path("*insertion_site.bigWig")
  )
  path("*_nrup_mqc.tsv"), emit: nrup_mqc
  path("*_rup_ratio_mqc.tsv"), emit: rup_ratio_mqc

  script:
  """
  total=\$(wc -l < ${input_bed})
  # Number of reads under peaks
  nrup=\$(bedtools intersect -iobuf 200M -a ${input_bed} -b ${input_peak} -f 0.5 -u | wc -l)
  # Fraction of reads in peaks
  rup_ratio=\$(echo "scale=4; \$nrup / \$total" | bc -l)

  echo -e "number_of_reads_under_peaks\\t\$nrup" > ${meta.id}_nrup_mqc.tsv
  echo -e "reads_under_peaks_ratio\\t\$rup_ratio" > ${meta.id}_rup_ratio_mqc.tsv

  # Add insertion site bigwig
  awk '{
    mid = int((\$2 + \$3) / 2);
    start = (\$6 == "+") ? mid : mid - 1;
    end = (\$6 == "+") ? mid + 1 : mid;
    print \$1, start, end, 1;
  }' OFS="\\t" ${input_bed} \
    | sort -k1,1 -k2,2n \
    | uniq -c \
    | awk '{print \$2, \$3, \$4, \$1}' OFS="\\t" \
    > ${meta.id}_insertion_site.bedGraph

  bedGraphToBigWig ${meta.id}_insertion_site.bedGraph ${nochrM_chrom_size} ${meta.id}_insertion_site.bigWig
  """
}
 
