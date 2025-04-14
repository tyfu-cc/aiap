#!/usr/bin/env nextflow

process NORMALIZE_BEDGRAPH {
  publishDir "${params.outdir}/${meta.id}/normalized_results", mode: "copy"

  input:
  tuple(
    val(meta),
    path(input_report),
    path(input_bedGraph)
  )
  path(nochrM_chrom_size)
  path(black_list)

  output:
  tuple(
    val(meta),
    path("${meta.id}_bl_removed*.bigWig"),
    path("${meta.id}_normalized_per_10M*.open.bedGraph"),
    path("${meta.id}_normalized_per_10M*.bigWig")
  )

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  # Add a new bigwig file without black list
  bedtools intersect -iobuf 200M -a ${input_bedGraph} -b ${black_list} -v > bl_removed.bedGraph
  bedGraphToBigWig bl_removed.bedGraph ${nochrM_chrom_size} ${meta.id}_bl_removed${suffix}.bigWig

  # Normalization
  norm=\$(grep 'non-redundant' ${input_report} | awk '{print \$6}')
  factor=\$(echo "scale=3; \$norm / 10000000" | bc -l)
  awk -v factor=\$factor '{print \$1, \$2, \$3, \$4 / factor}' OFS='\\t' bl_removed.bedGraph \
    > ${meta.id}_normalized_per_10M${suffix}.open.bedGraph
  bedGraphToBigWig ${meta.id}_normalized_per_10M${suffix}.open.bedGraph ${nochrM_chrom_size} \
    ${meta.id}_normalized_per_10M${suffix}.bigWig
  """
} 
