#!/usr/bin/env nextflow

process S3_2_NORMALIZE {
  publishDir "${params.outdir}-${meta.id}/3_2_normalize", mode: "copy"

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
    path("step3_2_rmbl*.bigWig"),
    path("step3_2_normalized_per_10M*.open.bedGraph"),
    path("step3_2_normalized_per_10M*.bigWig")
  )

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  # Add a new bigwig file without black list
  bedtools intersect -iobuf 200M -a ${input_bedGraph} -b ${black_list} -v > rmbl.bedGraph
  bedGraphToBigWig rmbl.bedGraph ${nochrM_chrom_size} step3_2_rmbl${suffix}.bigWig

  # Normalization
  norm=\$(grep 'non-redundant' ${input_report} | awk '{print \$6}')
  factor=\$(echo "scale=3; \$norm / 10000000" | bc -l)
  awk -v factor=\$factor '{print \$1, \$2, \$3, \$4 / factor}' OFS='\\t' rmbl.bedGraph \
    > step3_2_normalized_per_10M${suffix}.open.bedGraph
  bedGraphToBigWig step3_2_normalized_per_10M${suffix}.open.bedGraph ${nochrM_chrom_size} \
    step3_2_normalized_per_10M${suffix}.bigWig
  """
} 
