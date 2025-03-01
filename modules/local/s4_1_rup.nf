#!/usr/bin/env nextflow

process S4_1_RUP {
  publishDir "${params.outdir}-${meta.id}/4_1_rup", mode: "copy"

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
    path("step4_1_insertion_site.bedGraph"),
    path("step4_1_insertion_site.bigWig")
  )
  path("data_collection_s4_1.json"), emit: json

  script:
  """
  total=\$(wc -l < ${input_bed})
  # Number of reads under peaks
  nrup=\$(bedtools intersect -iobuf 200M -a ${input_bed} -b ${input_peak} -f 0.5 -u | wc -l)
  # Fraction of reads in peaks
  ratio=\$(echo "scale=4; \$nrup / \$total" | bc -l)

  jq --argjson nrup \$nrup --argjson rup_ratio \$ratio \
    '.number_of_reads_under_peak = \$nrup | .rup_ratio = \$rup_ratio' -n > data_collection_s4_1.json

  # 4.1.2
  # Add insertion site bigwig
  # echo "[AIAP Step 4.1.2 \$(date +"%H:%M:%S")] START: Add insertion site bigwig" 2>&1 | tee -a \$LOG
  awk '{
    mid = int((\$2 + \$3) / 2);
    start = (\$6 == "+") ? mid : mid - 1;
    end = (\$6 == "+") ? mid + 1 : mid;
    print \$1, start, end, 1;
  }' OFS="\\t" ${input_bed} \
    | sort -k1,1 -k2,2n \
    | uniq -c \
    | awk '{print \$2, \$3, \$4, \$1}' OFS="\\t" \
    > step4_1_insertion_site.bedGraph

  bedGraphToBigWig step4_1_insertion_site.bedGraph ${nochrM_chrom_size} step4_1_insertion_site.bigWig
  """
}
 
