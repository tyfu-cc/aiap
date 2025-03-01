#!/usr/bin/env nextflow

process S3_1_QA {
  publishDir "${params.outdir}-${meta.id}/3_1_methylQA", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(aligned_reads)
  )
  path(chrom_size)
  path(data_collection_s1_1)

  output:
  tuple(
    val(meta),
    path("step3_1_methylQA*"),
    path("step3_1_useful_reads*.result"),
    path("step3_1_insertion_distri*.result")
  )
  path("refined_chrom_size.txt"), emit: refined
  path("nochrM_chrom_size.txt"), emit: nochrM
  path("*.png"), emit: plot, optional: true
  path("data_collection_s3_1*.json"), emit: json

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  awk 'length(\$1) > 1 && length(\$1) < 6' OFS="\\t" ${chrom_size} > refined_chrom_size.txt
  grep -v '^chrM' refined_chrom_size.txt > nochrM_chrom_size.txt

  methylQA atac -X ${params.methylqa_cutoff} -o step3_1_methylQA${suffix} \
    nochrM_chrom_size.txt ${aligned_reads}
    
  # Mapping Status
  map_mapped=\$(grep '^mappable reads' step3_1_methylQA${suffix}.report | awk '{print \$NF}')
  map_uniq=\$(grep '^uniquely mapped reads' step3_1_methylQA${suffix}.report | awk '{print \$NF}')
  map_effect=\$(grep '^non-redundant uniquely mapped reads' step3_1_methylQA${suffix}.report | awk '{print \$NF}')
  raw_reads=\$(jq '.total' ${data_collection_s1_1})
  mapped_ratio=\$(echo "scale=2; \$map_mapped / \$raw_reads" | bc -l)
  effect_ratio=\$(echo "scale=2; \$map_effect / \$raw_reads" | bc -l)
  nodup_ratio=\$(echo "scale=3; \$map_effect / \$map_uniq" | bc -l)
  after_dup=\$(echo "scale=3; 1 - \$nodup_ratio" | bc -l)

  jq --argjson mapped \$map_mapped --argjson mapped_ratio \$mapped_ratio \
    --argjson uniq_mapped \$map_uniq --argjson non_redundant_uniq_mapped \$map_effect \
    --argjson effect_ratio \$effect_ratio --argjson after_align_dup \$after_dup \
    '.mapped = \$mapped | .mapped_ratio = \$mapped_ratio \
    | .uniq_mapped = \$uniq_mapped | .non_redundant_uniq_mapped = \$non_redundant_uniq_mapped \
    | .effect_ratio = \$effect_ratio | .after_align_dup = \$after_align_dup' -n > data_collection_s3_1${suffix}.json
  
  useful=\$(grep '^non-redundant uniquely mapped reads' step3_1_methylQA${suffix}.report | awk '{print \$NF}')
  single_end=\$(wc -l step3_1_methylQA${suffix}.open.bed | awk '{print \$1}')
  uf_ratio=\$(echo "scale=3; \$useful / \$raw_reads" | bc -l)

  {
    echo -e "file\\ttotal\\tuseful\\tuseful_ratio\\tsingle_end"
    echo -e "${meta.id}\\t\$raw_reads\\t\$useful\\t\$uf_ratio\\t\$single_end"
  } > step3_1_useful_reads${suffix}.result

  # Process insertion distribution
  sort -n step3_1_methylQA${suffix}.insertdistro \
    | uniq -c \
    | awk '{print \$2, \$1}' > step3_1_insertion_distri${suffix}.result

  # plots.R s3_1
  """
}
