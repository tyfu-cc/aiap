#!/usr/bin/env nextflow

process METHYLQA {
  publishDir "${params.outdir}/${meta.id}/methylQA_results", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(aligned_reads)
  )
  path(chrom_size)
  path(trimlog)
  path(dedup_pct)

  output:
  tuple(
    val(meta),
    path("${meta.id}*_methylQA*")
  )
  path("refined_chrom_size.txt"), emit: refined
  path("nochrM_chrom_size.txt"), emit: nochrM
  path("*mapping_status_mqc.tsv"), emit: mapping_status_mqc
  path("*dup_pct_mqc.tsv"), emit: dup_pct_mqc
  path("*insert_size_dist.tsv"), emit: insert_size
  path("*insert_size_dist_mqc.tsv"), emit: insert_size_mqc

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  awk 'length(\$1) > 1 && length(\$1) < 6' OFS="\\t" ${chrom_size} > refined_chrom_size.txt
  grep -v '^chrM' refined_chrom_size.txt > nochrM_chrom_size.txt

  methylQA atac -X ${params.methylqa_cutoff} -o ${meta.id}${suffix}_methylQA \
    nochrM_chrom_size.txt ${aligned_reads}

  # methylQA density -o ${meta.id}_methylQA_density${suffix} \
  #   nochrM_chrom_size.txt ${aligned_reads}
 
  # Mapping Status
  if [[ ${meta.paired} ]]; then
    raw_reads=\$(awk '/Total read pairs processed:/ {gsub(",", ""); print \$NF}' ${trimlog})
  else
    raw_reads=\$(awk '/Total reads processed:/ {gsub(",", ""); print \$NF}' ${trimlog})
  fi
  mappable=\$(grep '^mappable reads' ${meta.id}${suffix}_methylQA.report | awk '{print \$NF}')
  uniq_mapped=\$(grep '^uniquely mapped reads' ${meta.id}${suffix}_methylQA.report | awk '{print \$NF}')
  non_redundant_uniq_mapped=\$(grep '^non-redundant uniquely mapped reads' ${meta.id}${suffix}_methylQA.report | awk '{print \$NF}')
  useful_single_ends=\$(wc -l ${meta.id}${suffix}_methylQA.open.bed | awk '{print \$1}')

  before_aln_dedup=\$(awk 'NR > 1 {s += \$2} END {print s * 0.01 / (NR-1)}' ${dedup_pct})
  before_aln_dup=\$(python3 -c "print(f'{1 - \$before_aln_dedup:.4f}')")
  nodup_ratio=\$(echo "scale=4; \$non_redundant_uniq_mapped / \$uniq_mapped" | bc -l)
  after_aln_dup=\$(echo "scale=4; 1 - \$nodup_ratio" | bc -l)
  after_aln_dup=\$(python3 -c "print(f'{1 - \$non_redundant_uniq_mapped / \$uniq_mapped:.4f}')")
  {
    echo -e "before_alignment_duplicates_percentage\\t\$before_aln_dup"
    echo -e "after_alignment_duplicates_percentage\\t\$after_aln_dup"
  } > ${meta.id}${suffix}_dup_pct_mqc.tsv

  # Process insertion distribution
  sort -n ${meta.id}${suffix}_methylQA.insertdistro \
    | uniq -c \
    | awk 'BEGIN {OFS="\t"} {print \$2, \$1}' > ${meta.id}${suffix}_insert_size_dist.tsv
  density_estimation.R ${meta.id}${suffix}_insert_size_dist.tsv ${meta.id}${suffix}_insert_size_dist_mqc.tsv

  {
    echo -e "raw_reads\\t\$raw_reads"
    echo -e "mapped_reads\\t\$mappable"
    echo -e "uniquely_mapped_reads\\t\$uniq_mapped"
    echo -e "non_redundant_uniquely_mapped_reads\\t\$non_redundant_uniq_mapped"
    echo -e "useful_single_ends\\t\$useful_single_ends"
  } > ${meta.id}${suffix}_mapping_status_mqc.tsv
  """
}
