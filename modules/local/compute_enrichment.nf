#!/usr/bin/env nextflow

process COMPUTE_ENRICHMENT {
  publishDir "${params.outdir}/${meta.id}/enrichment_results", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(input_bed), 
    path(input_peak)
  )
  path(coding_promoter)

  output:
  tuple(
    val(meta),
    path("*sub10M_enrichment.tsv"),
    path("*promoter_enrichment.tsv")
  )
  path("*enrichment_scores_mqc.tsv"), emit: enrich_score_mqc

  script:
  """
  total=\$(wc -l < ${input_bed})
  num_of_useful_reads=\$(( \$total > 10000000 ? \$total : 10000000 ))

  # Subsample 10M reads or use all reads if fewer
  subsample_size=\$(( \$num_of_useful_reads > 10000000 ? 10000000 : \$num_of_useful_reads))
  shuf ${input_bed} | head -\$subsample_size > temp.open.bed

  # Peak calling with MACS
  macs3 callpeak -t temp.open.bed -g ${params.macs_genome} -q 0.01 -n temp_peak --keep-dup 1000 --nomodel --shift 0 --extsize 150
  peak_length=\$(awk '{s += \$3 - \$2 + 1} END {print s}' temp_peak_peaks.narrowPeak)

  # Count reads in peaks
  nrup=\$(bedtools intersect -iobuf 200M -a temp.open.bed -b temp_peak_peaks.narrowPeak -f 0.5 | wc -l)
  rup=\$(python3 -c "print(\$nrup / 1e7)")

  # Compute enrichment score
  upper=\$(python3 -c "print(\$nrup / \$peak_length + 1e7 / ${params.genome_size})")
  lower=\$(python3 -c "print((\$num_of_useful_reads + 1e7) / (${params.genome_size} - \$peak_length))")
  sub10M_enrich_score=\$(python3 -c "print(f'{\$upper / \$lower:.4f}')")

  {
    echo -e "name\\tnrup\\trup\\tcoverage\\tenrichment"
    echo -e "${meta.id}\\t\$nrup\\t\$rup\\t\$peak_length\\t\$sub10M_enrich_score"
  } > ${meta.id}_sub10M_enrichment.tsv


  # Coding promoter enrichment
  # enrichment = (reads in promoter / promoter length) / (total reads / genome size)

  # Find peaks overlapping coding promoters
  bedtools intersect -iobuf 200M -a ${input_peak} -b ${coding_promoter} -u > promoter_peak.bed

  # Count reads in promoter peaks
  reads_in_promoter=\$(bedtools intersect -iobuf 200M -a ${input_bed} -b promoter_peak.bed -f 0.5 -u | wc -l | awk '{print \$1}')
  promoter_number=\$(bedtools intersect -iobuf 200M -a ${coding_promoter} -b promoter_peak.bed -F 0.5 -u | wc -l | awk '{print \$1}')

  # Assume 2000 bp per promoter
  promoter_length=\$(echo "\$promoter_number * 2000 + 0.001" | bc -l)
  
  # Compute background read density
  denominator=\$(python3 -c "print(f'{\$total / ${params.genome_size}}')")

  # Compute promoter enrichment score
  promoter_enrich_score=\$(python3 -c "print(f'{\$reads_in_promoter / \$promoter_length / \$denominator:.4f}')")
  
  # Write results
  {
    echo -e "name\\ttotal_reads\\tpromoter_number\\treads_in_promoter\\tenrichment_ratio"
    echo -e "${meta.id}\\t\$total\\t\$promoter_number\\t\$reads_in_promoter\\t\$promoter_enrich_score"
  } > ${meta.id}_promoter_enrichment.tsv

  {
    echo -e "sub10M_enrichment_score\\t\$sub10M_enrich_score"
    echo -e "coding_promoter_enrichment_score\\t\$promoter_enrich_score" 
  } > ${meta.id}_enrichment_scores_mqc.tsv
  """
}
