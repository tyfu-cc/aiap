#!/usr/bin/env nextflow

process S4_2_ENRICH {
  publishDir "${params.outdir}-${meta.id}/4_2_enrich", mode: "copy"

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
    path("step4_2_sub10M_enrichment.result"),
    path("step4_2_enrichment_ratio_in_promoter.result")
  )
  path("data_collection_s4_2.json"), emit: json

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
  rupn=\$(bedtools intersect -iobuf 200M -a temp.open.bed -b temp_peak_peaks.narrowPeak -f 0.5 | wc -l)

  # Compute enrichment score
  upper=\$(python3 -c "print(\$rupn / \$peak_length + 1e7 / ${params.genome_size})")
  lower=\$(python3 -c "print((\$num_of_useful_reads + 1e7) / (${params.genome_size} - \$peak_length))")
  enrichment=\$(python3 -c "print(\$upper / \$lower)")
  
  rup=\$(python3 -c "print(\$rupn / 1e7)")

  {
    echo -e "name\\trupn\\trup\\tcoverage\\tenrichment"
    echo -e "${meta.id}\\t\$rupn\\t\$rup\\t\$peak_length\\t\$enrichment"
  } > step4_2_sub10M_enrichment.result


  # 4.2.2
  # Coding promoter enrichment
  # coding enrichment = (reads in promoter / promoter length)  /  (total reads / genome size)

  # Compute background read density
  denominator=\$(echo "scale=10; \$total / ${params.genome_size}" | bc -l)

  # Find peaks overlapping coding promoters
  bedtools intersect -iobuf 200M -a ${input_peak} -b ${coding_promoter} -u > promoter_peak.bed

  # Count reads in promoter peaks
  reads_in_promoter=\$(bedtools intersect -iobuf 200M -a ${input_bed} -b promoter_peak.bed -f 0.5 -u | wc -l | awk '{print \$1}')
  promoter_number=\$(bedtools intersect -iobuf 200M -a ${coding_promoter} -b promoter_peak.bed -F 0.5 -u | wc -l | awk '{print \$1}')

  # Assume 2000 bp per promoter
  promoter_length=\$(echo "\$promoter_number * 2000 + 0.001" | bc -l)

  # Compute promoter enrichment score
  enrichment_ratio=\$(echo "scale=3; \$reads_in_promoter / \$promoter_length / \$denominator" | bc -l)
  
  # Write results
  {
    echo -e "name\\ttotal_reads\\tpromoter_number\\treads_in_promoter\\tenrichment_ratio"
    echo -e "${meta.id}\\t\$total\\t\$promoter_number\\t\$reads_in_promoter\\t\$enrichment_ratio"
  } > step4_2_enrichment_ratio_in_promoter.result

  jq --argjson enrichment \$enrichment --argjson enrichment_ratio \$enrichment_ratio \
    '.sub10M_enrichment = \$enrichment | .coding_enrichment = \$enrichment_ratio' -n > data_collection_s4_2.json
  # plots.R s4_2
  """
}
