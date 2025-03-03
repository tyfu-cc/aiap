#!/usr/bin/env nextflow

process COMPUTE_SATURATION {
  publishDir "${params.outdir}/${meta.id}/saturation_analysis_results", mode: "copy"

  input:
   tuple(
    val(meta), 
    path(input_bed),
    path(input_peak)
  )

  output:
  path("*saturation.tsv")
  path("*saturation_mqc.tsv"), emit: saturation_mqc

  script:
  """
  # Subsample from the BED file obtained from peak calling on original sample
  total=\$(wc -l < ${input_bed})
  for number in 5 10 20 30 40 50 60 70 80 90; do
    sample_ratio=\$(( \$total * \$number / 100 ))
    shuf ${input_bed} | head -\$sample_ratio > trimmed_bl_removed_sample\$number.open.bed
  done

  # Peak calling for the subsamples
  for file in trimmed_bl_removed_sample*.open.bed; do
    macs3 callpeak -t \$file \
      -g ${params.macs_genome} \
      -q 0.01 \
      -n peakcall_\$file \
      --keep-dup 1000 --nomodel --shift 0 --extsize 150
  done
  
  # Summarise the results
  echo -e "5\\n\$(seq 10 10 100)" > saturation_points.txt
  for file in *open.bed; do
    read_num=\$(wc -l < \$file)
    echo "\$(bc -l <<< "scale=2; \$read_num / 1000000")"
  done | sort -n > saturation_reads.txt

  total_region=\$(awk '{s += \$3 - \$2 +1} END {print s}' ${input_peak})

  # Generate saturation peaks and ratios data
  for number in 5 10 20 30 40 50 60 70 80 90; do
    file=peakcall_trimmed_bl_removed_sample\$number.open.bed_peaks.narrowPeak
    npeak=\$(wc -l < \$file)
    echo "\$npeak" >> saturation_peaks.txt

    peak_region=\$(bedtools intersect -iobuf 200M -a \$file -b ${input_peak} \
      | awk '{s += \$3 - \$2 + 1} END {print s}')
    if [[ -z \$peak_region ]]; then
      peak_region=0
    fi
    echo "\$(bc -l <<< "scale=4; \$peak_region / \$total_region")" >> saturation_ratios.txt
  done

  echo \$(wc -l < ${input_peak}) >> saturation_peaks.txt
  echo 1 >> saturation_ratios.txt
  paste saturation_points.txt saturation_reads.txt saturation_peaks.txt saturation_ratios.txt > temp.txt

  {
    echo -e "ptc\\tread\\tpeak\\tratio\\tmarker"
    awk -v marker=${params.marker} '{print \$0, marker}' OFS="\t" temp.txt
  } > ${meta.id}_saturation.tsv

  tail -n +2 ${meta.id}_saturation.tsv | cut -f1,4 > ${meta.id}_saturation_mqc.tsv
  """
}

