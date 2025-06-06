#!/usr/bin/env nextflow

process COMPUTE_BACKGROUND {
  publishDir "${params.outdir}/${meta.id}/background_results", mode: "copy"
  
  input:
  tuple(
    val(meta),
    path(input_bed),
    path(input_peak)
  )
  path(promoter_file)
  path(nochrM_chrom_size)

  output:
  tuple(
    val(meta),
    path("${meta.id}_background.tsv"),
    path("${meta.id}_dichoto_bg.tsv")
  )
  path("bin_${meta.id}.txt"), optional: true
  path("*promoter_nonpromoter_peaks_mqc.tsv"), emit: pnp_peaks_mqc
  path("*promoter_nonpromoter_reads_mqc.tsv"), emit: pnp_reads_mqc
  path("*background_0_3777_mqc.tsv"), emit: background_0_3777_mqc

  script:
  """
  # signal part
  bedtools intersect -iobuf 200M -a ${input_peak} -b ${promoter_file} -u \
    | awk '{print \$1"\\t"\$2"\\t"\$3"\\t""1""\\t"\$9}' > promoter.narrowPeak
  bedtools intersect -iobuf 200M -a ${input_peak} -b ${promoter_file} -v \
    | awk '{print \$1"\\t"\$2"\\t"\$3"\\t""0""\\t"\$9}' > non_promoter.narrowPeak

  # Count peaks and reads in promoter and non-promoter regions
  npeaks_promoter=\$(wc -l < promoter.narrowPeak)
  nreads_promoter=\$(bedtools intersect -iobuf 200M -a ${input_bed} -b promoter.narrowPeak -u -f 0.50 | wc -l)
  npeaks_non_promoter=\$(wc -l < non_promoter.narrowPeak)
  nreads_non_promoter=\$(bedtools intersect -iobuf 200M -a ${input_bed} -b non_promoter.narrowPeak -u -f 0.50 | wc -l)

  {
    echo -e "num_of_peaks_in_promoter\\t\$npeaks_promoter"
    echo -e "num_of_peaks_in_non_promoter\\t\$npeaks_non_promoter"
  } > ${meta.id}_promoter_nonpromoter_peaks_mqc.tsv
  {
    echo -e "num_of_reads_in_promoter\\t\$nreads_promoter"
    echo -e "num_of_reads_in_non_promoter\\t\$nreads_non_promoter"
  } > ${meta.id}_promoter_nonpromoter_reads_mqc.tsv

  cat promoter.narrowPeak non-promoter.narrowPeak \
    | sort -k5 -n -r > top10k.narrowPeak

  lines=\$(wc -l < top10k.narrowPeak)
  if (( \$lines > 100 )); then
    promoter_bin.py top10k.narrowPeak bin_${meta.id}.txt
  else
    echo "Warning: total peak is fewer than 100, promoter bin step would be skipped. At least 100 peaks are required."
  fi

  # Background noise estimation
  # chr.peak will be generated by random_chr.py
  awk 'length(\$1) > 1 && length(\$1) < 6' OFS="\\t" ${nochrM_chrom_size} > temp.txt
  random_chr.py temp.txt chr.peak

  size=\$(wc -l < ${input_bed})

  awk '{
    x = int((\$2 + \$3) / 2)
    if (x - 100000 < 0) x = 100000
    print \$1, x - 100000, x + 100000, \$4
  }' OFS="\\t" ${input_peak} > temp.txt

  bedtools intersect -iobuf 200M -a chr.peak -b temp.txt -v | shuf - | head -50000 | sort -k1,1V -k2,2n > background
  bedtools intersect -iobuf 200M -a ${input_bed} -b background -u -f 0.5 | sort -k1,1V -k2,2n > temp.txt
  # read.txt will be generated by rpkm_bin
  rpkm_bin.py background temp.txt \$size ${meta.id}_background.tsv

  # Calculate background distribution
  bg_total=\$(wc -l < ${meta.id}_background.tsv)
  bg_half_thres=\$(awk '\$6 <= 0.188 {print \$0}' ${meta.id}_background.tsv | wc -l)
  bg_less=\$(awk '\$6 <= 0.377 {print \$0}' ${meta.id}_background.tsv | wc -l)
  bg_more=\$(awk '\$6 > 0.377 {print \$0}' ${meta.id}_background.tsv | wc -l)

  echo "bg total: \$bg_total"
  echo "3 bg: \$bg_half_thres, \$bg_less, \$bg_more"

  ra_half_thres=\$(echo "scale=2; \$bg_half_thres / \$bg_total" | bc -l)
  ra_less=\$(echo "scale=2; \$bg_less / \$bg_total" | bc -l)
  ra_more=\$(echo "scale=2; \$bg_more / \$bg_total" | bc -l)

  echo -e "\$ra_half_thres\\t\$ra_less\\t\$ra_more" > ${meta.id}_dichoto_bg.tsv

  echo -e "percentage_of_background_RPKM_larger_than_0.3777\\t\$ra_more" > ${meta.id}_background_0_3777_mqc.tsv
  """
}
