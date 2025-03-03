#!/usr/bin/env nextflow

process FASTQC {
  publishDir "${params.outdir}/${meta.id}/fastqc_results", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(trimmed)
  )

  output:
  path("*.zip"), emit: zip
  path("*.html"), emit: html
  path("*dedup_pct_mqc.tsv"), emit: dedup_pct_mqc
  path("*duplication_summary*.result")
  path("*fastqc_summary*.result")

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  fastqc -t ${params.threads} ${trimmed} -o .
  
  for zip in *.zip; do
    unzip -o \$zip
  done

  # Collect fastqc data for output
  echo -e "filename\\tdeduplication_percentage\\tmarker" > ${meta.id}${suffix}_dedup_pct_mqc.tsv
  end=1
  for dir in *fastqc/; do
    [[ -d \$dir ]] || continue

    # Append deduplication percentage
    out_value=\$(grep 'Total Deduplicated Percentage' \$dir/fastqc_data.txt | awk '{print \$4}')
    echo -e "${meta.id}_\$end\\t\$out_value\\t${params.marker}" >> ${meta.id}${suffix}_dedup_pct_mqc.tsv

    # Generate duplication summary
    echo -e "item\\t${meta.id}_\$end\\t${meta.id}_\$end" \
      > "${meta.id}${suffix}_duplication_summary_\$end.result"
    grep 'Sequence Duplication Levels' -A 15 \$dir/fastqc_data.txt \
      >> "${meta.id}${suffix}_duplication_summary_\$end.result"

    # Generate FastQC summary
    {
      echo -e "${meta.id}_\$end\\tfastqc_test"
      awk -F "\t" '{print \$1, \$2}' OFS="\\t" "\$dir/summary.txt"
    } > "${meta.id}${suffix}_fastqc_summary_\$end.result"

    (( end++ ))
  done

  # Get PE data R1 R2 deduplication difference percentage
  if [[ ${meta.paired} ]]; then
    per1=\$(tail -n 2 "${meta.id}${suffix}_dedup_pct.tsv" | awk '{print \$2}' | sed -n '1p')
    per2=\$(tail -n 2 "${meta.id}${suffix}_dedup_pct.tsv" | awk '{print \$2}' | sed -n '2p')
    dif=\$(echo "scale=2; (\$per1 - \$per2) * 200 / (\$per1 + \$per2)" | bc -l)
  else
    dif=0
  fi
  """
}
