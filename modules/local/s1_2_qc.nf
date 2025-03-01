#!/usr/bin/env nextflow

process S1_2_QC {
  publishDir "${params.outdir}-${meta.id}/1_2_qc", mode: "copy"

  input:
  tuple(
    val(meta), 
    path(trimmed)
  )

  output:
  path("*.zip"), emit: zip
  path("*.html"), emit: html, optional: true
  path("step1_2_1_dedup_percentage*.result")
  path("step1_2_1_duplication_summary_${meta.id}*.result")
  path("step1_2_1_fastqc_summary_${meta.id}*.result")
  path("data_collection_s1_2*.json"), emit: json

  script:
  def suffix = meta.cc ? (meta.case ? "_case" : "_ctrl") : ""
  """
  fastqc -t ${params.threads} ${trimmed} -o .
  
  for zip in *.zip; do
    unzip -o \$zip
  done

  # Collect fastqc data for output
  echo -e "filename\\tdeduplication_percentage\\tmarker" > step1_2_1_dedup_percentage.result
  end=1
  for dir in *fastqc/; do
    [[ -d \$dir ]] || continue

    # Append deduplication percentage
    out_value=\$(grep 'Total Deduplicated Percentage' \$dir/fastqc_data.txt | awk '{print \$4}')
    echo -e "${meta.id}_\$end\\t\$out_value\\t${params.marker}" >> "step1_2_1_dedup_percentage${suffix}.result"

    # Generate duplication summary
    echo -e "item\\t${meta.id}_\$end\\t${meta.id}_\$end" \
      > "step1_2_1_duplication_summary_${meta.id}_\$end.result"
    grep 'Sequence Duplication Levels' -A 15 \$dir/fastqc_data.txt \
      >> "step1_2_1_duplication_summary_${meta.id}${suffix}_\$end.result"

    # Generate FastQC summary
    {
      echo -e "${meta.id}_\$end\\tfastqc_test"
      awk -F "\t" '{print \$1, \$2}' OFS="\\t" "\$dir/summary.txt"
    } > "step1_2_1_fastqc_summary_${meta.id}${suffix}_\$end.result"

    (( end++ ))
  done

  before_dedup=\$(awk 'NR > 1 {s += \$2} END {print s * 0.01 / (NR-1)}' \
    "step1_2_1_dedup_percentage${suffix}.result")
  before_dup=\$(echo "scale=2; 1 - \$before_dedup" | bc -l)

  jq --argjson before_dup \$before_dup \
    '.fastqc_dup = \$before_dup' -n > "data_collection_s1_2${suffix}.json"

  # Get PE data R1 R2 deduplication difference percentage
  if [[ ${meta.paired} ]]; then
    per1=\$(tail -n 2 "step1_2_1_dedup_percentage${suffix}.result" | awk '{print \$2}' | sed -n '1p')
    per2=\$(tail -n 2 "step1_2_1_dedup_percentage${suffix}.result" | awk '{print \$2}' | sed -n '2p')
    dif=\$(echo "scale=2; (\$per1 - \$per2) * 200 / (\$per1 + \$per2)" | bc -l)
  else
    dif=0
  fi
  """
}
