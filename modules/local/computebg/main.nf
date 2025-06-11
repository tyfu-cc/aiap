process COMPUTE_BG {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fc325951871d402a00bdf9d0e712a5b81b8e0cb3:38034b9703d6561a40bcaf2f1ec16f8b158fde97-0' :
        'biocontainers/mulled-v2-fc325951871d402a00bdf9d0e712a5b81b8e0cb3:38034b9703d6561a40bcaf2f1ec16f8b158fde97-0' }"

    input:
    tuple val(meta), path(bed), path(peak)
    path promoter_file
    path chrom_sizes

    output:
    tuple val(meta), path("*background.txt"), emit: txt
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Peaks that are at promoters
    bedtools intersect -a ${peak} -b ${promoter_file} -u \\
        | awk '{print \$1"\\t"\$2"\\t"\$3"\\t""1""\\t"\$9}' > peaks.at.promoters.bed
    # Peaks that aren't at promoters
    bedtools intersect -a ${peak} -b ${promoter_file} -v \\
        | awk '{print \$1"\\t"\$2"\\t"\$3"\\t""0""\\t"\$9}' > peaks.at.non.promoters.bed

    NUM_OF_PEAKS_AT_PROMOTERS=\$(wc -l < peaks.at.promoters.bed)
    NUM_OF_READS_AT_PROMOTERS=\$(bedtools intersect -a ${bed} -b peaks.at.promoters.bed -u -f 0.5 | wc -l)
    NUM_OF_PEAKS_AT_NON_PROMOTERS=\$(wc -l < peaks.at.non.promoters.bed)
    NUM_OF_READS_AT_NONO_PROMOTERS=\$(bedtools intersect -a ${bed} -b peaks.at.non.promoters.bed -u -f 0.5 | wc -l)

    cat peaks.at.promoters.bed peaks.at.non.promoters.bed \\
        | sort -k5 -n -r > merged.bed

    LINES=\$(wc -l < merged.bed)
    if (( \$LINES > 100 )); then
        promoter_bin.py merged.bed ${prefix}.bin.txt
    else
        echo "Warning: total number of peaks is fewer than 100, promoter bin step would be skipped. At least 100 peaks are required."
    fi

    # Generate random genomic regions as background
    random_chr.py ${chrom_sizes} chr.bg

    TOTAL_READS=\$(wc -l < ${bed})

    # Create 200kb windows centered on peaks (used to exclude peak-rich regions from background)
    awk '{
        x = int((\$2 + \$3) / 2)
        if (x - 100000 < 0) x = 100000
        print \$1, x - 100000, x + 100000, \$4
    }' OFS="\\t" ${peak} > temp.txt

    # Get 50,000 random background regions not overlapping any peak regions
    bedtools intersect -a chr.bg -b temp.txt -v | shuf - | head -50000 | sort -V -k1,1 -k2,2n > background.txt

    # Extracts reads overlapping background regions.
    bedtools intersect -a ${bed} -b background.txt -u -f 0.5 | sort -V -k1,1 -k2,2n > temp.txt

    # Calculates RPKM values in background regions
    rpkm_bin.py background.txt temp.txt \$TOTAL_READS ${prefix}.background.rpkm.txt

    # Calculate background distribution
    NUM_OF_LOWER_BG=\$(awk '\$6 <= 0.188 {print \$0}' ${prefix}.background.rpkm.txt | wc -l)
    NUM_OF_LOW_BG=\$(awk '\$6 <= 0.377 {print \$0}' ${prefix}.background.rpkm.txt | wc -l)
    NUM_OF_HIGH_BG=\$(awk '\$6 > 0.377 {print \$0}' ${prefix}.background.rpkm.txt | wc -l)

    awk -v id="${prefix}" -v n="\$NUM_OF_HIGH_BG" 'BEGIN { OFS="\t"; print id, n / 50000 }' > ${prefix}.background.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
