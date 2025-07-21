process COMPUTE_SUBEN {
    tag "${meta.id}"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' :
        'biocontainers/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' }"

    input:
    tuple val(meta), path(bed), path(peak)
    val macs_gsize

    output:
    tuple val(meta), path("*.suben.txt"), emit: txt
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Subsample 10M reads or use all reads if fewer
    NUM_OF_USEFUL_READS=\$(wc -l < ${bed})
    SUBSAMPLE_SIZE=\$(( \$NUM_OF_USEFUL_READS > 10000000 ? 10000000 : \$NUM_OF_USEFUL_READS))
    shuf ${bed} | head -\$SUBSAMPLE_SIZE > subsample.bed

    # Peak calling with MACS
    macs2 \\
        callpeak \\
        -t subsample.bed \\
        -g ${macs_gsize} \\
        -q 0.01 \\
        -n temp \\
        --keep-dup 1000 \\
        --nomodel \\
        --shift 0 \\
        --extsize 150
    TOTAL_PEAK_LENGTH=\$(awk '{s += \$3 - \$2 + 1} END {print s}' temp_peaks.narrowPeak)

    # Count reads in peaks
    NUM_OF_READS_UNDER_PEAKS=\$(bedtools intersect -a subsample.bed -b temp_peaks.narrowPeak -f 0.5 | wc -l)

    # Compute enrichment score
    UPPER=\$(awk \\
        -v n=\$NUM_OF_READS_UNDER_PEAKS \\
        -v l=\$TOTAL_PEAK_LENGTH \\
        -v g=${macs_gsize} \\
        'BEGIN {
            print n / l + 1e7 / g
        }')
    LOWER=\$(awk \\
        -v n=\$NUM_OF_USEFUL_READS \\
        -v l=\$TOTAL_PEAK_LENGTH \\
        -v g=${macs_gsize} \\
        'BEGIN {
            print (n + 1e7) / (g - l)
        }')
    SUBEN=\$(awk -v u=\$UPPER -v l=\$LOWER 'BEGIN { print u / l }')

    echo -e "${prefix}\\t\$SUBEN" > ${prefix}.suben.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}

