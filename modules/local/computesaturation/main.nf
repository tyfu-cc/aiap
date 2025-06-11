process COMPUTE_SATURATION {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' :
        'biocontainers/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' }"

    input:
    tuple val(meta), path(bed), path(peak)
    val macs_gsize

    output:
    path("*saturation.txt"), emit: saturation
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Subsample from the BED file obtained from peak calling on original sample
    TOTAL_READS=\$(wc -l < ${bed})
    for number in 5 10 20 30 40 50 60 70 80 90; do
        RATIO=\$(( \$TOTAL_READS * \$number / 100 ))
        shuf ${bed} | head -\$RATIO > ${prefix}.subsample.\$number.bed
    done

    # Peak calling for the subsamples
    for file in ${prefix}.subsample.*.bed; do
        macs2 \\
            callpeak \\
            --treatment \$file \\
            --gsize ${macs_gsize} \\
            --qvalue 0.01 \\
            --name \$file \\
            --keep-dup 1000 --nomodel --shift 0 --extsize 150
    done
    
    # Summarise the results
    echo -e "5\\n\$(seq 10 10 100)" > saturation.points

    for number in 5 10 20 30 40 50 60 70 80 90; do
        awk '{ print \$1 / 1000000 }' <(wc -l ${prefix}.subsample.\$number.bed)
    done | sort -n > saturation.reads

    TOTAL_LENGTH_OF_PEAKS=\$(awk '{s += \$3 - \$2 +1} END {print s}' ${peak})

    for number in 5 10 20 30 40 50 60 70 80 90; do
        file=${prefix}.subsample.\$number.bed_peaks.narrowPeak

        # Count number of peaks
        wc -l < \$file >> saturation.peaks

        # Calculate total overlapping region length
        length_of_overlapping_regions=\$(bedtools intersect -a \$file -b ${peak} \\
            | awk '{s += \$3 - \$2 + 1} END {print s + 0}')
        
        # Compute fraction and write
        awk -v a="\$length_of_overlapping_regions" \\
            -v b="\$TOTAL_LENGTH_OF_PEAKS" \\
            'BEGIN { print a / b }' >> saturation.ratios
    done

    awk '{ print \$1 / 1000000 }' <(wc -l < "${bed}") >> saturation.reads
    wc -l < ${peak} >> saturation.peaks
    echo 1 >> saturation.ratios
    
    {
        echo -e "ptc\\tread\\tpeak\\tratio"
        paste saturation.points saturation.reads saturation.peaks saturation.ratios
    } > "${prefix}.saturation.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}

