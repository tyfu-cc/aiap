process COMPUTE_PROEN {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bed), path(peak)
    path coding_promoter
    val macs_gsize

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TOTAL_READS=\$(wc -l < ${bed})

    # Peaks that overlap with coding promoters
    bedtools intersect -a ${peak} -b ${coding_promoter} -u > peaks.at.promoters.bed

    # Reads that overalp with both coding promoters and peaks
    READS_IN_PROMOTERS_N_PEAKS=\$(bedtools intersect -a ${bed} -b peaks.at.promoters.bed -f 0.5 -u | wc -l)

    # Promoters that are >=50% overlapped by the selected peaks
    NUM_OF_PROMOTERS=\$(bedtools intersect -a ${coding_promoter} -b peaks.at.promoters.bed -F 0.5 -u | wc -l)

    # Assume 2000 bp per promoter
    TOTAL_LENGTH_OF_PROMOTERS=\$(awk -v n="\$NUM_OF_PROMOTERS" 'BEGIN { print n * 2000 }')
  
    # Coding promoter enrichment
    # enrichment = (reads in promoter / promoter length) / (total reads / genome size)
    awk -v id="${meta.id}" \\
        -v r="\$READS_IN_PROMOTERS_N_PEAKS" \\
        -v l="\$TOTAL_LENGTH_OF_PROMOTERS" \\
        -v t="\$TOTAL_READS" \\
        -v g="${macs_gsize}" \\
        'BEGIN { print id "\t" (r / l) / (t / g) }' > ${meta.id}.ProEn.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
