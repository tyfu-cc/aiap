process COMPUTE_RUPR {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bed), path(peak)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # READS_IN_PEAKS=\$(intersectBed -a ${bed} -b ${peak} ${args} | awk -F '\t' '{sum += \$NF} END {print sum}')
    READS_IN_PEAKS=\$(bedtools intersect -a ${bed} -b ${peak} ${args} | wc -l)
    TOTAL_READS=\$(wc -l < ${bed})
    awk -v a="\$READS_IN_PEAKS" -v t="\$TOTAL_READS" -v OFS='\t' '{ print "'${prefix}'", a / t }' <<< "" > ${prefix}.RUPr.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
