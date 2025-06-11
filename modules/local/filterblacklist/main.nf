process FILTER_BLACKLIST {
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bed)
    path blacklist

    output:
    tuple val(meta), path("*.blacklist.filtered.bed"), emit: bed
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_out = "${bed.baseName}.blacklist.filtered.bed"
    if (blacklist) {
        """
        bedtools \\
            intersect \\
            -v -a ${bed} \\
            -b ${blacklist} \\
            > ${file_out}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        cp ${bed} ${file_out}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
}
