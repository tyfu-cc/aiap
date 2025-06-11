process GETCHROMSIZES {
    tag "${fasta}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("${fasta}.sizes")         , emit: sizes
    tuple val(meta), path ("${fasta}.filtered.sizes"), emit: filtered_sizes
    tuple val(meta), path ("*.fai")                  , emit: fai
    tuple val(meta), path ("*.gzi")                  , emit: gzi, optional: true
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    samtools faidx ${fasta}
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes

    awk 'length(\$1) > 1 && length(\$1) < 6' OFS="\\t" ${fasta}.sizes \\
        | grep -v '^chrM' - > ${fasta}.filtered.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai
    touch ${fasta}.sizes
    if [[ "${fasta.extension}" == "gz" ]]; then
        touch ${fasta}.gzi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
