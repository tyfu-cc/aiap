process BEDTOOLS_GENOMECOV {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta), path(report)
    path chrom_sizes

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.txt")     , emit: scale_factor
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pe = meta.single_end ? "" : "-pc"
    """
    # Normalization using CPM?
    SCALING_FACTOR=\$(grep "non-redundant uniquely mapped reads" ${report} | awk '{print 1000000 / \$NF}')
    echo \$SCALING_FACTOR > ${prefix}.scaling_factor.txt

    bedtools \\
        genomecov \\
        -i ${bed} \\
        -g ${chrom_sizes} \\
        -bg \\
        -scale \$SCALING_FACTOR \\
        ${pe} \\
        ${args} \\
        > tmp.bg

    bedtools sort -i tmp.bg > ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
