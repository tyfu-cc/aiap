process GTFUTILS_TOBED {
    tag "${gtf}"
    label "process_low"

    container "quay.io/biocontainers/ngsutils:0.5.9--py27_0"

    input:
    path gtf

    output:
    path "*.bed"       , emit: bed, optional: true
    path "versions.yml", emit: versions, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    gtfutils \\
        tobed \\
        ${args} \\
        ${gtf} > tmp.txt
    """

    stub:
    """
    touch ${gtf.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
