process DENSITY_ESTIMATION {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-dcfb6eba6adda57b9d4a990b9096cb040320914f:588e2c290fce5c5c11ef2340b6184370efd2c628-0' :
        'quay.io/biocontainers/mulled-v2-dcfb6eba6adda57b9d4a990b9096cb040320914f:588e2c290fce5c5c11ef2340b6184370efd2c628-0' }"

    input:
    tuple val(meta), path(insertdistro), path(peak)
    path(insertsize_header)
    path(peaklength_header)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sort -n ${insertdistro} \
        | uniq -c \
        | awk 'BEGIN {OFS="\t"} {print \$2, \$1}' \
        > ${prefix}.insertsize.dist.tsv

    density_estimation.R ${prefix}.insertsize.dist.tsv ${prefix}.insertsize.tsv
    cat ${insertsize_header} ${prefix}.insertsize.tsv > ${prefix}.insertsize_mqc.tsv

    awk '{print \$3 - \$2 + 1}' ${peak} \
        | sort -n \
        | uniq -c \
        | awk 'BEGIN {OFS="\t"} {print \$2, \$1}' \
        > ${prefix}.peaklen.dist.tsv

    density_estimation.R ${prefix}.peaklen.dist.tsv ${prefix}.peaklen.tsv
    cat ${peaklength_header} ${prefix}.peaklen.tsv > ${prefix}.peaklen_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript --version | sed -e "s/Rscript (R) version //g")
    END_VERSIONS
    """
}

