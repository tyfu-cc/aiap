process REFORMAT_MAPPING_METRICS {
    tag "${meta.id}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(bed), path(report)
    path(mapping_metrics_header)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = !meta.cc ? "${meta.id}" : "${meta.id}.${meta.case ? 'case' : 'ctrl'}"
    """
    USEFUL_SINGLE_ENDS=\$(wc -l < ${bed})

    TOTAL_READS=\$(grep -E "^total reads" ${report} | awk -F': ' '{print \$2}')

    MAPPABLE_READS=\$(grep -E "^mappable reads" ${report} | awk -F': ' '{print \$2}')

    UNIQUELY_MAPPED=\$(grep -E "^uniquely mapped reads" ${report} | awk -F': ' '{print \$2}')

    NON_REDUNDANT=\$(grep -E "^non-redundant uniquely mapped" ${report} | awk -F': ' '{print \$2}')

    {
        cat ${mapping_metrics_header}
        echo -e "Sample\\tUseful single ends\\tTotal reads\\tMappable reads\\tUniquely mapped reads\\tNon-redundant uniquely mapped reads"
        echo -e "${prefix}\\t\$USEFUL_SINGLE_ENDS\\t\$TOTAL_READS\\t\$MAPPABLE_READS\\t\$UNIQUELY_MAPPED\\t\$NON_REDUNDANT"
    } > ${prefix}.mapping_metrics_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
