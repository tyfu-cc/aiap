process REFORMAT_PEAK_CALLING_METRICS {
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta),
          path(peak),
          path(rupr),
          path(background),
          path(proen),
          path(suben),
          path(saturation)
    path(peak_calling_metrics_header)
    path(saturation_header)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = !meta.cc ? "${meta.id}" : "${meta.id}.case"
    """
    NUM_OF_PEAKS=\$(wc -l < ${peak})

    RUPR=\$(cut -f 2 ${rupr})
    
    BACKGROUND=\$(cut -f 2 ${background})

    PROEN=\$(cut -f 2 ${proen})

    SUBEN=\$(cut -f 2 ${suben})
    
    {
        cat ${peak_calling_metrics_header}
        echo -e "Sample\\t# Peaks called\\tRUPr\\tBackground\\tProEn\\tSubEn"
        echo -e "${prefix}\\t\$NUM_OF_PEAKS\\t\$RUPR\\t\$BACKGROUND\\t\$PROEN\\t\$SUBEN"
    } > ${prefix}.peak_calling_metrics_mqc.tsv

    if [ -s ${saturation} ]; then
        {
            cat ${saturation_header}
            tail -n +2 ${saturation} | cut -f1,4
        } > ${prefix}.saturation_mqc.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
