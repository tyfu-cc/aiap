process MACS2_CALLPEAK {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' :
        'biocontainers/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' }"

    input:
    tuple val(meta), path(ipbed), path(controlbed)
    val macs_gsize

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*.xls")                   , emit: xls
    path  "versions.yml"                             , emit: versions

    tuple val(meta), path("*.gappedPeak"), optional:true, emit: gapped
    tuple val(meta), path("*.bed")       , optional:true, emit: bed
    tuple val(meta), path("*.bdg")       , optional:true, emit: bdg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    def format = meta.single_end ? "BED" : "BEDPE"
    def control = controlbed ? "--control ${controlbed}" : ""
    if (args_list.contains("--format")) {
        def id = args_list.findIndexOf{ it == "--format" }
        format = args_list[id + 1]
        args_list.remove(id + 1)
        args_list.remove(id)
    }
    """
    macs2 \\
        callpeak \\
        ${args_list.join(" ")} \\
        --gsize ${macs_gsize} \\
        --format ${format} \\
        --name ${prefix} \\
        --treatment ${ipbed} \\
        ${control}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: \$(macs2 --version | sed -e "s/macs2 //g")
    END_VERSIONS
    """
}
