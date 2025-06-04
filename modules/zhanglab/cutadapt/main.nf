process CUTADAPT {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/1758869538eb8e658077cc14cd7a4e76fd9b6d73d3a68f85a70bf292e39e27c5/data' :
        'community.wave.seqera.io/library/cutadapt:5.0--991bbd2e184b7014' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trimmed.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')             , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = [
      "--no-trim",
      "-q 15,10",
      "-a AGATCGGAAGAGCACACGTCTGAAC"
    ]
    if (meta.single_end) {
      args += [
        "-m 36"
      ]
    } else {
      args += [
        "-A AGATCGGAAGAGCGTCGTGTAGGGA",
        "-m 36:36"
      ]
    }
    args = args.join(" ").trim()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trimmed.fastq.gz" : "-o ${prefix}_1.trimmed.fastq.gz -p ${prefix}_2.trimmed.fastq.gz"
    """
    cutadapt \\
        -Z \\
        --cores ${task.cpus} \\
        ${args} \\
        ${trimmed} \\
        ${reads} \\
        > ${prefix}.cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "${prefix}.trimmed.fastq.gz" : "${prefix}_1.trimmed.fastq.gz ${prefix}_2.trimmed.fastq.gz"
    """
    touch ${prefix}.cutadapt.log
    touch ${trimmed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
