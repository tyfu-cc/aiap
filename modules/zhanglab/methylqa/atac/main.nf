process METHYLQA_ATAC {
    tag "${meta.id}"
    label "process_medium"

    container "ghcr.io/tyfu-cc/methylqa:0.2.1"

    input:
    tuple val(meta), path(bam)
    path(chrom_sizes)

    output:
    tuple val(meta), path("*.open.bed")     , emit: bed
    tuple val(meta), path("*.open.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.bigWig")       , emit: bigwig
    tuple val(meta), path("*.insertdistro") , emit: insertdistro
    tuple val(meta), path("*.report")       , emit: report
    
    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    methylQA \\
        atac \\
        ${args} \\
        -o ${prefix} \\
        ${chrom_sizes} \\
        ${bam}
    """
}
