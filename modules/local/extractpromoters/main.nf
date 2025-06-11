process EXTRACT_PROMOTERS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path gtf

    output:
    path "promoters.bed"       , emit: promoters_bed
    path "coding_promoters.bed", emit: coding_promoters_bed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk '\$3 == "gene" {
        chr = \$1;
        strand = \$7;
        tss = (strand == "+") ? \$4 : \$5;
        start = tss - 1000;
        end = tss + 1000;
        if (start < 0) start = 0;

        if (length(chr) > 1 && length(chr) < 6 && chr != "chrM") {
            print chr, start, end >> "promoters.bed";
            if (\$0 ~ /gene_type "protein_coding"/)
                print chr, start, end >> "coding_promoters.bed";
        }
    }' OFS="\\t" ${gtf}
    """
}
