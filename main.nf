
include { PREPARE_GENOME          } from "./subworkflows/local/prepare_genome"
include { AIAP                    } from "./workflows/aiap"

workflow ZHANGLAB_AIAP {

    take:
    samplesheet // channel: samplesheet read in from --samplesheet


    main:
    ch_versions = Channel.empty()

    PREPARE_GENOME(
        params.fasta,
        params.gtf
    )
    ch_versions = ch_versions.mix( PREPARE_GENOME.out.versions )
    ch_versions.view { v -> "Versions: $v" }
    // PREPARE_GENOME.out.bwa_index.view()

    // WORKFLOW: Run pipeline
    AIAP(
        samplesheet,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.filtered_chrom_sizes,
        PREPARE_GENOME.out.promoters_bed,
        PREPARE_GENOME.out.coding_promoters_bed,
        PREPARE_GENOME.out.blacklist,
        PREPARE_GENOME.out.genome_filtered_bed,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.macs_gsize,
        ch_versions
    )


    emit:
    ch_versions
    // multiqc_report = AIAP.out.multiqc_report // channel: /path/to/multiqc_report.html

}


workflow {

    main:
    ZHANGLAB_AIAP(
        params.input
    )

}
