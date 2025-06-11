include { BEDTOOLS_GENOMECOV    } from "../../../modules/local/bedtoolsgenomecov"
include { UCSC_BEDGRAPHTOBIGWIG } from "../../../modules/zhanglab/ucsc/bedgraphtobigwig/main"

workflow NORMALIZE_BEDGRAPH {

    take:
    ch_bed         // channel: [ val(meta), path(bed) ]
    ch_report      // channel: [ val(meta), path(*.report) ]
    ch_chrom_sizes // channel: [ path(chrom_sizes) ]
    

    main:
    ch_versions = Channel.empty()

    // Create bedGraph coverage track
    BEDTOOLS_GENOMECOV(
        ch_bed,
        ch_report,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    // Create bigWig coverage tracks
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())


    emit:
    bedgraph     = BEDTOOLS_GENOMECOV.out.bedgraph      // channel: [ val(meta), [ bedgraph ] ]
    scale_factor = BEDTOOLS_GENOMECOV.out.scale_factor  // channel: [ val(meta), [ txt ] ]
    bigwig       = UCSC_BEDGRAPHTOBIGWIG.out.bigwig     // channel: [ val(meta), [ bigwig ] ]
    versions     = ch_versions                          // channel: [ versions.yml ]

}
