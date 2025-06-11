
include { GUNZIP as GUNZIP_FASTA            } from "../../../modules/zhanglab/gunzip"
include { BWA_INDEX                         } from "../../../modules/zhanglab/bwa/index"
include { KHMER_UNIQUEKMERS                 } from "../../../modules/zhanglab/khmer/uniquekmers"

include { GENOME_BLACKLIST_REGIONS          } from "../../../modules/local/genomeblacklistregions"
include { EXTRACT_PROMOTERS                 } from "../../../modules/local/extractpromoters"
include { GETCHROMSIZES                     } from "../../../modules/local/getchromsizes"


workflow PREPARE_GENOME {

    take:
    fasta //      file: /path/to/genome.fasta
    gtf   //      file: /path/to/genome.gtf


    main:
    ch_versions = Channel.empty()

    ch_fasta = Channel.empty()
    if (fasta.endsWith(".gz")) {
        ch_fasta = GUNZIP_FASTA( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map{ it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    }

    // Create chromosome sizes file
    GETCHROMSIZES( ch_fasta.map{ [ [:], it ] } )
    ch_fai                  = GETCHROMSIZES.out.fai.map{ it[1] }
    ch_chrom_sizes          = GETCHROMSIZES.out.sizes.map{ it[1] }
    ch_filtered_chrom_sizes = GETCHROMSIZES.out.filtered_sizes.map{ it[1] }
    ch_versions             = ch_versions.mix(GETCHROMSIZES.out.versions)

    // Get something
    ch_gtf = Channel.empty()
    if (gtf.endsWith(".gz")) {
        ch_gtf = GUNZIP_GTF ( [ [:], file(gtf, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = Channel.value(file(gtf, checkIfExists: true))
    }

    ch_promoters_bed = Channel.empty()
    ch_coding_promoters_bed = Channel.empty()
    EXTRACT_PROMOTERS( ch_gtf )
    ch_promoters_bed = EXTRACT_PROMOTERS.out.promoters_bed
    ch_coding_promoters_bed = EXTRACT_PROMOTERS.out.coding_promoters_bed

    ch_bwa_index = Channel.empty()
    if (params.bwa_index) {
        if (params.bwa_index.endsWith(".tar.gz")) {
            ch_bwa_index = UNTAR_BWA_INDEX( [ [:], params.bwa_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
        } else {
            ch_bwa_index = [ [:], file(params.bwa_index) ]
        }
    } else {
        ch_bwa_index = BWA_INDEX( ch_fasta.map { [ [:], it ] } ).index
        ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
    }

    ch_blacklist = Channel.empty()
    if (params.blacklist) {
        if (params.blacklist.endsWith(".gz")) {
            ch_blacklist = GUNZIP_BLACKLIST( [ [:], params.blacklist ] ).gunzip.map{ it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.value(file(params.blacklist))
        }
    }

    ch_filtered_genome_bed = Channel.empty()
    GENOME_BLACKLIST_REGIONS(
        ch_chrom_sizes,
        ch_blacklist.ifEmpty([]),
        params.mito_name ?: "",
        params.keep_mito ?: false,
    )
    ch_genome_filtered_bed = GENOME_BLACKLIST_REGIONS.out.bed
    ch_versions = ch_versions.mix(GENOME_BLACKLIST_REGIONS.out.versions)


    ch_macs_gsize = params.macs_gsize
    if (!params.macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map{ it.text.trim() }
        ch_versions   = ch_versions.mix(KHMER_UNIQUEKMERS.out.versions)
    }


    emit: 
    fasta                = ch_fasta                      //    path: genome.fasta
    // fai           = ch_fai                        //    path: genome.fai
    // gtf           = ch_gtf                        //    path: genome.gtf
    chrom_sizes          = ch_chrom_sizes
    filtered_chrom_sizes = ch_filtered_chrom_sizes
    promoters_bed        = ch_promoters_bed
    coding_promoters_bed = ch_coding_promoters_bed
    blacklist            = ch_blacklist
    genome_filtered_bed  = ch_genome_filtered_bed
    bwa_index            = ch_bwa_index
    macs_gsize           = ch_macs_gsize

}
