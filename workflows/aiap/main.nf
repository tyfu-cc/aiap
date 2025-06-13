// include { MULTIQC } from "${projectDir}/modules/local/multiqc.nf"

include { NORMALIZE_BEDGRAPH                } from "../../subworkflows/local/normalize_bedgraph"

include { CUTADAPT                          } from "../../modules/zhanglab/cutadapt"
include { FASTQC                            } from "../../modules/zhanglab/fastqc"
include { BWA_MEM                           } from "../../modules/zhanglab/bwa/mem"
include { PRESEQ_LCEXTRAP                   } from "../../modules/zhanglab/preseq/lcextrap"
include { METHYLQA_ATAC                     } from "../../modules/zhanglab/methylqa/atac"
include { MACS2_CALLPEAK                    } from "../../modules/zhanglab/macs2/callpeak"

include { FILTER_BLACKLIST                  } from "../../modules/local/filterblacklist"
include { BEDTOOLS_GENOMECOV                } from "../../modules/local/bedtoolsgenomecov"
include { COMPUTE_RUPR                      } from "../../modules/local/computerupr"
include { COMPUTE_PROEN                     } from "../../modules/local/computeproen"
include { COMPUTE_BG                        } from "../../modules/local/computebg"
include { COMPUTE_SATURATION                } from "../../modules/local/computesaturation"
include { MULTIQC                           } from "../../modules/local/multiqc"


workflow AIAP {

    take:
    ch_samplesheet      // channel: samplesheet read in from --samplesheet
    ch_fasta            // channel: path(genome.fasta)
    // ch_fai              // channel: path(genome.fai)
    // ch_gtf              // channel: path(genome.gtf)
    ch_chrom_sizes      // channel: path(genome.sizes)
    ch_filtered_chrom_sizes      // channel: path(genome.sizes)
    ch_promoters_bed
    ch_coding_promoters_bed
    ch_blacklist
    ch_genome_filtered_bed
    // ch_gene_bed         // channel: path(gene.bed)
    ch_bwa_index        // channel: path(bwa/index/)
    ch_macs_gsize       //   value: 
    ch_versions         // channel: [ path(versions.yml) ]
    

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet = Channel.fromPath(params.samplesheet)

    // Parse it line by line
    ch_reads = ch_samplesheet.splitCsv(header:true).flatMap {
        case_r1 = it["case_r1"]
        case_r2 = it["case_r2"]
        ctrl_r1 = it["ctrl_r1"]
        ctrl_r2 = it["ctrl_r2"]

        // Detect whether single-end or paired-end
        single_end = case_r2.toString() == "" ? true : false
        
        // Detect whether there's control
        cc = ctrl_r1.toString() == "" ? false : true

        // The "meta" map 
        meta = [id: it["id"], single_end: single_end, cc: cc]

        // Case reads
        case_reads = single_end ? [case_r1] : [case_r1, case_r2]

        // Control reads
        ctrl_reads = single_end ? [ctrl_r1] : [ctrl_r1, ctrl_r2]
        
        // Return a nested map, the first entry is the meta map, the second one is the read(s)
        cc ? [[meta + ["case": true], case_reads], [meta + ["case": false], ctrl_reads]] : [[meta, case_reads]]
    }
    ch_reads.view()

    CUTADAPT( ch_reads )
    ch_trimmed_fastq = CUTADAPT.out.reads

    FASTQC( ch_trimmed_fastq )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    BWA_MEM( ch_trimmed_fastq, ch_bwa_index )

    if (!params.skip_preseq) {
        PRESEQ_LCEXTRAP( BWA_MEM.out.bam )
    }

    METHYLQA_ATAC(
        BWA_MEM.out.bam,
        ch_filtered_chrom_sizes
    )
    ch_multiqc_files = ch_multiqc_files.mix(METHYLQA_ATAC.out.report.collect{ it[1] })

    FILTER_BLACKLIST(
        METHYLQA_ATAC.out.bed,
        ch_blacklist
    )

    NORMALIZE_BEDGRAPH(
        FILTER_BLACKLIST.out.bed,
        METHYLQA_ATAC.out.report,
        ch_filtered_chrom_sizes
    )
   
    // Create channel: [ meta, case.*.bed, (ctrl.*.bed) ]
    FILTER_BLACKLIST.out.bed
        .map{ meta, bed ->
            def id = meta.id
            def entry = meta.cc ? (meta.case ? [case: bed, ctrl: null, meta: meta] : [case: null, ctrl: bed, meta: meta])
                : [case:bed, ctrl: null, meta: meta]
            tuple(id, entry)
        }
        .groupTuple(by: 0)
        .map{ id, entries ->
            def meta = entries[0].meta.subMap(["id", "single_end", "cc"])
            def case_bed = entries.find{ it.case != null }?.case
            def ctrl_bed = entries.find{ it.ctrl != null }?.ctrl
            if (!meta.cc) {
                ctrl_bed = []
            }
            tuple(meta, case_bed, ctrl_bed)
        }
        .set{ ch_merged }
    ch_merged.view()

    MACS2_CALLPEAK(
        ch_merged,
        ch_macs_gsize
    )
    ch_multiqc_files = ch_multiqc_files.mix(MACS2_CALLPEAK.out.xls.collect{ it[1] })


    // Create channels: [ meta, case_bed, peak ]
    ch_merged
        .join(MACS2_CALLPEAK.out.peak, by: [0])
        .map{
            meta, case_bed, control_bed, peak ->
            [ meta, case_bed, peak ]
        }
        .set{ ch_bed_peak }
    ch_bed_peak.view()

    COMPUTE_RUPR( ch_bed_peak )

    COMPUTE_PROEN(
        ch_bed_peak,
        ch_coding_promoters_bed,
        ch_macs_gsize
    )

    COMPUTE_BG(
        ch_bed_peak,
        ch_promoters_bed,
        ch_filtered_chrom_sizes
    )

    COMPUTE_SATURATION(
        ch_bed_peak,
        ch_macs_gsize
    )


    // MODULE: MultiQC
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc/multiqc_config.yaml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    ch_multiqc_files.collect().view()
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

}
