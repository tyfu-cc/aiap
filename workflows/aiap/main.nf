
include { TRIMGALORE                        } from "../../modules/zhanglab/trimgalore"
include { FASTQC                            } from "../../modules/zhanglab/fastqc"
include { BWA_MEM                           } from "../../modules/zhanglab/bwa/mem"
include { PRESEQ_LCEXTRAP                   } from "../../modules/zhanglab/preseq/lcextrap"
include { METHYLQA_ATAC                     } from "../../modules/zhanglab/methylqa/atac"
include { MACS2_CALLPEAK                    } from "../../modules/zhanglab/macs2/callpeak"

include { FILTER_BLACKLIST                  } from "../../modules/local/filterblacklist"
include { BEDTOOLS_GENOMECOV                } from "../../modules/local/bedtoolsgenomecov"
include { COMPUTE_RUPR                      } from "../../modules/local/computerupr"
include { COMPUTE_PROEN                     } from "../../modules/local/computeproen"
include { COMPUTE_SUBEN                     } from "../../modules/local/computesuben"
include { COMPUTE_BG                        } from "../../modules/local/computebg"
include { COMPUTE_SATURATION                } from "../../modules/local/computesaturation"
include { DENSITY_ESTIMATION                } from "../../modules/local/densityestimation"
include { REFORMAT_MAPPING_METRICS          } from "../../modules/local/reformatmappingmetrics"
include { REFORMAT_PEAK_CALLING_METRICS     } from "../../modules/local/reformatpeakcallingmetrics"
include { MULTIQC                           } from "../../modules/local/multiqc"

include { NORMALIZE_BEDGRAPH                } from "../../subworkflows/local/normalize_bedgraph"


workflow AIAP {

    take:
    ch_samplesheet      // channel: samplesheet read in from --samplesheet
    ch_fasta            // channel: path(genome.fasta)
    ch_chrom_sizes      // channel: path(genome.sizes)
    ch_filtered_chrom_sizes      // channel: path(genome.sizes)
    ch_promoters_bed
    ch_coding_promoters_bed
    ch_blacklist
    ch_genome_filtered_bed
    ch_bwa_index        // channel: path(bwa/index/)
    ch_macs_gsize       //   value: 
    ch_versions         // channel: [ path(versions.yml) ]

    // ch_fai              // channel: path(genome.fai)
    // ch_gtf              // channel: path(genome.gtf)
    // ch_gene_bed         // channel: path(gene.bed)
    

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet = Channel.fromPath(params.input)

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
        cc ? [[meta + ["case": true], case_reads], [meta + ["case": false], ctrl_reads]] : [[meta + ["case": true], case_reads]]
    }
    // ch_reads.view { v -> "Initial inputs: $v" }

    // Trim the adapters
    TRIMGALORE( ch_reads )
    ch_trimmed_fastq = TRIMGALORE.out.reads
    ch_multiqc_files = ch_multiqc_files.mix( TRIMGALORE.out.log.collect{ it[1] } )
    ch_multiqc_files = ch_multiqc_files.mix( TRIMGALORE.out.zip.collect{ it[1] } )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Do the alignment
    BWA_MEM( ch_trimmed_fastq, ch_bwa_index )

    // TODO: this step may cause a lot of trouble
    if (!params.skip_preseq) {
        PRESEQ_LCEXTRAP( BWA_MEM.out.bam )
    }

    // Run methylQA
    METHYLQA_ATAC(
        BWA_MEM.out.bam,
        ch_filtered_chrom_sizes
    )

    // Reformat the results from methylQA for MultiQC
    REFORMAT_MAPPING_METRICS(
        METHYLQA_ATAC.out.bed.join(METHYLQA_ATAC.out.report),
        Channel.value(file("${projectDir}/assets/multiqc/methylqa_table_header.txt"))
    )
    ch_multiqc_files = ch_multiqc_files.mix(REFORMAT_MAPPING_METRICS.out.mqc.collect { it[1] })
    
    // Filter out the blacklist regions from the BED file
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
            // def meta = entries[0].meta.subMap(["id", "single_end", "cc"])
            def meta = entries[0].meta
            def case_bed = entries.find{ it.case != null }?.case
            def ctrl_bed = entries.find{ it.ctrl != null }?.ctrl
            if (!meta.cc) {
                ctrl_bed = []
            }
            tuple(meta, case_bed, ctrl_bed)
        }
        .set{ ch_merged }
    // ch_merged.view { v -> "Merged: $v" }

    // Do the peak calling
    MACS2_CALLPEAK(
        ch_merged,
        ch_macs_gsize
    )
    ch_peak_count_mqc = MACS2_CALLPEAK.out.peak

    // Create channels: [ meta, case_bed, peak ]
    ch_merged
        .join(MACS2_CALLPEAK.out.peak, by: [0])
        .map{
            meta, case_bed, control_bed, peak ->
            [ meta, case_bed, peak ]
        }
        .set{ ch_bed_peak }
    // ch_bed_peak.view { v -> "Bed and peak channel: $v" }

    // Calculate reads under peaks ratio
    COMPUTE_RUPR( ch_bed_peak )

    COMPUTE_PROEN(
        ch_bed_peak,
        ch_coding_promoters_bed,
        ch_macs_gsize
    )

    COMPUTE_SUBEN(
        ch_bed_peak,
        ch_macs_gsize
    )

    COMPUTE_BG(
        ch_bed_peak,
        ch_promoters_bed,
        ch_filtered_chrom_sizes
    )

    ch_saturation = Channel.empty()
    if (!params.skip_saturation) {
        COMPUTE_SATURATION(
            ch_bed_peak,
            ch_macs_gsize
        )
    }
    ch_saturation = COMPUTE_SATURATION.out.txt

    // This is for plotting the peak length distribution
    DENSITY_ESTIMATION(
        METHYLQA_ATAC.out.insertdistro
            .join(MACS2_CALLPEAK.out.peak),
        Channel.value(file("${projectDir}/assets/multiqc/insertsize_header.txt")),
        Channel.value(file("${projectDir}/assets/multiqc/peaklength_header.txt"))
    )
    ch_multiqc_files = ch_multiqc_files.mix (
        DENSITY_ESTIMATION.out.mqc.map{ meta, files -> files }.flatten()
    )

    // Reformt the data so that they can be recognized by MultiQC
    REFORMAT_PEAK_CALLING_METRICS (
        MACS2_CALLPEAK.out.peak
            .join(COMPUTE_RUPR.out.txt)
            .join(COMPUTE_BG.out.txt)
            .join(COMPUTE_PROEN.out.txt)
            .join(COMPUTE_SUBEN.out.txt)
            .join(ch_saturation),
        Channel.value(file("${projectDir}/assets/multiqc/peak_calling_qc_header.txt")),
        Channel.value(file("${projectDir}/assets/multiqc/saturation_header.txt"))
    )
    ch_multiqc_files = ch_multiqc_files.mix (
        REFORMAT_PEAK_CALLING_METRICS.out.mqc.map{ meta, files -> files }.flatten()
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

    // ch_multiqc_files.collect().view { v -> "MultiQC inputs: $v }
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

}
