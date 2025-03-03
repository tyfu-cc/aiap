
include { CUTADAPT_TRIM } from "${projectDir}/modules/local/cutadapt.nf"
include { FASTQC } from "${projectDir}/modules/local/fastqc.nf"
include { BWA_ALIGN } from "${projectDir}/modules/local/bwa_align.nf"
include { PRESEQ } from "${projectDir}/modules/local/preseq.nf"
include { METHYLQA } from "${projectDir}/modules/local/methylqa.nf"
include { NORMALIZE_BEDGRAPH } from "${projectDir}/modules/local/normalize_bedgraph.nf"
include { MACS_PEAK_CALLING } from "${projectDir}/modules/local/macs_peak_calling.nf"
include { COMPUTE_RUP } from "${projectDir}/modules/local/compute_rup.nf"
include { COMPUTE_ENRICHMENT } from "${projectDir}/modules/local/compute_enrichment.nf"
include { COMPUTE_SATURATION } from "${projectDir}/modules/local/compute_saturation.nf"
include { COMPUTE_BACKGROUND } from "${projectDir}/modules/local/compute_background.nf"
include { MULTIQC } from "${projectDir}/modules/local/multiqc.nf"

workflow {
  // Channel for the samplesheet
  samplesheet_ch = Channel.fromPath(params.samplesheet)

  // Parse it line by line
  reads_ch = samplesheet_ch.splitCsv(header:true).flatMap {
    case_r1 = it["case_r1"]
    case_r2 = it["case_r2"]
    ctrl_r1 = it["ctrl_r1"]
    ctrl_r2 = it["ctrl_r2"]

    // Detect whether single-end or paired-end
    paired = case_r2.toString() == "" ? false : true
    
    // Detect whether there's control
    cc = ctrl_r1.toString() == "" ? false : true

    // The "meta" map 
    meta = [id: it["id"], paired: paired, cc: cc]

    // Case reads
    case_reads = paired ? [case_r1, case_r2] : [case_r1]

    // Control reads
    ctrl_reads = paired ? [ctrl_r1, ctrl_r2] : [ctrl_r1]
    
    // Return a nested map, the first entry is the meta map, the second one is the read(s)
    cc ? [[meta + ["case": true], case_reads], [meta + ["case": false], ctrl_reads]] : [[meta, case_reads]]
  }

  // reads_ch.view()

  CUTADAPT_TRIM(
    reads_ch,
    tuple(params.adapter_1, params.adapter_2)
  )

  FASTQC(
    CUTADAPT_TRIM.out[0],
  )

  ref_files = Channel.fromPath(params.bwa_ref).collect()
  BWA_ALIGN(
    CUTADAPT_TRIM.out[0],
    ref_files
  )

  ch_preseq_yield_mqc = Channel.empty()
  if (params.run_preseq) {
    PRESEQ(
      BWA_ALIGN.out  
    )
    ch_preseq_yield_mqc = PRESEQ.out
  }

  METHYLQA(
    BWA_ALIGN.out,
    file(params.chrom_size),
    CUTADAPT_TRIM.out.trimlog,
    FASTQC.out.dedup_pct_mqc
  )

  nochrM_chrom_size = METHYLQA.out.nochrM

  ch_s3_1_bedGraph = METHYLQA.out[0].map { 
    tuple(
      it[0],
      it[1..-1].flatten().find { it.name.endsWith("report") },
      it[1..-1].flatten().find { it.name.endsWith("open.bedGraph") }
    )
  }

  NORMALIZE_BEDGRAPH(
    ch_s3_1_bedGraph,
    nochrM_chrom_size,
    file(params.black_list)
  ) 

  grouped = METHYLQA.out[0].map {
    tuple(
      it[0].id, 
      it[0], 
      it[1..-1].flatten().find { it.name.endsWith("open.bed") }
    )
  }.groupTuple(by: 0).map {
    meta = it[1].first().subMap( ["id", "paired", "cc"] )
    tuple(meta, it[2])
  }

  // grouped.view()

  // Channel.of(
  //   [[id: "test", paired: true, cc: true, case: true], '/path/to/region1_chr1.vcf'],
  //   [[id: "test", paired: true, cc: true, case: false], '/path/to/region2_chr2.vcf']
  // ).map {
  //   [it[0].id, it[0], it[1]]
  // }.groupTuple(by: 0).map {
  //   meta = it[1].first().subMap( ["id", "paired", "cc"] )
  //   [meta, it[2]]
  // }.view()
  
  MACS_PEAK_CALLING(
    grouped,
    file(params.black_list)
  )

  ch_bed_and_peak = MACS_PEAK_CALLING.out[0].map {
    tuple(
      it[0],
      it[1..-1].flatten().find { it.name.endsWith("open.bed") },
      it[1..-1].flatten().find { it.name.endsWith("narrowPeak") }
    )
  }

  COMPUTE_RUP(
    ch_bed_and_peak,
    nochrM_chrom_size
  )

  COMPUTE_ENRICHMENT(
    ch_bed_and_peak,
    file(params.coding_promoter)
  )
  
  ch_saturation_result_mqc = Channel.empty()
  if (params.run_saturation) {
    COMPUTE_SATURATION(
      ch_bed_and_peak,
    )
    ch_saturation_result_mqc = COMPUTE_SATURATION.out.saturation_mqc
  }

  COMPUTE_BACKGROUND(
    ch_bed_and_peak,
    file(params.promoter_file),
    nochrM_chrom_size
  )

  // Channel.fromPath("${projectDir}-test-data/*", type: "file")
  //   .collect()
  //   .view()

  //MACS_PEAK_CALLING.out[0].map { 
  //  file = Channel.fromPath("${projectDir}-${it[0].id}/*", type: "file")
  //  tuple(it[0], file.view())
  //}.view()

  // CUTADAPT_TRIM.out[0].flatten().filter( Map ).view()
  // S4_6_SUMMARIZE(
  //   // Meta
  //   S1_1_TRIM.out[0].flatten().filter( Map ),

  //   // Trimlog
  //   S1_1_TRIM.out.trimlog,

  //   // Outputs from methylQA
  //   METHYLQA.out[0].map { it[1..-1].flatten().collect() },

  //   // Results from peak calling
  //   MACS_PEAK_CALLING.out[0].map { it[1..-1].flatten().collect() },


  //   // JSON files
  //   S1_1_TRIM.out.json,
  //   METHYLQA.out.json,
  //   S4_1_RUP.out.json,
  //   S4_2_ENRICH.out.json,
  //   S4_4_BACKGROUND.out.json
  // )

  ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yaml", checkIfExists: true)
  MULTIQC(
    ch_multiqc_config,
    CUTADAPT_TRIM.out.trimlog.mix(
      FASTQC.out.zip.flatten(),
      ch_preseq_yield_mqc,
      METHYLQA.out.insert_size_mqc,
      METHYLQA.out.mapping_status_mqc,
      METHYLQA.out.dup_pct_mqc,
      MACS_PEAK_CALLING.out[0].map { it[1..-1].flatten() },
      MACS_PEAK_CALLING.out.peak_len_dist_mqc,
      COMPUTE_RUP.out.nrup_mqc,
      COMPUTE_RUP.out.rup_ratio_mqc,
      COMPUTE_ENRICHMENT.out.enrich_score_mqc,
      ch_saturation_result_mqc,
      COMPUTE_BACKGROUND.out.pnp_peaks_mqc,
      COMPUTE_BACKGROUND.out.pnp_reads_mqc,
      COMPUTE_BACKGROUND.out.background_0_3777_mqc
    ).collect()
  )
}
