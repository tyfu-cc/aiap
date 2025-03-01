
include { S1_1_TRIM } from "${projectDir}/modules/local/s1_1_trim.nf"
include { S1_2_QC } from "${projectDir}/modules/local/s1_2_qc.nf"
include { S2_1_ALIGN } from "${projectDir}/modules/local/s2_1_align.nf"
include { S2_2_COMPLEXITY} from "${projectDir}/modules/local/s2_2_complexity.nf"
include { S3_1_QA } from "${projectDir}/modules/local/s3_1_qa.nf"
include { S3_2_NORMALIZE } from "${projectDir}/modules/local/s3_2_normalize.nf"
include { S3_3_PEAKS } from "${projectDir}/modules/local/s3_3_peaks.nf"
include { S4_1_RUP } from "${projectDir}/modules/local/s4_1_rup.nf"
include { S4_2_ENRICH } from "${projectDir}/modules/local/s4_2_enrich.nf"
include { S4_3_SATURATION } from "${projectDir}/modules/local/s4_3_saturation.nf"
include { S4_4_BACKGROUND } from "${projectDir}/modules/local/s4_4_background.nf"
include { S4_6_SUMMARIZE } from "${projectDir}/modules/local/s4_6_summarize.nf"

workflow {
  // Channel for the samplesheet
  samplesheet_ch = Channel.fromPath(params.samples)

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
    
    // We return a nested map, the first entry is the meta map, the second one is the read(s)
    cc ? [[meta + ["case": true], case_reads], [meta + ["case": false], ctrl_reads]] : [[meta, case_reads]]
  }

  reads_ch.view()

  S1_1_TRIM(
    reads_ch,
    tuple(params.adapter_1, params.adapter_2)
  )

  S1_2_QC(
    S1_1_TRIM.out[0],
  )

  ref_files = Channel.fromPath(params.bwa_ref).collect()
  // ref_files.view()
  S2_1_ALIGN(
    S1_1_TRIM.out[0],
    ref_files
  )

  if (params.run_preseq) {
    S2_2_COMPLEXITY(
      S2_1_ALIGN.out  
    )
  }

  S3_1_QA(
    S2_1_ALIGN.out,
    file(params.chrom_size),
    S1_1_TRIM.out.json
  )

  nochrM_chrom_size = S3_1_QA.out.nochrM

  ch_s3_1_bedGraph = S3_1_QA.out[0].map { 
    tuple(
      it[0],
      it[1..-1].flatten().find { it.name.endsWith("report") },
      it[1..-1].flatten().find { it.name.endsWith("open.bedGraph") }
    )
  }

  S3_2_NORMALIZE(
    ch_s3_1_bedGraph,
    nochrM_chrom_size,
    file(params.black_list)
  ) 

  grouped = S3_1_QA.out[0].map {
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
  
  S3_3_PEAKS(
    grouped,
    file(params.black_list)
  )

  ch_s3_3_bed_and_peak = S3_3_PEAKS.out[0].map {
    tuple(
      it[0],
      it[1..-1].flatten().find { it.name.endsWith("open.bed") },
      it[1..-1].flatten().find { it.name.endsWith("narrowPeak") }
    )
  }

  S4_1_RUP(
    ch_s3_3_bed_and_peak,
    nochrM_chrom_size
  )

  S4_2_ENRICH(
    ch_s3_3_bed_and_peak,
    file(params.coding_promoter)
  )

  if (params.run_saturation) {
    S4_3_SATURATION(
      ch_s3_3_bed_and_peak,
    )
  }

  S4_4_BACKGROUND(
    ch_s3_3_bed_and_peak,
    file(params.promoter_file),
    nochrM_chrom_size
  )

  // Channel.fromPath("${projectDir}-test-data/*", type: "file")
  //   .collect()
  //   .view()

  //S3_3_PEAKS.out[0].map { 
  //  file = Channel.fromPath("${projectDir}-${it[0].id}/*", type: "file")
  //  tuple(it[0], file.view())
  //}.view()

  S1_1_TRIM.out[0].flatten().filter( Map ).view()
  S4_6_SUMMARIZE(
    // Meta
    S1_1_TRIM.out[0].flatten().filter( Map ),

    // Trimlog
    S1_1_TRIM.out.trimlog,

    // Outputs from methylQA
    S3_1_QA.out[0].map { it[1..-1].flatten().collect() },

    // Results from peak calling
    S3_3_PEAKS.out[0].map { it[1..-1].flatten().collect() },


    // JSON files
    S1_1_TRIM.out.json,
    S3_1_QA.out.json,
    S4_1_RUP.out.json,
    S4_2_ENRICH.out.json,
    S4_4_BACKGROUND.out.json
  )
}
