nextflow.enable.dsl = 2

// ============================================================
// main.nf - LoReRNA v1.1.2
//
// Long-read RNA analysis:
// samplesheet -> MISER -> sort/index -> IsoQuant -> oarfish reference
//             -> Oarfish quantification -> fishpond/swish analysis
//
// Samplesheet format (3 columns, always):
//   sample_name,condition,bam
//   Control1,Control,/path/to/Control1.bam
//   KD1,KD,/path/to/KD1.bam
//
// Analysis design is controlled by flags, not samplesheet columns:
//   --swish_paired true   paired Wilcoxon test (row order defines pairing)
//   --swish_batch  true   batch-corrected test
// ============================================================

include { MISER                     } from './modules/local/miser'
include { SAMTOOLS_SORT_INDEX       } from './modules/local/samtools_sort_index'
include { ISOQUANT                  } from './modules/local/isoquant'
include { PREPARE_OARFISH_REFERENCE } from './modules/local/prepare_oarfish_ref'
include { OARFISH                   } from './modules/local/oarfish'
include { SWISH                     } from './modules/local/swish'

workflow LORERNA {

    take:
    ch_samplesheet   // tuple val(meta), path(bam)
                     // meta.id        — sample_name
                     // meta.condition — condition

    main:

    ch_genome_fasta        = Channel.value( file(params.reference_fasta, checkIfExists: true) )
    ch_annotation_bed      = Channel.value( file(params.reference_bed12, checkIfExists: true) )
    ch_isoquant_annotation = Channel.value( file(params.reference_db ?: params.reference_gtf, checkIfExists: true) )
    ch_swish_script        = Channel.value( file("${projectDir}/bin/swish_analysis.R", checkIfExists: true) )

    MISER (
        ch_samplesheet,
        ch_genome_fasta,
        ch_annotation_bed
    )

    SAMTOOLS_SORT_INDEX (
        MISER.out.bam
    )

    ISOQUANT (
        SAMTOOLS_SORT_INDEX.out.bam_bai,
        ch_genome_fasta,
        ch_isoquant_annotation
    )

    ch_transcript_models = ISOQUANT.out.transcript_model_gtf
        .map { meta, gtf -> gtf }
        .collect()

    ch_sorted_bams_collected = SAMTOOLS_SORT_INDEX.out.bam_bai
        .map { meta, bam, bai -> tuple(meta, bam) }
        .collect()
        .map { pairs ->
            tuple(
                pairs.collect { it[0] },
                pairs.collect { it[1] }
            )
        }

    PREPARE_OARFISH_REFERENCE (
        ch_transcript_models,
        ch_sorted_bams_collected,
        ch_genome_fasta
    )

    ch_merged_fa = PREPARE_OARFISH_REFERENCE.out.merged_fa.first()

    ch_clean_bams = PREPARE_OARFISH_REFERENCE.out.clean_bams
        .flatMap { metas, clean_bams ->
            def clean_bam_list = clean_bams instanceof List ? clean_bams : [clean_bams]
            def clean_bams_by_sample = clean_bam_list.collectEntries { bam ->
                def sample_id = bam.getName().replaceFirst(/\.clean\.bam$/, '')
                [(sample_id): bam]
            }
            metas.collect { meta ->
                def clean_bam = clean_bams_by_sample[meta.id]
                if (clean_bam == null) {
                    error "PREPARE_OARFISH_REFERENCE did not produce a clean BAM for sample '${meta.id}'"
                }
                tuple(meta, clean_bam)
            }
        }

    OARFISH (
        ch_clean_bams,
        ch_merged_fa
    )

    ch_quant_files = OARFISH.out.quant
        .map { meta, quant -> quant }
        .collect()

    ch_infrep_files = OARFISH.out.infreps
        .map { meta, infrep -> infrep }
        .collect()

    // ── Conditions CSV ────────────────────────────────────────────────────
    // Two columns: sample_name, condition.
    // sort: true → lexicographic ordering → reproducible column order
    // across parallel runs regardless of emission order.
    ch_conditions_csv = OARFISH.out.quant
        .map { meta, quant ->
            "${meta.id},${meta.condition}\n"
        }
        .collectFile(
            name: 'conditions.csv',
            seed: "sample_name,condition\n",
            sort: true
        )

    SWISH (
        ch_quant_files,
        ch_infrep_files,
        PREPARE_OARFISH_REFERENCE.out.tx2gene,
        ch_conditions_csv,
        ch_swish_script
    )

    emit:
    miser_bam         = MISER.out.bam
    sorted_bam        = SAMTOOLS_SORT_INDEX.out.bam_bai
    transcript_models = ISOQUANT.out.transcript_model_gtf
    transcript_counts = ISOQUANT.out.transcript_counts
    gene_counts       = ISOQUANT.out.gene_counts
    merged_fa         = PREPARE_OARFISH_REFERENCE.out.merged_fa
    tx2gene           = PREPARE_OARFISH_REFERENCE.out.tx2gene
    oarfish_quant     = OARFISH.out.quant
    oarfish_infreps   = OARFISH.out.infreps
    oarfish_meta      = OARFISH.out.meta_info
    swish_results     = SWISH.out.csv_results
    swish_plots       = SWISH.out.plots
    swish_log         = SWISH.out.log
}

workflow {

    // ── Required parameter validation ─────────────────────────────────────
    ['samplesheet', 'reference_fasta', 'reference_bed12'].each { p ->
        if (!params[p]) error "Missing required parameter: --${p}"
    }
    if (!params.reference_gtf && !params.reference_db) {
        error "Missing required IsoQuant annotation: provide either --reference_gtf or --reference_db"
    }
    if (params.reference_gtf && params.reference_db) {
        error "Provide only one of --reference_gtf or --reference_db, not both"
    }
    if (params.swish_paired && params.swish_batch) {
        error "--swish_paired and --swish_batch are mutually exclusive in fishpond — choose one"
    }

    // ── Samplesheet parsing ───────────────────────────────────────────────
    // Exactly three columns: sample_name, condition, bam.
    // Extra columns are silently ignored.
    ch_samplesheet = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->

            def sample_id = row['sample_name']
            def condition = row['condition']
            def bam_path  = row['bam']

            if (!sample_id) {
                error "Samplesheet is missing 'sample_name' column.\nFound: ${row.keySet().join(', ')}"
            }
            if (!condition) {
                error "Samplesheet row for '${sample_id}' is missing 'condition' column.\nFound: ${row.keySet().join(', ')}"
            }
            if (!bam_path) {
                error "Samplesheet row for '${sample_id}' is missing 'bam' column.\nFound: ${row.keySet().join(', ')}"
            }
            if (!(sample_id ==~ /[A-Za-z0-9._-]+/)) {
                error "sample_name '${sample_id}' contains unsupported characters — use letters, numbers, dots, underscores or hyphens"
            }

            def meta = [
                id:        sample_id.toString(),
                condition: condition.toString()
            ]

            tuple(meta, file(bam_path, checkIfExists: true))
        }

    LORERNA(ch_samplesheet)
}
