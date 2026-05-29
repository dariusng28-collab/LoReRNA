nextflow.enable.dsl = 2

// ============================================================
// main.nf — LoReRNA v1.0.0
//
// Execution DAG (→ = feeds into):
//
//   samplesheet
//       │
//       ▼
//     MISER ──────────────────────────────────────────────────┐
//       │                                                      │ .miser.bam
//       ▼                                                      │
//   MISER_QC (per-sample)          ┌──────────────────────────▼──────────┐
//       │                          │       SAMTOOLS_SORT_INDEX            │
//       ▼                          └───┬──────────┬────────────┬──────────┘
//   MISER_QC_MERGE                     │          │            │
//                                      ▼          ▼            ▼
//                               SAMTOOLS   CLEAN_BAM      ISOQUANT
//                                _STATS    (parallel)    (per-sample)
//                                  │           │              │
//                                  │           │         ┌────┴──────────┐
//                                  │           │    ISOQUANT_QC    transcript
//                                  │           │                   _models
//                                  │           │                       │
//                                  │           │              PREPARE_OARFISH
//                                  │           │               _REFERENCE
//                                  │           │              (GTF+FASTA+tx2gene)
//                                  │           │                       │
//                                  │           └──────► OARFISH ◄──────┘
//                                  │                       │
//                                  │                    SWISH
//                                  │                  ┌────┴─────┐
//                                  │            SWISH_PLOTS      │
//                                  │                             │
//                                  └─────────────────► MULTIQC ◄─┘
//
// Samplesheet columns:
//   required : sample_name, condition, bam
//   optional : pair   (activates paired Wilcoxon in swish)
//              batch  (activates batch-aware swish)
// ============================================================

include { MISER                     } from './modules/local/miser'
include { MISER_QC                  } from './modules/local/miser_qc'
include { MISER_QC_MERGE            } from './modules/local/miser_qc_merge'
include { SAMTOOLS_SORT_INDEX       } from './modules/local/samtools_sort_index'
include { SAMTOOLS_STATS            } from './modules/local/samtools_stats'
include { CLEAN_BAM                 } from './modules/local/clean_bam'
include { ISOQUANT                  } from './modules/local/isoquant'
include { ISOQUANT_QC               } from './modules/local/isoquant_qc'
include { PREPARE_OARFISH_REFERENCE } from './modules/local/prepare_oarfish_ref'
include { OARFISH                   } from './modules/local/oarfish'
include { SWISH                     } from './modules/local/swish'
include { SWISH_PLOTS               } from './modules/local/swish_plots'
include { MULTIQC                   } from './modules/local/multiqc'

workflow LORERNA {

    take:
    ch_samplesheet

    main:

    // ── Reference channels ────────────────────────────────────────────────
    ch_genome_fasta        = Channel.value(file(params.reference_fasta,  checkIfExists: true))
    ch_annotation_bed      = Channel.value(file(params.reference_bed12,  checkIfExists: true))
    ch_isoquant_annotation = Channel.value(file(params.reference_db ?: params.reference_gtf, checkIfExists: true))

    // ── Staged bin scripts ────────────────────────────────────────────────
    ch_swish_script        = Channel.value(file("${projectDir}/bin/swish_analysis.R",   checkIfExists: true))
    ch_swish_plots_script  = Channel.value(file("${projectDir}/bin/swish_plots.R",      checkIfExists: true))
    ch_miser_qc_script     = Channel.value(file("${projectDir}/bin/miser_qc_summary.py", checkIfExists: true))
    ch_isoquant_qc_script  = Channel.value(file("${projectDir}/bin/isoquant_qc.py",     checkIfExists: true))
    ch_multiqc_config      = Channel.value(file("${projectDir}/containers/multiqc_config.yml", checkIfExists: true))

    // ── MISER ─────────────────────────────────────────────────────────────
    MISER(ch_samplesheet, ch_genome_fasta, ch_annotation_bed)

    // ── MISER QC ──────────────────────────────────────────────────────────
    MISER_QC(MISER.out.micro_exon_bed, ch_miser_qc_script)

    MISER_QC_MERGE(
        MISER_QC.out.metrics.map { meta, tsv -> tsv }.collect()
    )

    // ── Sort + Index ──────────────────────────────────────────────────────
    SAMTOOLS_SORT_INDEX(MISER.out.bam)

    // ── Samtools stats (for MultiQC) — runs in parallel ───────────────────
    SAMTOOLS_STATS(SAMTOOLS_SORT_INDEX.out.bam_bai)

    // ── BAM cleaning — per-sample, parallel, starts before IsoQuant ──────
    CLEAN_BAM(SAMTOOLS_SORT_INDEX.out.bam_bai)

    // ── IsoQuant ──────────────────────────────────────────────────────────
    ISOQUANT(SAMTOOLS_SORT_INDEX.out.bam_bai, ch_genome_fasta, ch_isoquant_annotation)

    // ── IsoQuant QC — per-sample, parallel (for MultiQC) ─────────────────
    ISOQUANT_QC(ISOQUANT.out.read_assignments, ch_isoquant_qc_script)

    // ── Prepare oarfish reference (collected — waits for all IsoQuant) ────
    ch_transcript_models = ISOQUANT.out.transcript_model_gtf
        .map { meta, gtf -> gtf }
        .collect()

    PREPARE_OARFISH_REFERENCE(ch_transcript_models, ch_genome_fasta)

    ch_merged_fa = PREPARE_OARFISH_REFERENCE.out.merged_fa.first()

    // ── Oarfish quantification ────────────────────────────────────────────
    OARFISH(CLEAN_BAM.out.clean_bam, ch_merged_fa)

    // ── Conditions CSV ────────────────────────────────────────────────────
    // FIX v1.0.0: header was "sample_name" — must be "sample_id" for R stopifnot
    ch_conditions_csv = OARFISH.out.quant
        .map { meta, quant ->
            def pair_val  = (meta.pair  != null && meta.pair  != '') ? meta.pair  : ''
            def batch_val = (meta.batch != null && meta.batch != '') ? meta.batch : ''
            "${meta.id},${meta.condition},${pair_val},${batch_val}\n"
        }
        .collectFile(
            name: 'conditions.csv',
            seed: "sample_id,condition,pair,batch\n",
            sort: true
        )

    ch_quant_files  = OARFISH.out.quant.map   { meta, q  -> q  }.collect()
    ch_infrep_files = OARFISH.out.infreps.map { meta, ir -> ir }.collect()

    // ── SWISH ─────────────────────────────────────────────────────────────
    SWISH(
        ch_quant_files,
        ch_infrep_files,
        PREPARE_OARFISH_REFERENCE.out.tx2gene,
        ch_conditions_csv,
        ch_swish_script
    )

    // ── Publication plots (ggplot2) ───────────────────────────────────────
    SWISH_PLOTS(
        SWISH.out.csv_results.flatten().collect(),
        ch_quant_files,
        PREPARE_OARFISH_REFERENCE.out.tx2gene,
        ch_conditions_csv,
        ch_swish_plots_script
    )

    // ── MultiQC — collect all QC outputs ─────────────────────────────────
    ch_multiqc_inputs = Channel.empty()
        .mix( SAMTOOLS_STATS.out.flagstat.map { meta, f -> f } )
        .mix( SAMTOOLS_STATS.out.idxstats.map { meta, f -> f } )
        .mix( SAMTOOLS_STATS.out.stats.map    { meta, f -> f } )
        .mix( ISOQUANT_QC.out.mqc_tsv.map     { meta, f -> f } )
        .mix( MISER_QC_MERGE.out.merged_metrics )
        .mix( SWISH.out.log.flatten() )
        .collect()

    MULTIQC(ch_multiqc_inputs, ch_multiqc_config)

    emit:
    miser_bam         = MISER.out.bam
    sorted_bam        = SAMTOOLS_SORT_INDEX.out.bam_bai
    clean_bam         = CLEAN_BAM.out.clean_bam
    transcript_models = ISOQUANT.out.transcript_model_gtf
    merged_fa         = PREPARE_OARFISH_REFERENCE.out.merged_fa
    tx2gene           = PREPARE_OARFISH_REFERENCE.out.tx2gene
    oarfish_quant     = OARFISH.out.quant
    swish_results     = SWISH.out.csv_results
    swish_plots       = SWISH_PLOTS.out.summary
    multiqc_report    = MULTIQC.out.report
}

workflow {

    // ── Required parameter validation ─────────────────────────────────────
    ['samplesheet', 'reference_fasta', 'reference_bed12'].each { p ->
        if (!params[p]) error "Missing required parameter: --${p}"
    }
    if (!params.reference_gtf && !params.reference_db)
        error "Provide either --reference_gtf or --reference_db (not both)"
    if (params.reference_gtf && params.reference_db)
        error "Provide only one of --reference_gtf or --reference_db"

    // ── Samplesheet parsing ───────────────────────────────────────────────
    ch_samplesheet = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row['sample_name'] ?: row['sample']
            def condition = row['condition']
            def bam_path  = row['bam']
            def pair_val  = row['pair']  ?: null
            def batch_val = row['batch'] ?: null

            if (!sample_id) error "Samplesheet missing 'sample_name' column. Found: ${row.keySet().join(', ')}"
            if (!condition) error "Row '${sample_id}' missing 'condition' column"
            if (!bam_path)  error "Row '${sample_id}' missing 'bam' column"
            if (!(sample_id ==~ /[A-Za-z0-9._-]+/))
                error "sample_name '${sample_id}' has unsupported characters — use letters, numbers, . _ -"
            if (pair_val && batch_val)
                log.warn "Sample '${sample_id}' has both pair and batch — pair takes precedence"

            tuple(
                [ id: sample_id, condition: condition,
                  pair: pair_val?.toString() ?: null,
                  batch: batch_val?.toString() ?: null ],
                file(bam_path, checkIfExists: true)
            )
        }

    LORERNA(ch_samplesheet)
}
