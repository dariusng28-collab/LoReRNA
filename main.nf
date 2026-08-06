nextflow.enable.dsl = 2

// ============================================================
// main.nf — LoReRNA entry point
//
// Validates parameters and the samplesheet, builds the input channel,
// then hands off to the LORERNA workflow in workflows/lorerna.nf.
//
// Samplesheet columns:
//   required : sample_name, condition, bam
//   optional : pair   (activates paired Wilcoxon in swish)
//              batch  (activates batch-aware swish)
// ============================================================

include { LORERNA } from './workflows/lorerna'

// Parameter + samplesheet validation (nextflow_schema.json, assets/schema_input.json)
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {

    // ── Parameter validation ──────────────────────────────────────────────
    // Types, required params, enums (seq_tech, filter groups, publish_mode) and
    // the samplesheet columns are all validated against nextflow_schema.json /
    // assets/schema_input.json by nf-schema.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Cross-parameter rule that JSON Schema expresses poorly: exactly one of
    // --reference_gtf / --reference_db must be given.
    if (!params.reference_gtf && !params.reference_db)
        error "Provide either --reference_gtf or --reference_db (not both)"
    if (params.reference_gtf && params.reference_db)
        error "Provide only one of --reference_gtf or --reference_db"

    // ── Container build scratch ───────────────────────────────────────────
    // Singularity builds an image by unpacking OCI layers into a scratch
    // directory, resolved as SINGULARITY_TMPDIR -> TMPDIR -> /tmp. Schedulers
    // routinely point TMPDIR at node-local scratch, which exists on compute
    // nodes but not on the login node where `nextflow run` performs the pull.
    // Left to fail on its own this surfaces minutes later as a Singularity
    // error about a "build parent dir", which says nothing about the cause.
    // Checked here so it fails immediately, with something actionable.
    if (workflow.containerEngine == 'singularity') {
        def sing_tmp = System.getenv('SINGULARITY_TMPDIR')
        def plain_tmp = System.getenv('TMPDIR')
        def scratch = sing_tmp ?: plain_tmp ?: '/tmp'
        def source = sing_tmp ? 'SINGULARITY_TMPDIR' : (plain_tmp ? 'TMPDIR' : 'the default')
        if (!file(scratch).exists()) {
            error """
            |Singularity image build scratch does not exist: ${scratch}
            |  (taken from ${source})
            |
            |Image layers are unpacked there before the .img is written, so every
            |container pull will fail before any analysis starts. This usually means
            |the cluster points TMPDIR at node-local scratch that exists on compute
            |nodes but not on this one.
            |
            |Fix: export SINGULARITY_TMPDIR to an existing directory with a few GB
            |free, then re-run. For example:
            |
            |    export SINGULARITY_TMPDIR=/tmp
            """.stripMargin()
        }
    }

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
