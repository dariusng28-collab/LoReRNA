// ============================================================
// MISER_QC
//
// Parses the MisER missed_small.bed event table and produces
// three per-sample outputs:
//
//   <sample>.rescued_microexons.tsv  — unique rescued exon catalog
//   <sample>.rescued_microexons.bed  — IGV-compatible BED6 track
//   <sample>.rescue_metrics.tsv      — per-sample QC metrics
//
// Runs immediately after MISER, one job per sample.
// Lightweight — Python + pandas only, no alignment steps.
//
// The script is called through python3 so the run does not depend on
// executable bits surviving checkout/staging.
// ============================================================

process MISER_QC {

    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/01_miser_qc/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(micro_exon_bed)

    output:
    tuple val(meta), path("*.rescued_microexons.tsv"), emit: rescued_exons_tsv
    tuple val(meta), path("*.rescued_microexons.bed"), emit: rescued_exons_bed
    tuple val(meta), path("*.rescue_metrics.tsv"),     emit: metrics

    script:
    def min_support_arg = params.miser_qc_min_support > 1 ? "--min-support ${params.miser_qc_min_support}" : ''
    def all_events_arg  = params.miser_qc_all_events       ? '--all-events'                                 : ''
    """
    python3 "${projectDir}/bin/miser_qc_summary.py" \\
        --input  "${micro_exon_bed}" \\
        --outdir . \\
        ${min_support_arg} \\
        ${all_events_arg}
    """
}
