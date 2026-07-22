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
// Changelog v1.0.0:
//   miser_qc_summary.py is now staged via an input channel rather
//   than accessed through ${projectDir}/bin/. With singularity.autoMounts
//   = false and explicit bind mounts, projectDir may not be visible
//   inside the container on some HPC configurations. Staging via channel
//   ensures the script is always available in the work directory.
//   The call is updated to use the staged path: python3 ${script}
// ============================================================

process MISER_QC {

    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/01_miser_qc/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(micro_exon_bed)
    path script   // miser_qc_summary.py — staged from bin/ via channel in main.nf

    output:
    tuple val(meta), path("*.rescued_microexons.tsv"), emit: rescued_exons_tsv
    tuple val(meta), path("*.rescued_microexons.bed"), emit: rescued_exons_bed
    tuple val(meta), path("*.rescue_metrics.tsv"),     emit: metrics

    script:
    def min_support_arg = params.miser_qc_min_support > 1 ? "--min-support ${params.miser_qc_min_support}" : ''
    def all_events_arg  = params.miser_qc_all_events       ? '--all-events'                                 : ''
    """
    python3 "${script}" \\
        --input  "${micro_exon_bed}" \\
        --outdir . \\
        ${min_support_arg} \\
        ${all_events_arg}
    """

    stub:
    """
    touch ${meta.id}.rescued_microexons.tsv ${meta.id}.rescued_microexons.bed ${meta.id}.rescue_metrics.tsv
    """
}
