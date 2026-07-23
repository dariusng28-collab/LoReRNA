// ============================================================
// ISOQUANT_QC
//
// Parses IsoQuant's read_assignments.tsv.gz to produce a
// per-sample assignment-type breakdown:
//
//   unique                  — reads assigned to exactly one isoform
//   unique_minor_difference — unique with minor alignment differences
//   ambiguous               — consistent with multiple isoforms
//   inconsistent            — mapped but inconsistent with annotation
//   not_assigned            — not assigned (includes noninformative)
//
// Outputs two files per sample:
//   {sample}_isoquant_qc.tsv   — human-readable stats TSV
//   {sample}_isoquant_mqc.tsv  — MultiQC custom bargraph input
//
// Runs in parallel with CLEAN_BAM immediately after ISOQUANT.
// The *_mqc.tsv files are collected by MULTIQC at the end.
// ============================================================

process ISOQUANT_QC {

    tag  "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(read_assignments_gz)
    path  script     // isoquant_qc.py staged via Channel.value() in main.nf

    output:
    tuple val(meta), path("${meta.id}_isoquant_qc.tsv"),  emit: qc_tsv
    tuple val(meta), path("${meta.id}_isoquant_mqc.tsv"), emit: mqc_tsv

    script:
    """
    set -euo pipefail

    python3 "${script}" \\
        --input  "${read_assignments_gz}" \\
        --sample "${meta.id}" \\
        --outdir .
    """

    stub:
    """
    touch ${meta.id}_isoquant_qc.tsv ${meta.id}_isoquant_mqc.tsv
    """
}
