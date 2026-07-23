// ============================================================
// MISER_QC_MERGE
//
// Collects per-sample rescue_metrics.tsv files from MISER_QC
// and merges them into a single cohort-level summary TSV.
//
// Runs once after all MISER_QC jobs complete (.collect()).
//
// Outputs:
//   all_samples_rescue_metrics_mqc.tsv  — human-readable merged table (published)
//   miser_rescue_mqc.json               — MultiQC multi-dataset bargraph
//                                         (per-metric switcher, one scale per view)
//
// The TSV is published for users; MultiQC consumes the JSON (see main.nf), so
// the report shows one metric at a time on its own scale instead of a single
// bar crammed with 13 differently-scaled columns.
// ============================================================

process MISER_QC_MERGE {

    tag "all_samples"
    label 'process_single'

    publishDir "${params.outdir}/01_miser_qc", mode: params.publish_mode

    input:
    path metrics_files   // collected: all *.rescue_metrics.tsv files
    path script          // miser_rescue_multiqc.py — staged from bin/ via channel

    output:
    path "all_samples_rescue_metrics_mqc.tsv", emit: merged_metrics
    path "miser_rescue_mqc.json",              emit: multiqc_json

    script:
    """
    set -euo pipefail

    # Merge: keep the header from the first file, skip it in all subsequent files.
    # FNR==1 && NR!=1: first line of current file but not first line overall.
    awk 'FNR==1 && NR!=1 {next} {print}' *.rescue_metrics.tsv \\
        > all_samples_rescue_metrics_mqc.tsv

    # Build the MultiQC per-metric switcher bargraph from the merged table.
    python3 "${script}" \\
        --input  all_samples_rescue_metrics_mqc.tsv \\
        --output miser_rescue_mqc.json

    echo "Merged \$(ls *.rescue_metrics.tsv | wc -l) samples; wrote miser_rescue_mqc.json"
    echo "Rows (excluding header): \$(tail -n +2 all_samples_rescue_metrics_mqc.tsv | wc -l)"
    """

    stub:
    """
    touch all_samples_rescue_metrics_mqc.tsv miser_rescue_mqc.json
    """
}
