// ============================================================
// MISER_QC_MERGE
//
// Collects per-sample rescue_metrics.tsv files from MISER_QC
// and merges them into a single cohort-level summary TSV.
//
// Runs once after all MISER_QC jobs complete (.collect()).
//
// Output:
//   all_samples_rescue_metrics.tsv
//
// Format: one header row + one data row per sample, tab-separated.
// Ready for direct import into R, Excel, or MultiQC (future).
// ============================================================

process MISER_QC_MERGE {

    tag "all_samples"
    label 'process_single'

    publishDir "${params.outdir}/01_miser_qc", mode: params.publish_mode

    input:
    path metrics_files   // collected: all *.rescue_metrics.tsv files

    output:
    path "all_samples_rescue_metrics.tsv", emit: merged_metrics

    script:
    """
    # Keep header from first file only, skip it in all subsequent files.
    # FNR==1 && NR!=1: first line of current file but not first line overall.
    awk 'FNR==1 && NR!=1 {next} {print}' *.rescue_metrics.tsv \\
        > all_samples_rescue_metrics.tsv

    echo "Merged \$(ls *.rescue_metrics.tsv | wc -l) samples into all_samples_rescue_metrics.tsv"
    echo "Rows (excluding header): \$(tail -n +2 all_samples_rescue_metrics.tsv | wc -l)"
    """
}
