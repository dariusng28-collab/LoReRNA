// ============================================================
// SWISH — fishpond DTE / DTU / DGE analysis
//
// Runs after ALL per-sample OARFISH quantifications complete.
// Collects every *.quant + *.infreps.pq, then calls
// bin/swish_analysis.R to run three analyses:
//
//   DTE — Differential Transcript Expression
//   DTU — Differential Transcript Usage (isoform switching)
//   DGE — Differential Gene Expression (tx-level aggregation)
//
// Input staging
// ─────────────
// Nextflow stages collected path files flat in the work dir.
// The shell block reconstructs per-sample subdirs for swish_analysis.R:
//   quant_input/
//     <sample_id>/<sample_id>.quant
//     <sample_id>/<sample_id>.infreps.pq
//
// Conditions CSV (built in main.nf from meta):
//   sample_id,condition,pair,batch
//
// Paired / batch design is auto-detected by swish_analysis.R from
// non-empty 'pair' or 'batch' columns. No flags needed here.
//
// Changelog v1.2.0:
//   Removed dead pair_arg / batch_arg Groovy variables — these were
//   always empty strings and passed as blank arguments to Rscript.
//   swish_analysis.R auto-detects design from conditions CSV columns.
// ============================================================

process SWISH {

    tag "all_samples"
    label 'process_high'   // 8 CPUs, 32 GB — overridden by site config

    publishDir "${params.outdir}/06_swish", mode: params.publish_mode

    input:
    path quant_files        // collected: all *.quant files (staged flat)
    path infrep_files       // collected: all *.infreps.pq files (staged flat)
    path tx2gene            // tx2gene.tsv from PREPARE_OARFISH_REFERENCE
    path conditions_csv     // sample_id,condition[,pair,batch] CSV
    path swish_script       // swish_analysis.R staged from bin/ via channel

    output:
    path "results/*.csv",  emit: csv_results
    path "plots/*.pdf",    emit: plots
    path "logs/*.log",     emit: log

    script:
    """
    set -euo pipefail

    echo "========================================"
    echo " SWISH: fishpond DTE / DTU / DGE"
    echo " Conditions CSV : ${conditions_csv}"
    echo " tx2gene        : ${tx2gene}"
    echo " Condition A    : ${params.swish_condition_a}"
    echo " Condition B    : ${params.swish_condition_b}"
    echo " Design         : auto-detected from conditions CSV columns"
    echo " alpha          : ${params.swish_alpha}"
    echo " min_count      : ${params.swish_min_count}"
    echo " min_n          : ${params.swish_min_n}"
    echo " nperms         : ${params.swish_nperms}"
    echo "========================================"

    # ── Step 1: Reconstruct per-sample directory structure ────────────────
    mkdir -p quant_input

    for f in *.quant; do
        sample="\${f%.quant}"
        mkdir -p "quant_input/\${sample}"
        mv "\${f}" "quant_input/\${sample}/\${f}"
    done

    for f in *.infreps.pq; do
        sample="\${f%.infreps.pq}"
        mkdir -p "quant_input/\${sample}"
        mv "\${f}" "quant_input/\${sample}/\${f}"
    done

    echo "quant_input structure:"
    ls quant_input/

    # ── Step 2: Create output directories ────────────────────────────────
    mkdir -p results plots logs

    # ── Step 3: Run swish_analysis.R ─────────────────────────────────────
    # Paired/batch design is auto-detected from the conditions CSV:
    # a non-empty 'pair' column activates paired Wilcoxon;
    # a non-empty 'batch' column activates batch-aware swish.
    ${params.rscript_exe} "${swish_script}" \\
        --quant_dir   quant_input \\
        --tx2gene     "${tx2gene}" \\
        --conditions  "${conditions_csv}" \\
        --condition_a "${params.swish_condition_a}" \\
        --condition_b "${params.swish_condition_b}" \\
        --results_dir results \\
        --plots_dir   plots \\
        --logs_dir    logs \\
        --min_count   ${params.swish_min_count} \\
        --min_n       ${params.swish_min_n} \\
        --nperms      ${params.swish_nperms} \\
        --alpha       ${params.swish_alpha}

    # ── Step 4: Verify all expected outputs ──────────────────────────────
    for csv in results/DTE_full_results.csv \\
               results/DTU_full_results.csv \\
               results/DGE_full_results.csv \\
               results/DTE_all_significant.csv \\
               results/DTU_all_significant.csv \\
               results/DGE_all_significant.csv; do
        if [ ! -f "\${csv}" ]; then
            echo "ERROR: Expected output not found: \${csv}"
            exit 1
        fi
    done

    for pdf in plots/DTE_plots.pdf plots/DTU_plots.pdf plots/DGE_plots.pdf; do
        if [ ! -f "\${pdf}" ]; then
            echo "ERROR: Expected plot not found: \${pdf}"
            exit 1
        fi
    done

    echo "========================================"
    echo " SWISH complete"
    echo "  CSVs : \$(ls results/*.csv | wc -l) files"
    echo "  PDFs : \$(ls plots/*.pdf   | wc -l) files"
    echo "  Log  : \$(ls logs/*.log 2>/dev/null | head -1)"
    echo "========================================"
    """
}
