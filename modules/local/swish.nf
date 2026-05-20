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
// Nextflow stages collected path files flat in the work dir:
//   Control1.quant, Control2.quant, ..., STM2457_72h.quant
//   Control1.infreps.pq, ...
//
// The shell block reconstructs the per-sample subdirectory
// structure that swish_analysis.R expects before calling R:
//   quant_input/
//     Control1/Control1.quant + Control1.infreps.pq
//     Control2/Control2.quant + Control2.infreps.pq
//     ...
//
// Conditions CSV (built in main.nf from meta):
//   sample_id,condition,pair,batch
//   Control1,Control,,
//   Knockdown1,KD,,
//   STM2457_72h,STM,,
//
// pair / batch columns are empty strings when not used; the R
// Paired/batch design is auto-detected inside swish_analysis.R from the
// conditions CSV column content — no flags needed.
//
// Outputs
// ───────
//   results/    9 CSVs — all_sig / upregulated / downregulated
//               for DTE, DTU (increased/decreased usage), DGE
//   plots/      3 PDFs — DTE, DTU, DGE diagnostic plots
//   logs/       1 timestamped run log from R
//
// Resource label: process_high (nextflow.config: 8 CPUs, 32 GB, 48 h | UCL: 4 CPUs, 32 GB, 48 h)
//   R is single-threaded; 32 GB is needed for 7-sample ×
//   100-bootstrap infRep matrices on a human transcriptome.
//   Memory doubles on retry via check_max() in nextflow.config.
// ============================================================

process SWISH {

    tag "all_samples"
    label 'process_high'   // nextflow.config: 8 CPUs, 32 GB | UCL: 4 CPUs, 32 GB

    publishDir "${params.outdir}/06_swish", mode: params.publish_mode

    input:
    path quant_files        // collected: all *.quant files (staged flat)
    path infrep_files       // collected: all *.infreps.pq files (staged flat)
    path tx2gene            // tx2gene.tsv from PREPARE_OARFISH_REFERENCE
    path conditions_csv     // sample_id,condition,pair,batch CSV — generated in main.nf

    path swish_script

    output:
    path "results/*.csv",  emit: csv_results
    path "plots/*.pdf",    emit: plots
    path "logs/*.log",     emit: log

    script:
    // No flags needed — swish_analysis.R auto-detects paired/batch design
    // from the conditions.csv content (non-empty 'pair'/'batch' columns).
    def pair_arg  = ''
    def batch_arg = ''
    """
    set -euo pipefail

    echo "========================================"
    echo " SWISH: fishpond DTE / DTU / DGE"
    echo " Conditions CSV : ${conditions_csv}"
    echo " tx2gene        : ${tx2gene}"
    echo " Condition A    : ${params.swish_condition_a}"
    echo " Condition B    : ${params.swish_condition_b}"
    echo " Paired design  : auto-detected from 'pair' column in conditions CSV"
    echo " Batch variable : auto-detected from 'batch' column in conditions CSV"
    echo " alpha          : ${params.swish_alpha}"
    echo " min_count      : ${params.swish_min_count}"
    echo " min_n          : ${params.swish_min_n}"
    echo " nperms         : ${params.swish_nperms}"
    echo "========================================"

    # ── Step 1: Reconstruct per-sample directory structure ────────────────
    # swish_analysis.R expects:
    #   quant_input/<sample_id>/<sample_id>.quant
    #   quant_input/<sample_id>/<sample_id>.infreps.pq
    # Nextflow stages collected files flat; rebuild subdirs here.

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
    # Paired/batch design auto-detected from conditions CSV column content.
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
        --alpha       ${params.swish_alpha} \\
        ${pair_arg} \\
        ${batch_arg}

    # ── Step 4: Verify all expected outputs ──────────────────────────────
    echo "Verifying outputs..."

    for csv in results/DTE_all_significant.csv \\
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
