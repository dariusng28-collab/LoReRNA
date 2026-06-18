// ============================================================
// SWISH_PLOTS
//
// Generates publication-quality figures from swish_analysis.R results.
// Runs after SWISH completes, using:
//   - Full-results CSVs (*_full_results.csv) for MA and volcano plots
//     (shows all tested features, not just significant)
//   - Per-sample oarfish quant files for expression heatmaps
//   - tx2gene for gene name annotation
//   - conditions.csv for sample→condition mapping and heatmap annotation
//
// Output PDFs (one per analysis type + one summary):
//   DTE_publication_plots.pdf — p-value hist, MA, volcano, heatmap, direction bar
//   DTU_publication_plots.pdf — p-value hist, MA, volcano, heatmap, scatter
//   DGE_publication_plots.pdf — p-value hist, MA, volcano, heatmap, direction bar
//   summary_panel.pdf         — waterfall: up/down counts for DTE/DTU/DGE
//
// Packages required in lorerna container (environment.yml):
//   r-ggplot2, r-ggrepel, r-pheatmap, r-patchwork, r-viridis, r-rcolorbrewer
// ============================================================

process SWISH_PLOTS {

    tag "all_samples"
    label 'process_medium'

    publishDir "${params.outdir}/06_swish/publication_plots", saveAs: { it.replaceFirst("pub_plots/", "") }, mode: params.publish_mode

    input:
    path csv_results           // collected: all *.csv from SWISH results dir
    path quant_files           // collected: all *.quant files (flat-staged)
    path tx2gene
    path conditions_csv
    path script                // swish_plots.R staged via Channel.value()

    output:
    path "pub_plots/DTE_publication_plots.pdf", emit: dte_plots
    path "pub_plots/DTU_publication_plots.pdf", emit: dtu_plots
    path "pub_plots/DGE_publication_plots.pdf", emit: dge_plots
    path "pub_plots/summary_panel.pdf",         emit: summary

    script:
    """
    set -euo pipefail

    echo "========================================"
    echo " SWISH_PLOTS"
    echo " condition_a : ${params.swish_condition_a}"
    echo " condition_b : ${params.swish_condition_b}"
    echo " alpha       : ${params.swish_alpha}"
    echo " top_n       : ${params.swish_top_n}"
    echo "========================================"

    # ── Reconstruct results/ and quant_input/ from flat-staged files ─────────
    mkdir -p results quant_input pub_plots

    for f in *.csv; do
        [ "\${f}" = "conditions.csv" ] && continue  # keep in work dir for R
        mv "\${f}" results/
    done

    for f in *.quant; do
        sample="\${f%.quant}"
        mkdir -p "quant_input/\${sample}"
        mv "\${f}" "quant_input/\${sample}/\${f}"
    done

    ${params.rscript_exe} "${script}" \\
        --results_dir  results \\
        --quant_dir    quant_input \\
        --tx2gene      "${tx2gene}" \\
        --conditions   "${conditions_csv}" \\
        --condition_a  "${params.swish_condition_a}" \\
        --condition_b  "${params.swish_condition_b}" \\
        --plots_dir    pub_plots \\
        --alpha        ${params.swish_alpha} \\
        --top_n        ${params.swish_top_n}
    """
}
