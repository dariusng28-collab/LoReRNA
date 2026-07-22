// ============================================================
// MULTIQC
//
// Runs once after all per-sample QC processes complete.
// Collects outputs from:
//   SAMTOOLS_STATS  → *_flagstat.txt, *_idxstats.txt, *_stats.txt
//   ISOQUANT_QC     → *_isoquant_mqc.tsv
//   MISER_QC_MERGE  → all_samples_rescue_metrics_mqc.tsv
//   SWISH           → logs/*_mqc.tsv  (summary stats for custom content)
//
// MultiQC auto-detects samtools files by filename pattern.
// IsoQuant, MisER, and SWISH files are parsed via multiqc_config.yml
// custom_data sections.
//
// Outputs:
//   multiqc_report.html      — main interactive HTML report
//   multiqc_report_data/     — underlying data files (TSV, JSON)
//   multiqc_plots/           — static plot images (PNG)
// ============================================================

process MULTIQC {

    tag "all_samples"
    label 'process_medium'

    publishDir "${params.outdir}/07_multiqc", mode: params.publish_mode

    input:
    path multiqc_files    // collected: all QC files to aggregate
    path multiqc_config   // multiqc_config.yml

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_report_data/",        emit: data
    path "multiqc_plots/",       emit: plots, optional: true

    script:
    def title   = "LoReRNA QC — ${params.swish_condition_b} vs ${params.swish_condition_a}"
    def comment = "Pipeline v${workflow.manifest.version} | Run: ${workflow.runName}"
    """
    set -euo pipefail

    echo "========================================"
    echo " MultiQC: aggregating pipeline QC"
    echo " Files staged: \$(ls | wc -l)"
    echo "========================================"

    ls -1

    multiqc . \\
        --config   "${multiqc_config}" \\
        --title    "${title}" \\
        --comment  "${comment}" \\
        --filename multiqc_report.html \\
        --outdir   . \\
        --force \\
        --data-dir \\
        --cl-config 'log_filesize_limit: 10000000'

    if [ ! -f "multiqc_report.html" ]; then
        echo "ERROR: multiqc_report.html not produced"
        exit 1
    fi

    echo "MultiQC report: \$(ls -lh multiqc_report.html | awk '{print \$5}')"
    echo "SUCCESS: multiqc_report.html"
    """

    stub:
    """
    mkdir -p multiqc_report_data multiqc_plots
    touch multiqc_report.html
    """
}
