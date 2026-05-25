// ============================================================
// SAMTOOLS_STATS
//
// Runs three samtools QC commands on the sorted BAM from
// SAMTOOLS_SORT_INDEX. Outputs are named for MultiQC auto-detection:
//
//   *_flagstat.txt  — alignment summary (reads mapped, paired, etc.)
//   *_idxstats.txt  — per-chromosome read counts
//   *_stats.txt     — detailed alignment statistics
//
// MultiQC natively parses all three without any custom config.
// The samtools section in the MultiQC report shows alignment rates,
// coverage summary, and duplicate rates for each sample side-by-side.
//
// This process runs in parallel with ISOQUANT and CLEAN_BAM immediately
// after SAMTOOLS_SORT_INDEX — it has no downstream dependencies so it
// does not add to the critical path.
// ============================================================

process SAMTOOLS_STATS {

    tag  "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/07_multiqc/samtools", mode: params.publish_mode

    input:
    tuple val(meta), path(sorted_bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_flagstat.txt"),  emit: flagstat
    tuple val(meta), path("${meta.id}_idxstats.txt"),  emit: idxstats
    tuple val(meta), path("${meta.id}_stats.txt"),     emit: stats

    script:
    """
    set -euo pipefail

    echo "========================================"
    echo " samtools QC: ${meta.id}"
    echo " Input BAM  : ${sorted_bam}"
    echo " Threads    : ${task.cpus}"
    echo "========================================"

    samtools flagstat -@ ${task.cpus} "${sorted_bam}" > "${meta.id}_flagstat.txt"
    samtools idxstats "${sorted_bam}"                  > "${meta.id}_idxstats.txt"
    samtools stats    -@ ${task.cpus} "${sorted_bam}" > "${meta.id}_stats.txt"

    # Quick summary to stdout
    echo "Flagstat summary:"
    grep 'mapped (' "${meta.id}_flagstat.txt" | head -2 || true
    """
}
