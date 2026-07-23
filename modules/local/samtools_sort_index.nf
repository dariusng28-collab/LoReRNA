// ============================================================
// SAMTOOLS_SORT_INDEX
//
// Sorts and indexes the MisER output BAM. The sorted BAM feeds
// both ISOQUANT and CLEAN_BAM downstream.
//
// Memory management:
//   samtools sort -m controls per-thread in-memory sort budget.
//   Previously hardcoded to 3G regardless of label; now computed
//   dynamically from task.memory and task.cpus so the same module
//   works correctly across different site resource configs.
//
//   Formula: floor(task.memory_MB * 0.80 / task.cpus) MB per thread
//   Uses 80% of allocated RAM to leave headroom for I/O buffers and
//   the samtools index step that follows.
//   Minimum: 512 MB (avoids degenerate configs with very low memory).
//
// Changelog v1.0.0:
//   Replaced hardcoded -m 3G with dynamic per-thread memory
//   calculation using shell arithmetic on Nextflow task variables.
// ============================================================

process SAMTOOLS_SORT_INDEX {

    tag "${meta.id}"
    label 'process_sort'    // cpus=4, memory=16 GB — see nextflow.config

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.miser.sorted.bam"), path("${meta.id}.miser.sorted.bam.bai"), emit: bam_bai

    script:
    def sample_id   = meta.id
    def total_mb    = task.memory.toMega() as int
    def n_threads   = task.cpus as int
    // 80% of total RAM shared across threads; minimum 512 MB per thread.
    // Shell integer arithmetic: (total_mb * 4/5) / n_threads, floored.
    """
    set -euo pipefail

    echo "========================================"
    echo " Sort + Index: ${sample_id}"
    echo " Input BAM  : ${bam}"
    echo " Threads    : ${n_threads}"
    echo " Memory     : ${total_mb} MB allocated"
    echo "========================================"

    # Per-thread sort memory: 80% of allocated RAM ÷ threads, floor to nearest MB.
    # Minimum 512 M so degenerate (low-memory) configs don't produce a bad -m flag.
    MEM_PER_THREAD_MB=\$(( ${total_mb} * 4 / 5 / ${n_threads} ))
    [ "\${MEM_PER_THREAD_MB}" -lt 512 ] && MEM_PER_THREAD_MB=512
    echo "Sort memory : \${MEM_PER_THREAD_MB} M per thread (${n_threads} threads)"

    echo "[1/2] Sorting BAM"
    samtools sort \\
        -@ ${n_threads} \\
        -m "\${MEM_PER_THREAD_MB}M" \\
        -o "${sample_id}.miser.sorted.bam" \\
        "${bam}"

    echo "[2/2] Indexing BAM"
    samtools index "${sample_id}.miser.sorted.bam"

    echo "SUCCESS: ${sample_id}.miser.sorted.bam + .bai"
    echo "Output size: \$(ls -lh ${sample_id}.miser.sorted.bam | awk '{print \$5}')"
    """

    stub:
    """
    touch ${meta.id}.miser.sorted.bam ${meta.id}.miser.sorted.bam.bai
    """
}
