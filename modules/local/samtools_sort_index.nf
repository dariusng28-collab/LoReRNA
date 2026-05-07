// ============================================================
// SAMTOOLS_SORT_INDEX
//
// Replicates the pre-IsoQuant sort/index loop exactly:
//   samtools sort -@ <cpus> -m <mem_per_thread> -o <sample>.sorted.bam <input.bam>
//   samtools index <sample>.sorted.bam
//
// In Nextflow the miserbam.txt param file is not needed — the
// sorted BAM channel feeds directly into ISOQUANT via channel
// chaining. See main.nf for the MISER | SAMTOOLS_SORT_INDEX |
// ISOQUANT chain.
//
// Outputs both the sorted BAM and its .bai index as a single
// tuple so ISOQUANT receives them together and Nextflow can
// verify both exist before the next process starts.
//
// Changelog v1.1.1 (fix #6):
//   Label changed from process_medium (4 GB) to process_sort (16 GB).
//   ONT direct RNA BAMs routinely reach 10–50 GB on disk; samtools
//   sort with only 4 GB of address space exhausts its in-memory sort
//   budget, generates many temp files, and slows dramatically.
//   process_sort matches the -m 3G per-thread budget (4 threads × 3 G
//   = 12 G active + headroom for I/O buffers within the 16 GB limit).
// ============================================================

process SAMTOOLS_SORT_INDEX {

    tag "${meta.id}"
    label 'process_sort'    // cpus=4, memory=16 GB — see nextflow.config

    publishDir "${params.outdir}/02_sorted_bam", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.miser.sorted.bam"), path("${meta.id}.miser.sorted.bam.bai"), emit: bam_bai

    script:
    def sample_id = meta.id
    """
    set -euo pipefail

    echo "========================================"
    echo " Sort + Index: ${sample_id}"
    echo " Input BAM  : ${bam}"
    echo " Threads    : ${task.cpus}"
    echo "========================================"

    # -m 3G: allow each sort thread 3 GB of RAM for in-memory sort buffers.
    # 4 threads × 3 G = 12 G peak active + I/O headroom stays within 16 G label.
    echo "[1/2] Sorting BAM"
    samtools sort \\
        -@ ${task.cpus} \\
        -m 3G \\
        -o "${sample_id}.miser.sorted.bam" \\
        "${bam}"

    echo "[2/2] Indexing BAM"
    samtools index "${sample_id}.miser.sorted.bam"

    echo "SUCCESS: ${sample_id}.miser.sorted.bam + .bai"
    """
}
