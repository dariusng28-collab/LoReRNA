// ============================================================
// CLEAN_BAM — strip duplicate ts tags and filter secondary/supplementary reads
//
// Previously a sequential for-loop inside PREPARE_OARFISH_REFERENCE.
// Running it as a per-sample process means all samples are cleaned in
// parallel and the process can start as soon as SAMTOOLS_SORT_INDEX
// finishes — it no longer has to wait for IsoQuant + gffcompare.
//
// Why -F 2304:
//   Exclude secondary alignments (flag 256) and supplementary alignments
//   (flag 2048). These are artefacts of the genome aligner; oarfish
//   re-aligns in uBAM mode and does not use them.
//
// Why -x ts:
//   MisER adds a ts (strand) tag on top of the ts tag already written by
//   minimap2. oarfish's strict BAM parser rejects duplicate tags with
//   InvalidData(DuplicateTag(Tag("ts"))). Removing the tag is safe
//   because oarfish does not use ts — it re-aligns internally.
//
// Changelog v1.2.0:
//   Extracted from PREPARE_OARFISH_REFERENCE step 5 into its own process.
//   Adds -@ threads to both samtools commands.
// ============================================================

process CLEAN_BAM {

    tag "${meta.id}"
    label 'process_sort'    // cpus=4, memory=16 GB — samtools view is I/O-bound

    publishDir "${params.outdir}/02_sorted_bam", mode: params.publish_mode,
        pattern: "*.clean.bam"

    input:
    tuple val(meta), path(sorted_bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.clean.bam"), emit: clean_bam

    script:
    def sample_id = meta.id
    """
    set -euo pipefail

    echo "========================================"
    echo " CLEAN_BAM: ${sample_id}"
    echo " Input     : ${sorted_bam}"
    echo " Threads   : ${task.cpus}"
    echo " Filter    : -F 2304 (exclude secondary + supplementary)"
    echo " Tag strip : -x ts   (remove duplicate MisER ts tag)"
    echo "========================================"

    if ! samtools quickcheck "${sorted_bam}"; then
        echo "ERROR: Input BAM failed integrity check: ${sorted_bam}"
        exit 1
    fi

    # -F 2304 = exclude secondary (256) + supplementary (2048)
    # -x ts   = remove ts strand tag (written twice: minimap2 + MisER)
    # -@ splits threads between the two samtools processes in the pipe
    samtools view \\
        -h -F 2304 -x ts \\
        -@ ${task.cpus} \\
        "${sorted_bam}" \\
    | samtools view \\
        -b \\
        -@ 1 \\
        -o "${sample_id}.clean.bam"

    if [ ! -f "${sample_id}.clean.bam" ]; then
        echo "ERROR: Output clean BAM not produced"
        exit 1
    fi

    if ! samtools quickcheck "${sample_id}.clean.bam"; then
        echo "ERROR: Output clean BAM failed integrity check"
        exit 1
    fi

    READS_IN=\$(samtools view -c "${sorted_bam}")
    READS_OUT=\$(samtools view -c "${sample_id}.clean.bam")
    echo "Reads in  : \${READS_IN}"
    echo "Reads out : \${READS_OUT} (after -F 2304 filter)"
    echo "Clean BAM : \$(ls -lh ${sample_id}.clean.bam | awk '{print \$5}')"
    echo "SUCCESS: ${sample_id}.clean.bam"
    """
}
