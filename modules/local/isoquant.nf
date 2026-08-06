// ============================================================
// ISOQUANT
//
// Replicates the isoquant bash script exactly:
//   isoquant \
//     --reference  <genome.fa>          (REF)
//     --bam        <sample.sorted.bam>  (sorted_bam from SAMTOOLS_SORT_INDEX)
//     --genedb     <annotation.gtf|annotation.gtf.db>
//     --sqanti_output
//     --data_type  nanopore
//     --prefix     <sample_id>
//     -o           <sample_id>
//
// IMPORTANT:
//   IsoQuant's --genedb option accepts either a GTF/GFF annotation or a
//   pre-built gffutils database. This module receives whichever one the user
//   chose via --reference_gtf or --reference_db.
//
// Output layout:
//   With --prefix <sample_id> and -o <sample_id>, IsoQuant v3 writes all
//   results to:
//     <sample_id>/<sample_id>/<sample_id>.transcript_models.gtf
//     <sample_id>/<sample_id>/<sample_id>.transcript_counts.tsv
//     etc.
//   The 'OUT/OUT.*' layout only applies when --prefix is omitted.
//   publishDir then copies the whole directory to results/03_isoquant/.
//
//   Key emitted files for downstream steps:
//     *.transcript_counts.tsv   — transcript-level counts
//     *.gene_counts.tsv         — gene-level counts
//     *.read_assignments.tsv.gz — per-read transcript assignments
//
// Changelog v1.0.0 (fix #1 — critical):
//   Corrected transcript_models.gtf path from
//     ${sample_id}/OUT/OUT.transcript_models.gtf   (default --prefix layout)
//   to
//     ${sample_id}/${sample_id}/${sample_id}.transcript_models.gtf  (--prefix layout, IsoQuant v3 nested output)
//   The old path caused every ISOQUANT job to fail at the cp step.
// ============================================================

process ISOQUANT {

    tag "${meta.id}"
    label 'process_very_high'   // cpus=4, memory=64 GB — see conf/sge.config

    input:
    tuple val(meta), path(sorted_bam), path(bai)   // from SAMTOOLS_SORT_INDEX
    path  genome_fasta                              // reference_fasta param
    path  isoquant_annotation                       // reference_gtf or reference_db param

    output:
    tuple val(meta), path("${meta.id}/${meta.id}/*.transcript_counts.tsv"),    emit: transcript_counts
    tuple val(meta), path("${meta.id}/${meta.id}/*.gene_counts.tsv"),          emit: gene_counts
    tuple val(meta), path("${meta.id}/${meta.id}/*.read_assignments.tsv.gz"),  emit: read_assignments
    tuple val(meta), path("${meta.id}.transcript_models.gtf"),      emit: transcript_model_gtf

    script:
    def sample_id  = meta.id
    def extra_args = params.isoquant_args ?: ''
    def complete_genedb_arg = params.isoquant_complete_genedb ? '--complete_genedb' : ''
    """
    set -euo pipefail

    echo "========================================"
    echo " IsoQuant: ${sample_id}"
    echo " BAM      : ${sorted_bam}"
    echo " Genome   : ${genome_fasta}"
    echo " Annotation: ${isoquant_annotation}"
    echo " Threads  : ${task.cpus}"
    echo "========================================"

    # Validate sorted BAM and its index are present
    if [ ! -f "${sorted_bam}" ]; then
        echo "ERROR: sorted BAM not found: ${sorted_bam}"
        exit 1
    fi
    if [ ! -f "${bai}" ]; then
        echo "ERROR: BAM index not found: ${bai}"
        exit 1
    fi

    # Validate BAM integrity
    echo "Checking BAM integrity..."
    samtools quickcheck "${sorted_bam}"

    # --prefix ${sample_id} sets the output file basename.
    # -o ${sample_id}        sets the output directory.
    # With both set, IsoQuant v3 writes:
    #   ${sample_id}/${sample_id}/${sample_id}.transcript_models.gtf
    #   ${sample_id}/${sample_id}/${sample_id}.transcript_counts.tsv
    #   etc.
    # Do NOT omit --prefix — the default 'OUT/OUT.*' layout is not used here.
    ${params.isoquant_exe} \\
        --reference  "${genome_fasta}" \\
        --bam        "${sorted_bam}" \\
        --genedb     "${isoquant_annotation}" \\
        ${complete_genedb_arg} \\
        --sqanti_output \\
        --data_type  ${params.isoquant_data_type} \\
        --prefix     ${sample_id} \\
        --threads    ${task.cpus} \\
        -o           ${sample_id} \\
        ${extra_args}

    # With --prefix, IsoQuant v3 writes to ${sample_id}/${sample_id}/${sample_id}.transcript_models.gtf
    # (nested: -o sets outer dir, --prefix sets inner dir AND file basename)
    # (not OUT/OUT.transcript_models.gtf — that is the no-prefix default layout).
    if [ ! -f "${sample_id}/${sample_id}/${sample_id}.transcript_models.gtf" ]; then
        echo "ERROR: IsoQuant transcript model GTF not found."
        echo "Expected: ${sample_id}/${sample_id}/${sample_id}.transcript_models.gtf"
        echo "Output directory contents:"
        ls -lh "${sample_id}/" || true
        exit 1
    fi
    cp "${sample_id}/${sample_id}/${sample_id}.transcript_models.gtf" "${sample_id}.transcript_models.gtf"

    echo "SUCCESS: IsoQuant finished for ${sample_id}"
    """

    stub:
    """
    mkdir -p ${meta.id}/${meta.id}
    touch ${meta.id}/${meta.id}/${meta.id}.transcript_counts.tsv \\
          ${meta.id}/${meta.id}/${meta.id}.gene_counts.tsv \\
          ${meta.id}/${meta.id}/${meta.id}.read_assignments.tsv.gz \\
          ${meta.id}.transcript_models.gtf
    """
}
