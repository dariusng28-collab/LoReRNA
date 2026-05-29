// ============================================================
// OARFISH
//
// Replicates oarfish_quant.sh exactly (per sample):
//   oarfish \\
//     --threads    <task.cpus>
//     --reads      <sample.clean.bam>      (uBAM mode)
//     --reference  <merged_transcripts.fa>
//     --seq-tech   <params.seq_tech>       (ont-drna for direct RNA)
//     -o           <sample_id>
//     --filter-group no-filters
//     --model-coverage
//     --num-bootstraps <params.oarfish_num_bootstraps>
//     --verbose
//
// Why --reads (uBAM mode) not --alignments:
//   The clean BAM is a genome-aligned coordinate-sorted BAM from MisER.
//   oarfish does not support spliced genome alignments in --alignments mode.
//   In --reads mode, oarfish extracts raw read sequences from the BAM,
//   ignores genome coordinates entirely, and re-aligns internally to the
//   merged transcriptome FASTA using its built-in minimap2-rs engine.
//   This is the same approach validated in issue #54 of the oarfish repo.
//
// Why --reference not --index:
//   oarfish v0.6.5 does not support --index for reading a prebuilt index.
//   Only --index-out (write) exists in this version. Each sample rebuilds
//   the minimap2 index internally from the FASTA at runtime (~10 minutes).
//   This is a per-run cost with no effect on quantification correctness.
//
// Why --filter-group no-filters:
//   Multi-mapping reads are handled probabilistically by oarfish's EM
//   algorithm. Pre-filtering discards genuine reads that map to multiple
//   isoforms — exactly the reads most informative for isoform switching.
//   EM assigns fractional counts based on coverage model and prior.
//
// Why --model-coverage:
//   Direct RNA-seq reads have non-uniform positional coverage across
//   transcripts (3' bias). Coverage modelling corrects for this during
//   EM — improves quantification accuracy for long reads.
//
// Why --num-bootstraps 100:
//   fishpond/swish DTU testing requires inferential replicates to estimate
//   per-transcript count uncertainty. 100 bootstraps stored in .infreps.pq
//   provides stable variance estimates for the 3v3 design.
//
// Outputs per sample:
//   *.quant         — transcript-level point estimate counts
//   *.infreps.pq    — 100 bootstrap replicates in Parquet format (~43 MB)
//   *.meta_info.json — run metadata: valid_best_aln, discard_supp, params
//   *.ambig_info.tsv — unique vs ambiguously mapped read breakdown
// ============================================================

process OARFISH {

    tag "${meta.id}"
    label 'process_very_high'   // cpus=4, memory=64 GB, time=72 h — see conf/sge.config

    publishDir "${params.outdir}/05_oarfish/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(clean_bam)   // from CLEAN_BAM.out.clean_bam
    path  merged_fa                    // from PREPARE_OARFISH_REFERENCE.out.merged_fa

    output:
    tuple val(meta), path("${meta.id}.quant"),          emit: quant
    tuple val(meta), path("${meta.id}.infreps.pq"),     emit: infreps
    tuple val(meta), path("${meta.id}.meta_info.json"), emit: meta_info
    tuple val(meta), path("${meta.id}.ambig_info.tsv"), emit: ambig_info

    script:
    def sample_id   = meta.id
    def extra_args  = params.oarfish_args ?: ''
    """
    set -euo pipefail

    echo "========================================"
    echo " OARFISH: ${sample_id}"
    echo " BAM      : ${clean_bam}"
    echo " Reference: ${merged_fa}"
    echo " seq-tech : ${params.seq_tech}"
    echo " Bootstraps: ${params.oarfish_num_bootstraps}"
    echo " Threads  : ${task.cpus}"
    echo "========================================"

    # ── Verify inputs ────────────────────────────────────────────────────
    if [ ! -f "${clean_bam}" ]; then
        echo "ERROR: Clean BAM not found: ${clean_bam}"
        exit 1
    fi

    if [ ! -f "${merged_fa}" ]; then
        echo "ERROR: Reference FASTA not found: ${merged_fa}"
        exit 1
    fi

    # ── Run oarfish ──────────────────────────────────────────────────────
    # --reads      : uBAM mode — extracts sequences, ignores genome coords
    # --reference  : merged IsoQuant transcriptome (gffcompare-deduplicated)
    # --seq-tech   : ont-drna for direct RNA-seq (controls minimap2 preset)
    # --filter-group no-filters : EM handles multi-mappers probabilistically
    # --model-coverage          : corrects 3' bias in direct RNA coverage
    # --num-bootstraps          : inferential replicates for swish DTU

    ${params.oarfish_exe} \\
        --threads          ${task.cpus} \\
        --reads            "${clean_bam}" \\
        --reference        "${merged_fa}" \\
        --seq-tech         ${params.seq_tech} \\
        -o                 "${sample_id}" \\
        --filter-group     ${params.oarfish_filter_group} \\
        --num-bootstraps   ${params.oarfish_num_bootstraps} \\
        ${extra_args} \\
        --verbose

    # ── Verify all expected outputs exist ────────────────────────────────
    for f in "${sample_id}.quant" "${sample_id}.infreps.pq" \\
              "${sample_id}.meta_info.json" "${sample_id}.ambig_info.tsv"; do
        if [ ! -f "\${f}" ]; then
            echo "ERROR: Expected output not found: \${f}"
            exit 1
        fi
    done

    # ── Log key stats from meta_info ─────────────────────────────────────
    VALID=\$(grep valid_best_aln "${sample_id}.meta_info.json" | awk '{print \$2}' | tr -d ',')
    SUPP=\$(grep  discard_supp   "${sample_id}.meta_info.json" | awk '{print \$2}' | tr -d ',')
    echo ""
    echo "========================================"
    echo " OARFISH complete: ${sample_id}"
    echo "  valid_best_aln : \${VALID}"
    echo "  discard_supp   : \${SUPP} (secondary/supplementary — expected)"
    echo "  quant lines    : \$(wc -l < ${sample_id}.quant)"
    echo "  infreps size   : \$(ls -lh ${sample_id}.infreps.pq | awk '{print \$5}')"
    echo "========================================"
    """
}
