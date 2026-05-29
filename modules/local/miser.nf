process MISER {

    tag "${meta.id}"
    label 'process_miser'

    publishDir "${params.outdir}/01_miser", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path  genome_fasta
    path  annotation_bed12

    output:
    tuple val(meta), path("${meta.id}.miser.bam"), emit: bam
    tuple val(meta), path("${meta.id}.missed_small.bed"), emit: micro_exon_bed
    path  "${meta.id}.miser.log",                   emit: log

    script:
    def sample_id    = meta.id
    def scratch_root = params.miser_scratch_root ?: '/scratch0'
    """
    {
    set -euo pipefail

    echo "============================================"
    echo " MisER: ${sample_id}"
    echo " Input BAM : ${bam}"
    echo " Genome    : ${genome_fasta}"
    echo " Annotation: ${annotation_bed12}"
    echo " Threads   : ${task.cpus}"
    echo "============================================"

    # ── Step 1: Validate input BAM ─────────────────────────────────────────
    echo "[1/7] samtools quickcheck on input BAM"
    samtools quickcheck "${bam}" || { echo "ERROR: BAM failed quickcheck"; exit 1; }

    # ── Step 2: Set up scratch ─────────────────────────────────────────────
    echo "[2/7] Setting up scratch"
    SCRATCH_BASE="${scratch_root}/\${USER:-nfuser}"
    mkdir -p "\${SCRATCH_BASE}"
    tmp_dir=\$(mktemp -d "\${SCRATCH_BASE}/miser_${sample_id}_XXXXXX")
    echo "Scratch directory : \${tmp_dir}"
    df -h "\${tmp_dir}"
    echo "Input BAM size    : \$(ls -lh ${bam} | awk '{print \$5}')"

    # ── Step 3: Check scratch space >= 40 G ────────────────────────────────
    echo "[3/7] Checking scratch space"
    avail_gb=\$(df -BG "\${tmp_dir}" | awk 'NR==2 {gsub(/G/,""); print \$4}')
    if [ "\${avail_gb}" -lt 40 ]; then
        echo "ERROR: Not enough scratch space (\${avail_gb} G free, need >= 40 G)"
        rm -rf "\${tmp_dir}"
        exit 1
    fi
    echo "Scratch space OK  : \${avail_gb} G available"

    # ── Step 4: Run MisER ──────────────────────────────────────────────────
    echo "[4/7] Running MisER"
    set +e
    ${params.miser_exe} \\
        -c ${task.cpus} \\
        ${params.miser_args} \\
        "${bam}" \\
        "${genome_fasta}" \\
        "${annotation_bed12}" \\
        "${sample_id}.missed_small.bed" \\
        -o "\${tmp_dir}/${sample_id}.miser.bam"
    miser_exit=\$?
    set -e
    echo "MisER exit code: \${miser_exit}"

    # ── Step 5: Check MisER exit code ──────────────────────────────────────
    echo "[5/7] Checking MisER exit code"
    if [ "\${miser_exit}" -ne 0 ]; then
        echo "ERROR: MisER failed (exit \${miser_exit})"
        ls -lh "\${tmp_dir}" || true
        rm -rf "\${tmp_dir}"
        exit "\${miser_exit}"
    fi

    # ── Step 6: Check output BAM exists ────────────────────────────────────
    echo "[6/7] Checking output BAM exists"
    if [ ! -f "\${tmp_dir}/${sample_id}.miser.bam" ]; then
        echo "ERROR: Output BAM missing"
        ls -lh "\${tmp_dir}" || true
        rm -rf "\${tmp_dir}"
        exit 1
    fi
    echo "Output BAM size: \$(ls -lh \${tmp_dir}/${sample_id}.miser.bam | awk '{print \$5}')"

    # ── Step 7: Validate output BAM ────────────────────────────────────────
    echo "[7/7] samtools quickcheck on output BAM"
    samtools quickcheck "\${tmp_dir}/${sample_id}.miser.bam" || {
        echo "ERROR: Output BAM failed integrity check"
        rm -rf "\${tmp_dir}"
        exit 1
    }

    mv "\${tmp_dir}/${sample_id}.miser.bam" "${sample_id}.miser.bam"
    rm -rf "\${tmp_dir}"
    echo "SUCCESS: ${sample_id}.miser.bam"

    } 2>&1 | tee "${sample_id}.miser.log"
    """
}
