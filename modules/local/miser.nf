process MISER {

    tag "${meta.id}"
    label 'process_miser'

    publishDir "${params.outdir}/01_miser", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path  genome_fasta
    path  annotation_bed12

    output:
    tuple val(meta), path("${meta.id}.miser.bam"),        emit: bam
    tuple val(meta), path("${meta.id}.missed_small.bed"), emit: micro_exon_bed
    path  "${meta.id}.miser.log",                         emit: log

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

    # ── Step 3: Check scratch space ─────────────────────────────────────────
    # Need space for: input BAM copy + output BAM + working files
    # Rule: need at least 2.5x the input BAM size free
    BAM_SIZE_GB=\$(du -BG "${bam}" | awk '{gsub(/G/,""); print \$1}')
    NEEDED_GB=\$(( BAM_SIZE_GB * 3 ))
    [ "\${NEEDED_GB}" -lt 40 ] && NEEDED_GB=40
    echo "[3/7] Scratch space check (need \${NEEDED_GB} G for \${BAM_SIZE_GB} G BAM)"
    avail_gb=\$(df -BG "\${tmp_dir}" | awk 'NR==2 {gsub(/G/,""); print \$4}')
    if [ "\${avail_gb}" -lt "\${NEEDED_GB}" ]; then
        echo "ERROR: Not enough scratch (\${avail_gb} G free, need \${NEEDED_GB} G)"
        rm -rf "\${tmp_dir}"
        exit 1
    fi
    echo "Scratch space OK  : \${avail_gb} G available"

    # ── Step 4: Copy BAM to local scratch ───────────────────────────────────
    # Root cause of previous timeouts: Nextflow stages BAM as a symlink to
    # /SAN. MisER does heavy random-access reads on it, saturating NFS I/O.
    # Copying to local scratch first matches the behaviour of the bash script
    # that ran successfully and removes the /SAN dependency during compute.
    echo "[4/7] Copying BAM to local scratch (removes NFS dependency)"
    cp "${bam}" "\${tmp_dir}/${sample_id}.input.bam"
    echo "BAM copy complete : \$(ls -lh \${tmp_dir}/${sample_id}.input.bam | awk '{print \$5}')"
    echo "Indexing local BAM..."
    samtools index "\${tmp_dir}/${sample_id}.input.bam"
    echo "Index done"

    # ── Step 5: Run MisER on local copy ────────────────────────────────────
    echo "[5/7] Running MisER on local scratch BAM"
    set +e
    ${params.miser_exe} \\
        -c ${task.cpus} \\
        ${params.miser_args} \\
        "\${tmp_dir}/${sample_id}.input.bam" \\
        "${genome_fasta}" \\
        "${annotation_bed12}" \\
        "${sample_id}.missed_small.bed" \\
        -o "\${tmp_dir}/${sample_id}.miser.bam"
    miser_exit=\$?
    set -e
    echo "MisER exit code: \${miser_exit}"

    # ── Step 6: Check MisER exit code ──────────────────────────────────────
    echo "[6/7] Checking MisER exit code"
    if [ "\${miser_exit}" -ne 0 ]; then
        echo "ERROR: MisER failed (exit \${miser_exit})"
        ls -lh "\${tmp_dir}" || true
        rm -rf "\${tmp_dir}"
        exit "\${miser_exit}"
    fi

    # ── Step 7: Validate and move output BAM ───────────────────────────────
    echo "[7/7] Validating output BAM"
    if [ ! -f "\${tmp_dir}/${sample_id}.miser.bam" ]; then
        echo "ERROR: Output BAM missing"
        ls -lh "\${tmp_dir}" || true
        rm -rf "\${tmp_dir}"
        exit 1
    fi
    samtools quickcheck "\${tmp_dir}/${sample_id}.miser.bam" || {
        echo "ERROR: Output BAM failed integrity check"
        rm -rf "\${tmp_dir}"
        exit 1
    }
    echo "Output BAM size: \$(ls -lh \${tmp_dir}/${sample_id}.miser.bam | awk '{print \$5}')"

    mv "\${tmp_dir}/${sample_id}.miser.bam" "${sample_id}.miser.bam"

    # Clean up scratch — input copy and index no longer needed
    rm -rf "\${tmp_dir}"
    echo "SUCCESS: ${sample_id}.miser.bam"

    } 2>&1 | tee "${sample_id}.miser.log"
    """

    stub:
    """
    touch ${meta.id}.miser.bam ${meta.id}.missed_small.bed ${meta.id}.miser.log
    """
}
