// ============================================================
// MISER — micro-exon rescue and splice-junction correction
//
// Replicates misertest_scratch.sh exactly:
//   Step 1  — samtools quickcheck on input BAM
//   Step 2  — node-local scratch at /scratch0/$USER via mktemp
//   Step 3  — 40 G free-space guard on scratch
//   Step 4  — MisER run (flags from params.miser_args)
//   Step 5  — explicit MisER exit-code check + scratch diagnostics
//   Step 6  — output BAM existence check
//   Step 7  — samtools quickcheck on output BAM
//   Cleanup — mv scratch → work dir, rm -rf scratch
//
// Changelog v1.1.1 (fix #2):
//   set -e is suspended around the MisER call with 'set +e … set -e'
//   so that miser_exit=$? is actually reachable. Previously, set -e
//   inside the { ... } group would abort immediately on a non-zero
//   exit, making step 5's check dead code — the failure was handled
//   by set -e but never logged. Now MisER failures are caught,
//   logged with scratch diagnostics, and the scratch dir is cleaned
//   up before the process exits non-zero.

process MISER {

    tag "${meta.id}"
    label 'process_miser'          // cpus=4, memory=64 GB, time=48 h — nextflow.config

    publishDir "${params.outdir}/01_miser", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path  genome_fasta
    path  annotation_bed12

    output:
    tuple val(meta), path("${meta.id}.miser.bam"), emit: bam
    path  "${meta.id}.missed_small.bed",            emit: micro_exon_bed
    path  "${meta.id}.miser.log",                   emit: log

    script:
    def sample_id = meta.id
    """
    # Synchronous log capture — no background subshell race.
    # '{ ... } 2>&1 | tee' is synchronous: bash waits for tee to drain
    # before exiting, so Nextflow always stages out a complete log.
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
    if ! samtools quickcheck "${bam}"; then
        echo "ERROR: Input BAM failed integrity check: ${bam}"
        exit 1
    fi

    # ── Step 2: Node-local scratch ─────────────────────────────────────────
    echo "[2/7] Setting up node-local scratch"
    SCRATCH_ROOT="${params.miser_scratch_root ?: ''}"
    if [ -z "\${SCRATCH_ROOT}" ]; then
        if [ -n "\${NXF_SCRATCH:-}" ] && [ -d "\${NXF_SCRATCH}" ]; then
            SCRATCH_ROOT="\${NXF_SCRATCH}"
        elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
            SCRATCH_ROOT="\${TMPDIR}"
        elif [ -d /scratch0 ]; then
            SCRATCH_ROOT="/scratch0/\${USER:-nf}"
        else
            SCRATCH_ROOT="\${PWD}/miser_scratch"
        fi
    fi
    mkdir -p "\${SCRATCH_ROOT}"
    tmp_dir=\$(mktemp -d "\${SCRATCH_ROOT%/}/miser_${sample_id}_XXXXXX")
    echo "Scratch directory : \${tmp_dir}"
    echo "Scratch space:"
    df -h "\${tmp_dir}"
    echo "Input BAM size    : \$(ls -lh ${bam} | awk '{print \$5}')"

    # ── Step 3: Guard — need >= 40 G free on scratch ────────────────────────
    echo "[3/7] Checking scratch space (need >= 40 G)"
    avail_gb=\$(df -BG "\${tmp_dir}" | awk 'NR==2 {print \$4}' | sed 's/G//')
    if [ "\${avail_gb}" -lt 40 ]; then
        echo "ERROR: Not enough scratch space (\${avail_gb} G free, need >= 40 G)"
        rm -rf "\${tmp_dir}"
        exit 1
    fi
    echo "Scratch space OK  : \${avail_gb} G available"

    # ── Step 4: Run MisER ──────────────────────────────────────────────────
    # Flags from params.miser_args (nextflow.config). -c from task.cpus (= 4).
    # Micro-exon BED written to NF work dir (not scratch) for direct staging.
    #
    # set +e / set -e: suspend errexit around MisER so miser_exit=$? is
    # reachable. Without this, set -e aborts the script before the
    # assignment and step 5 is dead code (fix #2).
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
        echo "Scratch contents at failure:"
        ls -lh "\${tmp_dir}" || true
        rm -rf "\${tmp_dir}"
        exit "\${miser_exit}"
    fi

    # ── Step 6: Check output BAM exists ────────────────────────────────────
    echo "[6/7] Checking output BAM exists"
    if [ ! -f "\${tmp_dir}/${sample_id}.miser.bam" ]; then
        echo "ERROR: Output BAM missing after MisER run"
        echo "Scratch contents:"
        ls -lh "\${tmp_dir}" || true
        rm -rf "\${tmp_dir}"
        exit 1
    fi
    echo "Output BAM size: \$(ls -lh \${tmp_dir}/${sample_id}.miser.bam | awk '{print \$5}')"

    # ── Step 7: Validate output BAM ────────────────────────────────────────
    echo "[7/7] samtools quickcheck on output BAM"
    if ! samtools quickcheck "\${tmp_dir}/${sample_id}.miser.bam"; then
        echo "ERROR: Output BAM failed integrity check"
        rm -rf "\${tmp_dir}"
        exit 1
    fi

    # ── Move from scratch → Nextflow work dir ──────────────────────────────
    mv "\${tmp_dir}/${sample_id}.miser.bam" "${sample_id}.miser.bam"

    # ── Cleanup scratch ─────────────────────────────────────────────────────
    rm -rf "\${tmp_dir}"
    echo "SUCCESS: ${sample_id}.miser.bam"

    } 2>&1 | tee "${sample_id}.miser.log"
    """
}
