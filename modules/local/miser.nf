// ============================================================
// MISER — micro-exon rescue and splice-junction correction
//
// Replicates misertest_scratch.sh exactly:
//   Step 1  — samtools quickcheck on input BAM
//   Step 2  — node-local scratch resolution
//   Step 3  — 40 G free-space guard on scratch
//   Step 4  — MisER run (flags from params.miser_args)
//   Step 5  — explicit MisER exit-code check + scratch diagnostics
//   Step 6  — output BAM existence check
//   Step 7  — samtools quickcheck on output BAM
//   Cleanup — mv scratch → work dir, rm -rf scratch
//
// Scratch resolution order (scheduler-agnostic, most-preferred first):
//   1. params.miser_scratch_root  — explicit override in site config
//   2. /scratch0/$USER            — SGE/UGE node-local (tscratch-allocated)
//   3. $SLURM_TMPDIR              — SLURM --tmp allocation (most HPC sites)
//   4. /scratch/$USER             — generic HPC scratch mount
//   5. /local_scratch/$USER       — alternative generic mount
//   6. $TMPDIR                    — PBS/Torque/LSF fallback
//   7. $PWD/miser_scratch         — shared filesystem last resort (always succeeds)
//
// Each candidate is only accepted if the directory EXISTS and has >= 40 G free.
// If a candidate exists but lacks space, the resolver moves to the next option
// rather than hard-failing — this is the key difference from v1.1.x.
//
// Changelog v1.2.0:
//   fix #8  — scratch fallback order corrected: /scratch0 is now tried before
//             $TMPDIR. Previously, scratch = true in the site config set
//             $TMPDIR → /tmp (a small partition), which was accepted before
//             /scratch0 was ever checked, causing the 40 G guard to always fail.
//   fix #8b — 40 G free-space check is now inside the resolver (_try_scratch),
//             not only in step 3. A candidate with insufficient space is skipped
//             automatically instead of aborting the whole job.
//
// Changelog v1.1.1 (fix #2):
//   set -e is suspended around the MisER call so miser_exit=$? is reachable.

process MISER {

    tag "${meta.id}"
    label 'process_miser'   // cpus=4, memory=64 GB, time=48 h — see nextflow.config

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
    # Resolve the best available scratch directory.
    # _try_scratch accepts a candidate path; returns 0 (success) only if
    # the directory EXISTS and has at least 40 G free — then sets SCRATCH_ROOT.
    # This means a full but present /scratch0 does not block the fallback chain.
    echo "[2/7] Setting up node-local scratch"

    SCRATCH_ROOT="${params.miser_scratch_root ?: ''}"

    if [ -z "\${SCRATCH_ROOT}" ]; then

        _try_scratch() {
            local candidate="\$1"
            [ -n "\${candidate}" ] || return 1
            [ -d "\${candidate}" ] || { mkdir -p "\${candidate}" 2>/dev/null || return 1; }
            local free_gb
            free_gb=\$(df -BG "\${candidate}" 2>/dev/null | awk 'NR==2{gsub(/G/,""); print \$4}')
            [ "\${free_gb:-0}" -ge 40 ] || return 1
            SCRATCH_ROOT="\${candidate}"
            return 0
        }

        # Preference order — most-local scratch first, shared filesystem last.
        # /scratch0 is tried first since tscratch=60G allocates there.
        ||  # SGE/UGE node-local (tscratch-allocated)
        _try_scratch "\${SLURM_TMPDIR:-}"               ||  # SLURM    (--tmp allocation)
        _try_scratch "/scratch/\${USER:-nf}"             ||  # generic HPC scratch
        _try_scratch "/local_scratch/\${USER:-nf}"       ||  # alternate generic
        _try_scratch "\${TMPDIR:-}"                      ||  # PBS/Torque/LSF $TMPDIR
        SCRATCH_ROOT="${PWD}/miser_scratch"                 # shared filesystem last resort

        echo "Scratch resolver  : \${SCRATCH_ROOT}"
    fi

    mkdir -p "\${SCRATCH_ROOT}"
    tmp_dir=\$(mktemp -d "\${SCRATCH_ROOT%/}/miser_${sample_id}_XXXXXX")
    echo "Scratch directory : \${tmp_dir}"
    echo "Scratch space:"
    df -h "\${tmp_dir}"
    echo "Input BAM size    : \$(ls -lh ${bam} | awk '{print \$5}')"

    # ── Step 3: Guard — verify scratch meets the 40 G floor ────────────────
    # The resolver already checks 40 G for scheduler-provided paths; this step
    # catches the PWD/NFS fallback and any externally set miser_scratch_root.
    echo "[3/7] Checking scratch space (need >= 40 G)"
    avail_gb=\$(df -BG "\${tmp_dir}" | awk 'NR==2 {gsub(/G/,""); print \$4}')
    if [ "\${avail_gb}" -lt 40 ]; then
        echo "ERROR: Not enough scratch space (\${avail_gb} G free, need >= 40 G)"
        echo "Resolved scratch root: \${SCRATCH_ROOT}"
        echo "To override, set params.miser_scratch_root in your site config."
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
        echo "Scratch contents at failure:"
        ls -lh "\${tmp_dir}" || true
        rm -rf "\${tmp_dir}"
        exit "\${miser_exit}"
    fi

    # ── Step 6: Check output BAM exists ────────────────────────────────────
    echo "[6/7] Checking output BAM exists"
    if [ ! -f "\${tmp_dir}/${sample_id}.miser.bam" ]; then
        echo "ERROR: Output BAM missing after MisER run"
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
