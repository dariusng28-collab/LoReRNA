#!/usr/bin/env bash
# =============================================================================
# validate_preflight.sh — LoReRNA v1.0.0 pre-flight checks
#
# Runs on the login node. No cluster jobs submitted. ~5 minutes.
# Checks every fix in-place by inspecting the .nf files directly.
#
# Usage:
#   cd /path/to/lorerna
#   bash test/validate_preflight.sh
#
# Exits 1 if any check fails.
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PASS=0; FAIL=0

_pass() { echo "  PASS  $1"; PASS=$(( PASS + 1 )); }
_fail() { echo "  FAIL  $1"; FAIL=$(( FAIL + 1 )); }
_sec()  { echo; echo "── $1 ──────────────────────────────────────────────────"; }

echo "============================================================"
echo " LoReRNA v1.0.0 pre-flight validation"
echo " Pipeline : ${PIPELINE_DIR}"
echo " Date     : $(date)"
echo "============================================================"

# ── 1. Conditions CSV header ─────────────────────────────────────────────────
_sec "Fix 1 — conditions CSV header (sample_id not sample_name)"

SEED=$(grep 'seed:' "${PIPELINE_DIR}/main.nf" | head -1)
if echo "${SEED}" | grep -q '"sample_id,condition'; then
    _pass "conditions CSV seed uses 'sample_id' — R stopifnot will pass"
else
    _fail "conditions CSV seed does NOT use 'sample_id' — SWISH will fail on every run"
    echo "        Found : ${SEED}"
    echo "        Needed: seed: \"sample_id,condition,pair,batch\""
fi

if grep -q 'sample_id' "${PIPELINE_DIR}/bin/swish_analysis.R" 2>/dev/null; then
    _pass "swish_analysis.R references 'sample_id' — consistent with CSV seed"
else
    _fail "swish_analysis.R does not reference 'sample_id' — check bin/swish_analysis.R"
fi

# ── 2. environment_miser.yml samtools version ─────────────────────────────────
_sec "Fix 2 — environment_miser.yml samtools version string"

SAMTOOLS_LINE=$(grep 'samtools' "${PIPELINE_DIR}/containers/environment_miser.yml" | head -1)
if echo "${SAMTOOLS_LINE}" | grep -qE 'samtools=[0-9]+\.[0-9]+'; then
    _pass "samtools version is valid: ${SAMTOOLS_LINE}"
else
    _fail "samtools version is invalid (was 'samtools=1.' — conda build fails): ${SAMTOOLS_LINE}"
fi

# ── 3. Scratch resolver priority ─────────────────────────────────────────────
_sec "Fix 3 — scratch resolver: /scratch0 tried before \$TMPDIR"

MISER="${PIPELINE_DIR}/modules/local/miser.nf"
if [ ! -f "${MISER}" ]; then
    _fail "modules/local/miser.nf not found"
else
    L0=$(grep -n 'scratch0' "${MISER}" | grep '_try_scratch' | head -1 | cut -d: -f1)
    LT=$(grep -n 'TMPDIR'   "${MISER}" | grep '_try_scratch' | head -1 | cut -d: -f1)
    if [ -z "${L0}" ]; then
        _fail "/scratch0 not referenced in miser.nf _try_scratch resolver"
    elif [ -z "${LT}" ]; then
        _fail "\$TMPDIR not referenced in miser.nf _try_scratch resolver"
    elif [ "${L0}" -lt "${LT}" ]; then
        _pass "/scratch0 (line ${L0}) tried before \$TMPDIR (line ${LT})"
    else
        _fail "/scratch0 (line ${L0}) tried AFTER \$TMPDIR (line ${LT}) — wrong priority"
    fi
    grep -q '_try_scratch()' "${MISER}" \
        && _pass "_try_scratch() exists — 40 G check is inside resolver" \
        || _fail "_try_scratch() not found — resolver not updated"
fi

# ── 4. CLEAN_BAM exists and is wired ─────────────────────────────────────────
_sec "Fix 4 — CLEAN_BAM process (parallel BAM cleaning)"

[ -f "${PIPELINE_DIR}/modules/local/clean_bam.nf" ] \
    && _pass "modules/local/clean_bam.nf exists" \
    || _fail "modules/local/clean_bam.nf MISSING"

grep -q 'CLEAN_BAM' "${PIPELINE_DIR}/main.nf" \
    && _pass "CLEAN_BAM imported and called in main.nf" \
    || _fail "CLEAN_BAM not referenced in main.nf"

grep -q 'samtools view.*2304' \
    "${PIPELINE_DIR}/modules/local/prepare_oarfish_ref.nf" 2>/dev/null \
    && _fail "prepare_oarfish_ref.nf still has -F 2304 cleaning — should be in clean_bam.nf only" \
    || _pass "prepare_oarfish_ref.nf has no BAM cleaning (correctly moved to clean_bam.nf)"

# ── 5. Dynamic sort memory ───────────────────────────────────────────────────
_sec "Fix 5 — samtools sort memory is dynamic (not hardcoded 3G)"

SORT="${PIPELINE_DIR}/modules/local/samtools_sort_index.nf"
if grep -q '\-m 3G' "${SORT}" 2>/dev/null; then
    _fail "samtools_sort_index.nf still hardcodes -m 3G"
else
    grep -qE 'MEM_PER_THREAD|toMega|memory.*cpus' "${SORT}" 2>/dev/null \
        && _pass "samtools_sort_index.nf uses dynamic per-thread memory" \
        || _fail "samtools_sort_index.nf: 3G removed but no dynamic calculation found"
fi

# ── 6. miser_qc.nf staged script ─────────────────────────────────────────────
_sec "Fix 6 — miser_qc.nf stages Python script via channel"

grep -q 'projectDir.*bin' "${PIPELINE_DIR}/modules/local/miser_qc.nf" 2>/dev/null \
    && _fail "miser_qc.nf still uses \${projectDir}/bin/ — Singularity may not see it" \
    || _pass "miser_qc.nf does not access \${projectDir}/bin/ directly"

grep -q 'path script' "${PIPELINE_DIR}/modules/local/miser_qc.nf" 2>/dev/null \
    && _pass "miser_qc.nf receives script via 'path script' channel input" \
    || _fail "miser_qc.nf: 'path script' input not found"

grep -q 'miser_qc_summary.py' "${PIPELINE_DIR}/main.nf" \
    && _pass "main.nf creates ch_miser_qc_script channel" \
    || _fail "main.nf does not stage miser_qc_summary.py"

# ── 7. pair/batch support ────────────────────────────────────────────────────
_sec "Fix 7 — pair/batch samplesheet columns flow to conditions CSV"

grep -q 'meta.pair'  "${PIPELINE_DIR}/main.nf" \
    && _pass "main.nf reads 'pair' column into meta" \
    || _fail "main.nf does not read 'pair' column"

grep -q 'meta.batch' "${PIPELINE_DIR}/main.nf" \
    && _pass "main.nf reads 'batch' column into meta" \
    || _fail "main.nf does not read 'batch' column"

grep 'seed:' "${PIPELINE_DIR}/main.nf" | head -1 | grep -q 'pair.*batch' \
    && _pass "conditions CSV seed includes pair and batch columns" \
    || _fail "conditions CSV seed does not include pair/batch columns"

# ── 8. Dead code removed from swish.nf ───────────────────────────────────────
_sec "Fix 8 — dead pair_arg/batch_arg removed from swish.nf"

grep -qE 'pair_arg|batch_arg' "${PIPELINE_DIR}/modules/local/swish.nf" 2>/dev/null \
    && _fail "swish.nf still contains dead pair_arg / batch_arg variables" \
    || _pass "swish.nf: dead variables removed"

# ── 9. New modules exist and routed ──────────────────────────────────────────
_sec "Fix 9 — new modules exist (samtools_stats, isoquant_qc, swish_plots, multiqc)"

for mod in samtools_stats isoquant_qc swish_plots multiqc; do
    [ -f "${PIPELINE_DIR}/modules/local/${mod}.nf" ] \
        && _pass "modules/local/${mod}.nf exists" \
        || _fail "modules/local/${mod}.nf MISSING"
done

grep -q 'MULTIQC' "${PIPELINE_DIR}/main.nf" \
    && _pass "MULTIQC wired in main.nf" \
    || _fail "MULTIQC not wired in main.nf"

grep -q 'SWISH_PLOTS' "${PIPELINE_DIR}/main.nf" \
    && _pass "SWISH_PLOTS wired in main.nf" \
    || _fail "SWISH_PLOTS not wired in main.nf"

# ── 10. Container routing complete ───────────────────────────────────────────
_sec "Fix 10 — nextflow.config routes all new processes to a container"

for proc in CLEAN_BAM SAMTOOLS_STATS ISOQUANT_QC SWISH_PLOTS MULTIQC; do
    grep -q "${proc}" "${PIPELINE_DIR}/nextflow.config" \
        && _pass "nextflow.config: ${proc} routed to container" \
        || _fail "nextflow.config: ${proc} NOT in container routing"
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo
echo "============================================================"
echo " Results:  PASS ${PASS}   FAIL ${FAIL}"
echo "============================================================"
[ "${FAIL}" -eq 0 ] \
    && echo " All checks passed. Safe to proceed to Stage 1 stub run." \
    && exit 0 \
    || echo " Fix all FAILs before running on the cluster." \
    && exit 1
