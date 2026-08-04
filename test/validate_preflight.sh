#!/usr/bin/env bash
# =============================================================================
# validate_preflight.sh — LoReRNA v1.0.0 pre-flight checks
#
# Runs on the login node. No cluster jobs submitted. ~5 minutes.
# Checks every fix in-place by inspecting the .nf and source files directly.
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

# The workflow body moved from main.nf to workflows/lorerna.nf in the nf-core
# layout refactor, leaving main.nf a thin entry point. Checks below search both,
# so they keep working wherever a given statement now lives.
# An array, not a string: the checkout path may contain spaces.
NF_SRC=("${PIPELINE_DIR}/main.nf" "${PIPELINE_DIR}/workflows/lorerna.nf")

_pass() { echo "  PASS  $1"; PASS=$(( PASS + 1 )); }
_fail() { echo "  FAIL  $1"; FAIL=$(( FAIL + 1 )); }
_sec()  { echo; echo "── $1 ──────────────────────────────────────────────────"; }

echo "============================================================"
echo " LoReRNA v1.0.0 pre-flight validation"
echo " Pipeline : ${PIPELINE_DIR}"
echo " Date     : $(date)"
echo "============================================================"

# ── 1. Conditions CSV header ─────────────────────────────────────────────────
_sec "1 — conditions CSV seed uses sample_id"

SEED=$(grep 'seed:' "${NF_SRC[@]}" 2>/dev/null | head -1 || true)
if echo "${SEED}" | grep -q '"sample_id,condition'; then
    _pass "conditions CSV seed uses 'sample_id'"
else
    _fail "conditions CSV seed does NOT use 'sample_id' — SWISH will fail"
    echo "        Found : ${SEED}"
fi

grep -q 'sample_id' "${PIPELINE_DIR}/bin/swish_analysis.R" \
    && _pass "swish_analysis.R references 'sample_id'" \
    || _fail "swish_analysis.R does not reference 'sample_id'"

# ── 2. samtools version valid in both environment files ──────────────────────
_sec "2 — samtools version strings are valid"

for yml in containers/environment.yml containers/environment_miser.yml; do
    LINE=$(grep 'samtools' "${PIPELINE_DIR}/${yml}" | head -1)
    if echo "${LINE}" | grep -qE 'samtools=[0-9]+\.[0-9]+'; then
        _pass "${yml}: ${LINE}"
    else
        _fail "${yml}: invalid samtools version — ${LINE}"
    fi
done

# ── 3. IsoQuant executable renamed (v1.0.0 fix) ─────────────────────────────
_sec "3 — IsoQuant 3.13.0 executable is 'isoquant' not 'isoquant.py'"

IQ_EXE=$(grep 'isoquant_exe' "${PIPELINE_DIR}/nextflow.config" | grep -v '//' | head -1)
if echo "${IQ_EXE}" | grep -q "'isoquant'"; then
    _pass "nextflow.config: isoquant_exe = 'isoquant'"
else
    _fail "nextflow.config: isoquant_exe is not 'isoquant' — IsoQuant 3.13.0 will not be found"
    echo "        Found : ${IQ_EXE}"
fi

grep -q "isoquant\.py" "${PIPELINE_DIR}/containers/build_containers.sh" \
    && _fail "build_containers.sh still calls isoquant.py in sanity check" \
    || _pass "build_containers.sh uses 'isoquant' in sanity check"

# ── 4. IsoQuant emit: outdir removed (v1.0.0 fix) ───────────────────────────
_sec "4 — isoquant.nf emit: outdir removed"

grep -q 'emit: outdir' "${PIPELINE_DIR}/modules/local/isoquant.nf" \
    && _fail "isoquant.nf still has 'emit: outdir' — triple-nested output dir will appear" \
    || _pass "isoquant.nf: emit: outdir removed"

grep -q 'ISOQUANT\.out\.outdir' "${NF_SRC[@]}" \
    && _fail "main.nf still references ISOQUANT.out.outdir" \
    || _pass "main.nf: no reference to ISOQUANT.out.outdir"

# ── 5. Container tags are 1.1.0 (v1.0.0 fix) ────────────────────────────────
_sec "5 — container tags updated to 1.1.0"

# Extract the actual image tag (docker:// contains '//', so don't filter on '//')
LORERNA_TAG=$(grep -oE "lorerna:[0-9]+\.[0-9]+\.[0-9]+" "${PIPELINE_DIR}/nextflow.config" | head -1)
MISER_TAG=$(grep -oE "miser:[0-9]+\.[0-9]+\.[0-9]+"     "${PIPELINE_DIR}/nextflow.config" | head -1)

[ "${LORERNA_TAG}" = "lorerna:1.1.0" ] \
    && _pass "container_lorerna: ${LORERNA_TAG}" \
    || _fail "container_lorerna is not 1.1.0: ${LORERNA_TAG:-<not found>}"

[ "${MISER_TAG}" = "miser:1.1.0" ] \
    && _pass "container_miser: ${MISER_TAG}" \
    || _fail "container_miser is not 1.1.0: ${MISER_TAG:-<not found>}"

# ── 6. swish publication plots publishDir fix (v1.0.0 fix) ──────────────────
_sec "6 — swish_plots.nf publishDir has saveAs to strip pub_plots/"

# publishDir moved out of the modules into conf/modules.config, so the saveAs
# that strips the pub_plots/ prefix now lives there.
grep -q 'saveAs.*replaceFirst.*pub_plots' "${PIPELINE_DIR}/conf/modules.config" \
    && _pass "conf/modules.config: SWISH_PLOTS saveAs strips pub_plots/ prefix" \
    || _fail "conf/modules.config: missing saveAs — PDFs will nest under publication_plots/pub_plots/"

grep -q 'conditions\.csv.*continue' "${PIPELINE_DIR}/modules/local/swish_plots.nf" \
    && _pass "swish_plots.nf: conditions.csv correctly excluded from results/ move" \
    || _fail "swish_plots.nf: conditions.csv not excluded — will be moved to results/ and lost"

# ── 7. top_labels() fix in swish_plots.R (v1.0.0 fix) ───────────────────────
_sec "7 — swish_plots.R top_labels() uses sig\$gene not sig\$label"

grep -q 'head(sig\$gene' "${PIPELINE_DIR}/bin/swish_plots.R" \
    && _pass "swish_plots.R: top_labels() returns sig\$gene" \
    || _fail "swish_plots.R: top_labels() may still use sig\$label — gene labels will be blank"

grep -q 'head(sig\$label' "${PIPELINE_DIR}/bin/swish_plots.R" \
    && _fail "swish_plots.R: sig\$label still present in top_labels()" \
    || _pass "swish_plots.R: sig\$label not found in top_labels()"

# ── 8. Multi-condition filter in swish scripts (v1.0.0 fix) ─────────────────
_sec "8 — condition filter in swish_analysis.R and swish_plots.R"

grep -q 'cond_tbl.*condition.*%in%.*cond_a.*cond_b' "${PIPELINE_DIR}/bin/swish_analysis.R" \
    && _pass "swish_analysis.R: condition filter present" \
    || _fail "swish_analysis.R: condition filter missing — 3-condition samplesheets will fail"

grep -q 'cond_tbl.*condition.*%in%.*cond_a.*cond_b' "${PIPELINE_DIR}/bin/swish_plots.R" \
    && _pass "swish_plots.R: condition filter present" \
    || _fail "swish_plots.R: condition filter missing — 3-condition samplesheets will produce wrong plots"

grep -q 'nrow(cond_tbl).*0' "${PIPELINE_DIR}/bin/swish_analysis.R" \
    && _pass "swish_analysis.R: empty-condition guard present" \
    || _fail "swish_analysis.R: no empty-condition guard — mistyped condition names give cryptic errors"

# ── 9. MultiQC output dir name (v1.0.0 fix) ─────────────────────────────────
_sec "9 — multiqc.nf output is multiqc_report_data/ not multiqc_data/"

grep -q 'multiqc_report_data' "${PIPELINE_DIR}/modules/local/multiqc.nf" \
    && _pass "multiqc.nf: emits multiqc_report_data/ (correct for --filename flag)" \
    || _fail "multiqc.nf: still emits multiqc_data/ — collect will fail with MultiQC 1.25.1"

grep -q '"multiqc_data/"' "${PIPELINE_DIR}/modules/local/multiqc.nf" \
    && _fail "multiqc.nf: old multiqc_data/ path still present" \
    || _pass "multiqc.nf: old multiqc_data/ path removed"

# ── 10. Miser scratch uses /scratch0 ────────────────────────────────────────
_sec "10 — miser.nf scratch resolver targets /scratch0"

grep -q 'scratch0' "${PIPELINE_DIR}/modules/local/miser.nf" \
    && _pass "miser.nf: /scratch0 referenced as scratch path" \
    || _fail "miser.nf: /scratch0 not found — MisER will not use node-local scratch"

grep -q 'miser_scratch_root' "${PIPELINE_DIR}/modules/local/miser.nf" \
    && _pass "miser.nf: params.miser_scratch_root override present" \
    || _fail "miser.nf: no params.miser_scratch_root — scratch path cannot be overridden per site"

# ── 11. CLEAN_BAM exists and wired ───────────────────────────────────────────
_sec "11 — CLEAN_BAM process exists and is wired in main.nf"

[ -f "${PIPELINE_DIR}/modules/local/clean_bam.nf" ] \
    && _pass "modules/local/clean_bam.nf exists" \
    || _fail "modules/local/clean_bam.nf MISSING"

grep -q 'CLEAN_BAM' "${NF_SRC[@]}" \
    && _pass "CLEAN_BAM wired in main.nf" \
    || _fail "CLEAN_BAM not referenced in main.nf"

# ── 12. Dynamic sort memory ──────────────────────────────────────────────────
_sec "12 — samtools_sort_index.nf uses dynamic per-thread memory"

grep -q 'MEM_PER_THREAD' "${PIPELINE_DIR}/modules/local/samtools_sort_index.nf" \
    && _pass "samtools_sort_index.nf: dynamic MEM_PER_THREAD calculation" \
    || _fail "samtools_sort_index.nf: MEM_PER_THREAD not found"

# Only flag -m 3G in an actual command line, not in the // changelog comment
grep -v '^[[:space:]]*//' "${PIPELINE_DIR}/modules/local/samtools_sort_index.nf" | grep -q '\-m 3G' \
    && _fail "samtools_sort_index.nf: hardcoded -m 3G still present in a command" \
    || _pass "samtools_sort_index.nf: hardcoded -m 3G removed"

# ── 13. pair/batch columns flow to conditions CSV ────────────────────────────
_sec "13 — pair/batch samplesheet columns flow to conditions CSV"

grep -q 'meta.pair'  "${NF_SRC[@]}" \
    && _pass "main.nf: 'pair' column read into meta" \
    || _fail "main.nf: 'pair' column not read"

grep -q 'meta.batch' "${NF_SRC[@]}" \
    && _pass "main.nf: 'batch' column read into meta" \
    || _fail "main.nf: 'batch' column not read"

grep 'seed:' "${NF_SRC[@]}" 2>/dev/null | head -1 | grep -q 'pair.*batch' \
    && _pass "conditions CSV seed includes pair and batch columns" \
    || _fail "conditions CSV seed missing pair/batch columns"

# ── 14. Container routing complete ───────────────────────────────────────────
_sec "14 — nextflow.config routes all processes to a container"

for proc in CLEAN_BAM SAMTOOLS_STATS ISOQUANT_QC SWISH_PLOTS MULTIQC; do
    grep -q "${proc}" "${PIPELINE_DIR}/nextflow.config" \
        && _pass "nextflow.config: ${proc} routed to container" \
        || _fail "nextflow.config: ${proc} NOT in container routing"
done

# ── 15. prepare_oarfish_ref.nf — no duplicate sort ───────────────────────────
_sec "15 — prepare_oarfish_ref.nf has no duplicate sort -u"

COUNT=$(grep -c 'sort -u tx2gene' "${PIPELINE_DIR}/modules/local/prepare_oarfish_ref.nf" || true)
if [ "${COUNT}" -eq 1 ]; then
    _pass "prepare_oarfish_ref.nf: exactly one sort -u tx2gene"
elif [ "${COUNT}" -gt 1 ]; then
    _fail "prepare_oarfish_ref.nf: duplicate sort -u tx2gene (${COUNT} occurrences)"
else
    _fail "prepare_oarfish_ref.nf: sort -u tx2gene not found at all"
fi

# ── 16. conditions.csv publishing (Issue 9) ─────────────────────────────────
_sec "16 — swish.nf publishes a real conditions.csv without self-copy"

grep -q "stageAs: 'conditions_in.csv'" "${PIPELINE_DIR}/modules/local/swish.nf" \
    && _pass "swish.nf: conditions_csv staged as conditions_in.csv (no self-copy on cp)" \
    || _fail "swish.nf: missing stageAs — 'cp conditions.csv conditions.csv' will abort under set -e"

grep -q 'emit: conditions_out' "${PIPELINE_DIR}/modules/local/swish.nf" \
    && _pass "swish.nf: conditions.csv emitted (real file, survives nextflow clean)" \
    || _fail "swish.nf: conditions_out not emitted — Issue 9 not applied"

# ── 17. MisER rescue metrics reach MultiQC (Issue 10a) ──────────────────────
_sec "17 — MisER merged metrics use _mqc.tsv suffix"

grep -q 'all_samples_rescue_metrics_mqc.tsv' "${PIPELINE_DIR}/modules/local/miser_qc_merge.nf" \
    && _pass "miser_qc_merge.nf: output is all_samples_rescue_metrics_mqc.tsv" \
    || _fail "miser_qc_merge.nf: missing _mqc suffix — MultiQC will ignore MisER metrics"

# The MisER section is now a self-describing custom-content file written by
# bin/miser_rescue_multiqc.py, so it deliberately has no custom_data / sp entry
# in multiqc_config.yml. Check the emitting module instead.
grep -q 'miser_rescue_mqc.json' "${PIPELINE_DIR}/modules/local/miser_qc_merge.nf" \
    && _pass "miser_qc_merge.nf: emits miser_rescue_mqc.json for MultiQC" \
    || _fail "miser_qc_merge.nf: no miser_rescue_mqc.json — MisER section will be absent from the report"

# ── 18. SWISH summary reaches MultiQC (Issue 10b) ───────────────────────────
_sec "18 — SWISH summary TSV wired into MultiQC"

grep -q 'swish_summary_mqc.tsv' "${PIPELINE_DIR}/bin/swish_analysis.R" \
    && _pass "swish_analysis.R: writes swish_summary_mqc.tsv" \
    || _fail "swish_analysis.R: does not write swish_summary_mqc.tsv"

grep -q 'emit: mqc_tsv' "${PIPELINE_DIR}/modules/local/swish.nf" \
    && _pass "swish.nf: mqc_tsv emitted" \
    || _fail "swish.nf: mqc_tsv not emitted"

grep -q 'SWISH.out.mqc_tsv' "${NF_SRC[@]}" \
    && _pass "main.nf: MultiQC consumes SWISH.out.mqc_tsv" \
    || _fail "main.nf: still mixing SWISH.out.log instead of mqc_tsv"

grep -q 'swish_summary' "${PIPELINE_DIR}/containers/multiqc_config.yml" \
    && _pass "multiqc_config.yml: swish_summary section present" \
    || _fail "multiqc_config.yml: swish_summary section missing"

# ── 19. Container images pinned to correct tool versions ────────────────────
_sec "19 — environment.yml tool versions match the 1.1.0 container"

for pin in 'isoquant=3.13.0' 'oarfish=0.9.4' 'samtools=1.23.1'; do
    grep -q "${pin}" "${PIPELINE_DIR}/containers/environment.yml" \
        && _pass "environment.yml pins ${pin}" \
        || _fail "environment.yml does NOT pin ${pin}"
done

# ── 20. Container runtime environment ────────────────────────────────────────
# Checks 1-19 inspect the pipeline's own files. This one checks the environment
# the pipeline is about to run in, which is where cluster users actually get
# stuck: Singularity builds an image by unpacking OCI layers into a scratch
# directory, resolving it as SINGULARITY_TMPDIR -> TMPDIR -> /tmp. Schedulers
# routinely point TMPDIR at node-local scratch, which exists on compute nodes
# but not on the login node where `nextflow run` performs the pull. The result
# is a Singularity-internal error about a "build parent dir" that says nothing
# about the actual cause.
_sec "20 — container runtime environment"

_ENGINE=""
for e in singularity apptainer docker; do
    command -v "${e}" >/dev/null 2>&1 && { _ENGINE="${e}"; break; }
done
if [ -n "${_ENGINE}" ]; then
    _pass "container engine found: ${_ENGINE}"
else
    _fail "no container engine on PATH (need singularity, apptainer or docker)"
fi

if [ "${_ENGINE}" = "singularity" ] || [ "${_ENGINE}" = "apptainer" ]; then
    # Resolve the build scratch the same way Singularity does.
    _TMP="${SINGULARITY_TMPDIR:-${TMPDIR:-/tmp}}"
    _SRC="SINGULARITY_TMPDIR"
    [ -n "${SINGULARITY_TMPDIR:-}" ] || _SRC="TMPDIR"
    [ -n "${SINGULARITY_TMPDIR:-}" ] || [ -n "${TMPDIR:-}" ] || _SRC="default"

    if [ ! -d "${_TMP}" ]; then
        _fail "image build scratch does not exist: ${_TMP}  (from ${_SRC})"
        echo "        Singularity cannot unpack image layers there, so every"
        echo "        container pull will fail before any analysis starts."
        echo "        Fix: export SINGULARITY_TMPDIR=/tmp   (or any existing"
        echo "        directory with a few GB free), then re-run."
    elif [ ! -w "${_TMP}" ]; then
        _fail "image build scratch is not writable: ${_TMP}  (from ${_SRC})"
        echo "        Fix: export SINGULARITY_TMPDIR to a writable directory."
    else
        _AVAIL=$(df -Pk "${_TMP}" 2>/dev/null | awk 'NR==2 {print int($4/1048576)}')
        if [ -n "${_AVAIL}" ] && [ "${_AVAIL}" -lt 5 ]; then
            _fail "image build scratch has only ${_AVAIL} G free: ${_TMP}"
            echo "        Unpacking the layers needs roughly 5 G. Point"
            echo "        SINGULARITY_TMPDIR at a filesystem with more room."
        else
            _pass "image build scratch usable: ${_TMP} (${_AVAIL:-?} G free, from ${_SRC})"
        fi
    fi
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo
echo "============================================================"
echo " Results:  PASS ${PASS}   FAIL ${FAIL}"
echo "============================================================"
if [ "${FAIL}" -eq 0 ]; then
    echo " All checks passed. Safe to proceed to Stage 1 stub run."
    exit 0
else
    echo " Fix all FAILs before running on the cluster."
    exit 1
fi
