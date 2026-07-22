#!/usr/bin/env bash
# =============================================================================
# container_cli_contract.sh
#
# Asserts that the tool CLIs and R packages LoReRNA depends on actually exist
# in the containers, with the EXACT flags the pipeline passes. This is the
# guard that would have caught the oarfish 0.9.4 break (--reference removed,
# renamed to --annotated) BEFORE it reached a run: it checks that
# `oarfish --help` still advertises every flag oarfish.nf uses.
#
# Run it in CI on every PR, and locally before bumping any tool version.
#
# Usage (Docker — CI default):
#   bash test/container_cli_contract.sh
#   bash test/container_cli_contract.sh dariusng28/lorerna:1.1.0 dariusng28/miser:1.1.0
#
# Usage (Singularity — e.g. on an HPC, point at local .img files):
#   ENGINE=singularity \
#     LORERNA_IMG=/path/dariusng28-lorerna-1.1.0.img \
#     MISER_IMG=/path/dariusng28-miser-1.1.0.img \
#     bash test/container_cli_contract.sh
# =============================================================================
set -uo pipefail

ENGINE="${ENGINE:-docker}"
LORERNA_IMG="${1:-${LORERNA_IMG:-dariusng28/lorerna:1.1.0}}"
MISER_IMG="${2:-${MISER_IMG:-dariusng28/miser:1.1.0}}"

run() {  # <image> <shell-command>
  case "${ENGINE}" in
    docker)      docker run --rm "$1" bash -lc "$2" ;;
    singularity) singularity exec "$1" bash -lc "$2" ;;
    *) echo "unknown ENGINE='${ENGINE}' (use docker|singularity)"; exit 2 ;;
  esac
}

PASS=0; FAIL=0
check() {  # <description> <image> <command that must exit 0>
  if run "$2" "$3" >/dev/null 2>&1; then
    printf '  PASS  %s\n' "$1"; PASS=$(( PASS + 1 ))
  else
    printf '  FAIL  %s\n' "$1"; FAIL=$(( FAIL + 1 ))
  fi
}

echo "============================================================"
echo " Container CLI + R-package contract"
echo " engine  : ${ENGINE}"
echo " lorerna : ${LORERNA_IMG}"
echo " miser   : ${MISER_IMG}"
echo "============================================================"

# ── lorerna: executables the pipeline invokes ───────────────────────────────
echo; echo "── lorerna: executables present ──"
for tool in samtools isoquant oarfish gffcompare gffread multiqc Rscript; do
    check "exe present: ${tool}" "${LORERNA_IMG}" "command -v ${tool}"
done

# ── oarfish flag contract (modules/local/oarfish.nf) ────────────────────────
# THE regression guard. oarfish.nf passes exactly these flags in --reads mode.
echo; echo "── oarfish flag contract ──"
for flag in --reads --annotated --seq-tech --filter-group --num-bootstraps --model-coverage --threads --output; do
    check "oarfish --help advertises ${flag}" "${LORERNA_IMG}" \
          "oarfish --help 2>&1 | grep -q -- '${flag}'"
done

# ── isoquant flag contract (modules/local/isoquant.nf) ──────────────────────
echo; echo "── isoquant flag contract ──"
check "isoquant --version runs" "${LORERNA_IMG}" "isoquant --version"
for flag in --reference --bam --genedb --sqanti_output --data_type --prefix --threads; do
    check "isoquant --help advertises ${flag}" "${LORERNA_IMG}" \
          "isoquant --help 2>&1 | grep -q -- '${flag}'"
done

# ── R packages every bin/*.R script loads ───────────────────────────────────
echo; echo "── R package contract ──"
check "R packages load (fishpond, data.table, arrow, ggplot2, ggrepel, pheatmap, patchwork, scales, SummarizedExperiment)" \
      "${LORERNA_IMG}" \
      "Rscript -e 'pkgs<-c(\"fishpond\",\"data.table\",\"arrow\",\"ggplot2\",\"ggrepel\",\"pheatmap\",\"patchwork\",\"scales\",\"SummarizedExperiment\"); invisible(lapply(pkgs, library, character.only=TRUE))'"

# ── miser container (modules/local/miser.nf, miser_qc.nf) ───────────────────
echo; echo "── miser image ──"
check "exe present: samtools"    "${MISER_IMG}" "command -v samtools"
check "python: import MisER"     "${MISER_IMG}" "python3 -c 'import MisER'"
check "python: import pysam"     "${MISER_IMG}" "python3 -c 'import pysam'"
check "python: import pandas"    "${MISER_IMG}" "python3 -c 'import pandas'"

# ── summary ─────────────────────────────────────────────────────────────────
echo
echo "============================================================"
echo " contract results:  PASS ${PASS}   FAIL ${FAIL}"
echo "============================================================"
if [ "${FAIL}" -ne 0 ]; then
    echo " A tool CLI or package the pipeline depends on has changed."
    echo " Fix the module/env before shipping — this is a real break."
    exit 1
fi
echo " All contracts satisfied."
