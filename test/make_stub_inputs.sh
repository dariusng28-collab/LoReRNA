#!/usr/bin/env bash
# =============================================================================
# make_stub_inputs.sh — create placeholder inputs for `-stub-run` wiring tests.
#
# `-stub-run` never reads these files (stub blocks just touch outputs), but
# main.nf validates the paths exist via checkIfExists. Empty files satisfy that.
#
#   bash test/make_stub_inputs.sh
#   nextflow run main.nf -profile test -stub-run
# =============================================================================
set -euo pipefail

DATA="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data"
mkdir -p "${DATA}"

: > "${DATA}/genome.fa"
: > "${DATA}/annotation.gtf"
: > "${DATA}/annotation.bed12"
: > "${DATA}/WT1.bam"
: > "${DATA}/KD1.bam"

cat > "${DATA}/samplesheet_stub.csv" <<'CSV'
sample_name,condition,bam
WT1,WT,test/data/WT1.bam
KD1,KD,test/data/KD1.bam
CSV

echo "Placeholder stub inputs written to ${DATA}"
ls -la "${DATA}"
