# LoReRNA v1.0.0 — Validation Runbook

Four stages in order of cost. Stop if any stage fails — there is no point
spending queue time when a cheaper stage has already found the problem.

---

## Stage 0 — Pre-flight (login node, ~5 min, no jobs)

Checks every code fix by inspecting `.nf` files directly.
No cluster access required.

```bash
cd /path/to/lorerna
bash test/validate_preflight.sh
```

**Expected:** `PASS: 10  FAIL: 0`

Fix every FAIL before proceeding.

---

## Stage 1 — Stub run (login node, ~2 min, no jobs)

Runs the full pipeline DAG with stub outputs — each process returns
immediately with empty files. Validates channel wiring, parameter
resolution, and that `LORERNA:CLEAN_BAM` is present in the DAG.

```bash
# 1. Fill in two real BAM paths in test/samplesheet_2sample.csv
# 2. Run stub:
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  -c conf/validate.config \
  --samplesheet   test/samplesheet_2sample.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  -stub-run
```

**What to check in the execution summary:**

```
LORERNA:MISER                      1
LORERNA:MISER_QC                   1
LORERNA:MISER_QC_MERGE             1
LORERNA:SAMTOOLS_SORT_INDEX        1
LORERNA:SAMTOOLS_STATS             1   ← new in v1.0.0
LORERNA:CLEAN_BAM                  1   ← new in v1.0.0 (was sequential for-loop)
LORERNA:ISOQUANT                   1
LORERNA:ISOQUANT_QC                1   ← new in v1.0.0
LORERNA:PREPARE_OARFISH_REFERENCE  1
LORERNA:OARFISH                    1
LORERNA:SWISH                      1
LORERNA:SWISH_PLOTS                1   ← new in v1.0.0
LORERNA:MULTIQC                    1   ← new in v1.0.0
```

If `CLEAN_BAM` or `MULTIQC` are absent, the include or channel wiring in
`main.nf` is broken. Fix before proceeding.

---

## Stage 2 — 2-sample cluster run (~6–8 h)

The minimum run that exercises every process including SWISH.
`conf/validate.config` reduces bootstraps (10) and permutations (10)
so oarfish and swish finish in ~1–2 h instead of ~4–6 h.

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  -c conf/validate.config \
  --samplesheet   test/samplesheet_2sample.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD
```

### Checks after completion

**MisER scratch — was the critical failure mode before v1.0.0:**
```bash
grep 'Scratch resolver' lorerna_validation/01_miser/*/miser.log
# Must show:  Scratch resolver  : /scratch0/...
# Must NOT:   Scratch resolver  : /tmp/...
```

**CLEAN_BAM — new parallel process:**
```bash
ls -lh lorerna_validation/02_sorted_bam/*.clean.bam
# Expected: one .clean.bam per sample, non-zero size
```

**Conditions CSV header — the other critical fix:**
```bash
head -1 lorerna_validation/06_swish/conditions.csv
# Must be:  sample_id,condition,pair,batch
# Must NOT: sample_name,condition
```

**SWISH produced results — confirms R stopifnot did not fire:**
```bash
ls lorerna_validation/06_swish/results/*.csv
# Expected: DTE_full_results.csv, DTE_all_significant.csv,
#           DTU_*, DGE_* equivalents
```

**Publication plots generated:**
```bash
ls lorerna_validation/06_swish/publication_plots/*.pdf
# Expected: DTE_publication_plots.pdf, DTU_publication_plots.pdf,
#           DGE_publication_plots.pdf, summary_panel.pdf
```

**MultiQC report:**
```bash
ls lorerna_validation/07_multiqc/multiqc_report.html
# Open in a browser and check:
#   - samtools section populated (alignment rates per sample)
#   - IsoQuant assignment section populated
#   - MisER QC section populated
```

---

## Stage 3 — Paired design validation (optional, ~5 min, login node)

Confirms the pair/batch samplesheet columns flow correctly through
to the conditions CSV without a cluster run.

```bash
# Fill in paired BAM paths in test/samplesheet_paired.csv, then:
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  -c conf/validate.config \
  --samplesheet   test/samplesheet_paired.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  -stub-run

# Then check the conditions CSV for pair column content:
find work -name 'conditions.csv' -newer main.nf | head -1 | xargs cat
# Expected header: sample_id,condition,pair,batch
# Expected data  : pair column populated with 1, 2, ...
```

---

## Stage 4 — Full cohort (production)

Only run this after all previous stages pass.

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  --samplesheet   samplesheet.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  --outdir /path/to/final/results \
  -resume
```

The `-resume` flag reuses any work completed in previous runs.
Nextflow matches cached results by hashing inputs — if inputs have
not changed, nothing is recomputed.

---

## Failure reference

| Failure | Most likely cause | Where to look |
|---------|------------------|---------------|
| MISER scratch on `/tmp` | `scratch = false` missing for MISER | `conf/sge.config` `withName: 'MISER'` block |
| MISER scratch < 40 G | Scheduler did not allocate local scratch | `tscratch=60G` (SGE) or `--tmp=61440` (SLURM) in MISER `clusterOptions` |
| `/scratch0` not visible in container | Bind mount missing | `containerOptions = "-B /scratch0"` in `withName: 'MISER'` |
| SWISH: `sample_id not found` | Wrong CSV header | `seed:` line in `main.nf` must be `"sample_id,condition,pair,batch\n"` |
| CLEAN_BAM not in DAG | Process not wired | Check `include { CLEAN_BAM }` and `CLEAN_BAM(...)` in `main.nf` |
| SWISH_PLOTS: missing R package | Wrong container tag | Ensure `container_lorerna = 'docker://your-hub/lorerna:1.0.0'` |
| MultiQC: empty sections | File naming mismatch | Check `sp:` patterns in `containers/multiqc_config.yml` vs actual filenames |
| Memory error on retry | Flat memory — no scaling | Use `memory = { 64.GB * task.attempt }` not `memory = 64.GB` |
