# LoReRNA v1.0.0 — Push Guide

Target repository: https://github.com/dariusng28-collab/LoReRNA

---

## Overview

You already have the repo. The plan is:
1. Clone it locally
2. Remove old files, copy in new files
3. Commit and push as v1.0.0
4. Build and push containers
5. Validate on the cluster

Work through every step in order.

---

## Part 1 — Local machine: build containers first

Build the containers before touching git — if the build fails you
want to know before committing anything.

### Step 1 — Log in to Docker Hub

```bash
docker login
# Enter your Docker Hub username and password
```

Expected: `Login Succeeded`

### Step 2 — Open nextflow.config and confirm your Docker Hub username

```bash
grep 'container_' nextflow.config
```

Should show:
```
container_lorerna = 'docker://dariusng28/lorerna:1.0.0'
container_miser   = 'docker://dariusng28/miser:1.0.0'
```

If your Docker Hub username is different, update those two lines and also
update `DOCKER_USER` on line ~10 of `containers/build_containers.sh`.

### Step 3 — Build both images (no push yet)

```bash
bash containers/build_containers.sh
```

Wait for the sanity checks to complete. Expected at the end:

```
 lorerna:
samtools 1.21 ...
MultiQC v1.25.1

 miser:
samtools 1.21 ...
MisER OK

 To push to Docker Hub:
   bash containers/build_containers.sh --push --version 1.0.0
```

**If the build fails:** fix the error before continuing.
Most common cause: Docker daemon not running — start Docker Desktop.

### Step 4 — Push containers to Docker Hub

```bash
bash containers/build_containers.sh --push
```

Expected:
```
 Pushed to Docker Hub:
   https://hub.docker.com/r/dariusng28/lorerna/tags
   https://hub.docker.com/r/dariusng28/miser/tags
```

**Verify before continuing:** open hub.docker.com in your browser,
sign in, confirm both `lorerna:1.0.0` and `miser:1.0.0` tags exist.

---

## Part 2 — Local machine: push to GitHub

### Step 5 — Clone the existing repo

```bash
cd ~
git clone https://github.com/dariusng28-collab/LoReRNA.git
cd LoReRNA
```

### Step 6 — Remove all old files

```bash
git rm -rf .
```

Expected: a list of deleted files. The working directory is now empty
(except the hidden `.git/` folder which must stay).

Confirm:
```bash
ls -la
# Should show only .git/ — nothing else
```

### Step 7 — Copy all new files in

Replace `/path/to/delivered/lorerna` with wherever you saved the
files from this session.

```bash
DELIVERED=/path/to/delivered/lorerna

# Root files
cp ${DELIVERED}/main.nf             .
cp ${DELIVERED}/nextflow.config     .
cp ${DELIVERED}/CHANGES.md          .
cp ${DELIVERED}/LICENSE             .
cp ${DELIVERED}/README.md           .
cp ${DELIVERED}/SETUP.md            .
cp ${DELIVERED}/.gitignore          .

# conf/
mkdir -p conf
cp ${DELIVERED}/conf/sge.config        conf/
cp ${DELIVERED}/conf/slurm.config      conf/
cp ${DELIVERED}/conf/lsf.config        conf/
cp ${DELIVERED}/conf/validate.config   conf/

# modules/local/
mkdir -p modules/local
cp ${DELIVERED}/modules/local/miser.nf                 modules/local/
cp ${DELIVERED}/modules/local/miser_qc.nf              modules/local/
cp ${DELIVERED}/modules/local/samtools_sort_index.nf   modules/local/
cp ${DELIVERED}/modules/local/samtools_stats.nf        modules/local/
cp ${DELIVERED}/modules/local/clean_bam.nf             modules/local/
cp ${DELIVERED}/modules/local/isoquant_qc.nf           modules/local/
cp ${DELIVERED}/modules/local/prepare_oarfish_ref.nf   modules/local/
cp ${DELIVERED}/modules/local/swish.nf                 modules/local/
cp ${DELIVERED}/modules/local/swish_plots.nf           modules/local/
cp ${DELIVERED}/modules/local/multiqc.nf               modules/local/

# bin/
mkdir -p bin
cp ${DELIVERED}/bin/isoquant_qc.py    bin/
cp ${DELIVERED}/bin/swish_plots.R     bin/

# containers/
mkdir -p containers
cp ${DELIVERED}/containers/environment.yml         containers/
cp ${DELIVERED}/containers/environment_miser.yml   containers/
cp ${DELIVERED}/containers/multiqc_config.yml      containers/
cp ${DELIVERED}/containers/build_containers.sh     containers/

# test/
mkdir -p test
cp ${DELIVERED}/test/validate_preflight.sh     test/
cp ${DELIVERED}/test/samplesheet_2sample.csv   test/
cp ${DELIVERED}/test/samplesheet_paired.csv    test/
cp ${DELIVERED}/test/VALIDATION_RUNBOOK.md     test/
```

### Step 8 — Copy the 7 unchanged files from your existing pipeline

These files did not change and are not in the delivered set.
Copy them from your original pipeline directory.

```bash
ORIGINAL=/path/to/your/original/pipeline

cp ${ORIGINAL}/modules/local/miser_qc_merge.nf   modules/local/
cp ${ORIGINAL}/modules/local/isoquant.nf          modules/local/
cp ${ORIGINAL}/modules/local/oarfish.nf           modules/local/
cp ${ORIGINAL}/bin/miser_qc_summary.py            bin/
cp ${ORIGINAL}/bin/swish_analysis.R               bin/
cp ${ORIGINAL}/containers/Dockerfile              containers/
cp ${ORIGINAL}/containers/Dockerfile.miser        containers/
```

### Step 9 — Update your cluster config paths

Open `conf/sge.config` and update the three placeholder lines:

```bash
nano conf/sge.config
```

Change:
```groovy
singularity.runOptions = "-B ${HOME},${PWD},/path/to/shared/filesystem"
singularity.cacheDir   = "/path/to/singularity/cache"
queue    = 'your.queue'
```

To your actual values. Also update `conf/validate.config` if the
`outdir` path needs adjusting for your cluster.

### Step 10 — Verify the file count

```bash
find . -type f | grep -v '.git' | sort
```

You should see exactly **32 files**:

```bash
find . -type f | grep -v '.git' | wc -l
# Expected: 32
```

If the count is wrong, run this to find what is missing:

```bash
for f in \
  main.nf nextflow.config CHANGES.md LICENSE README.md SETUP.md .gitignore \
  conf/sge.config conf/slurm.config conf/lsf.config conf/validate.config \
  modules/local/miser.nf modules/local/miser_qc.nf modules/local/miser_qc_merge.nf \
  modules/local/samtools_sort_index.nf modules/local/samtools_stats.nf \
  modules/local/clean_bam.nf modules/local/isoquant.nf modules/local/isoquant_qc.nf \
  modules/local/prepare_oarfish_ref.nf modules/local/oarfish.nf \
  modules/local/swish.nf modules/local/swish_plots.nf modules/local/multiqc.nf \
  bin/miser_qc_summary.py bin/swish_analysis.R bin/isoquant_qc.py bin/swish_plots.R \
  containers/Dockerfile containers/Dockerfile.miser \
  containers/environment.yml containers/environment_miser.yml \
  containers/multiqc_config.yml containers/build_containers.sh \
  test/validate_preflight.sh test/samplesheet_2sample.csv \
  test/samplesheet_paired.csv test/VALIDATION_RUNBOOK.md
do
  [ -f "$f" ] && echo "  OK  $f" || echo "  MISSING  $f"
done
```

### Step 11 — Stage everything

```bash
git add .
```

Review what is staged:

```bash
git status
```

Expected: 32 files under "Changes to be committed".

Make sure these are **not** present:
```
work/      results/     *.sif      *.img
.nextflow/ *.bam        *.fa       *.db
```

If any appear, they are not covered by `.gitignore` — run
`git rm --cached <filename>` for each one.

### Step 12 — Commit

```bash
git commit -m "v1.0.0 — complete rewrite, initial public release

Long-read RNA-seq: MisER → IsoQuant → Oarfish → swish (fishpond)

New modules: SAMTOOLS_STATS, CLEAN_BAM, ISOQUANT_QC, SWISH_PLOTS, MULTIQC
New configs: sge.config, slurm.config, lsf.config (scheduler-agnostic)
New scripts: isoquant_qc.py, swish_plots.R (publication figures)
New tests:   validate_preflight.sh, VALIDATION_RUNBOOK.md

Key fixes:
- conditions CSV header sample_name → sample_id (SWISH was always failing)
- environment_miser.yml samtools=1. → samtools=1.21
- miser.nf scratch resolver: /scratch0 tried before \$TMPDIR
- samtools_sort_index.nf: dynamic -m (was hardcoded 3G)
- clean_bam.nf: parallel per-sample (was sequential for-loop in collected process)
- pair/batch samplesheet columns now flow to conditions CSV
- all bin/ scripts staged via channel (Singularity autoMounts=false safe)"
```

Expected:
```
[main xxxxxxx] v1.0.0 — complete rewrite, initial public release
 32 files changed, XXXX insertions(+), XX deletions(-)
```

### Step 13 — Push

```bash
git push origin main
```

Expected:
```
Enumerating objects: XX, done.
...
main -> main
```

Open https://github.com/dariusng28-collab/LoReRNA in your browser.
You should see all 32 files and the README rendered on the front page.

### Step 14 — Tag the release

```bash
git tag -a v1.0.0 -m "LoReRNA v1.0.0 — initial public release"
git push origin v1.0.0
```

On GitHub:
1. Go to the repo → **Releases** → **Draft a new release**
2. Select tag **v1.0.0** from the dropdown
3. Title: `LoReRNA v1.0.0`
4. Copy the contents of `CHANGES.md` into the description
5. Click **Publish release**

---

## Part 3 — HPC cluster: validate

### Step 15 — Clone the repo on the cluster

```bash
git clone https://github.com/dariusng28-collab/LoReRNA.git
cd LoReRNA
```

### Step 16 — Check Nextflow version

```bash
nextflow -version
# Must be >= 23.10.0
```

### Step 17 — Run pre-flight checks

```bash
bash test/validate_preflight.sh
```

Expected final lines:
```
 Pre-flight results
   PASS : 10
   FAIL : 0
 All checks passed. Safe to proceed to Stage 1: stub run.
```

**Do not submit any jobs if this shows any FAIL.**

### Step 18 — Fill in the test samplesheet

```bash
nano test/samplesheet_2sample.csv
```

Replace both `/path/to/your/...` lines with absolute paths to
one Control BAM and one KD BAM from your existing data.

### Step 19 — Stub run (no jobs, ~2 min)

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  -c conf/validate.config \
  --samplesheet test/samplesheet_2sample.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  -stub-run
```

**Must see** `LORERNA:CLEAN_BAM` in the process list.

### Step 20 — 2-sample validation run (~6–8 h)

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  -c conf/validate.config \
  --samplesheet test/samplesheet_2sample.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD
```

**Four checks when it completes:**

```bash
# 1. MisER used node-local scratch
grep 'Scratch resolver' lorerna_validation/01_miser/*/miser.log
# Good: /scratch0/... or /scratch/... or $SLURM_TMPDIR

# 2. Correct conditions CSV header
head -1 lorerna_validation/06_swish/conditions.csv
# Must be: sample_id,condition,pair,batch

# 3. SWISH produced result CSVs
ls lorerna_validation/06_swish/results/*.csv

# 4. MultiQC report exists
ls lorerna_validation/07_multiqc/multiqc_report.html
```

### Step 21 — Full cohort

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  --samplesheet /path/to/full/samplesheet.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  --outdir /path/to/final/results \
  -resume
```

---

## Troubleshooting

| Step | Failure | Fix |
|------|---------|-----|
| 3–4 | Docker build fails | Check Docker daemon is running; check internet for conda downloads |
| 4 | Push denied | Run `docker login`; confirm repo exists on Docker Hub |
| 6 | `git rm -rf .` leaves files | Run `git status` and manually `git rm` anything remaining |
| 10 | File count < 32 | Run the checklist loop above to find which file is missing |
| 11 | Unexpected files in `git status` | Run `git rm --cached <file>` for each; verify `.gitignore` is staged |
| 17 | Any FAIL in pre-flight | Fix the specific check; re-run until PASS: 10 |
| 19 | CLEAN_BAM absent | Check `include { CLEAN_BAM }` and `CLEAN_BAM(...)` in main.nf |
| 20 | MisER scratch on /tmp | Add `scratch = false`, `containerOptions = "-B /scratch0"`, `tscratch=60G` to `withName: 'MISER'` in conf/sge.config |
| 20 | SWISH `sample_id not found` | Check `seed:` line in main.nf — must be `"sample_id,condition,pair,batch\n"` |
| 20 | SWISH_PLOTS R package missing | Container `lorerna:1.0.0` not pulled; check Docker Hub tag exists; check `container_lorerna` in nextflow.config |
| 20 | MultiQC report missing | Check `work/` for the MULTIQC process log |
