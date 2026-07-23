# LoReRNA: Usage

## Introduction

LoReRNA analyses long-read (ONT direct-RNA / cDNA, PacBio) RNA-seq from
aligned BAMs through to differential transcript and gene results:

```
MisER → samtools sort/index → CLEAN_BAM + IsoQuant → gffcompare/gffread
      → oarfish → swish (fishpond) → publication plots + MultiQC
```

Input is **genome-aligned BAM** (e.g. minimap2 output), not FASTQ.

## Samplesheet

Provide a CSV with `--samplesheet`. It is validated against
[`assets/schema_input.json`](../assets/schema_input.json) at launch.

| Column | Required | Description |
|---|---|---|
| `sample_name` | yes | Sample ID. Letters, numbers, `.`, `_`, `-` only (no spaces). |
| `condition` | yes | Group label. Must match `--swish_condition_a` / `--swish_condition_b`. |
| `bam` | yes | Path to the genome-aligned BAM. |
| `pair` | no | Pairing ID — activates a **paired** Wilcoxon test in swish. |
| `batch` | no | Batch ID — activates **batch-aware** swish. |

Minimal example:

```csv
sample_name,condition,bam
Control1,WT,/data/Control1.sorted.bam
Control2,WT,/data/Control2.sorted.bam
Knockdown1,KD,/data/Knockdown1.sorted.bam
Knockdown2,KD,/data/Knockdown2.sorted.bam
```

Paired design (e.g. matched donors):

```csv
sample_name,condition,bam,pair
Ctrl_donor1,WT,/data/Ctrl_donor1.bam,1
KD_donor1,KD,/data/KD_donor1.bam,1
Ctrl_donor2,WT,/data/Ctrl_donor2.bam,2
KD_donor2,KD,/data/KD_donor2.bam,2
```

If both `pair` and `batch` are set for a sample, `pair` takes precedence
(the pipeline logs a warning).

## Reference files

| Parameter | Required | Description |
|---|---|---|
| `--reference_fasta` | yes | Genome FASTA (e.g. GRCh38 primary assembly). |
| `--reference_bed12` | yes | BED12 annotation, used by MisER. |
| `--reference_gtf` | one of | Annotation GTF/GFF for IsoQuant. |
| `--reference_db` | one of | Pre-built gffutils database for IsoQuant. |

Supply **exactly one** of `--reference_gtf` or `--reference_db`. Using a
pre-built `.db` avoids IsoQuant rebuilding the database on every run.

## Running the pipeline

```bash
nextflow run main.nf \
  -profile singularity \
  --samplesheet samplesheet.csv \
  --reference_fasta /refs/GRCh38.primary_assembly.genome.fa \
  --reference_db    /refs/gencode.v42.gffutils.db \
  --reference_bed12 /refs/gencode.v42.annotation.bed12 \
  --swish_condition_a WT \
  --swish_condition_b KD \
  --outdir results
```

`--swish_condition_a` is the **reference** (denominator) and
`--swish_condition_b` the **comparison** (numerator), so log2FC is
reported as `B / A`.

### Profiles

| Profile | Effect |
|---|---|
| `singularity` | Run all processes in Singularity containers (recommended on HPC). |
| `docker` | Run in Docker containers. |
| `conda` | Build environments from `containers/environment*.yml` instead of containers. |
| `test` | Small resource caps and reduced bootstraps/permutations. |
| `standard` | Local executor, no scheduler. |

### Cluster / site configuration

Templates are provided for common schedulers; stack one with `-c`:

```bash
nextflow run main.nf -profile singularity -c conf/sge.config ...
```

`conf/sge.config`, `conf/slurm.config` and `conf/lsf.config` are **templates**
— update the queue/partition, account, `singularity.cacheDir` and bind mounts
for your site. Keep site-specific settings in your own config file rather than
committing them.

> **MisER scratch.** MisER copies the input BAM to node-local scratch to avoid
> saturating shared storage. The default path is `/scratch0` (SGE convention).
> On SLURM/LSF set `--miser_scratch_root` to your scheduler's local scratch
> (e.g. `$TMPDIR`), and make sure that path is bound into the container.

## Multi-condition experiments

LoReRNA compares **exactly two conditions per run**. For three or more
conditions, run it twice and reuse the cache.

**Run 1 — all samples, compare WT vs cond1.** Include every sample: the extra
condition still contributes to the shared transcriptome reference, and is
filtered out of the statistical comparison automatically.

```bash
nextflow run main.nf --samplesheet samplesheet_all.csv \
  --swish_condition_a WT --swish_condition_b cond1 \
  --outdir results_WT_vs_cond1
```

**Run 2 — compare WT vs cond2.** With `-resume`, MisER, IsoQuant and oarfish
results are reused; only swish re-runs.

```bash
nextflow run main.nf --samplesheet samplesheet_all.csv \
  --swish_condition_a WT --swish_condition_b cond2 \
  --outdir results_WT_vs_cond2 -resume
```

## Resuming

```bash
nextflow run main.nf ... -resume
```

Use a **single dash** (`-resume`). `--resume` is parsed as a pipeline
parameter and silently starts a fresh run. To resume a specific earlier run:

```bash
nextflow log                      # list runs, find the session ID
nextflow run main.nf ... -resume <session-id>
```

Editing a module invalidates only that process and its dependents; upstream
results stay cached.

## Key parameters

| Parameter | Default | Notes |
|---|---|---|
| `--seq_tech` | `ont-drna` | One of `ont-drna`, `ont-cdna`, `pac-bio`, `pac-bio-hifi`. Sets the oarfish/minimap2 preset. |
| `--oarfish_num_bootstraps` | `100` | Inferential replicates; required for swish DTU. |
| `--swish_min_count` | `10` | Minimum count to retain a feature. |
| `--swish_min_n` | `3` | Minimum samples meeting `min_count`. |
| `--swish_nperms` | `100` | Permutations for the swish test. |
| `--swish_alpha` | `0.05` | FDR threshold. |
| `--publish_mode` | `copy` | Use `copy` for production so results survive `nextflow clean`. |

Full list with types and defaults: `nextflow run main.nf --help`.

## Resource limits

Per-process requests are capped by `--max_cpus`, `--max_memory` and
`--max_time`. Processes retry with increased memory on OOM-type exit codes
(up to `maxRetries`).

MisER is the most demanding step: it needs roughly **3× the input BAM size**
in node-local scratch (minimum 40 GB) and can run for many hours on large
BAMs. Request walltime accordingly.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `Missing required parameter` at launch | A required param is unset; check the message from `validateParameters()`. |
| swish errors about conditions | `--swish_condition_a/b` don't match values in the samplesheet `condition` column. |
| MisER fails on scratch space | Not enough node-local scratch; request more or set `--miser_scratch_root`. |
| `command not found` inside a container | Container PATH stripped; ensure the site config sets `beforeScript` or uses the provided images. |
| Everything re-runs on `-resume` | You used `--resume`, changed a module, or the work directory was cleaned. |
