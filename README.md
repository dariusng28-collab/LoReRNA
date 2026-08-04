# LoReRNA

**Long-read RNA-seq differential analysis pipeline**

LoReRNA is a [Nextflow](https://www.nextflow.io/) DSL2 pipeline for end-to-end
analysis of long-read direct RNA-seq data. It performs micro-exon rescue,
transcript assembly, transcript-level quantification with inferential replicates,
and differential transcript/gene expression and usage analysis.

```
Input BAMs (genome-aligned, e.g. from minimap2)
    │
    ▼
 MisER            micro-exon rescue + splice-junction correction
    │
    ▼
 samtools         sort + index  ──► samtools flagstat/stats (→ MultiQC)
    │
    ├──► CLEAN_BAM (parallel)   remove secondary reads + duplicate ts tags
    │        │
    ├──► IsoQuant (parallel)    per-sample transcript assembly
    │        │
    │    IsoQuant QC ──────────────────────────────────────► MultiQC
    │        │
    │    PREPARE_OARFISH_REFERENCE
    │    (gffcompare + gffread + tx2gene)
    │        │
    └────────┴──► Oarfish       transcript quantification with inferential replicates
                      │
                   swish (fishpond)   DTE · DTU · DGE
                      │
                      ├──► publication plots   MA · volcano · heatmap · summary
                      │
                      └──► MultiQC             unified QC report
```

---

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Samplesheet format](#samplesheet-format)
- [Parameters](#parameters)
- [Output structure](#output-structure)
- [Exploring the results](#exploring-the-results)
- [Site configuration](#site-configuration)
- [Building containers](#building-containers)
- [Validation](#validation)
- [Troubleshooting](#troubleshooting)

---

## Requirements

| Dependency | Version | Notes |
|------------|---------|-------|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.10.0 | Requires Java 11+ |
| [Singularity](https://sylabs.io/) / [Apptainer](https://apptainer.org/) | ≥ 3.8 | Recommended on HPC |
| [Docker](https://www.docker.com/) | any | For local runs |
| conda / mamba | any | Alternative to containers |

All bioinformatics tools (MisER, IsoQuant, oarfish, samtools, R/fishpond, MultiQC)
are packaged in containers — no manual tool installation needed.

---

## Installation

```bash
# Clone
git clone https://github.com/your-org/lorerna.git
cd lorerna

# Install Nextflow (skip if already installed)
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

# Verify
nextflow -version   # should be >= 23.10.0
```

Containers are pulled automatically on first run. To pre-pull onto an HPC cluster:

```bash
singularity pull docker://your-dockerhub/lorerna:1.1.0
singularity pull docker://your-dockerhub/miser:1.1.0
```

---

## Quick start

### Minimal run (Docker, local machine)

```bash
nextflow run main.nf \
  -profile docker \
  --samplesheet samplesheet.csv \
  --reference_fasta /path/to/GRCh38.genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  --outdir results
```

### HPC cluster (SGE/SLURM/LSF)

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \       # or slurm.config / lsf.config
  --samplesheet samplesheet.csv \
  --reference_fasta /path/to/GRCh38.genome.fa \
  --reference_db    /path/to/annotation.gffutils.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  --outdir /path/to/output
```

### Resume a failed or interrupted run

```bash
nextflow run main.nf [same options] -resume
```

Nextflow hashes all inputs. Any process whose inputs are unchanged reuses its
cached result — no recomputation.

---

## Samplesheet format

A CSV with a header row. Three columns are required; two are optional.

```
sample_name,condition,bam[,pair,batch]
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample_name` | ✓ | Unique sample ID. Alphanumeric plus `.` `_` `-` only. |
| `condition` | ✓ | Condition label. Must match `--swish_condition_a` or `--swish_condition_b`. |
| `bam` | ✓ | Absolute path to genome-aligned BAM. Must pass `samtools quickcheck`. |
| `pair` | optional | Pairing key (integer or string). Activates paired Wilcoxon swish. Samples sharing the same value across conditions are treated as matched pairs. |
| `batch` | optional | Batch covariate. Activates batch-aware swish. Mutually exclusive with `pair` — if both are present, `pair` takes precedence. |

### Standard two-condition

```csv
sample_name,condition,bam
Control_1,Control,/data/Control_1.sorted.bam
Control_2,Control,/data/Control_2.sorted.bam
Control_3,Control,/data/Control_3.sorted.bam
KD_1,KD,/data/KD_1.sorted.bam
KD_2,KD,/data/KD_2.sorted.bam
KD_3,KD,/data/KD_3.sorted.bam
```

### Paired design (e.g. same donor treated two ways)

```csv
sample_name,condition,bam,pair
Donor1_Control,Control,/data/Donor1_Control.bam,1
Donor1_KD,KD,/data/Donor1_KD.bam,1
Donor2_Control,Control,/data/Donor2_Control.bam,2
Donor2_KD,KD,/data/Donor2_KD.bam,2
```

### Batch correction (e.g. samples run on different days)

```csv
sample_name,condition,bam,batch
Control_1,Control,/data/Control_1.bam,batch1
KD_1,KD,/data/KD_1.bam,batch1
Control_2,Control,/data/Control_2.bam,batch2
KD_2,KD,/data/KD_2.bam,batch2
```

Paired and batch designs are **auto-detected** from non-empty samplesheet columns.
No additional pipeline flags are needed.

---

## Multi-condition experiments

LoReRNA compares exactly two conditions per run (`--swish_condition_a` vs `--swish_condition_b`).
For three or more conditions, run the pipeline once per pairwise comparison.
The `swish_analysis.R` script automatically filters the conditions CSV to only the two
conditions being compared, so you can safely pass a samplesheet containing all samples.

**Example: CTRL vs cond1 and CTRL vs cond2**

Run 1 — all samples in samplesheet, compare CTRL vs cond1:
```bash
nextflow run main.nf \
  --samplesheet samplesheet_all.csv \
  --swish_condition_a CTRL \
  --swish_condition_b cond1 \
  --outdir results_CTRL_vs_cond1
```

Run 2 — same full samplesheet, compare CTRL vs cond2 (uses `-resume` to reuse CTRL Oarfish cache):
```bash
nextflow run main.nf \
  --samplesheet samplesheet_all.csv \
  --swish_condition_a CTRL \
  --swish_condition_b cond2 \
  --outdir results_CTRL_vs_cond2 \
  -resume
```

With `-resume`, all CTRL sample steps (MisER → IsoQuant → Oarfish) are reused from the
first run's cache (~4 h total vs ~48 h from scratch for a 15 G BAM cohort).

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--samplesheet` | Path to samplesheet CSV |
| `--reference_fasta` | Genome FASTA (e.g. GRCh38 primary assembly) |
| `--reference_bed12` | Gene annotation in BED12 format (required by MisER) |
| `--reference_db` | Pre-built IsoQuant gffutils DB — provide this **or** `--reference_gtf`, not both |
| `--reference_gtf` | GTF annotation file — used by IsoQuant when no gffutils DB is available |
| `--swish_condition_a` | Reference condition (log2FC denominator) |
| `--swish_condition_b` | Test condition (log2FC numerator — positive = higher in `b`) |

### Commonly changed

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--swish_alpha` | `0.05` | FDR threshold |
| `--swish_min_count` | `10` | Minimum count filter before swish |
| `--swish_min_n` | `3` | Minimum samples per condition |
| `--swish_nperms` | `100` | Permutations for swish |
| `--swish_top_n` | `30` | Top-N features shown in heatmaps |
| `--oarfish_num_bootstraps` | `100` | Inferential replicates (minimum 50 recommended for DTU) |
| `--isoquant_complete_genedb` | `true` | Use `--complete_genedb` — correct for GENCODE/Ensembl |
| `--miser_scratch_root` | `null` | Override scratch path for MisER (auto-resolved if null) |
| `--publish_mode` | `copy` | Nextflow publishDir mode (`copy`, `symlink`, `link`) |

### Resource caps (set per site in conf/)

| Parameter | Default |
|-----------|---------|
| `--max_cpus` | `16` |
| `--max_memory` | `128.GB` |
| `--max_time` | `72.h` |

---

## Output structure

```
results/
├── 01_miser/
│   └── <sample>/
│       ├── <sample>.miser.bam           MisER-corrected BAM
│       ├── <sample>.missed_small.bed    micro-exon event table
│       └── <sample>.miser.log           full run log with scratch diagnostics
│
├── 01_miser_qc/
│   └── <sample>/
│       ├── <sample>.rescued_microexons.tsv
│       ├── <sample>.rescued_microexons.bed
│       └── <sample>.rescue_metrics.tsv
│
├── 02_sorted_bam/
│   ├── <sample>.miser.sorted.bam        sorted + indexed BAM (feeds IsoQuant)
│   ├── <sample>.miser.sorted.bam.bai
│   └── <sample>.clean.bam               filtered BAM (feeds Oarfish)
│
├── 03_isoquant/
│   └── <sample>/
│       ├── <sample>.transcript_counts.tsv
│       ├── <sample>.gene_counts.tsv
│       ├── <sample>.transcript_models.gtf
│       └── <sample>.read_assignments.tsv.gz
│
├── 04_oarfish_reference/
│   ├── merged_expressed_transcripts.fa  merged transcriptome FASTA
│   ├── tx2gene.tsv                      transcript → gene name mapping
│   ├── gffcmp.combined.gtf
│   └── gffcmp.stats
│
├── 05_oarfish/
│   └── <sample>/
│       ├── <sample>.quant               transcript counts + TPM
│       └── <sample>.infreps.pq          inferential replicates (parquet)
│
├── 06_swish/
│   ├── results/
│   │   ├── DTE_full_results.csv         all tested transcripts (use for custom plots)
│   │   ├── DTE_all_significant.csv      FDR < alpha
│   │   ├── DTE_upregulated.csv
│   │   ├── DTE_downregulated.csv
│   │   ├── DTU_*                        differential transcript usage
│   │   └── DGE_*                        differential gene expression
│   ├── plots/                           base-R QC plots (p-value histograms, MA, infRep)
│   ├── publication_plots/
│   │   ├── DTE_publication_plots.pdf    MA · volcano · heatmap · direction bar
│   │   ├── DTU_publication_plots.pdf    proportion shift · scatter · heatmap
│   │   ├── DGE_publication_plots.pdf    gene MA · volcano · gene heatmap
│   │   └── summary_panel.pdf           up/down counts across DTE/DTU/DGE
│   └── logs/
│       └── swish_analysis.log
│
├── 07_multiqc/
│   ├── multiqc_report.html              unified QC report — open in a browser
│   ├── multiqc_data/                    underlying TSV/JSON
│   ├── multiqc_plots/                   static PNG exports
│   ├── samtools/                        per-sample flagstat, idxstats, stats
│   └── isoquant_qc/                     per-sample read assignment breakdown
│
└── pipeline_info/
    ├── timeline.html                    per-process wallclock times
    ├── report.html                      CPU/memory/I/O resource usage
    ├── dag.html                         execution graph
    └── trace.txt                        raw per-task trace
```

---

## Exploring the results

`06_swish/publication_plots/` covers the top `--swish_top_n` features (default
30) and writes DTE, DTU and DGE to separate PDFs. For the full tables, or to
inspect one gene across all three analyses at once, `lorerna-explorer/`
contains a Shiny application that reads the CSVs in `06_swish/results/`:

A hosted copy runs at
**<https://darius28.shinyapps.io/lorerna-explorer/>** and needs nothing
installed.

To run it yourself:

```bash
Rscript -e 'shiny::runApp("lorerna-explorer")'
```

Or without a local checkout:

```bash
Rscript -e 'shiny::runGitHub("LoReRNA", "dariusng28-collab", ref = "release/container-v1.1.0", subdir = "lorerna-explorer")'
```

Run locally it reads local files and nothing leaves your machine; a hosted copy
processes uploads on its host's servers. Requirements, input format and known
limitations are documented in
[`lorerna-explorer/README.md`](lorerna-explorer/README.md).

---

## Site configuration

The pipeline separates scheduler-specific settings into `conf/` files.
All analysis logic lives in `modules/local/` and is scheduler-agnostic.

### Provided configs

| File | Scheduler | Notes |
|------|-----------|-------|
| `conf/sge.config` | SGE / UGE | Update queue, bind paths, cacheDir |
| `conf/slurm.config` | SLURM | Update partition, account |
| `conf/lsf.config` | LSF / bsub | Update queue, account, cacheDir |

### Adapting to a new cluster

1. Copy the config closest to your scheduler:
   ```bash
   cp conf/slurm.config conf/mycluster.config
   ```

2. Edit exactly three things:
   ```groovy
   queue    = 'your_partition'
   clusterOptions = "--account=your_account ..."

   // In withName: 'MISER' — change the local scratch resource flag:
   clusterOptions = "... --tmp=61440"      // SLURM: 60 G scratch in MB
   // or for SGE:
   clusterOptions = "... ,tscratch=60G"
   ```

3. Run with `-c conf/mycluster.config`.

### MisER local scratch

MisER requires ≥ 40 G of fast local scratch for BAM I/O. The pipeline
auto-resolves the best available path in this priority order:

1. `params.miser_scratch_root` — explicit override in your site config
2. `/scratch0/$USER` — SGE/UGE node-local (tscratch-allocated)
3. `$SLURM_TMPDIR` — SLURM `--tmp` allocation
4. `/scratch/$USER` — generic HPC scratch mount
5. `/local_scratch/$USER` — alternative HPC scratch mount
6. `$TMPDIR` — PBS/Torque/LSF fallback
7. `$PWD/miser_scratch` — shared filesystem last resort (slower, always succeeds)

Each candidate is only accepted if it **exists** and has **≥ 40 G free**.
A full-but-present path falls through to the next candidate automatically.

---

## Building containers

Containers are hosted on Docker Hub.
Build on a **local machine with Docker** — not on HPC nodes.

```bash
cd lorerna

# Build and sanity-check (no push)
bash containers/build_containers.sh

# Build and push to Docker Hub
bash containers/build_containers.sh --push

# Rebuild only the MisER image
bash containers/build_containers.sh --push --miser-only

# Specify a version tag
bash containers/build_containers.sh --push --version 1.1.0
```

After pushing, update `nextflow.config`:
```groovy
container_lorerna = 'docker://your-dockerhub/lorerna:1.1.0'
container_miser   = 'docker://your-dockerhub/miser:1.1.0'
```

### Container contents

| Image | Key tools |
|-------|-----------|
| `lorerna:1.1.0` | samtools 1.23.1 · IsoQuant 3.13.0 · oarfish 0.9.4 · R 4.4.2 + fishpond 2.12.0 · ggplot2 · pheatmap · MultiQC 1.25.1 |
| `miser:1.1.0` | samtools 1.23.1 · MisER 1.0.0 · pysam · parasail |

---

## Validation

Before running a full cohort, work through the four stages in
`test/VALIDATION_RUNBOOK.md`. The stages are designed so each one
catches a different class of issue at the lowest possible cost.

### Stage 0 — Pre-flight (login node, ~5 min, no jobs submitted)

```bash
bash test/validate_preflight.sh
# Expected output: PASS: 10  FAIL: 0
```

Inspects every changed `.nf` file and checks all fixes in-place.
Only move to Stage 1 when this returns zero failures.

### Stage 1 — Stub run (login node, ~2 min, no jobs submitted)

```bash
# Fill in real BAM paths in test/samplesheet_2sample.csv first
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  -c conf/validate.config \
  --samplesheet test/samplesheet_2sample.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  -stub-run
```

Check `LORERNA:CLEAN_BAM` appears in the execution summary.

### Stage 2 — 2-sample cluster run (~6–8 h)

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  -c conf/validate.config \
  --samplesheet test/samplesheet_2sample.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD
```

Key checks after completion:

```bash
# MisER used node-local scratch (not /tmp)
grep 'Scratch resolver' lorerna_validation/01_miser/*/miser.log

# conditions CSV has the correct header
head -1 lorerna_validation/06_swish/conditions.csv
# Expected: sample_id,condition,pair,batch

# SWISH produced result CSVs
ls lorerna_validation/06_swish/results/*.csv

# MultiQC report was generated
ls lorerna_validation/07_multiqc/multiqc_report.html
```

### Stage 3 — Full cohort

```bash
nextflow run main.nf \
  -profile singularity \
  -c conf/sge.config \
  --samplesheet samplesheet.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_db    /path/to/annotation.db \
  --reference_bed12 /path/to/annotation.bed12 \
  --swish_condition_a Control \
  --swish_condition_b KD \
  --outdir /path/to/final/results
```

---

## Troubleshooting

### MisER: "Not enough scratch space"

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `Scratch resolver: /tmp/...` | `scratch = false` missing for MISER | Add `withName: 'MISER' { scratch = false }` to your site config |
| `/scratch0` not visible in container | Bind mount missing | Add `containerOptions = "-B /scratch0"` in `withName: 'MISER'` |
| Scratch not allocated by scheduler | Resource flag missing | Add `tscratch=60G` (SGE) or `--tmp=61440` (SLURM) to MISER `clusterOptions` |
| All local options full | Busy node | Set `params.miser_scratch_root = '/path/to/large/filesystem'` |

### SWISH: `sample_id not found` (R stopifnot error)

The conditions CSV header is wrong. This is the bug fixed in v1.0.0.
Verify `main.nf` has:
```groovy
seed: "sample_id,condition,pair,batch\n"
```

### SWISH_PLOTS: R package not found

Ensure you are using the `lorerna:1.1.0` container which includes
`r-ggplot2`, `r-ggrepel`, `r-pheatmap`, `r-patchwork`, `r-viridis`.
Earlier container versions do not include these packages.

### IsoQuant QC: input file not found

Verify the glob pattern for `read_assignments.tsv.gz` in `isoquant.nf` matches
what IsoQuant actually produces in your version:
```bash
ls work/*/*/*.tsv.gz | head -5
```

### Memory errors on retry

Ensure your site config uses the `task.attempt` multiplier:
```groovy
memory = { 64.GB * task.attempt }   // doubles each retry: 64→128→192
# NOT:
memory = 64.GB                       // flat — never scales up
```

---

## Citation

If you use LoReRNA, please cite the tools it wraps:

- **MisER** — [github.com/comprna/MisER](https://github.com/comprna/MisER)
- **IsoQuant** — Prjibelski et al. *Nature Biotechnology* 2023
- **oarfish** — [github.com/COMBINE-lab/oarfish](https://github.com/COMBINE-lab/oarfish)
- **fishpond / swish** — Zhu et al. *NAR* 2019; Zhu et al. *Genome Biology* 2021
- **MultiQC** — Ewels et al. *Bioinformatics* 2016
- **Nextflow** — Di Tommaso et al. *Nature Biotechnology* 2017

---

## Contributing

Issues and pull requests are welcome.
When adding a new site config, only edit the scheduler-specific block —
do not modify `modules/local/` or `main.nf` for site-specific reasons.
