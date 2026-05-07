# LoReRNA

**Long-read RNA-seq analysis pipeline: micro-exon rescue → transcript assembly → isoform quantification → differential analysis**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://conda.io/miniconda.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

**LoReRNA** is a [Nextflow](https://www.nextflow.io/) DSL2 pipeline for end-to-end analysis of Oxford Nanopore direct RNA-seq data. It takes genome-aligned BAMs as input and produces differential transcript expression (DTE), differential transcript usage (DTU), and differential gene expression (DGE) results.

The pipeline integrates:

1. **[MisER](https://github.com/vpc-ccg/miser)** — rescues micro-exons missed by standard long-read aligners
2. **[IsoQuant](https://github.com/ablab/IsoQuant)** — builds per-sample transcript models from corrected alignments
3. **[gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)** — merges all sample transcript models into a single non-redundant reference
4. **[oarfish](https://github.com/COMBINE-lab/oarfish)** — quantifies transcript abundance with bootstrap inferential replicates
5. **[fishpond/swish](https://bioconductor.org/packages/fishpond/)** — non-parametric differential analysis accounting for mapping uncertainty

## Pipeline summary

```
Genome-aligned BAMs (per sample)
         │
         ▼
  ┌─────────────┐
  │    MisER    │  Micro-exon rescue and splice-junction correction
  └──────┬──────┘
         │
         ▼
  ┌─────────────────────┐
  │  SAMtools sort+index │  Coordinate-sort for IsoQuant input
  └──────────┬──────────┘
             │
             ▼
  ┌─────────────┐
  │  IsoQuant   │  Per-sample transcript model building
  └──────┬──────┘
         │ transcript_models.gtf (× N samples)
         ▼
  ┌──────────────────────────────────┐
  │  Prepare oarfish reference       │  Runs once after all samples complete
  │  gffcompare → gffread → tx2gene  │  Builds shared transcriptome FASTA
  └────────────────┬─────────────────┘
                   │
                   ▼
  ┌─────────────┐
  │   oarfish   │  Per-sample quantification (100 bootstraps)
  └──────┬──────┘
         │ .quant + .infreps.pq (× N samples)
         ▼
  ┌─────────────┐
  │    swish    │  DTE · DTU · DGE via fishpond
  └──────┬──────┘
         │
         ▼
  results/06_swish/
    DTE / DTU / DGE  →  CSV + PDF per analysis
```

## Quick start

1. **Install Nextflow** (≥ 23.10.0):

   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. **Install a container runtime** — Docker, Singularity, or conda.

3. **Run the test profile** to verify your setup:

   ```bash
   # Docker
   nextflow run dariusng28-collab/LoReRNA -profile test,docker

   # Singularity
   nextflow run dariusng28-collab/LoReRNA -profile test,singularity
   ```

4. **Run with your own data:**

   ```bash
   nextflow run dariusng28-collab/LoReRNA \
       -profile docker \
       --samplesheet samplesheet.csv \
       --reference_fasta /path/to/GRCh38.fa \
       --reference_gtf   /path/to/gencode.v45.gtf \
       --reference_bed12 /path/to/gencode.v45.bed12 \
       --outdir results
   ```

## Samplesheet

Create a CSV with the following columns:

```csv
sample_name,condition,bam
Control1,Control,/data/Control1.bam
Control2,Control,/data/Control2.bam
Control3,Control,/data/Control3.bam
KD1,KD,/data/KD1.bam
KD2,KD,/data/KD2.bam
KD3,KD,/data/KD3.bam
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample_name` | ✅ | Unique sample identifier. Used as the filename prefix throughout the pipeline. |
| `condition` | ✅ | Group label. Must match `--swish_condition_a` and `--swish_condition_b` exactly (case-sensitive). |
| `bam` | ✅ | Absolute path to a genome-aligned BAM (minimap2 spliced alignment). |
| `pair` | Optional | Integer pair label (1, 2, …). Add this column to activate a **paired** design — no flags needed, auto-detected. |
| `batch` | Optional | Batch label. Add this column to activate a **batch-corrected** design — auto-detected. |

### Paired design

Add a `pair` column linking one sample per condition into a pair:

```csv
sample_name,condition,bam,pair
Control1,Control,/data/Control1.bam,1
Control2,Control,/data/Control2.bam,2
Control3,Control,/data/Control3.bam,3
KD1,KD,/data/KD1.bam,1
KD2,KD,/data/KD2.bam,2
KD3,KD,/data/KD3.bam,3
```

### Multi-condition experiment

Include all conditions in one samplesheet. Run the pipeline once per comparison — `-resume` caches all quantification steps so only the R analysis reruns:

```csv
sample_name,condition,bam
Control1,Control,/data/Control1.bam
Control2,Control,/data/Control2.bam
Control3,Control,/data/Control3.bam
KD1,KD,/data/KD1.bam
KD2,KD,/data/KD2.bam
KD3,KD,/data/KD3.bam
STM1,STM,/data/STM1.bam
STM2,STM,/data/STM2.bam
STM3,STM,/data/STM3.bam
```

```bash
# KD vs Control
nextflow run dariusng28-collab/LoReRNA -profile docker \
    --samplesheet samplesheet.csv ... \
    --swish_condition_a Control --swish_condition_b KD \
    --outdir results_KD

# STM vs Control — quantification steps are cached
nextflow run dariusng28-collab/LoReRNA -profile docker \
    --samplesheet samplesheet.csv ... \
    --swish_condition_a Control --swish_condition_b STM \
    --outdir results_STM -resume
```

## Running the pipeline

### Minimal run

```bash
nextflow run dariusng28-collab/LoReRNA \
    -profile docker \
    --samplesheet     samplesheet.csv \
    --reference_fasta /path/to/genome.fa \
    --reference_gtf   /path/to/annotation.gtf \
    --reference_bed12 /path/to/annotation.bed12 \
    --outdir          results
```

### Using a params file (recommended for regular use)

```bash
nextflow run dariusng28-collab/LoReRNA \
    -profile singularity \
    -params-file params/my_experiment.yml
```

`params/my_experiment.yml`:
```yaml
samplesheet:     "/data/samplesheet.csv"
reference_fasta: "/data/GRCh38.fa"
reference_gtf:   "/data/gencode.v45.gtf"
reference_bed12: "/data/gencode.v45.bed12"
outdir:          "/data/results"

isoquant_complete_genedb: true
swish_condition_a: "Control"
swish_condition_b: "KD"
```

### HPC clusters

Site-specific configuration is handled via `-c` rather than `-profile`, following the [nf-core site config convention](https://nf-co.re/docs/usage/configuration). A site config sets the executor, queue, memory format, and container bind mounts — the pipeline itself does not need to change.

```bash
nextflow run dariusng28-collab/LoReRNA \
    -profile singularity \
    -c conf/your_cluster.config \
    --samplesheet samplesheet.csv \
    ... 
```

A template site config is provided at `conf/cluster_template.config`. See [docs/site_config.md](docs/site_config.md) for instructions on adapting it to your cluster.

### Resume a failed run

Nextflow caches completed steps. Pass `-resume` to skip them on re-run:

```bash
nextflow run dariusng28-collab/LoReRNA -profile docker ... -resume
```

## Parameters

### Input / output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `null` | Path to samplesheet CSV **(required)** |
| `--reference_fasta` | `null` | Genome FASTA **(required)** |
| `--reference_gtf` | `null` | GTF/GFF annotation — provide exactly one of `--reference_gtf` or `--reference_db` **(required)** |
| `--reference_db` | `null` | Pre-built gffutils DB — faster for repeated runs; provide exactly one of `--reference_gtf` or `--reference_db` |
| `--reference_bed12` | `null` | BED12 annotation for MisER **(required)** |
| `--outdir` | `results` | Output directory |

### Sequencing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--seq_tech` | `ont-drna` | Sequencing technology preset for oarfish/minimap2. Use `ont-cdna` for cDNA libraries. |

### IsoQuant

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--isoquant_complete_genedb` | `false` | Recommended `true` for complete annotations (GENCODE, Ensembl) |
| `--isoquant_data_type` | `nanopore` | IsoQuant data type flag |
| `--isoquant_args` | `--count_exons` | Additional IsoQuant arguments |

### oarfish

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--oarfish_num_bootstraps` | `100` | Bootstrap replicates for swish DTU — 100 required for stable estimates |
| `--oarfish_filter_group` | `no-filters` | Multi-mapper handling — EM assigns fractional counts probabilistically |
| `--oarfish_args` | `--model-coverage` | Additional oarfish arguments — `--model-coverage` corrects 3′ positional bias |

### Differential analysis (swish)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--swish_condition_a` | `Control` | Reference condition — log2FC denominator |
| `--swish_condition_b` | `KD` | Test condition — log2FC numerator. log2FC > 0 means higher in condition_b |
| `--swish_min_count` | `10` | Minimum mean count for `labelKeep()` filter |
| `--swish_min_n` | `3` | Minimum number of samples passing count filter |
| `--swish_nperms` | `100` | Permutations for Wilcoxon rank-sum test |
| `--swish_alpha` | `0.05` | FDR threshold for result CSV output |

### Resource caps

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `72.h` | Maximum walltime per process |

## Outputs

```
results/
├── 01_miser/
│   └── <sample>/
│       ├── <sample>.miser.bam
│       ├── <sample>.missed_small.bed
│       └── <sample>.miser.log
│
├── 02_sorted_bam/
│   └── <sample>/
│       ├── <sample>.miser.sorted.bam
│       └── <sample>.miser.sorted.bam.bai
│
├── 03_isoquant/
│   └── <sample>/
│       ├── <sample>.transcript_counts.tsv
│       ├── <sample>.gene_counts.tsv
│       └── <sample>.read_assignments.tsv.gz
│
├── 04_oarfish_reference/
│   ├── merged_expressed_transcripts.fa
│   ├── tx2gene.tsv
│   ├── gffcmp.combined.gtf
│   └── <sample>.clean.bam (× N)
│
├── 05_oarfish/
│   └── <sample>/
│       ├── <sample>.quant
│       ├── <sample>.infreps.pq
│       ├── <sample>.meta_info.json
│       └── <sample>.ambig_info.tsv
│
├── 06_swish/
│   ├── results/
│   │   ├── DTE_all_significant.csv
│   │   ├── DTE_upregulated.csv
│   │   ├── DTE_downregulated.csv
│   │   ├── DTU_all_significant.csv
│   │   ├── DTU_increased_usage.csv
│   │   ├── DTU_decreased_usage.csv
│   │   ├── DGE_all_significant.csv
│   │   ├── DGE_upregulated.csv
│   │   └── DGE_downregulated.csv
│   ├── plots/
│   │   ├── DTE_plots.pdf
│   │   ├── DTU_plots.pdf
│   │   └── DGE_plots.pdf
│   └── logs/
│       └── swish_run_YYYYMMDD_HHMMSS.log
│
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── dag.html
    └── trace.txt
```

### Result CSV columns

**DTE / DGE**

| Column | Description |
|--------|-------------|
| `transcript_id` / `gene_id` | Feature identifier |
| `gene_id` | Parent gene (DTE only) |
| `log10mean` | Mean log10 scaled count |
| `log2FC` | Median fold change over bootstraps (`condition_b / condition_a`) |
| `stat` | Wilcoxon rank-sum statistic |
| `pvalue` | Permutation-based p-value |
| `qvalue` | BH-adjusted q-value |
| `meanInfRV` | Mean inferential relative variance (DTE only) |

**DTU**

Same as DTE but `log2FC` is the fold change of the isoform **proportion** [0,1], not raw counts. A gene can have one isoform with increased and another with decreased usage simultaneously — use `gene_id` to group them.

## Statistical background

**Why swish?**
swish uses a non-parametric Wilcoxon rank-sum test over bootstrap replicates, propagating mapping uncertainty (inferential variance) into the test. This is important for isoform-level analysis where many reads align to multiple transcripts — conventional tools that ignore this uncertainty produce inflated false discovery rates.

**Why 100 bootstraps?**
oarfish generates bootstrap quantification replicates by sampling from its EM uncertainty distribution. swish treats these as inferential replicates. 100 per sample provides stable variance estimates for balanced designs (≥ 3 samples per condition).

**log2FC direction**
`log2FC = log2(condition_b / condition_a)`. Positive values indicate higher expression or usage in `condition_b`. Set `--swish_condition_a` and `--swish_condition_b` to match your experimental design.

**q-value banding**
Volcano plots show horizontal bands — many features sharing the same q-value. This is expected with swish's plug-in permutation approach and is not a bug. Within a q-value tier, features are ranked by `|log2FC|`.

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `InvalidData(DuplicateTag(Tag("ts")))` in oarfish log | ts tag not stripped from MisER BAM | Check PREPARE_OARFISH_REFERENCE completed successfully |
| `More than 5% of transcripts unmatched in tx2gene` | tx2gene built from a different reference | Ensure tx2gene and oarfish reference come from the same gffcompare run |
| SWISH OOM (exit 137) | infRep matrices exceed available memory | Increase `process_high` memory, or reduce `--oarfish_num_bootstraps` |
| IsoQuant GTF not found | Wrong output path | Ensure `--isoquant_complete_genedb` is set correctly for your annotation |
| MisER fails with exit 1 | Insufficient scratch space | Set `--miser_scratch_root` to a path with ≥ 40 GB free |
| Samples queued indefinitely on HPC | Per-core memory request exceeds node capacity | Lower `--max_memory` in your site config or reduce process label memory |

## Credits

LoReRNA was developed by the [VYP Lab](https://github.com/dariusng28-collab), UCL.

Pipeline developed by Darius.

## Citation

If you use LoReRNA in your research, please cite the following tools:

- **MisER**: [citation pending]
- **IsoQuant**: Prjibelski et al., *Nature Biotechnology* (2023). https://doi.org/10.1038/s41587-022-01555-3
- **oarfish**: Srivastava et al., *bioRxiv* (2024). https://doi.org/10.1101/2024.02.12.580051
- **fishpond/swish**: Zhu et al., *Nucleic Acids Research* (2019). https://doi.org/10.1093/nar/gkz622
- **Nextflow**: Di Tommaso et al., *Nature Biotechnology* (2017). https://doi.org/10.1038/nbt.3820

## Licence

LoReRNA is released under the [MIT licence](LICENSE).
