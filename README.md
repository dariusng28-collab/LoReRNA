# Long-read RNA-seq Nextflow scaffold

This repository contains a `Nextflow DSL2` scaffold for a long-read RNA-seq workflow built around the study logic you described:

1. `MisER` corrects splice-junction artifacts and rescues missed micro-exons from sorted long-read genome BAMs.
2. `IsoQuant` reconstructs transcript models and quantifies transcript abundance from the corrected BAMs.
3. `Oarfish` quantifies transcript abundance with bootstrap uncertainty on an IsoQuant-derived transcriptome reference.
4. `swish` tests differential transcript expression and differential transcript usage using inferential replicates.

## Why the workflow is wired this way

Your step table is the right high-level story, but two tool-level details affect the concrete pipeline design:

- `Oarfish` quantifies against a transcriptome reference, not a genome-spliced BAM. In this scaffold, the corrected `MisER` BAM is converted back to FASTQ, and `Oarfish` is run in read-based mode against a transcriptome FASTA generated from `IsoQuant` transcript models.
- `swish` requires inferential replicates for valid uncertainty-aware testing. Because `IsoQuant` counts alone do not provide those replicates, this scaffold feeds `swish` from the `Oarfish` quantification layer rather than the raw `IsoQuant` count table.

That keeps the biological order you want, while making the statistical handoff explicit.

## Workflow graph

`sorted BAM -> MisER -> .miser.bam -> IsoQuant -> transcript_models.gtf + counts`

`transcript_models.gtf + genome FASTA -> transcriptome FASTA`

`.miser.bam -> FASTQ -> Oarfish -> per-sample quant + bootstrap inferential replicates`

`Oarfish quants + metadata + tx2gene -> swish -> DTE + DTU result tables`

## Expected inputs

Required parameters:

- `--samplesheet`
- `--reference_fasta`
- `--reference_gtf`

Optional but supported:

- `--reference_bed12`

The samplesheet must contain:

- `sample`
- `condition`
- `bam`

An example is provided in [assets/samplesheet.csv](/C:/Users/dariu/OneDrive/Documents/New project/assets/samplesheet.csv).

If `--reference_bed12` is not supplied, the pipeline derives a transcript-level BED12 annotation from the input GTF before running `MisER`, because `MisER` expects BED12 rather than GTF.

## Main outputs

The pipeline publishes results under `results/` by default:

- `results/01_miser/*.miser.bam`
- `results/02_isoquant/*.transcript_models.gtf`
- `results/02_isoquant/combined_transcript_counts.tsv`
- `results/02_isoquant/*.exon_grouped_counts.tsv` when exon counting is enabled
- `results/03_oarfish_reference/isoquant.transcriptome.fa`
- `results/03_oarfish_reference/tx2gene.tsv`
- `results/05_oarfish/*.quant`
- `results/05_oarfish/*.infreps.pq`
- `results/06_oarfish_merged/oarfish.transcript_counts.tsv`
- `results/06_oarfish_merged/oarfish.isoform_usage.tsv`
- `results/06_oarfish_merged/oarfish.infreps_manifest.tsv`
- `results/07_swish/swish.dte.tsv`
- `results/07_swish/swish.dtu.tsv`

## Important caveat about PSI tables

This scaffold produces an isoform-usage matrix from `Oarfish` transcript counts, which is directly usable for transcript-usage analysis in `swish`.

It does not currently emit event-level PSI tables for cassette exons, alternative 5'/3' splice sites, or retained introns, because native `Oarfish` output is transcript quantification rather than an event catalog. If event-level PSI is a hard requirement, the clean extension would be to add an explicit event-definition layer after `IsoQuant`, using either:

- `IsoQuant --count_exons` outputs plus a custom event summariser
- a dedicated event-based splicing tool alongside `Oarfish`

## Container support

A unified Docker image definition is provided at [containers/Dockerfile](/C:/Users/dariu/OneDrive/Documents/New project/containers/Dockerfile), with its dependency spec at [containers/environment.yml](/C:/Users/dariu/OneDrive/Documents/New project/containers/environment.yml).

Build it locally with:

```bash
docker build -t longread-rnaseq-miser-isoquant-oarfish-swish:0.1.0 -f containers/Dockerfile .
```

Then run the workflow with the `docker` profile:

```bash
nextflow run main.nf \
  -profile docker \
  --samplesheet assets/samplesheet.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_gtf /path/to/annotation.gtf \
  --outdir results
```

The default container image name is set in [nextflow.config](/C:/Users/dariu/OneDrive/Documents/New project/nextflow.config) as `longread-rnaseq-miser-isoquant-oarfish-swish:0.1.0`. Override it with `--container_image` if you push the image to a registry.

An `apptainer` profile is also defined for HPC environments that prefer Apptainer/Singularity-style execution.

## Tool assumptions

If you do not use the container, the following executables must be available in your runtime environment:

- `miser`
- `isoquant.py`
- `gffread`
- `samtools`
- `oarfish`
- `python`
- `Rscript`

The `swish` script additionally expects R packages such as `fishpond`, `SummarizedExperiment`, `arrow`, `readr`, `dplyr`, and `optparse`.

## Running without containers

```bash
nextflow run main.nf \
  --samplesheet assets/samplesheet.csv \
  --reference_fasta /path/to/genome.fa \
  --reference_gtf /path/to/annotation.gtf \
  --outdir results
```

## Parameter notes

Useful defaults are defined in [nextflow.config](/C:/Users/dariu/OneDrive/Documents/New project/nextflow.config).

Parameters you are most likely to adjust first:

- `--seq_tech` for `Oarfish` (`ont-cdna`, `ont-drna`, `pac-bio`, `pac-bio-hifi`)
- `--isoquant_data_type` for `IsoQuant` (`nanopore`, `pacbio_ccs`, etc.)
- `--oarfish_num_bootstraps` for `swish` compatibility
- `--miser_args`, `--isoquant_args`, and `--oarfish_args` for site-specific tuning

## Source notes

The scaffold logic reflects the tool interfaces described in:

- [Oarfish documentation](https://combine-lab.github.io/oarfish/)
- [IsoQuant documentation](https://ablab.github.io/IsoQuant/)
- [Swish vignette](https://bioconductor.org/packages/devel/bioc/vignettes/fishpond/inst/doc/swish.html)

The exact `MisER` command-line invocation varies across installations and packaging, so the `MISER` process is left as a thin wrapper around a configurable executable and argument string.
