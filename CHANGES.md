# LoReRNA — Changelog

## Unreleased — nf-core-quality hardening

Structural and quality work on top of v1.0.0. **No change to the analysis**:
every tool invocation is byte-identical, which was verified by resuming a
completed run and observing all 13 processes report `cached`, and by diffing
the published output tree (189 files, no differences).

### Repository layout

- `main.nf` reduced to a thin entry point; the workflow moved to
  `workflows/lorerna.nf`.
- Process resources split into `conf/base.config`; container routing and all
  13 `publishDir` declarations centralised in `conf/modules.config`.
- Example site configs for SGE, SLURM and LSF; `conf/ucl.config` untracked.

### Validation

- `nextflow_schema.json` and `assets/schema_input.json` added; parameters and
  the samplesheet are validated by nf-schema at launch.
- `.nf-core.yml` added. Checks specific to pipelines hosted under the nf-core
  organisation are disabled with a stated reason for each; all substantive
  checks remain enabled.

### Resources

- Migrated from a hand-written `check_max()` and `params.max_*` to native
  `process.resourceLimits`. Raises the required Nextflow version to 24.04.
  `conf/base.config` deliberately sets no limit — the ceiling is a
  site-specific fact set by the `-c` config.

### Provenance

- New `MISER_VERSION` and `SOFTWARE_VERSIONS` processes report the tool
  versions inside each of the two containers; the result appears as a
  "Software Versions" section in the MultiQC report and is published to
  `pipeline_info/software_versions.yml`.
- One dump per image rather than nf-core's per-process `versions.yml`: LoReRNA
  runs every process from one of two images, so per-process files would repeat
  the same versions while invalidating the task hash of all 13 modules.

### Testing

- `test/container_cli_contract.sh` asserts the tool CLIs and R packages match
  the flags the modules actually pass — the guard that would have caught
  oarfish renaming `--reference` to `--annotated`.
- Every module has a `stub:` block, and CI runs the whole DAG with
  `-profile stub -stub-run` to check channel wiring without tools or data.

### Results browser

- `lorerna-explorer/`, a Shiny application for the swish result tables,
  showing one gene across DGE, DTE and DTU together — the interpretive case
  the capped publication PDFs cannot cover. Hosted at
  <https://darius28.shinyapps.io/lorerna-explorer/>.

---

## v1.0.0 — Initial public release

This is the first public release of LoReRNA. It represents a complete,
production-validated pipeline for long-read direct RNA-seq analysis,
incorporating all fixes and scalability improvements developed during
internal testing.

---

### Bug fixes

**1. SWISH always failed — conditions CSV column name mismatch**
- `main.nf`: `collectFile` seed was `"sample_name,condition\n"`
- `bin/swish_analysis.R` checks `stopifnot(c("sample_id","condition") %in% colnames(...))`
- Fixed: seed is now `"sample_id,condition,pair,batch\n"`
- Impact: SWISH produced a cryptic R error on every run regardless of inputs

**2. MisER conda environment could not be built**
- `containers/environment_miser.yml`: `samtools=1.` is not a valid conda version specifier
- Fixed: `samtools=1.23.1`
- Impact: fresh conda installs of the MisER environment always failed

---

### Scalability

**3. MisER node-local scratch**
- `modules/local/miser.nf`: MisER does heavy random-access I/O on the input BAM. Running it against the NFS-staged symlink saturated shared I/O and caused timeouts.
- Fixed: the input BAM is copied to node-local scratch first (`/scratch0/$USER/...` by default, overridable per-site with `--miser_scratch_root`), and MisER runs against the local copy.
- A space guard requires at least 3× the input BAM size (minimum 40 G) of free scratch and fails fast with a clear message if the node cannot provide it.

**4. Site config memory scaling on retry**
- `conf/sge.config` and all site configs: labels used flat memory values (e.g. `memory = 64.GB`) with no `task.attempt` multiplier
- Fixed: all labels now use `memory = { N.GB * task.attempt }`
- Impact: OOM jobs retried with identical resources and failed again

**5. MISER NF work dir consumed local scratch before MisER ran**
- `conf/sge.config`: MISER inherited the global `scratch = true`
- Nextflow staged the work dir into `/tmp` before MisER started, exhausting the same partition MisER would use
- Fixed: `withName: 'MISER' { scratch = false }` — NF work dir stays on the shared filesystem; full local scratch is available to MisER

**6. samtools sort memory hardcoded**
- `modules/local/samtools_sort_index.nf`: `-m 3G` regardless of label or site config
- Fixed: `MEM_PER_THREAD_MB=$(( total_mb * 4 / 5 / n_threads ))` computed at runtime; minimum 512 MB

**7. pair/batch samplesheet columns silently ignored**
- `main.nf`: `meta` only carried `id` and `condition`; pair/batch columns were never read
- `bin/swish_analysis.R` had full paired/batch auto-detection but never received the data
- Fixed: samplesheet parser reads `pair` and `batch`; conditions CSV includes both columns; paired and batch swish designs now work

**8. BAM cleaning was sequential — critical path bottleneck**
- `modules/local/prepare_oarfish_ref.nf`: all-sample BAM cleaning ran as a sequential for-loop inside a `.collect()` process, blocked behind IsoQuant + gffcompare
- Fixed: extracted as `modules/local/clean_bam.nf` — a per-sample parallel process that starts as soon as `SAMTOOLS_SORT_INDEX` finishes, overlapping with IsoQuant + gffcompare
- Impact: for a 10-sample cohort, this removes ~3 hours from wall-clock time

**9. miser_qc.nf Python script not visible inside Singularity**
- `modules/local/miser_qc.nf` accessed `${projectDir}/bin/miser_qc_summary.py` directly
- With `singularity.autoMounts = false`, the project directory is not guaranteed to be mounted inside the container
- Fixed: script is staged via `path script` channel input — always present in the work directory

---

### New modules

**10. `modules/local/samtools_stats.nf`**
Runs `samtools flagstat`, `idxstats`, and `stats` on each sorted BAM.
Outputs are named for MultiQC auto-detection (`*_flagstat.txt`, etc.).

**11. `modules/local/isoquant_qc.nf` + `bin/isoquant_qc.py`**
Parses IsoQuant `read_assignments.tsv.gz` per sample to produce a
read-assignment-type breakdown (unique / ambiguous / inconsistent /
not_assigned). Outputs a MultiQC-compatible TSV per sample.

**12. `modules/local/swish_plots.nf` + `bin/swish_plots.R`**
Generates publication-quality ggplot2 figures from swish results:
- MA plots with ggrepel-labelled top hits
- Volcano plots
- pheatmap heatmaps (top-N transcripts/genes, Z-score log1p TPM)
- Summary panel (up/down counts across DTE/DTU/DGE)

**13. `modules/local/multiqc.nf`**
Aggregates all per-sample QC outputs (samtools, IsoQuant, MisER)
into a single interactive HTML report via MultiQC.

---

### New infrastructure

| File | Purpose |
|------|---------|
| `conf/sge.config` | SGE/UGE site config template |
| `conf/slurm.config` | SLURM site config template |
| `conf/lsf.config` | LSF/bsub site config template |
| `conf/validate.config` | Fast 2-sample validation config |
| `containers/build_containers.sh` | Local Docker build + push script |
| `containers/multiqc_config.yml` | MultiQC custom data sections |
| `containers/environment.yml` | Added multiqc, r-ggplot2, r-ggrepel, r-pheatmap, r-patchwork, r-viridis, r-rcolorbrewer |
| `test/validate_preflight.sh` | Login-node pre-flight script (static checks, no cluster jobs) |
| `test/VALIDATION_RUNBOOK.md` | 4-stage validation guide |
| `test/samplesheet_2sample.csv` | Minimal validation samplesheet template |
| `test/samplesheet_paired.csv` | Paired design samplesheet template |
| `.gitignore` | Excludes work dirs, SIFs, results, BAMs |
| `LICENSE` | MIT licence |
