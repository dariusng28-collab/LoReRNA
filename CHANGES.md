# LoReRNA — Changelog

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
- Fixed: `samtools=1.21`
- Impact: fresh conda installs of the MisER environment always failed

---

### Scalability

**3. Scratch resolver order**
- `modules/local/miser.nf`: `$TMPDIR` was tried before `/scratch0`
- With `scratch = true` in the site config, Nextflow sets `$TMPDIR → /tmp/JOBID/` (a small partition) which was accepted before `/scratch0` was ever checked
- Fixed: `/scratch0` is now the highest-priority automatic candidate
- The 40 G guard is now inside `_try_scratch()` so a full-but-present path falls through to the next candidate rather than aborting the job

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
| `test/validate_preflight.sh` | 10-check login-node pre-flight script |
| `test/VALIDATION_RUNBOOK.md` | 4-stage validation guide |
| `test/samplesheet_2sample.csv` | Minimal validation samplesheet template |
| `test/samplesheet_paired.csv` | Paired design samplesheet template |
| `.gitignore` | Excludes work dirs, SIFs, results, BAMs |
| `LICENSE` | MIT licence |
