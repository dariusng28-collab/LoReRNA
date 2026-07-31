# LoReRNA: Output

All results are written under `--outdir`, in numbered directories that follow
the order of the pipeline.

## Directory overview

| Directory | Produced by | Contents |
|---|---|---|
| `01_miser/` | MISER | Splice-corrected BAMs and micro-exon events |
| `01_miser_qc/` | MISER_QC, MISER_QC_MERGE | Per-sample and cohort micro-exon rescue QC |
| `02_sorted_bam/` | SAMTOOLS_SORT_INDEX, CLEAN_BAM | Sorted/indexed and cleaned BAMs |
| `03_isoquant/` | ISOQUANT | Per-sample transcript models and counts |
| `04_oarfish_reference/` | PREPARE_OARFISH_REFERENCE | Merged transcriptome and tx2gene map |
| `05_oarfish/` | OARFISH | Per-sample transcript quantification |
| `06_swish/` | SWISH, SWISH_PLOTS | Differential results and figures |
| `07_multiqc/` | MULTIQC, SAMTOOLS_STATS, ISOQUANT_QC | Aggregated QC report |
| `pipeline_info/` | Nextflow | Execution reports and trace |

---

## 01_miser/

MisER corrects misaligned splice junctions and rescues micro-exons.

- `<sample>.miser.bam` — splice-corrected alignments
- `<sample>.missed_small.bed` — rescued micro-exon events (read-level)
- `<sample>.miser.log` — full MisER run log

## 01_miser_qc/

- `<sample>/<sample>.rescued_microexons.tsv` — unique rescued exon catalogue
  with support counts and realignment scores
- `<sample>/<sample>.rescued_microexons.bed` — BED6 track for IGV/UCSC
- `<sample>/<sample>.rescue_metrics.tsv` — per-sample QC metrics
- `all_samples_rescue_metrics_mqc.tsv` — cohort table (one row per sample)
- `miser_rescue_mqc.json` — MultiQC input for the micro-exon rescue section

Key metrics: `pass_rate_pct` (fraction of read events passing MisER's
realignment filter), `unique_microexons`, and the single- vs multi-read
support split (a proxy for confidence).

## 02_sorted_bam/

- `<sample>.miser.sorted.bam` (+ `.bai`) — coordinate-sorted, indexed
- `<sample>.clean.bam` — secondary/supplementary alignments removed
  (`-F 2304`) and the duplicate `ts` tag stripped, ready for oarfish

## 03_isoquant/

One directory per sample:

- `<sample>.transcript_counts.tsv` — transcript-level counts
- `<sample>.gene_counts.tsv` — gene-level counts
- `<sample>.read_assignments.tsv.gz` — per-read assignment detail
- `<sample>.transcript_models.gtf` — discovered transcript models

## 04_oarfish_reference/

Built once from all samples' transcript models:

- `gffcmp.combined.gtf` — gffcompare structural merge across samples
- `gffcmp.stats` — merge statistics
- `merged_expressed_transcripts.fa` — transcriptome FASTA (gffread)
- `tx2gene.tsv` — transcript → gene mapping

Transcripts with no annotated gene fall back to their own transcript ID as
the gene name, so novel transcripts still appear in gene-level results.

## 05_oarfish/

One directory per sample:

- `<sample>.quant` — transcript counts (`tname`, `len`, `num_reads`)
- `<sample>.infreps.pq` — bootstrap replicates (Parquet), used by swish
- `<sample>.meta_info.json` — run metadata and mapping statistics
- `<sample>.ambig_info.tsv` — unique vs ambiguous read breakdown

## 06_swish/

Differential analysis via fishpond/swish. Three analyses are run:

- **DTE** — differential transcript expression
- **DTU** — differential transcript usage (isoform proportion shifts)
- **DGE** — differential gene expression (transcript counts aggregated)

```
06_swish/
├── results/
│   ├── {DTE,DTU,DGE}_full_results.csv      all tested features
│   ├── {DTE,DTU,DGE}_all_significant.csv   FDR < alpha
│   └── {DTE,DGE}_{up,down}regulated.csv
├── plots/                                   diagnostic PDFs
├── logs/                                    run log + MultiQC summary
├── publication_plots/                       see below
└── conditions.csv                           sample → condition mapping used
```

`publication_plots/` contains ggplot2 figures:

- `DTE_publication_plots.pdf`, `DTU_publication_plots.pdf`,
  `DGE_publication_plots.pdf` — p-value histogram, MA plot, volcano,
  expression heatmap and direction summary
- `summary_panel.pdf` — up/down counts across DTE, DTU and DGE

Result columns include `log2FC` (as `condition_b / condition_a`), `pvalue`,
`qvalue` and `log10mean`.

### Browsing these results interactively

The heatmaps are capped by `--swish_top_n`, and DTE, DTU and DGE are written to
separate PDFs. To explore the full tables, or to see a single gene across all
three analyses at once, `lorerna-explorer/` contains a Shiny application that
reads the CSVs in `06_swish/results/`:

```bash
Rscript -e 'shiny::runApp("lorerna-explorer")'
```

It runs locally and reads local files. See
[`lorerna-explorer/README.md`](../lorerna-explorer/README.md) for requirements
and usage.

## 07_multiqc/

- `multiqc_report.html` — the aggregated QC report (start here)
- `multiqc_report_data/` — parsed data behind the report
- `samtools/` — flagstat, idxstats and stats per sample
- `isoquant_qc/` — per-sample read-assignment breakdowns

The report covers alignment statistics (samtools), IsoQuant read-assignment
composition, MisER micro-exon rescue, and a swish results summary.

## pipeline_info/

Nextflow execution records: `timeline.html`, `report.html`, `dag.html` and
`trace.txt`. `trace.txt` is the quickest way to see per-process runtime,
memory and exit status — useful when tuning resources.
