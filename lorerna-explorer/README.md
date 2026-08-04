# LoReRNA Explorer

A Shiny application for browsing the swish differential-expression tables that
LoReRNA writes to `results/06_swish/results/`. It loads the DGE, DTE and DTU
results for a contrast and allows navigation between them, with a per-gene
panel showing all three analyses together.

This is a browser for one pipeline's output, not a general
differential-expression toolkit. For broader interactive exploration of
RNA-seq results, see [iSEE](https://bioconductor.org/packages/iSEE/),
[IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/)
and [Glimma](https://bioconductor.org/packages/Glimma/).

## Introduction

LoReRNA already produces publication figures and a MultiQC report. This
application covers what those do not:

- Heatmaps are capped by `--swish_top_n` (default 30) out of tens of thousands
  of tested features. The remainder exists only in the CSVs.
- DGE, DTE and DTU are written to three separate PDFs, and nothing joins them.
- MultiQC reports whether the run succeeded, not what it found.

The per-gene cross-analysis panel is the part not otherwise available in the
pipeline output. A gene can show no change in total expression while its
isoforms shift markedly in both abundance and relative usage — the pattern that
differential transcript usage exists to detect, and one associated with cryptic
splicing in TDP-43 proteinopathy. That signature is invisible in a gene-level
volcano plot and requires DGE, DTE and DTU side by side to recognise.

## Requirements

R 4.1 or later, and:

| Package | Minimum |
| --- | --- |
| shiny | 1.7.0 |
| bslib | 0.6.0 |
| dplyr | 1.1.0 |
| tidyr | 1.2.0 |
| ggplot2 | 3.4.0 |
| plotly | 4.10.0 |
| DT | 0.20 |
| readr | 2.0.0 |

`bslib >= 0.6` is a hard requirement: `page_sidebar()`, `accordion()`,
`layout_columns()` and `full_screen` are not available before that release. A
guard at the top of `app.R` runs before any `library()` call, so an
unsatisfied dependency produces an actionable message rather than a missing
function error.

```r
install.packages(c("shiny", "bslib", "dplyr", "tidyr", "ggplot2", "plotly", "DT", "readr"))
```

## Usage

A hosted copy runs at
**<https://darius28.shinyapps.io/lorerna-explorer/>** and needs nothing
installed. Note that uploads to it are processed on a third-party server — see
the note below before using it for data you may not share externally.

From a local checkout:

```bash
Rscript -e 'shiny::runApp("lorerna-explorer")'
```

Directly from GitHub, without cloning:

```bash
Rscript -e 'shiny::runGitHub("LoReRNA", "dariusng28-collab", ref = "release/container-v1.1.0", subdir = "lorerna-explorer")'
```

The `ref` argument is required while the application resides on the release
branch, and can be dropped once that branch is merged to `main`.

Results are then loaded through the sidebar.

Files are read by whichever machine runs the application. Run locally, nothing
leaves your machine. Run from a hosted copy, uploads are processed on that
host's servers and held for the session. Do not upload data you are not
permitted to share with a third-party host — run it locally instead.

## Input

Any of the swish result tables. Both the current pipeline naming
(`DTE_full_results.csv`) and the earlier `*_all_results.csv` convention are
accepted: the analysis is determined from the `DGE`, `DTE` or `DTU` substring
in the filename, and can be overridden from the sidebar if the assignment is
wrong.

Required columns, after case-insensitive alias resolution: a gene identifier,
`log2FC`, `pvalue` and `qvalue`. All other columns are optional and are used
when present. The alias table covers common spellings from other differential
expression tools — `padj`, `FDR`, `adj.P.Val`, `logFC`, `log2FoldChange` and
others — so the application is not restricted to swish output. To support a
further tool, extend the alias lists in `ALIASES` rather than adding
special cases.

Files may be supplied two ways: through the sidebar, or by placing them in
`lorerna-explorer/data/`, which is read at startup. That directory is
gitignored.

Where several files claim the same analysis, the table with the most rows is
used, on the basis that full results supersede significant-only subsets. Every
loaded file is listed in the sidebar with its assigned analysis, row count and
the reason for the assignment, including files that were skipped and why.

`sig` and `direction` columns, where present, are deliberately ignored. Status
is recomputed from the thresholds set in the interface; retaining a pre-computed
verdict alongside would allow two sources of truth to disagree on screen.

## Interface behaviour

Two properties of swish output shape the defaults.

**Permutation floor.** `--swish_nperms` (default 100) imposes a hard lower
bound on any achievable p-value. Significant features accumulate at that bound,
so the q-value threshold has limited discriminating power and the volcano
plot's y-axis collapses to a line of points above a cloud. The MA plot is
therefore the default view, and the sidebar computes and reports the floor for
the loaded table. Increasing `--swish_nperms` and re-running is the only way to
obtain finer resolution; where that has been done, the volcano becomes the more
informative default.

**Effect size scale.** DTU fold changes are computed on isoform *fractions*
rather than expression, and are consequently an order of magnitude smaller than
DGE or DTE fold changes on the same data. A fixed `|log2FC|` default would
return nothing for DTU while being far too permissive for DGE, so thresholds
are rescaled from the loaded table whenever the active analysis changes.

## Performance

The result tables are large enough that a naive implementation becomes
unresponsive. Two measures prevent this and should be preserved:

- The overview plot is constructed from a purpose-built four-column frame, never
  from the full results table. Passing a complete table to `plot_ly()`
  serialises every column of every row into the page; on a transcript-level
  table this is several megabytes of JSON per redraw.
- Above 15,000, non-significant points are randomly subsampled for display.
  Significant points are never subsampled, and tables and downloads always use
  the complete table. The caption beneath the plot states how many points are
  hidden. That caption should not be removed: a thinned scatter plot that does
  not declare itself is misleading.

Threshold inputs are debounced at 400 ms, and reactive inputs are frozen while
being reprogrammed after a load, so the reactive chain evaluates once per
change rather than once per keystroke or twice per analysis switch.

## Verifying a change

The data-loading functions have no Shiny dependency, and can be sourced
directly out of `app.R` for testing without a refactor:

```r
src   <- readLines("lorerna-explorer/app.R")
start <- grep("^ANALYSES <- ", src)[1]
stop  <- grep("^options\\(shiny\\.maxRequestSize", src)[1] - 1L
eval(parse(text = paste(src[start:stop], collapse = "\n")), envir = globalenv())

store <- list()
for (f in list.files("path/to/06_swish/results", full.names = TRUE)) {
  store <- ingest(f, basename(f), store)
}
store_manifest(store)
d <- build_datasets(store)
```

From there, apply a q-value threshold and count passing features, affected
genes and reciprocal genes per analysis. Record those counts against a dataset
you control, and treat them as a regression baseline: a change that alters them
has changed behaviour.

In the running application:

1. Load all three tables together. The analysis selector should offer DGE, DTE
   and DTU.
2. Check the value boxes against the baseline counts for each analysis.
3. Open "Loaded files". Every supplied file should appear with its analysis,
   row count and assignment reason. Nothing should be silently absent.
4. Select a gene known to show isoform switching and confirm the detail panel
   renders all three analyses.
5. Switch between analyses and confirm thresholds rescale and the selected gene
   persists.
6. Download both CSVs and confirm the row counts match the value boxes.

## Limitations

- **One contrast at a time.** LoReRNA compares exactly two conditions per
  invocation, so an experiment with more than two conditions produces separate
  output directories that this application does not join. This is the principal
  functional gap.
- Sample-level counts are not read, so the gene panel shows summary statistics
  rather than the underlying points.
- No outbound links to Ensembl or other transcript annotation resources.
- No automated test suite in the repository. The loading functions are written
  without a Shiny dependency specifically so that they can be covered by
  `testthat`, and the reactive graph is exercisable with `shiny::testServer()`.

The server logic has been checked headlessly against real pipeline output:
file ingestion and analysis assignment, precedence when two files claim the
same analysis, manual reassignment, threshold rescaling between analyses,
subsampling behaviour, the cross-analysis gene panel including genes absent
from one or more tables, alias resolution for limma and DESeq2 column
spellings, tables with no expression column, and unreadable input.

Not covered by those checks, and therefore unverified: the rendered ggplot in
the gene panel, including faceting when only one analysis is loaded and the
monospace axis font, which is device-dependent; selection by clicking a point
in the overview plot or a row in either table, which depends on browser-side
events; and any platform other than Windows.

## Credits

Written for the [LoReRNA](https://github.com/dariusng28-collab/LoReRNA)
pipeline.

## Citations

Differential testing is performed by the pipeline, not by this application. If
you use the results, cite swish and fishpond as listed in the pipeline's
[`CITATIONS.md`](../CITATIONS.md).

This application is built with
[shiny](https://cran.r-project.org/package=shiny),
[bslib](https://cran.r-project.org/package=bslib),
[plotly](https://cran.r-project.org/package=plotly),
[DT](https://cran.r-project.org/package=DT) and the
[tidyverse](https://www.tidyverse.org/).
