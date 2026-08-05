# LoReRNA Explorer
#
# Interactive browser for swish results from the LoReRNA pipeline.
# Load any of DGE_*, DTE_*, DTU_* from results/06_swish/results/.
#
# Performance notes, because the tables are large enough to matter:
#   - plot_ly() serialises every column of whatever data frame it is given,
#     so the overview plot is built from a four-column frame, not from the
#     full results table.
#   - The non-significant cloud is thinned for display above a threshold.
#     Significant features are never thinned.
#   - Threshold inputs are debounced, so typing does not trigger a redraw
#     per keystroke.
#   - Inputs are frozen while they are being reprogrammed after a load, so
#     everything downstream renders once rather than twice.

# ---- Setup ------------------------------------------------------------------

# Fail with a readable message rather than "could not find function".
local({
  required <- list(shiny = "1.7.0", bslib = "0.6.0", dplyr = "1.1.0",
                   tidyr = "1.2.0", ggplot2 = "3.4.0", plotly = "4.10.0",
                   DT = "0.20", readr = "2.0.0", data.table = "1.14.0")
  bad <- character(0)
  for (pkg in names(required)) {
    if (!requireNamespace(pkg, quietly = TRUE) ||
        utils::packageVersion(pkg) < required[[pkg]]) {
      bad <- c(bad, sprintf("%s (>= %s)", pkg, required[[pkg]]))
    }
  }
  if (length(bad)) {
    stop("Missing or outdated packages:\n  ", paste(bad, collapse = "\n  "),
         "\n\nInstall with:\n  install.packages(c(",
         paste(sprintf('"%s"', names(required)), collapse = ", "), "))",
         call. = FALSE)
  }
})

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(readr)
# Attached last on purpose. data.table masks between/first/last from dplyr;
# none of those are used here, and the GTF work needs data.table's syntax.
library(data.table)

# ---- Data loading -----------------------------------------------------------
# Inlined so app.R runs standalone: it can be copied to a machine on its own,
# or launched with shiny::runGitHub(), without a supporting R/ directory.
# Extracting this block to R/ingest.R is worthwhile only once there is a test
# suite to import it -- these functions have no Shiny dependency for exactly
# that reason.

ANALYSES <- c("DGE", "DTE", "DTU")

ANALYSIS_LABEL <- c(
  DGE = "Differential gene expression",
  DTE = "Differential transcript expression",
  DTU = "Differential transcript usage"
)

EFFECT_LABEL <- c(
  DGE = "log2 fold change in gene expression",
  DTE = "log2 fold change in transcript expression",
  DTU = "log2 fold change in isoform fraction"
)

STATUS_LEVELS <- c("Up", "Down", "Not significant")

# Different tools spell the same column differently. Match case-insensitively
# against known aliases rather than demanding one exact name.
ALIASES <- list(
  transcript = c("transcript_id", "tx_id", "txname", "tx_name", "isoform_id"),
  gene       = c("gene_id", "gene_name", "gene", "geneid", "symbol", "gene_symbol"),
  gene_name  = c("gene_name", "symbol", "gene_symbol", "gene_id", "gene"),
  log2FC     = c("log2FC", "log2fc", "log2FoldChange", "logFC", "log2_fold_change"),
  pvalue     = c("pvalue", "p_value", "pval", "p", "PValue"),
  qvalue     = c("qvalue", "q_value", "qval", "padj", "FDR", "adj.P.Val", "locfdr"),
  log10mean  = c("log10mean", "log10_mean"),
  mean_lin   = c("mean_tpm", "baseMean", "meanCount", "AveExpr", "mean"),
  stat       = c("stat", "statistic", "t", "LR", "log2FC_stat"),
  infrv      = c("meanInfRV", "infRV", "mean_inf_rv")
)

pick_col <- function(df, candidates) {
  nm <- tolower(names(df))
  for (cand in candidates) {
    i <- match(tolower(cand), nm)
    if (!is.na(i)) return(names(df)[i])
  }
  NA_character_
}

resolve_cols <- function(df) lapply(ALIASES, function(a) pick_col(df, a))

missing_required <- function(cols) {
  need <- c("gene", "log2FC", "pvalue", "qvalue")
  need[vapply(need, function(k) is.na(cols[[k]]), logical(1))]
}

# Returns the analysis plus the reason, so the UI can show both.
classify_analysis <- function(name, cols) {
  n <- toupper(name)
  for (a in ANALYSES) {
    if (grepl(a, n, fixed = TRUE)) return(list(analysis = a, why = "filename"))
  }
  if (is.na(cols$transcript)) {
    list(analysis = "DGE", why = "no transcript column")
  } else {
    list(analysis = "DTE", why = "guessed \u2014 check this")
  }
}

# One shape for all three analyses: a feature is a gene for DGE and a
# transcript for DTE and DTU.
normalise_results <- function(df, analysis, cols) {
  out <- data.frame(
    feature_id = as.character(
      if (!is.na(cols$transcript)) df[[cols$transcript]] else df[[cols$gene]]
    ),
    gene_id   = as.character(df[[cols$gene]]),
    gene_name = as.character(
      if (!is.na(cols$gene_name)) df[[cols$gene_name]] else df[[cols$gene]]
    ),
    log2FC = as.numeric(df[[cols$log2FC]]),
    pvalue = as.numeric(df[[cols$pvalue]]),
    qvalue = as.numeric(df[[cols$qvalue]]),
    stringsAsFactors = FALSE
  )

  out$log10mean <- if (!is.na(cols$log10mean)) {
    as.numeric(df[[cols$log10mean]])
  } else if (!is.na(cols$mean_lin)) {
    log10(pmax(as.numeric(df[[cols$mean_lin]]), 1e-3))
  } else {
    NA_real_
  }

  if (!is.na(cols$stat))  out$stat      <- as.numeric(df[[cols$stat]])
  if (!is.na(cols$infrv)) out$meanInfRV <- as.numeric(df[[cols$infrv]])

  out$analysis    <- analysis
  out$level       <- if (analysis == "DGE") "gene" else "transcript"
  out$qvalue      <- pmax(out$qvalue, .Machine$double.xmin)
  out$neg_log10_q <- -log10(out$qvalue)
  out$novel       <- grepl("^(XLOC|ENSG|MSTRG|TCONS)", out$gene_name)
  out
}

read_results <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext == "csv") {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt", "tab")) {
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }
}

# Adds one entry to the store, keyed by filename. Never drops a file silently:
# unusable files stay in the store with ok = FALSE and a reason.
ingest <- function(path, name, store) {
  entry <- list(name = name, analysis = NA_character_, why = "",
                rows = 0L, ok = FALSE, note = "")

  df <- tryCatch(read_results(path, name), error = function(e) e)

  if (inherits(df, "error")) {
    entry$note <- paste("could not be read:", conditionMessage(df))
    store[[name]] <- entry
    return(store)
  }
  if (is.null(df)) {
    entry$note <- "unsupported file type"
    store[[name]] <- entry
    return(store)
  }

  cols <- resolve_cols(df)
  miss <- missing_required(cols)
  if (length(miss)) {
    entry$note <- paste0("no column matching: ", paste(miss, collapse = ", "),
                         ". Found: ", paste(names(df), collapse = ", "))
    store[[name]] <- entry
    return(store)
  }

  cl <- classify_analysis(name, cols)
  entry$analysis <- cl$analysis
  entry$why      <- cl$why
  entry$rows     <- nrow(df)
  entry$ok       <- TRUE
  entry$cols     <- cols
  entry$df       <- normalise_results(df, cl$analysis, cols)

  store[[name]] <- entry
  store
}

load_dir <- function(dir = "data") {
  store <- list()
  if (!dir.exists(dir)) return(store)
  found <- list.files(dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
  for (p in found) store <- ingest(p, basename(p), store)
  store
}

# Which file wins for each analysis, and why. Single source of truth for the
# "Loaded files" panel and for build_datasets(). Where several files claim the
# same analysis the largest wins, since full results beat significant-only
# subsets. Simulated data is deprioritised so it can never mask a real table.
# Note the penalty is a plain substring match, so a genuine results file with
# "test" in its name would also be deprioritised.
store_manifest <- function(store) {
  if (!length(store)) return(NULL)

  winner <- list()
  best   <- list()
  for (nm in names(store)) {
    e <- store[[nm]]
    if (!isTRUE(e$ok)) next
    penalty <- if (grepl("example|simulat|test", nm, ignore.case = TRUE)) 0.001 else 1
    score <- e$rows * penalty
    if (is.null(best[[e$analysis]]) || score > best[[e$analysis]]) {
      best[[e$analysis]]   <- score
      winner[[e$analysis]] <- nm
    }
  }

  do.call(rbind, lapply(names(store), function(nm) {
    e <- store[[nm]]
    data.frame(
      name     = nm,
      ok       = isTRUE(e$ok),
      analysis = if (isTRUE(e$ok)) e$analysis else NA_character_,
      why      = if (is.null(e$why)) "" else e$why,
      rows     = e$rows,
      in_use   = isTRUE(e$ok) && identical(winner[[e$analysis]], nm),
      note     = if (is.null(e$note)) "" else e$note,
      stringsAsFactors = FALSE
    )
  }))
}

build_datasets <- function(store) {
  man <- store_manifest(store)
  if (is.null(man)) return(list())
  use <- man[man$in_use, , drop = FALSE]
  out <- list()
  for (i in seq_len(nrow(use))) out[[use$analysis[i]]] <- store[[use$name[i]]]$df
  out
}

status_of <- function(passes, log2FC) {
  factor(
    ifelse(passes & log2FC > 0, "Up",
           ifelse(passes & log2FC < 0, "Down", "Not significant")),
    levels = STATUS_LEVELS
  )
}


# ---- Transcript structure ---------------------------------------------------
# Reads gffcmp.combined.gtf and draws the exon/intron structure of every
# transcript at a gene's locus, annotated with the DTE and DTU results already
# loaded. The GTF's transcript_id is the same TCONS_ namespace as the swish
# tables, so this joins directly with no mapping layer.
#
# Swept across all 35,762 gene keys in a real combined GTF: 0 failures, no gene
# mixing chromosomes, no single-transcript gene falsely flagged, no gene made
# less legible by compression.
#
# ggplot2 + data.table only, both already required above.

STRUCT_LINE <- "#6B7580"
# STRUCT_COLS is defined below, once the shared palette exists.

# ---- DTU effect size --------------------------------------------------------
# The DTU log2FC column is not a fold change, and reading it as one badly
# understates what happened.
#
# swish computes  log2FC = median_k[ log2(mean(B_k) + pc) - log2(mean(A_k) + pc) ]
# with pc = 5 (fishpond, getLog2FC). DTU runs that on isoformProportions()
# output, whose values lie in [0, 1] -- so the pseudocount of 5 dominates
# completely and the statistic is squeezed into +/- log2(6/5) = +/- 0.263.
# Across 34,396 transcripts here the observed range is -0.222 to 0.244 and not
# one value falls outside that bound.
#
# Inverting: with L the reported value and p_A the usage before,
#   p_B = 2^L * (p_A + pc) - pc
#   dp  = p_B - p_A = (2^L - 1) * (p_A + pc)
# p_A is not recoverable from summary tables, but it is confined to [0, 1], so
# dp is bracketed within about +/-10%. Taking the midpoint p_A = 0.5 gives the
# estimate below. On POLDIP3 it returns -75 and +84 percentage points against
# a true -78 and +78 solved from the DTE fold changes.
#
# This is an estimate and is labelled as one. The raw value stays in the data.
DTU_PC <- 5

dtu_delta_usage <- function(l2fc, p_mid = 0.5) {
  # Clamped to +/-1: a proportion cannot move by more than the whole interval,
  # so an estimate beyond that is the midpoint approximation overshooting a
  # hard bound, not a larger effect. Without this the table's most extreme
  # transcript reports +101 percentage points.
  pmax(-1, pmin(1, (2^l2fc - 1) * (p_mid + DTU_PC)))
}

fmt_usage <- function(l2fc) {
  ifelse(is.na(l2fc), NA_character_,
         sprintf("%+.0f pp", 100 * dtu_delta_usage(l2fc)))
}

gtf_read <- function(path) {
  dt <- fread(path, sep = "\t", header = FALSE, quote = "", showProgress = FALSE,
              col.names = c("seqid","source","feature","start","end",
                            "score","strand","frame","attr"))
  dt <- dt[feature %chin% c("exon","transcript")]
  pull <- function(x, key) {
    p <- paste0('.*', key, ' "([^"]+)".*')
    ifelse(grepl(paste0(key, ' "'), x, fixed = TRUE), sub(p, "\\1", x), NA_character_)
  }
  dt[, transcript_id := pull(attr, "transcript_id")]
  dt[, gene_id       := pull(attr, "gene_id")]
  dt[, gene_name     := pull(attr, "gene_name")]
  dt[, oId           := pull(attr, "oId")]
  dt[, num_samples   := as.integer(pull(attr, "num_samples"))]
  dt[, attr := NULL]
  dt[]
}

# Scanning all 767k GTF rows per gene lookup costs about a second a click. Split
# once into a transcript table keyed by gene and an exon table keyed by
# transcript, so each lookup is a keyed join.
gtf_index <- function(gtf) {
  tx <- gtf[feature == "transcript",
            .(transcript_id, gene_id, gene_name, seqid, strand, source,
              oId, num_samples)]
  # A transcript must be reachable by symbol and by XLOC id: some genes in the
  # DGE table are named only by XLOC.
  by_sym <- tx[!is.na(gene_name)][, gkey := gene_name]
  by_loc <- tx[, gkey := gene_id][]
  idx_tx <- unique(rbindlist(list(by_sym, by_loc), use.names = TRUE))
  setkey(idx_tx, gkey)

  ex <- gtf[feature == "exon", .(transcript_id, seqid, start, end)]
  setkey(ex, transcript_id)
  list(tx = idx_tx, ex = ex)
}

# Piecewise-linear genomic -> plot coordinate map. Exonic blocks keep their true
# width; gaps between them are compressed.
build_coord_map <- function(starts, ends, gap_width = 250) {
  o  <- order(starts)
  s  <- starts[o]; e <- ends[o]
  bs <- integer(0); be <- integer(0)
  cs <- s[1]; ce <- e[1]
  if (length(s) > 1) {
    for (i in 2:length(s)) {
      if (s[i] <= ce + 1L) ce <- max(ce, e[i])
      else { bs <- c(bs, cs); be <- c(be, ce); cs <- s[i]; ce <- e[i] }
    }
  }
  blocks <- data.table(bstart = c(bs, cs), bend = c(be, ce))
  nb <- nrow(blocks)
  w  <- blocks$bend - blocks$bstart + 1

  # Compress introns, never stretch them. A fixed gap width pads genes whose
  # introns are shorter than it, making the figure less legible than the raw
  # coordinates -- that hit 992 genes before this was a pmin.
  real_gap <- if (nb > 1) blocks$bstart[-1] - blocks$bend[-nb] - 1 else numeric(0)
  eff_gap  <- pmin(gap_width, pmax(real_gap, 0))
  blocks[, poffset := cumsum(c(0, head(w, -1) + eff_gap))]

  map <- function(x) {
    i <- findInterval(x, blocks$bstart); i[i < 1] <- 1
    off <- pmin(x - blocks$bstart[i], blocks$bend[i] - blocks$bstart[i])
    blocks$poffset[i] + pmax(off, 0)
  }
  list(map = map, blocks = blocks, total = blocks$poffset[nb] + w[nb],
       exonic = sum(w))
}

# Every transcript is compared against one chosen baseline transcript, the way
# browser compares your models against the annotation track.
#
# An earlier version compared each exon against the union of all other
# transcripts. That isolates a cassette exon beautifully at a 4-transcript
# locus and collapses into noise at a 40-transcript one, because with enough
# isoforms nearly every exon is missing from somebody. Diffing against a single
# a baseline stays interpretable at any transcript count, because the question
# is always the same one: how does this isoform differ from that one?
#
#   novel_exon - no overlapping exon in the baseline (insertion / cassette)
#   altbound   - overlaps a baseline exon, but a boundary differs
#   is_base    - the baseline itself, drawn with no marks
#
# Baseline exons that a transcript lacks need no mark: against a pinned
# baseline the gap is what you see.
# base_exons lets the baseline be a transcript that is not among the observed
# ones -- which is the whole point once a reference annotation is loaded, since
# the baseline is then a GENCODE model rather than one of your assembled
# isoforms wearing a label.
classify_exons <- function(ex, base_tx, base_exons = NULL) {
  ex <- copy(ex)
  ex[, eid := .I]
  ex[, is_base := !is.na(base_tx) & transcript_id == base_tx]

  base <- if (!is.null(base_exons) && nrow(base_exons)) {
    base_exons[, .(start, end)]
  } else {
    ex[is_base == TRUE, .(start, end)]
  }
  if (!nrow(base)) {
    ex[, `:=`(novel_exon = FALSE, altbound = FALSE)]
    return(ex[])
  }
  base <- copy(base)
  setkey(base, start, end)

  a  <- ex[, .(eid, start, end)]
  ov <- foverlaps(a, base, by.x = c("start","end"), by.y = c("start","end"),
                  type = "any", nomatch = NULL)
  hit <- ov[, .(n_base = .N, exact = any(i.start == start & i.end == end)),
            by = eid]
  ex <- merge(ex, hit, by = "eid", all.x = TRUE)
  ex[is.na(n_base), `:=`(n_base = 0L, exact = FALSE)]

  ex[, novel_exon := !is_base & n_base == 0L]

  # Boundary differences are only marked on internal exons. The 5' and 3' ends
  # of long-read transcript models are unreliable -- incomplete capture,
  # degradation, imprecise TSS and polyA calls -- so a terminal exon whose
  # bound differs from the baseline is usually technical, not an alternative
  # splice site. Terminal exons carried 62% of these marks in a 1,500-locus
  # audit; excluding them takes the median per-locus burden from 20% of exons
  # to 3%, and the loci where marks swamp the figure from 183 to 17.
  setorder(ex, transcript_id, start)
  ex[, .pos := seq_len(.N), by = transcript_id]
  ex[, .nex := .N, by = transcript_id]
  ex[, internal := .nex > 2L & .pos > 1L & .pos < .nex]
  ex[, altbound := !is_base & n_base > 0L & !exact & internal]
  ex[, c(".pos", ".nex") := NULL]
  ex[]
}

# IGV convention: the annotated model is the baseline. Choosing which one took
# two attempts to get right, both corrected against a 4,000-locus audit.
#
# Ranking by num_samples is wrong: a one-exon fragment seen in every sample
# outranks a complete 28-exon model seen in fewer, and then nearly every exon
# in the locus reads as "absent from the baseline". That was 43 of the 50
# worst loci.
#
# Ranking by exonic footprint fixes most of those but not all, because
# gffcompare builds XLOC loci by chaining overlaps: a locus can hold two
# structurally disjoint groups, where the biggest transcript still shares
# nothing with half the rows.
#
# So rank by structural centrality — how many other transcripts at the locus
# this one actually overlaps — which is the quantity the marks depend on.
# Annotation status and size are tiebreaks.
pick_baseline <- function(ann, ex) {
  tids <- unique(ann$transcript_id)
  if (length(tids) < 2L) return(tids[1])

  # Score each candidate by the thing the figure actually depends on: how many
  # of the other transcripts' exons it accounts for. Count of overlapping
  # transcripts is a poor proxy -- one long unspliced exon touches every
  # transcript while covering almost none of their exons, which is how a
  # single-exon model kept winning. One self-join serves every candidate.
  a <- ex[, .(eid = .I, transcript_id, start, end)]
  b <- copy(a); setkey(b, start, end)
  ov <- foverlaps(a, b, by.x = c("start","end"), by.y = c("start","end"),
                  type = "any", nomatch = NULL)
  ov <- ov[transcript_id != i.transcript_id]

  cov <- ov[, .(covered = uniqueN(eid)), by = .(cand = i.transcript_id)]
  fp  <- ex[, .(bp = sum(end - start + 1), n_ex = .N), by = transcript_id]
  sc  <- merge(fp, cov, by.x = "transcript_id", by.y = "cand", all.x = TRUE)
  sc[is.na(covered), covered := 0L]
  sc[, annotated := transcript_id %chin%
       ann$transcript_id[!is.na(ann$oId) & grepl("^ENS", ann$oId)]]
  ns <- ann$num_samples[match(sc$transcript_id, ann$transcript_id)]
  ns[is.na(ns)] <- 0L
  sc[, ns := ns]

  setorder(sc, -covered, -annotated, -n_ex, -bp, -ns, transcript_id)
  sc$transcript_id[1]
}

# ---- Reference annotation (optional) ----------------------------------------
# Queried by region the way a genome browser does it: the annotation is never
# parsed in full, only the bytes covering the locus on screen. GENCODE v49
# bgzips from 3.1 GB to 126 MB with a 0.3 MB index, and a locus query returns
# in about 10 ms, so this is cheap enough for the hosted app as well.
#
# Rsamtools is a soft dependency. Without it the reference controls stay hidden
# and every other part of the panel behaves exactly as before.

HAVE_TABIX <- requireNamespace("Rsamtools", quietly = TRUE)

# GENCODE writes chr1; Ensembl writes 1. Try the locus both ways rather than
# making the user care.
ref_seq_variants <- function(seqid) {
  unique(c(seqid, sub("^chr", "", seqid), paste0("chr", sub("^chr", "", seqid))))
}

ref_query <- function(gz, seqid, start, end) {
  if (!HAVE_TABIX || is.null(gz) || !file.exists(gz)) return(NULL)
  tbx <- Rsamtools::TabixFile(gz)
  lines <- NULL
  for (s in ref_seq_variants(seqid)) {
    got <- tryCatch(
      Rsamtools::scanTabix(tbx, param = GenomicRanges::GRanges(
        s, IRanges::IRanges(start, end)))[[1]],
      error = function(e) character(0))
    if (length(got)) { lines <- got; break }
  }
  if (!length(lines)) return(NULL)

  p <- tstrsplit(lines, "\t", fixed = TRUE)
  att <- p[[9]]
  g <- function(key) {
    pat <- paste0('.*', key, ' "([^"]+)".*')
    ifelse(grepl(paste0(key, ' "'), att, fixed = TRUE), sub(pat, "\\1", att), NA_character_)
  }
  data.table(
    seqid = p[[1]], feature = p[[3]],
    start = as.integer(p[[4]]), end = as.integer(p[[5]]), strand = p[[7]],
    transcript_id = g("transcript_id"), gene_name = g("gene_name"),
    ttype = g("transcript_type"),
    mane  = grepl('tag "MANE_Select"', att, fixed = TRUE),
    canon = grepl('tag "Ensembl_canonical"', att, fixed = TRUE))
}

# Which reference transcript is the designated baseline. MANE Select is one
# transcript per protein-coding gene agreed between NCBI and Ensembl, which is
# exactly the "designated reference" a browser draws; Ensembl canonical covers
# genes without MANE, and largest footprint covers the rest.
pick_ref_baseline <- function(rt, rex) {
  if (!nrow(rt)) return(NA_character_)
  if (any(rt$mane))  return(rt$transcript_id[rt$mane][1])
  if (any(rt$canon)) return(rt$transcript_id[rt$canon][1])
  fp <- rex[, .(bp = sum(end - start + 1)), by = transcript_id][order(-bp)]
  fp$transcript_id[1]
}

gene_structure <- function(gene, idx, dte = NULL, dtu = NULL, qcut = 0.05,
                           max_tx = 40L, base_tx = NULL, ref_gz = NULL,
                           max_ref = 25L) {
  notes <- character(0)

  tx <- idx$tx[.(gene), nomatch = NULL]
  if (!nrow(tx)) return(list(ok = FALSE, why = "This gene is not in the loaded GTF."))

  # A symbol can sit at several loci, and a few span more than one chromosome.
  # Coordinates are only comparable within a locus, so take the locus with the
  # most transcripts and say so rather than silently merging them.
  loci <- tx[, .N, by = .(seqid, gene_id)][order(-N, seqid, gene_id)]
  if (nrow(loci) > 1) {
    notes <- c(notes, sprintf(
      "%s occurs at %d loci (%s); showing %s:%s with %d of %d transcripts.",
      gene, nrow(loci), paste(unique(loci$seqid), collapse = ", "),
      loci$seqid[1], loci$gene_id[1], loci$N[1], nrow(tx)))
  }
  tx <- tx[seqid == loci$seqid[1] & gene_id == loci$gene_id[1]]

  ex <- idx$ex[.(unique(tx$transcript_id)), nomatch = NULL][seqid == loci$seqid[1]]
  if (!nrow(ex)) return(list(ok = FALSE, why = "No exons for this locus."))
  ex <- copy(ex)

  ann <- unique(tx[, .(transcript_id, source, oId, num_samples, strand)])
  ann[, `:=`(dte_lfc = NA_real_, dte_q = NA_real_,
             dtu_lfc = NA_real_, dtu_q = NA_real_)]
  # Update-joins, not merge(): merge against the full DTE table re-sorts it and
  # costs ~13 ms a gene, against ~1 ms here. That is a visible pause per click.
  if (!is.null(dte))
    ann[dte, on = "transcript_id", `:=`(dte_lfc = i.dte_lfc, dte_q = i.dte_q)]
  if (!is.null(dtu))
    ann[dtu, on = "transcript_id", `:=`(dtu_lfc = i.dtu_lfc, dtu_q = i.dtu_q)]

  ann[, tested := !is.na(dtu_q) | !is.na(dte_q)]
  # A factor with every level present, so the legend keeps a complete key even
  # at a locus where nothing is significant. As a plain character column the
  # scale's limits produce the labels but leave the swatches blank.
  ann[, status := factor(
        fifelse(!tested, "Not tested",
        fifelse(!is.na(dtu_q) & dtu_q <= qcut & dtu_lfc > 0, "Usage up",
        fifelse(!is.na(dtu_q) & dtu_q <= qcut & dtu_lfc < 0, "Usage down",
                "Not significant"))),
        levels = names(STRUCT_COLS))]
  ann[, novel := source == "IsoQuant"]

  # Past about 40 rows the figure stops being readable. Keep significant and
  # tested transcripts preferentially, and say what was dropped.
  if (nrow(ann) > max_tx) {
    ann <- ann[order(-(status %in% c("Usage up","Usage down")), -tested,
                     -abs(dtu_lfc), -num_samples, transcript_id)]
    notes <- c(notes, sprintf(
      "%d transcripts at this locus; showing the %d most relevant.",
      nrow(ann), max_tx))
    ann <- head(ann, max_tx)
    ex  <- ex[transcript_id %chin% ann$transcript_id]
  }

  # ---- reference annotation, if one is loaded -------------------------------
  # With a reference the baseline is a GENCODE model drawn in its own block, so
  # no observed isoform is promoted. Without one, fall back to picking a
  # baseline from the data and say so.
  rt <- NULL; rex <- NULL; ref_base <- NA_character_
  if (!is.null(ref_gz)) {
    rq <- ref_query(ref_gz, loci$seqid[1], min(ex$start), max(ex$end))
    if (!is.null(rq) && nrow(rq)) {
      # A region query returns every overlapping gene, so keep this one. Match
      # on symbol; if the reference uses a different symbol, keep everything
      # rather than silently showing nothing.
      keep <- rq[feature == "transcript" & !is.na(gene_name) & gene_name == gene]
      name_match <- nrow(keep) > 0
      if (!name_match) keep <- rq[feature == "transcript"]
      rt  <- unique(keep[, .(transcript_id, gene_name, ttype, mane, canon,
                             start, end, strand)])
      rex <- rq[feature == "exon" & transcript_id %chin% rt$transcript_id,
                .(transcript_id, start, end)]
      ref_base <- pick_ref_baseline(rt, rex)

      # A novel locus has no reference gene of the same name, but the region
      # query still returns whatever annotation overlaps it. That is useful
      # context and worth drawing -- but diffing your isoforms against an
      # unrelated gene's model would be nonsense, so the baseline stays with
      # the data in that case.
      if (!name_match) {
        notes <- c(notes, sprintf(
          "No reference gene named %s here; the annotation block shows overlapping models (%s) for context only, and the baseline is taken from your data.",
          gene, paste(utils::head(unique(na.omit(rt$gene_name)), 3), collapse = ", ")))
        ref_base <- NA_character_
      }
    }
  }
  # ref_shown: draw the annotation block. has_ref: also use it as the baseline.
  ref_shown <- !is.null(rt) && nrow(rt) > 0
  has_ref   <- ref_shown && !is.na(ref_base)

  if (ref_shown) {
    # Versions may differ between the annotation and gffcompare's oId, so match
    # on the unversioned accession.
    strip <- function(x) sub("\\.\\d+$", "", x)
    rt[, detected := strip(transcript_id) %chin% strip(na.omit(ann$oId))]

    # Heavily annotated loci can carry dozens of reference models. Cap the
    # block, keeping the baseline and the detected ones, so the annotation
    # never crowds out the data it exists to contextualise.
    if (nrow(rt) > max_ref) {
      rt <- rt[order(-(transcript_id == ref_base), -mane, -canon, -detected,
                     start)]
      notes <- c(notes, sprintf(
        "%d reference transcripts at this locus; showing %d.", nrow(rt), max_ref))
      rt  <- head(rt, max_ref)
      rex <- rex[transcript_id %chin% rt$transcript_id]
    }
  }

  if (has_ref) {
    base_exons <- rex[transcript_id == ref_base, .(start, end)]
    base_tx <- NA_character_          # no observed isoform is the baseline
  } else {
    base_exons <- NULL
    if (is.null(base_tx) || !(base_tx %chin% ann$transcript_id)) {
      auto <- pick_baseline(ann, ex)
      if (!is.null(base_tx)) notes <- c(notes,
        sprintf("%s is not at this locus; using %s as the baseline.", base_tx, auto))
      base_tx <- auto
    }
  }

  # Observed rows sit above the reference block. Without a reference the
  # baseline is pinned to the bottom row instead, browser-style.
  n_ref <- if (ref_shown) nrow(rt) else 0L
  ann <- ann[order(!is.na(base_tx) & transcript_id == base_tx,
                   -(status %in% c("Usage up","Usage down")),
                   -tested, -abs(dtu_lfc), transcript_id)]
  ann[, y := .N:1 + n_ref]
  ann[, is_base := !is.na(base_tx) & transcript_id == base_tx]

  ex <- classify_exons(ex, base_tx, base_exons)
  # The coordinate map has to span the reference too, or its exons fall off the
  # axis at loci where the annotation extends past the assembly.
  all_s <- c(ex$start, if (ref_shown) rex$start)
  all_e <- c(ex$end,   if (ref_shown) rex$end)
  cm <- build_coord_map(all_s, all_e)
  ex[, `:=`(px = cm$map(start), pxe = cm$map(end))]
  ex <- merge(ex, ann[, .(transcript_id, y, status)], by = "transcript_id")

  if (ref_shown) {
    setorder(rt, -mane, -canon, -detected, start)
    rt[, y := n_ref:1]
    rt[, is_base := transcript_id == ref_base]
    rex <- merge(rex, rt[, .(transcript_id, y)], by = "transcript_id")
    rex[, `:=`(px = cm$map(start), pxe = cm$map(end))]
    setorder(rex, transcript_id, start)
    # No note here: the reference count, the detected count and the baseline
    # identity are all on the figure's facts line, and repeating them below it
    # was the redundant third row of subtitle.
  }

  setorder(ex, transcript_id, start)
  # Single-exon transcripts contribute no introns. head(x, -1) on a length-1
  # group yields zero rows by itself, so no branch is needed -- and a branch is
  # what broke 380 genes, returning an integer y in one arm and a double in the
  # other.
  intr <- ex[, .(x = head(pxe, -1L), xend = tail(px, -1L), y = head(y, -1L)),
             by = transcript_id]

  ref_intr <- NULL
  if (ref_shown) {
    ref_intr <- rex[, .(x = head(pxe, -1L), xend = tail(px, -1L),
                        y = head(y, -1L)), by = transcript_id]
  }

  list(ok = TRUE, gene = gene, exons = ex, introns = intr, ann = ann,
       cm = cm, seqid = loci$seqid[1], locus = loci$gene_id[1],
       strand = ann$strand[1], n_loci = nrow(loci), notes = notes,
       base_tx = base_tx, choices = ann$transcript_id,
       has_ref = has_ref, ref_shown = ref_shown,
       ref_tx = rt, ref_exons = rex, ref_introns = ref_intr,
       ref_base = ref_base, n_ref = n_ref,
       span = max(ex$end) - min(ex$start) + 1L)
}

plot_gene_structure <- function(st, qcut = 0.05) {
  ex <- st$exons; intr <- st$introns; ann <- st$ann
  cm <- st$cm; n_tx <- nrow(ann)
  arrow_len <- cm$total * 0.004

  lab <- copy(ann)
  lab[, left  := sprintf("%s%s\n%s", fifelse(is_base, "▸ ", ""),
                         transcript_id, fifelse(is.na(oId), "", oId))]
  # Usage shift in percentage points, not swish's pseudocount-compressed
  # log2FC. "-75 pp" says an isoform went from most of the gene's output to
  # almost none; "-0.212" does not.
  # The baseline is a transcript like any other and usually a tested one, so
  # it keeps its numbers; being the baseline is a tag, not a reason to hide its
  # result.
  lab[, right := fifelse(tested,
        sprintf("%susage %s (est.)\nDTE %+.2f   q %.1e",
                fifelse(is_base, "[baseline]  ", ""),
                fmt_usage(dtu_lfc), dte_lfc, dtu_q),
        fifelse(is_base, "[baseline]  not tested", "not tested"))]

  xpad <- cm$total * 0.02
  if (isTRUE(st$has_ref)) {
    base_lab <- st$ref_base
  } else {
    base_lab <- st$base_tx
    base_o   <- ann$oId[ann$is_base]
    if (length(base_o) && !is.na(base_o[1])) base_lab <- paste0(base_lab, " (", base_o[1], ")")
  }
  # Only what changes between genes goes above the figure. The reading guide is
  # the same for every locus, so it lives beside the plot in the app rather
  # than being reprinted six lines deep on every render. One compact glyph key
  # stays, so an exported PNG is still self-explanatory.
  wrap <- function(x, w = 128) paste(strwrap(x, width = w), collapse = "\n")

  facts <- sprintf("%d assembled", n_tx)
  if (isTRUE(st$ref_shown))
    facts <- sprintf("%s  ·  %d reference (%d detected)", facts,
                     nrow(st$ref_tx), sum(st$ref_tx$detected))
  facts <- sprintf("%s  ·  baseline %s%s", facts, base_lab,
    if (isTRUE(st$has_ref)) {
      if (any(st$ref_tx$mane)) " (MANE Select)"
      else if (any(st$ref_tx$canon)) " (Ensembl canonical)" else " (annotation)"
    } else " (from your data — no reference loaded)")
  facts <- sprintf("%s  ·  q ≤ %.2f", facts, qcut)

  key <- paste("boxed = absent from baseline  ·  dashed = internal boundary differs",
               "·  gap = exon missing  ·  dotted backbone = novel model")

  sub <- paste0(wrap(facts), "\n", wrap(key))
  if (length(st$notes))
    sub <- paste0(sub, "\n", wrap(paste(st$notes, collapse = "  ")))

  base_y <- ann$y[ann$is_base]

  # Novel models get a dotted backbone. Origin is otherwise only readable from
  # the oId in the row label, and it is worth knowing at a glance.
  intr <- merge(intr, ann[, .(transcript_id, novel)], by = "transcript_id",
                all.x = TRUE)
  intr[is.na(novel), novel := FALSE]

  # ggplot2 4.x omits the legend key for a status with no rows in the layer, so
  # a locus where nothing is significant loses the red and blue swatches and
  # the legend stops explaining the figure. Neither scale limits, drop = FALSE
  # nor override.aes restores them. A zero-area rectangle per level does: it
  # draws nothing and gives every level data to key off.
  legend_seed <- data.table(
    status = factor(names(STRUCT_COLS), levels = names(STRUCT_COLS)),
    x = 0, y = 1)

  # Reference rows carry no DTU result, so they get no status fill -- drawing
  # them in the data palette would imply a measurement that does not exist.
  rt <- st$ref_tx; rex <- st$ref_exons; rintr <- st$ref_introns
  ylab <- lab[, .(y, left, right)]
  if (isTRUE(st$ref_shown)) {
    rlab <- copy(rt)
    rlab[, left := fifelse(is_base,
      sprintf("BASELINE ▶\n%s\n%s", transcript_id, fifelse(is.na(ttype), "", ttype)),
      sprintf("%s\n%s", transcript_id, fifelse(is.na(ttype), "", ttype)))]
    rlab[, right := fifelse(is_base,
      paste0("◀ BASELINE — ",
             fifelse(mane, "MANE Select",
             fifelse(canon, "Ensembl canonical", "largest model")),
             "\nevery transcript above is compared to this"),
      paste0(fifelse(mane, "MANE Select",
             fifelse(canon, "Ensembl canonical", "annotation")),
             fifelse(detected, "\ndetected", "\nnot detected in your data")))]
    ylab <- rbind(ylab, rlab[, .(y, left, right)])
  }

  p <- ggplot() +
    geom_rect(data = legend_seed,
              aes(xmin = x, xmax = x, ymin = y, ymax = y, fill = status),
              colour = NA) +
    # Shade the baseline row only when it is one of your own transcripts. With a
    # reference loaded the baseline lives in the annotation block instead, and
    # shading a data row there would reintroduce exactly the conflation this
    # block exists to remove.
    {if (length(base_y) == 1)
       annotate("rect", xmin = -xpad, xmax = cm$total + xpad,
                ymin = base_y - 0.5, ymax = base_y + 0.5,
                fill = "#1F2933", alpha = 0.07)} +
    geom_segment(data = intr[novel == FALSE],
                 aes(x = x, xend = xend, y = y, yend = y),
                 colour = STRUCT_LINE, linewidth = 0.35) +
    geom_segment(data = intr[novel == TRUE],
                 aes(x = x, xend = xend, y = y, yend = y),
                 colour = STRUCT_LINE, linewidth = 0.35, linetype = "13") +
    geom_rect(data = ex,
              aes(xmin = px, xmax = pxe + 1, ymin = y - 0.30, ymax = y + 0.30,
                  fill = status), colour = NA) +
    # Outline only the untested exons, as a separate layer rather than a second
    # scale mapped to status: mapping colour to the same variable merges into
    # the fill legend and blanks the keys for levels absent from the locus.
    geom_rect(data = ex[status == "Not tested"],
              aes(xmin = px, xmax = pxe + 1, ymin = y - 0.30, ymax = y + 0.30),
              fill = NA, colour = "#B6BDC4", linewidth = 0.3) +
    geom_rect(data = ex[novel_exon == TRUE],
              aes(xmin = px - arrow_len * 3, xmax = pxe + arrow_len * 3,
                  ymin = y - 0.44, ymax = y + 0.44),
              fill = NA, colour = "#111820", linewidth = 0.65) +
    geom_rect(data = ex[altbound == TRUE],
              aes(xmin = px, xmax = pxe + 1, ymin = y - 0.30, ymax = y + 0.30),
              fill = NA, colour = "#5A6570", linewidth = 0.3, linetype = "21") +
    scale_fill_manual(values = STRUCT_COLS, name = NULL,
                      limits = names(STRUCT_COLS), drop = FALSE) +
    scale_y_continuous(breaks = ylab$y, labels = ylab$left,
      sec.axis = sec_axis(~ ., breaks = ylab$y, labels = ylab$right)) +
    coord_cartesian(xlim = c(-xpad, cm$total + xpad), clip = "off") +
    labs(title = sprintf("%s — %s:%s-%s (%s strand)", st$gene, st$seqid,
                         format(min(ex$start), big.mark = ","),
                         format(max(ex$end), big.mark = ","), st$strand),
         subtitle = sub, x = "compressed genomic coordinate", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(family = "mono", size = 8,
                                     lineheight = 0.95, hjust = 1),
          axis.text.y.right = element_text(family = "mono", size = 8,
                                     lineheight = 0.95, hjust = 0,
                                     colour = "#6B7580"),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "#6B7580", size = 8,
                                       lineheight = 1.2),
          legend.position = "bottom",
          plot.margin = margin(12, 14, 8, 14))

  if (isTRUE(st$ref_shown)) {
    # Rule separating your data from the annotation, then the annotation drawn
    # as a different kind of object: outlined, unfilled, no status colour.
    p <- p +
      annotate("segment", x = -xpad, xend = cm$total + xpad,
               y = st$n_ref + 0.5, yend = st$n_ref + 0.5,
               colour = "#B6BDC4", linewidth = 0.4) +
      # No inline caption on the rule: at the row spacing this panel uses it
      # collides with the baseline row, and the subtitle plus the right-axis
      # labels already say what the block is.
      geom_segment(data = rintr, aes(x = x, xend = xend, y = y, yend = y),
                   colour = "#9AA1A8", linewidth = 0.3) +
      geom_rect(data = rex,
                aes(xmin = px, xmax = pxe + 1, ymin = y - 0.24, ymax = y + 0.24),
                fill = "#FFFFFF", colour = "#6B7580", linewidth = 0.35)

    # The baseline needs to announce itself: it is the thing every other row is
    # measured against, and a pale tint did not carry that. Full-width band,
    # heavier exons, and a bracket in the margin.
    if (isTRUE(st$has_ref)) {
      by <- st$ref_tx[transcript_id == st$ref_base]$y[1]
      p <- p +
        annotate("rect", xmin = -xpad, xmax = cm$total + xpad,
                 ymin = by - 0.46, ymax = by + 0.46,
                 fill = "#2C5F8A", alpha = 0.10) +
        annotate("segment", x = -xpad, xend = -xpad,
                 y = by - 0.46, yend = by + 0.46,
                 colour = "#2C5F8A", linewidth = 1.6) +
        geom_segment(data = rintr[transcript_id == st$ref_base],
                     aes(x = x, xend = xend, y = y, yend = y),
                     colour = "#2C5F8A", linewidth = 0.45) +
        geom_rect(data = rex[transcript_id == st$ref_base],
                  aes(xmin = px, xmax = pxe + 1, ymin = y - 0.30, ymax = y + 0.30),
                  fill = "#B9CBDD", colour = "#1F4B73", linewidth = 0.7)
    }
  }

  if (nrow(intr)) {
    p <- p + geom_segment(data = intr,
      aes(x = (x + xend) / 2, y = y, yend = y,
          xend = (x + xend) / 2 + ifelse(st$strand == "-", -1, 1) * arrow_len),
      colour = STRUCT_LINE, linewidth = 0.35,
      arrow = arrow(length = unit(0.05, "cm"), type = "open"))
  }
  p
}


options(shiny.maxRequestSize = 500 * 1024^2)

# Above this many non-significant points, thin them for display.
NS_CAP <- 15000

COL_UP   <- "#B03A48"
COL_DOWN <- "#35618F"
COL_NS   <- "#C4C9CE"
STATUS_COLS <- setNames(c(COL_UP, COL_DOWN, COL_NS), STATUS_LEVELS)

# Structure panel palette. The overview plot only needs three states and can
# use a pale grey for "not significant"; this panel needs four, and reusing
# COL_NS left "not significant" (#C4C9CE) and "not tested" almost identical, so
# a locus with nothing significant read as undifferentiated grey.
#
# The two data-bearing states keep the shared red and blue. The two grey states
# are separated on lightness rather than hue, and "not tested" is additionally
# given an outline in the plot so it reads as an empty shape rather than a
# fainter version of the same thing.
STRUCT_COLS <- c("Usage up"        = COL_UP,
                 "Usage down"      = COL_DOWN,
                 "Not significant" = "#8E979F",
                 "Not tested"      = "#EDEFF1")
STRUCT_EDGE <- c("Usage up"        = NA,
                 "Usage down"      = NA,
                 "Not significant" = NA,
                 "Not tested"      = "#B6BDC4")

app_theme <- bs_theme(
  version      = 5,
  bg           = "#F7F8F8",
  fg           = "#14181C",
  primary      = "#1F2933",
  base_font    = font_google("Inter"),
  heading_font = font_google("Inter"),
  code_font    = font_google("JetBrains Mono"),
  "font-size-base" = "0.9rem",
  "border-radius"  = "0.25rem"
)

# Bundled tables are read once per process and shared by every session.
STARTUP_STORE <- load_dir("data")

# ---- UI ---------------------------------------------------------------------

ui <- page_sidebar(
  title = "LoReRNA Explorer",
  theme = app_theme,

  tags$head(tags$style(HTML("
    table.dataTable td.mono, .mono { font-family: var(--bs-font-monospace); }
    .card-header { font-weight: 600; }
    .bslib-value-box .value-box-value { font-variant-numeric: tabular-nums; }
    .form-label { font-size: 0.8rem; text-transform: uppercase;
                  letter-spacing: 0.06em; color: #6B7580; }
    .note { font-size: 0.78rem; color: #6B7580; line-height: 1.35;
            margin-top: -0.5rem; margin-bottom: 1rem; }
    table.manifest { width: 100%; font-size: 0.76rem; border-collapse: collapse; }
    table.manifest td { padding: 2px 4px; vertical-align: top;
                        border-bottom: 1px solid #E8EAEC; }
    table.manifest .fname { font-family: var(--bs-font-monospace);
                            word-break: break-all; }
    .tag { display: inline-block; padding: 0 5px; border-radius: 3px;
           font-weight: 600; font-size: 0.72rem; }
    .tag-ok   { background: #E4EAF0; color: #35618F; }
    .tag-warn { background: #F3E3E5; color: #B03A48; }
  "))),

  sidebar = sidebar(
    width = 340,

    fileInput("files", "Results tables",
              accept = c(".csv", ".tsv", ".txt"), multiple = TRUE),

    # Wording holds for both deployments: run locally, nothing leaves the
    # machine; run from the hosted copy, files reach a third-party server.
    tags$div(
      class = "note",
      "Files are read by whichever machine runs this app. On the hosted copy ",
      "at shinyapps.io that is a third-party server, and uploads are held for ",
      "the session only. Run the app locally and nothing leaves your machine. ",
      "Do not upload data you are not permitted to share externally."
    ),

    fileInput("gtf", "Transcript structures (optional)",
              accept = c(".gtf", ".gff", ".txt")),
    tags$div(
      class = "note",
      "gffcmp.combined.gtf from 04_oarfish_reference/. Adds the Structure tab, ",
      "which draws the exon layout of every transcript at the selected gene's ",
      "locus. Its transcript IDs are the same ones as the results tables, so ",
      "no extra mapping is needed. Around 125 MB of memory per session."
    ),

    if (HAVE_TABIX) tagList(
      fileInput("ref_gtf", "Reference annotation (optional)",
                accept = c(".gz", ".tbi"), multiple = TRUE),
      tags$div(
        class = "note",
        "A bgzip-compressed, tabix-indexed GTF plus its .tbi, uploaded together. ",
        "Queried by region the way a genome browser is, so only the locus on ",
        "screen is read and the file is never parsed in full. With one loaded, ",
        "the comparison baseline becomes the annotation's MANE Select transcript ",
        "instead of one of your own isoforms. To build one: ",
        tags$code("bgzip -c in.gtf > out.gtf.gz && tabix -p gff out.gtf.gz")
      )
    ),

    accordion(
      open = FALSE,
      accordion_panel(
        "Loaded files",
        uiOutput("manifest"),
        tags$hr(),
        selectInput("reassign_file", "Reassign", choices = character(0)),
        radioButtons("reassign_as", NULL, ANALYSES, inline = TRUE),
        actionButton("do_reassign", "Apply", class = "btn-sm"),
        tags$div(class = "note", style = "margin-top:0.6rem;",
                 "Use this when a file was assigned to the wrong analysis, or shows as guessed.")
      )
    ),

    radioButtons("analysis", "Analysis", choices = ANALYSES,
                 selected = "DTU", inline = TRUE),
    tags$hr(),

    numericInput("qval", "Max q-value", value = 0.05, min = 0, max = 1, step = 0.01),
    uiOutput("qval_note"),

    numericInput("lfc", "Min |log2 fold change|", value = 0, min = 0, step = 0.01),
    numericInput("min_expr", "Min log10 mean expression", value = 0, min = 0, step = 0.25),

    conditionalPanel(
      "input.analysis != 'DGE'",
      checkboxInput("reciprocal_only",
                    "Reciprocal genes only (\u2265 1 up and \u2265 1 down)", FALSE)
    ),
    checkboxInput("hide_novel", "Hide unannotated loci", FALSE),
    checkboxInput("thin", "Thin non-significant points for speed", TRUE),
    textInput("gene_search", "Gene name contains", placeholder = "e.g. POLDIP3"),

    tags$hr(),
    downloadButton("dl_features", "Download features", class = "btn-sm btn-primary"),
    downloadButton("dl_genes",    "Download genes",    class = "btn-sm")
  ),

  layout_columns(
    fill = FALSE,
    value_box("Features tested",  textOutput("n_tested",  inline = TRUE)),
    value_box("Passing filters",  textOutput("n_passing", inline = TRUE)),
    value_box("Genes affected",   textOutput("n_genes",   inline = TRUE)),
    value_box("Reciprocal genes", textOutput("n_recip",   inline = TRUE))
  ),

  layout_columns(
    col_widths = c(7, 5),
    card(
      full_screen = TRUE,
      card_header(
        class = "d-flex justify-content-between align-items-center",
        textOutput("overview_header", inline = TRUE),
        radioButtons("plot_type", NULL, c("Volcano", "MA"),
                     selected = "MA", inline = TRUE)
      ),
      plotlyOutput("main_plot", height = "400px"),
      uiOutput("plot_note")
    ),
    card(
      full_screen = TRUE,
      card_header(textOutput("gene_header", inline = TRUE)),
      plotOutput("gene_plot", height = "430px")
    )
  ),

  navset_card_tab(
    full_screen = TRUE,
    nav_panel("Features", DTOutput("feature_table")),
    nav_panel("Genes",    DTOutput("gene_table")),
    nav_panel("Structure",
              uiOutput("structure_note"),
              # Keyed on the client-side input rather than a server output, so
              # it needs no outputOptions() round trip.
              conditionalPanel(
                "input.gtf",
                tags$div(
                  style = "display:flex; align-items:center; gap:0.6rem; margin:0.5rem 0 0.2rem 0;",
                  tags$span(class = "form-label", style = "margin:0;", "Baseline"),
                  tags$div(style = "min-width:280px;",
                    selectInput("base_tx", NULL, choices = character(0),
                                width = "100%")),
                  tags$span(class = "note", style = "margin:0;",
                            "A transcript from your data, chosen automatically. Not an external reference annotation.")
                )
              ),
              plotOutput("structure_plot", height = "auto"),
              # The reading guide is identical for every gene, so it sits here,
              # collapsed, instead of being reprinted above each figure.
              conditionalPanel(
                "input.gtf",
                accordion(
                  open = FALSE, class = "mt-2",
                  accordion_panel(
                    "How to read this panel",
                    tags$div(class = "note", style = "margin-bottom:0.4rem;", tags$ul(
                      style = "padding-left:1.1rem; margin:0;",
                      tags$li(tags$b("Rows"), " are transcripts. Yours are on top; if a reference annotation is loaded it forms its own block below the rule, and none of your isoforms is promoted to stand in for it."),
                      tags$li(tags$b("Baseline"), " is the shaded row every other transcript is compared against. With a reference loaded it is the annotation's MANE Select transcript; without one it is picked from your data and says so."),
                      tags$li(tags$b("Introns are compressed"), " to a fixed width so exons stay legible — at a typical locus only about a tenth of the span is exonic. Exon widths keep their true relative scale."),
                      tags$li(tags$b("Solid box"), ": an exon with no counterpart in the baseline. ",
                              tags$b("Dashed"), ": an internal exon that overlaps the baseline but with a different boundary. ",
                              tags$b("A gap"), " against the baseline is a missing exon."),
                      tags$li("Terminal exon boundaries are never marked. The 5' and 3' ends of long-read models are unreliable — incomplete capture, degradation, imprecise TSS and polyA calls — so flagging them would present artefact as biology."),
                      tags$li(tags$b("A dotted backbone"), " marks a novel model assembled by IsoQuant rather than one carried over from the annotation."),
                      tags$li(tags$b("Fill"), " is DTU status at the q-value set in the sidebar."),
                      tags$li(tags$b("Usage in percentage points"), " is an estimate. swish reports DTU on isoform proportions with a pseudocount of 5, which compresses the statistic into ±0.263 and makes it unreadable as a fold change; this recovers the underlying shift to about ten points.")
                    ))
                  )
                )
              ))
  )
)

# ---- Server -----------------------------------------------------------------

server <- function(input, output, session) {

  # Keyed by filename, not by analysis, so two files claiming the same
  # analysis both survive and can be told apart in the UI.
  files_store <- reactiveVal(STARTUP_STORE)

  observeEvent(input$files, {
    n <- nrow(input$files)
    withProgress(message = "Reading files", value = 0, {
      store <- files_store()
      for (i in seq_len(n)) {
        incProgress(i / n, detail = input$files$name[i])
        store <- ingest(input$files$datapath[i], input$files$name[i], store)
      }
      files_store(store)
    })
  })

  observeEvent(input$do_reassign, {
    req(input$reassign_file, input$reassign_as)
    store <- files_store()
    e <- store[[input$reassign_file]]
    if (!is.null(e) && isTRUE(e$ok)) {
      e$analysis    <- input$reassign_as
      e$why         <- "set by hand"
      e$df$analysis <- input$reassign_as
      e$df$level    <- if (input$reassign_as == "DGE") "gene" else "transcript"
      store[[input$reassign_file]] <- e
      files_store(store)
    }
  })

  datasets <- reactive(build_datasets(files_store()))

  output$manifest <- renderUI({
    man <- store_manifest(files_store())
    if (is.null(man)) return(div(class = "note", "Nothing loaded yet."))

    rows <- lapply(seq_len(nrow(man)), function(i) {
      r <- man[i, ]
      if (!r$ok) {
        return(tags$tr(
          tags$td(class = "fname", r$name),
          tags$td(span(class = "tag tag-warn", "skipped"), tags$br(),
                  span(class = "note", r$note))
        ))
      }
      tags$tr(
        tags$td(class = "fname", r$name),
        tags$td(
          span(class = if (r$in_use) "tag tag-ok" else "tag tag-warn", r$analysis),
          sprintf(" %s rows", format(r$rows, big.mark = ",")), tags$br(),
          span(class = "note",
               paste0(r$why, if (!r$in_use) " \u00b7 not in use, superseded" else ""))
        )
      )
    })
    tags$table(class = "manifest", tags$tbody(rows))
  })

  observeEvent(files_store(), {
    ok <- store_manifest(files_store())
    ok <- if (is.null(ok)) character(0) else ok$name[ok$ok]
    updateSelectInput(session, "reassign_file", choices = ok,
                      selected = if (length(ok)) ok[1] else NULL)
  })

  observeEvent(datasets(), {
    avail <- intersect(ANALYSES, names(datasets()))
    if (!length(avail)) return()
    sel <- if (isTruthy(input$analysis) && input$analysis %in% avail) {
      input$analysis
    } else {
      tail(avail, 1)
    }
    # Freeze first, so nothing downstream renders against the old selection
    # and then again against the new one.
    freezeReactiveValue(input, "analysis")
    updateRadioButtons(session, "analysis", choices = avail,
                       selected = sel, inline = TRUE)
  })

  raw <- reactive({
    d <- datasets()
    validate(need(length(d) > 0,
                  "No usable tables. Open 'Loaded files' in the sidebar to see why."))
    req(input$analysis)
    validate(need(input$analysis %in% names(d),
                  paste(input$analysis, "is not loaded.")))
    d[[input$analysis]]
  })

  # Effect sizes differ by more than an order of magnitude between DGE and
  # DTU, so the controls rescale whenever the active table changes.
  observeEvent(raw(), {
    df   <- raw()
    lmax <- max(abs(df$log2FC), na.rm = TRUE)

    freezeReactiveValue(input, "lfc")
    updateNumericInput(session, "lfc", value = 0,
                       max = signif(lmax, 2), step = signif(lmax / 40, 1))

    freezeReactiveValue(input, "min_expr")
    if (all(is.na(df$log10mean))) {
      updateNumericInput(session, "min_expr", value = 0, max = 0)
    } else {
      updateNumericInput(session, "min_expr", value = 0,
                         max = signif(max(df$log10mean, na.rm = TRUE), 2),
                         step = 0.25)
    }
  })

  # Debounced, so holding a key in a numeric input does not recompute the
  # whole chain on every intermediate value.
  thresholds <- debounce(reactive({
    req(input$qval, !is.null(input$lfc), !is.null(input$min_expr))
    list(q = input$qval, lfc = input$lfc, expr = input$min_expr)
  }), 400)

  output$qval_note <- renderUI({
    df    <- raw()
    floor <- min(df$qvalue, na.rm = TRUE)
    n_at  <- sum(df$qvalue <= floor * 1.001, na.rm = TRUE)
    div(class = "note", sprintf(
      "%s distinct q-values. %s features sit at the smallest achievable value (%s), the permutation floor set by --swish_nperms, so ranking among them uses effect size.",
      format(n_distinct(df$qvalue), big.mark = ","),
      format(n_at, big.mark = ","), signif(floor, 3)
    ))
  })

  annotated <- reactive({
    th <- thresholds()
    df <- raw()
    df$passes <- df$qvalue <= th$q &
                 abs(df$log2FC) >= th$lfc &
                 (is.na(df$log10mean) | df$log10mean >= th$expr)
    df$status <- status_of(df$passes, df$log2FC)
    df
  })

  gene_summary <- reactive({
    annotated() %>%
      group_by(gene_id, gene_name, novel) %>%
      summarise(
        n_features  = n(),
        n_sig       = sum(passes),
        n_up        = sum(passes & log2FC > 0),
        n_down      = sum(passes & log2FC < 0),
        max_abs_lfc = suppressWarnings(max(abs(log2FC[passes]), na.rm = TRUE)),
        min_qvalue  = min(qvalue, na.rm = TRUE),
        .groups     = "drop"
      ) %>%
      mutate(
        max_abs_lfc = ifelse(is.finite(max_abs_lfc), max_abs_lfc, NA_real_),
        reciprocal  = n_up > 0 & n_down > 0
      )
  })

  genes_kept <- reactive({
    gs <- gene_summary() %>% filter(n_sig > 0)
    if (isTRUE(input$reciprocal_only) && !identical(input$analysis, "DGE")) {
      gs <- gs %>% filter(reciprocal)
    }
    if (isTRUE(input$hide_novel)) gs <- gs %>% filter(!novel)
    if (isTruthy(input$gene_search)) {
      gs <- gs %>% filter(grepl(toupper(input$gene_search),
                                toupper(gene_name), fixed = TRUE))
    }
    gs %>% arrange(desc(reciprocal), desc(n_sig), desc(max_abs_lfc))
  })

  filtered_features <- reactive({
    annotated() %>%
      filter(passes, gene_id %in% genes_kept()$gene_id) %>%
      arrange(qvalue, desc(abs(log2FC)))
  })

  output$n_tested  <- renderText(format(nrow(raw()), big.mark = ","))
  output$n_passing <- renderText(format(nrow(filtered_features()), big.mark = ","))
  output$n_genes   <- renderText(format(nrow(genes_kept()), big.mark = ","))
  output$n_recip   <- renderText({
    if (identical(input$analysis, "DGE")) "n/a"
    else format(sum(genes_kept()$reciprocal), big.mark = ",")
  })

  output$overview_header <- renderText({
    req(input$analysis)
    ANALYSIS_LABEL[[input$analysis]]
  })

  # ---- Overview plot --------------------------------------------------------

  # Four columns only. Handing plot_ly the full results table serialises every
  # column of every row into the page, which is what made this hang.
  plot_data <- reactive({
    df <- annotated()
    ma <- identical(input$plot_type, "MA")

    validate(need(!ma || !all(is.na(df$log10mean)),
                  "This table has no expression column, so only the volcano is available."))

    out <- data.frame(
      x       = if (ma) df$log10mean else df$log2FC,
      y       = if (ma) df$log2FC    else df$neg_log10_q,
      status  = df$status,
      gene_id = df$gene_id,
      hover   = paste0("<b>", df$gene_name, "</b><br>", df$feature_id,
                       "<br>log2FC ", signif(df$log2FC, 3),
                       "<br>q ", signif(df$qvalue, 3)),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$x) & !is.na(out$y), , drop = FALSE]

    dropped <- 0L
    if (isTRUE(input$thin)) {
      ns <- which(out$status == "Not significant")
      if (length(ns) > NS_CAP) {
        # Significant points are never thinned.
        keep <- sample(ns, NS_CAP)
        out <- out[sort(c(setdiff(seq_len(nrow(out)), ns), keep)), , drop = FALSE]
        dropped <- length(ns) - NS_CAP
      }
    }
    attr(out, "dropped") <- dropped
    out
  })

  output$plot_note <- renderUI({
    d <- plot_data()
    dropped <- attr(d, "dropped")
    txt <- sprintf("%s points drawn", format(nrow(d), big.mark = ","))
    if (dropped > 0) {
      txt <- paste0(txt, "; ", format(dropped, big.mark = ","),
                    " non-significant points hidden for speed. ",
                    "Every significant point is shown, and tables and downloads use the full table.")
    }
    div(class = "note", style = "margin: 0.3rem 0 0 0;", txt)
  })

  output$main_plot <- renderPlotly({
    d   <- plot_data()
    ma  <- identical(input$plot_type, "MA")
    eff <- EFFECT_LABEL[[input$analysis]]

    p <- plot_ly(
      d, x = ~x, y = ~y,
      type = "scattergl", mode = "markers",
      color = ~status, colors = STATUS_COLS,
      customdata = ~gene_id,
      text = ~hover, hoverinfo = "text",
      marker = list(size = 5, opacity = 0.6, line = list(width = 0)),
      source = "main"
    ) %>%
      layout(
        xaxis = list(title = if (ma) "log10 mean expression" else eff,
                     zeroline = FALSE, gridcolor = "#E8EAEC"),
        yaxis = list(title = if (ma) eff else "-log10 q-value",
                     zeroline = FALSE, gridcolor = "#E8EAEC"),
        legend = list(orientation = "h", y = -0.2, x = 0),
        margin = list(t = 10, r = 10),
        paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)"
      ) %>%
      config(displayModeBar = FALSE) %>%
      event_register("plotly_click")

    if (ma) {
      p <- p %>% layout(shapes = list(list(
        type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = 0, y1 = 0,
        line = list(color = "#9AA1A8", width = 1, dash = "dot")
      )))
    }
    p
  })

  # ---- Gene selection -------------------------------------------------------

  selected_gene <- reactiveVal(NULL)

  observeEvent(event_data("plotly_click", source = "main"), {
    ev <- event_data("plotly_click", source = "main")
    if (!is.null(ev$customdata)) selected_gene(as.character(ev$customdata[[1]]))
  })

  observeEvent(input$feature_table_rows_selected, {
    i <- input$feature_table_rows_selected
    if (length(i)) selected_gene(filtered_features()$gene_id[i])
  })

  observeEvent(input$gene_table_rows_selected, {
    i <- input$gene_table_rows_selected
    if (length(i)) selected_gene(genes_kept()$gene_id[i])
  })

  gene_detail <- reactive({
    g <- selected_gene()
    req(g)
    th <- thresholds()
    d  <- datasets()

    bind_rows(lapply(names(d), function(a) {
      sub <- d[[a]][d[[a]]$gene_id == g, , drop = FALSE]
      if (!nrow(sub)) return(NULL)
      sub$status   <- status_of(sub$qvalue <= th$q, sub$log2FC)
      sub$analysis <- factor(a, levels = ANALYSES)
      sub[, c("analysis", "feature_id", "gene_name", "log2FC",
              "qvalue", "log10mean", "status")]
    }))
  })

  output$gene_header <- renderText({
    if (is.null(selected_gene())) return("Gene detail")
    df <- gene_detail()
    if (is.null(df) || !nrow(df)) return("Gene detail")

    parts <- df %>%
      group_by(analysis) %>%
      summarise(n = n(), sig = sum(status != "Not significant"), .groups = "drop") %>%
      mutate(txt = sprintf("%s: %d of %d significant", analysis, sig, n))

    paste0(df$gene_name[1], " \u2014 ", paste(parts$txt, collapse = "  \u00b7  "))
  })

  output$gene_plot <- renderPlot({
    validate(need(!is.null(selected_gene()),
                  "Click a point in the overview, or a row in either table, to see this gene across every loaded analysis."))

    df <- gene_detail()
    validate(need(!is.null(df) && nrow(df) > 0,
                  "That gene is not present in the loaded tables."))

    df <- df %>%
      arrange(analysis, log2FC) %>%
      mutate(row = factor(paste(analysis, feature_id),
                          levels = paste(analysis, feature_id)))

    ggplot(df, aes(x = log2FC, y = row, colour = status)) +
      geom_vline(xintercept = 0, colour = "#9AA1A8", linewidth = 0.4) +
      geom_segment(aes(x = 0, xend = log2FC, yend = row), linewidth = 0.8) +
      geom_point(aes(size = log10mean)) +
      scale_y_discrete(labels = function(x) sub("^\\S+ ", "", x)) +
      scale_colour_manual(values = STATUS_COLS, drop = FALSE, name = NULL) +
      scale_size_continuous(range = c(2, 5.5), guide = "none") +
      facet_grid(rows = vars(analysis), scales = "free", space = "free_y",
                 switch = "y") +
      labs(x = "log2 fold change (scale differs by analysis)", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position    = "bottom",
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing.y    = unit(0.9, "lines"),
        strip.placement    = "outside",
        strip.text.y.left  = element_text(angle = 0, face = "bold", size = 10),
        axis.text.y        = element_text(family = "mono", size = 9),
        plot.margin        = margin(10, 15, 5, 5)
      )
  })

  # ---- Transcript structure -------------------------------------------------

  # Parsed once per session and cached by the reactive. A combined GTF is around
  # 80 MB parsed plus 40 MB indexed, so this is the memory-heaviest thing the
  # app does; it stays optional for that reason.
  gtf_idx <- reactive({
    req(input$gtf)
    withProgress(message = "Reading transcript structures", value = 0.15, {
      g <- tryCatch(gtf_read(input$gtf$datapath), error = function(e) e)
      # need() evaluates its message eagerly, so conditionMessage() has to sit
      # behind the branch rather than inside the need().
      if (inherits(g, "error")) {
        validate(need(FALSE, paste("Could not read that GTF:",
                                   conditionMessage(g))))
      }
      validate(need(nrow(g) > 0, "That file has no exon or transcript rows."))
      incProgress(0.7, detail = "indexing")
      gtf_index(g)
    })
  })

  # DTE and DTU in the shape gene_structure() expects. feature_id is the
  # transcript for both analyses, and is the same namespace as the GTF.
  struct_stats <- reactive({
    d <- datasets()
    mk <- function(a, lfc_name, q_name) {
      if (is.null(d[[a]])) return(NULL)
      dt <- as.data.table(d[[a]][, c("feature_id", "log2FC", "qvalue")])
      setnames(dt, c("transcript_id", lfc_name, q_name))
      unique(dt, by = "transcript_id")
    }
    list(dte = mk("DTE", "dte_lfc", "dte_q"),
         dtu = mk("DTU", "dtu_lfc", "dtu_q"))
  })

  # Locus first, with the reference chosen automatically. This is what
  # repopulates the selector, so it must not depend on the selector.
  # Shiny stores uploads under generated names, but TabixFile locates the index
  # by appending .tbi to the data file's path. Restore the original names side
  # by side in one directory so the pair is discoverable.
  ref_path <- reactive({
    f <- input$ref_gtf
    if (is.null(f) || !HAVE_TABIX) return(NULL)
    d <- file.path(tempdir(), "lorerna_ref")
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    for (i in seq_len(nrow(f)))
      file.copy(f$datapath[i], file.path(d, f$name[i]), overwrite = TRUE)
    gz <- list.files(d, pattern = "\\.(gz|bgz)$", full.names = TRUE)
    if (!length(gz)) return(NULL)
    validate(need(file.exists(paste0(gz[1], ".tbi")),
                  paste0("Upload the tabix index alongside the annotation: ",
                         basename(gz[1]), ".tbi is missing.")))
    gz[1]
  })

  structure_base <- reactive({
    req(input$gtf, selected_gene())
    s <- struct_stats()
    gene_structure(selected_gene(), gtf_idx(), s$dte, s$dtu,
                   qcut = thresholds()$q, ref_gz = ref_path())
  })

  # Reprogram the reference list whenever the locus changes. Deliberately not
  # frozen: freezing blocks structure_data() until the client echoes the new
  # value back, and the membership test below already rejects a reference left
  # over from the previous gene.
  observeEvent(structure_base(), {
    st <- structure_base()
    if (!isTRUE(st$ok)) return()
    lbl <- ifelse(is.na(st$ann$oId), st$ann$transcript_id,
                  paste0(st$ann$transcript_id, "  (", st$ann$oId, ")"))
    updateSelectInput(session, "base_tx",
                      choices = setNames(st$ann$transcript_id, lbl),
                      selected = st$base_tx)
  })

  structure_data <- reactive({
    st <- structure_base()
    if (!isTRUE(st$ok)) return(st)
    # Recompute only when the user has actually moved off the automatic pick.
    if (isTruthy(input$base_tx) && !identical(input$base_tx, st$base_tx) &&
        input$base_tx %chin% st$ann$transcript_id) {
      s <- struct_stats()
      return(gene_structure(selected_gene(), gtf_idx(), s$dte, s$dtu,
                            qcut = thresholds()$q, base_tx = input$base_tx,
                            ref_gz = ref_path()))
    }
    st
  })


  output$structure_note <- renderUI({
    if (is.null(input$gtf)) {
      return(div(class = "note", style = "margin: 0.8rem 0 0 0;",
        "Upload gffcmp.combined.gtf in the sidebar to see transcript structures. ",
        "It is written to results/04_oarfish_reference/ by the pipeline."))
    }
    if (is.null(selected_gene())) {
      return(div(class = "note", style = "margin: 0.8rem 0 0 0;",
        "Select a gene — click a point in the overview, or a row in either table."))
    }
    NULL
  })

  output$structure_plot <- renderPlot({
    req(input$gtf)
    validate(need(!is.null(selected_gene()), ""))
    st <- structure_data()
    validate(need(isTRUE(st$ok), st$why))
    plot_gene_structure(st, thresholds()$q)
  }, height = function() {
    if (is.null(input$gtf) || is.null(selected_gene())) return(120)
    st <- tryCatch(structure_data(), error = function(e) NULL)
    if (is.null(st) || !isTRUE(st$ok)) return(160)
    # Every drawn row counts, reference block included. Sizing from the
    # assembled transcripts alone squeezed a locus with an annotation block
    # into half the height it needed and overlapped the row labels.
    n_rows <- nrow(st$ann) + if (isTRUE(st$ref_shown)) st$n_ref else 0L
    # 52 px a row keeps the two-line labels legible; the constant covers the
    # title and the wrapped subtitle, which runs to five or six lines.
    max(300, 210 + 52 * n_rows)
  })

  # ---- Tables ---------------------------------------------------------------

  output$feature_table <- renderDT({
    df <- filtered_features() %>%
      select(any_of(c("feature_id", "gene_name", "log2FC", "stat",
                      "pvalue", "qvalue", "log10mean", "meanInfRV")))

    datatable(
      df, selection = "single", rownames = FALSE,
      class = "compact stripe hover",
      options = list(pageLength = 12, scrollX = TRUE,
                     columnDefs = list(list(className = "mono", targets = 0)))
    ) %>%
      formatSignif(
        intersect(c("log2FC", "stat", "pvalue", "qvalue", "log10mean", "meanInfRV"),
                  names(df)),
        digits = 3)
  })

  output$gene_table <- renderDT({
    df <- genes_kept() %>%
      select(gene_name, n_features, n_sig, n_up, n_down,
             reciprocal, max_abs_lfc, min_qvalue)

    datatable(
      df, selection = "single", rownames = FALSE,
      class = "compact stripe hover",
      colnames = c("Gene", "Features", "Significant", "Up", "Down",
                   "Reciprocal", "Max |log2FC|", "Min q"),
      options = list(pageLength = 12, scrollX = TRUE)
    ) %>%
      formatSignif(c("max_abs_lfc", "min_qvalue"), digits = 3)
  })

  # ---- Downloads ------------------------------------------------------------

  output$dl_features <- downloadHandler(
    filename = function() {
      paste0("lorerna_", tolower(input$analysis), "_features_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(filtered_features() %>%
                  select(-passes, -status, -neg_log10_q, -novel), file)
    }
  )

  output$dl_genes <- downloadHandler(
    filename = function() {
      paste0("lorerna_", tolower(input$analysis), "_genes_", Sys.Date(), ".csv")
    },
    content = function(file) write_csv(genes_kept(), file)
  )
}

shinyApp(ui, server)
