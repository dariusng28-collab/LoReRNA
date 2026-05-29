#!/usr/bin/env Rscript
# =============================================================================
# swish_analysis.R — fishpond DTE / DTU / DGE from oarfish quantifications
#
# Called by modules/local/swish.nf after all OARFISH jobs complete.
# All paths and parameters are passed as CLI arguments so the script
# is also testable standalone outside Nextflow.
#
# Usage:
#   swish_analysis.R \
#     --quant_dir   <path>    # base dir: quant_dir/sample_id/sample_id.quant
#     --tx2gene     <file>    # tx2gene.tsv (transcript_id <TAB> gene_id)
#     --conditions  <file>    # CSV: sample_id,condition[,pair,batch] (header required)
#     --condition_a <str>     # reference condition label   [default: Control]
#     --condition_b <str>     # test condition label        [default: KD]
#     --results_dir <path>    # output directory for CSVs  [default: results]
#     --plots_dir   <path>    # output directory for PDFs  [default: plots]
#     --logs_dir    <path>    # output directory for log   [default: logs]
#     --min_count   <int>     # labelKeep minCount          [default: 10]
#     --min_n       <int>     # labelKeep minN              [default: 3]
#     --nperms      <int>     # swish() nperms              [default: 100]
#     --alpha       <float>   # FDR threshold for CSV output [default: 0.05]
#     (no flags needed for paired/batch — auto-detected from conditions CSV columns)
#
# Changelog v1.1.1:
#   fix #3  — paired/batch design auto-detected from conditions CSV column content.
#             No flags needed — add 'pair'/'batch' column to samplesheet and the
#             pipeline activates the appropriate design automatically.
#   fix #5  — se_tx freed immediately after se_tx_dtu is created (before
#             isoformProportions), reducing peak memory during DTU analysis
#             by the full size of one transcript-level SE (~1/3 of peak RAM)
# =============================================================================

# ── 0. Parse CLI arguments ────────────────────────────────────────────────────

get_arg <- function(args, key, default = NULL) {
  idx <- which(args == paste0("--", key))
  if (length(idx) == 0L) return(default)
  if (idx[1L] + 1L > length(args)) return(default)
  args[idx[1L] + 1L]
}

args        <- commandArgs(trailingOnly = TRUE)
quant_dir   <- get_arg(args, "quant_dir")
tx2gene_f   <- get_arg(args, "tx2gene")
cond_csv    <- get_arg(args, "conditions")
cond_a      <- get_arg(args, "condition_a", "Control")
cond_b      <- get_arg(args, "condition_b", "KD")
results_dir <- get_arg(args, "results_dir", "results")
plots_dir   <- get_arg(args, "plots_dir",   "plots")
logs_dir    <- get_arg(args, "logs_dir",    "logs")
min_count   <- as.integer(get_arg(args, "min_count", "10"))
min_n       <- as.integer(get_arg(args, "min_n",     "3"))
nperms      <- as.integer(get_arg(args, "nperms",    "100"))
alpha       <- as.numeric(get_arg(args, "alpha",     "0.05"))
# Optional: paired and batch design support (fix #3)
# paired/batch design is auto-detected below from conditions CSV column content.
# No CLI flags needed.

stopifnot("--quant_dir is required"  = !is.null(quant_dir))
stopifnot("--tx2gene is required"    = !is.null(tx2gene_f))
stopifnot("--conditions is required" = !is.null(cond_csv))

# ── 1. Session log ────────────────────────────────────────────────────────────

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir,    recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir,
                      paste0("swish_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("=================================================\n")
cat("swish_analysis.R\n")
cat("Run started  :", format(Sys.time()), "\n")
cat("quant_dir    :", quant_dir,   "\n")
cat("tx2gene      :", tx2gene_f,   "\n")
cat("conditions   :", cond_csv,    "\n")
cat("condition_a  :", cond_a,      "\n")
cat("condition_b  :", cond_b,      "\n")
cat("alpha        :", alpha,       "\n")
cat("min_count    :", min_count,   "\n")
cat("min_n        :", min_n,       "\n")
cat("nperms       :", nperms,      "\n")
cat("paired/batch : auto-detected from conditions CSV\n")
cat("=================================================\n\n")

# ── 2. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(SummarizedExperiment)
  library(fishpond)
})

# ── 3. Load conditions CSV ────────────────────────────────────────────────────
# Required columns: sample_id, condition
# Optional columns: pair, batch  (auto-detected from column presence and content)

cond_tbl <- fread(cond_csv, header = TRUE)
stopifnot(
  "conditions CSV must have columns: sample_id, condition" =
    all(c("sample_id", "condition") %in% colnames(cond_tbl))
)

samples    <- cond_tbl$sample_id
conditions <- cond_tbl$condition

# Auto-detect paired design — 'pair' column present and non-empty
pair_col <- NULL
if ("pair" %in% colnames(cond_tbl)) {
  pairs_vec <- cond_tbl[["pair"]]
  pairs_vec[pairs_vec == ""] <- NA_character_
  if (!all(is.na(pairs_vec))) {
    pair_col <- "pair"
    cat("Paired design detected (non-empty 'pair' column)\n")
    print(table(pairs_vec))
    cat("\n")
  }
}

# Auto-detect batch variable — 'batch' column present and non-empty
batch_col <- NULL
if ("batch" %in% colnames(cond_tbl)) {
  batch_vec <- cond_tbl[["batch"]]
  batch_vec[batch_vec == ""] <- NA_character_
  if (!all(is.na(batch_vec))) {
    if (!is.null(pair_col)) {
      warning("Both 'pair' and 'batch' columns are present. ",
              "pair and batch are mutually exclusive in most fishpond versions. ",
              "Using paired design and ignoring batch.")
    } else {
      batch_col <- "batch"
      cat("Batch variable detected (non-empty 'batch' column)\n")
      print(table(batch_vec))
      cat("\n")
    }
  }
}

cat("Samples loaded from conditions CSV:", length(samples), "\n")
cat("Condition breakdown:\n")
print(table(conditions))
cat("\n")

stopifnot(
  "condition_a not found in conditions CSV" = cond_a %in% conditions,
  "condition_b not found in conditions CSV" = cond_b %in% conditions
)

# ── 4. Load counts + effective lengths ───────────────────────────────────────

counts_list <- list()
length_list <- list()

for (s in samples) {
  quant_file <- file.path(quant_dir, s, paste0(s, ".quant"))
  stopifnot(
    paste0("quant file not found: ", quant_file) = file.exists(quant_file)
  )
  q <- fread(quant_file)
  counts_list[[s]] <- q$num_reads
  length_list[[s]] <- q$len
  if (!exists("tx_names")) tx_names <- q$tname
}

counts     <- do.call(cbind, counts_list)
eff_length <- do.call(cbind, length_list)

rownames(counts)     <- tx_names
rownames(eff_length) <- tx_names
colnames(counts)     <- samples
colnames(eff_length) <- samples

rm(counts_list, length_list, q)
gc()

# ── 5. Load inferential replicates ───────────────────────────────────────────
# Build inf_assays directly — avoids 3D array intermediate (~1/3 memory saving)

inf_list <- lapply(samples, function(s) {
  inf_file <- file.path(quant_dir, s, paste0(s, ".infreps.pq"))
  stopifnot(
    paste0("infreps file not found: ", inf_file) = file.exists(inf_file)
  )
  as.matrix(read_parquet(inf_file))
})
names(inf_list) <- samples

stopifnot(
  "infRep row count does not match .quant transcript count — check oarfish output" =
    all(sapply(inf_list, nrow) == length(tx_names))
)
nboot_per_sample <- sapply(inf_list, ncol)
stopifnot(
  "Samples have differing numbers of inferential replicates — check oarfish output" =
    all(nboot_per_sample == nboot_per_sample[1])
)
nboot <- nboot_per_sample[1]
cat("Inferential replicates per sample:", nboot, "\n")

inf_assays <- lapply(seq_len(nboot), function(b) {
  mat <- do.call(cbind, lapply(inf_list, function(m) m[, b]))
  colnames(mat) <- samples
  mat
})
names(inf_assays) <- paste0("infRep", seq_len(nboot))

rm(inf_list, nboot_per_sample)
gc()

# ── 6. Load tx2gene + strip version suffixes ─────────────────────────────────

tx2gene <- fread(tx2gene_f, header = FALSE)
colnames(tx2gene) <- c("tx", "gene")

clean_tx       <- function(x) sub("\\..*$", "", x)
tx_names_clean <- clean_tx(tx_names)
tx2gene$tx     <- clean_tx(tx2gene$tx)

n_dup <- sum(duplicated(tx2gene$tx))
if (n_dup > 0) {
  warning(n_dup, " duplicate tx IDs in tx2gene after version stripping — keeping first occurrence")
  tx2gene <- tx2gene[!duplicated(tx2gene$tx), ]
}

m             <- match(tx_names_clean, tx2gene$tx)
n_unmatched   <- sum(is.na(m))
pct_unmatched <- n_unmatched / length(m) * 100

cat("Unmatched transcripts (no gene):", n_unmatched,
    sprintf("(%.1f%%)\n", pct_unmatched))

if (n_unmatched > 0) {
  warning(n_unmatched, " transcripts have no tx2gene entry — excluded from DGE and DTU")
}
if (pct_unmatched > 5) {
  stop("More than 5% of transcripts unmatched in tx2gene (",
       round(pct_unmatched, 1),
       "%) — verify tx2gene was built from the same reference used by oarfish")
}

tx2gene <- tx2gene[m, ]

# ── 7. Build transcript-level SE (all samples) ───────────────────────────────

se <- SummarizedExperiment(
  assays = c(
    list(counts = counts,
         length = eff_length),
    inf_assays
  )
)

rownames(se)     <- tx_names_clean
rowData(se)$gene <- tx2gene$gene

se$condition <- factor(conditions, levels = unique(conditions))

# Attach optional pair / batch columns to SE colData (fix #3)
if (!is.null(pair_col)) {
  se[[pair_col]] <- factor(cond_tbl[[pair_col]])
}
if (!is.null(batch_col)) {
  se[[batch_col]] <- factor(cond_tbl[[batch_col]])
}

cat("Note: analysis compares", cond_b, "vs", cond_a,
    "— other conditions loaded but excluded from tests\n\n")

rm(counts, eff_length, inf_assays, tx_names, tx_names_clean, tx2gene, m)
gc()

# ── 8. Helper: write analysis CSVs ───────────────────────────────────────────

write_sig_csvs <- function(sig_df, prefix, up_label, dn_label) {
  sig_up <- sig_df[sig_df$log2FC > 0, ]
  sig_dn <- sig_df[sig_df$log2FC < 0, ]

  write.csv(sig_df,  file.path(results_dir, paste0(prefix, "_all_significant.csv")), row.names = FALSE)
  write.csv(sig_up,  file.path(results_dir, paste0(prefix, "_", up_label, ".csv")),  row.names = FALSE)
  write.csv(sig_dn,  file.path(results_dir, paste0(prefix, "_", dn_label, ".csv")),  row.names = FALSE)

  cat(prefix, "significant:", nrow(sig_df), "\n")
  cat("  ", up_label, ":", nrow(sig_up), "\n")
  cat("  ", dn_label, ":", nrow(sig_dn), "\n")
  cat(prefix, "CSVs written to", results_dir, "\n\n")

  list(all = sig_df, up = sig_up, dn = sig_dn)
}

# =============================================================================
# 9. DTE — DIFFERENTIAL TRANSCRIPT EXPRESSION
# =============================================================================

se_tx <- se[, se$condition %in% c(cond_a, cond_b)]
se_tx$condition <- droplevels(se_tx$condition)
se_tx$condition <- relevel(se_tx$condition, ref = cond_a)

se_tx <- scaleInfReps(se_tx)
se_tx <- labelKeep(se_tx, minCount = min_count, minN = min_n)

cat("DTE — transcripts before labelKeep:", nrow(se_tx), "\n")
cat("DTE — transcripts retained        :", sum(mcols(se_tx)$keep),
    sprintf("(%.1f%%)\n", sum(mcols(se_tx)$keep) / nrow(se_tx) * 100))

se_tx <- se_tx[mcols(se_tx)$keep, ]

# computeInfRV must run BEFORE swish(); after swish() some fishpond versions
# reorganise the object and can no longer locate the infRep assays.
se_tx <- computeInfRV(se_tx)

set.seed(1)
# pair_col / batch_col are 'pair'/'batch' when flags are set, NULL otherwise.
y_tx <- swish(se_tx, x = "condition",
              pair  = pair_col,
              batch = batch_col,
              nperms = nperms)

dte_df <- data.frame(
  transcript_id = rownames(y_tx),
  gene_id       = rowData(y_tx)$gene,
  log10mean     = mcols(y_tx)$log10mean,
  log2FC        = mcols(y_tx)$log2FC,
  stat          = mcols(y_tx)$stat,
  pvalue        = mcols(y_tx)$pvalue,
  qvalue        = mcols(y_tx)$qvalue,
  meanInfRV     = mcols(y_tx)$meanInfRV,
  stringsAsFactors = FALSE
)

write.csv(dte_df[order(dte_df$qvalue, -abs(dte_df$log2FC)), ],
          file.path(results_dir, "DTE_full_results.csv"), row.names = FALSE)
cat("DTE full results written:", nrow(dte_df), "transcripts tested\n")
sig_dte <- dte_df[!is.na(dte_df$qvalue) & dte_df$qvalue < alpha, ]
sig_dte <- sig_dte[order(sig_dte$qvalue, -abs(sig_dte$log2FC)), ]

dte_out <- write_sig_csvs(sig_dte, "DTE", "upregulated", "downregulated")

# =============================================================================
# 10. DTE PLOTS — plots/DTE_plots.pdf  (6 pages)
# =============================================================================

sig_dte_idx <- !is.na(mcols(y_tx)$qvalue) & mcols(y_tx)$qvalue < alpha
hi_tx <- which(sig_dte_idx)[order(-mcols(y_tx)$log2FC[sig_dte_idx])]
lo_tx <- which(sig_dte_idx)[order( mcols(y_tx)$log2FC[sig_dte_idx])]

pdf(file.path(plots_dir, "DTE_plots.pdf"), width = 8, height = 6)

hist(mcols(y_tx)$pvalue,
     breaks = 40, col = "grey75", border = "white",
     main = paste0("DTE: P-value distribution\n(", cond_b, " vs ", cond_a, ")"),
     xlab = "P-value", ylab = "Frequency")

plotMASwish(y_tx, alpha = alpha,
            main = paste0("DTE: MA plot (FDR < ", alpha, " coloured)"))

if (any(sig_dte_idx)) {
  plotInfReps(y_tx, idx = hi_tx[1], x = "condition",
              main = paste0("DTE top up-regulated\n", rownames(y_tx)[hi_tx[1]]))
  plotInfReps(y_tx, idx = lo_tx[1], x = "condition",
              main = paste0("DTE top down-regulated\n", rownames(y_tx)[lo_tx[1]]))
}

# pmax(..., 1e-300): prevents log10(0) = Inf for features at the permutation
# distribution extreme (qvalue = 0 is possible with swish's plug-in approach).
plot(mcols(y_tx)$log2FC,
     -log10(pmax(mcols(y_tx)$qvalue, 1e-300)),
     col = ifelse(sig_dte_idx, "firebrick", "grey70"),
     pch = 20, cex = 0.4,
     main = paste0("DTE: Effect size vs significance\n",
                   "(horizontal banding is expected behaviour)"),
     xlab = paste0("log2 Fold Change (", cond_b, " / ", cond_a, ")"),
     ylab = "-log10(q-value)")
abline(v = 0, col = "navy", lty = 2)
legend("topright",
       legend = c(paste0("Sig (n=", sum(sig_dte_idx), ")"), "Not sig"),
       col = c("firebrick", "grey70"), pch = 20, bty = "n")

hist(log10(mcols(y_tx)$meanInfRV + 1e-6),
     breaks = 30, col = "steelblue", border = "white",
     main = "DTE: Inferential uncertainty (InfRV)\nHeavy right tail = high mapping ambiguity",
     xlab = "log10(mean InfRV)")

dev.off()
cat("DTE plots written to", file.path(plots_dir, "DTE_plots.pdf"), "\n\n")

rm(y_tx, dte_df, hi_tx, lo_tx, sig_dte_idx)
gc()

# =============================================================================
# 11. DTU — DIFFERENTIAL TRANSCRIPT USAGE
# =============================================================================

# Subset to transcripts with a valid gene assignment BEFORE freeing se_tx.
# se_tx is the only object that holds the filtered, condition-subset SE at
# this point. It must be subsetted here, then freed immediately to release
# the infRep assay memory before isoformProportions duplicates the object
# (fix #5: previously se_tx was held until the DGE rm() call on line ~431,
# keeping an extra copy of ~10-20 GB in memory through the DTU analysis).
se_tx_dtu <- se_tx[!is.na(rowData(se_tx)$gene), ]
rm(se_tx)   # free transcript SE now — isoformProportions will make a new copy
gc()

# isoformProportions: scaledTPM → per-gene isoform proportions [0,1].
# All infRep assays are converted. Single-isoform genes are dropped.
# Must run after scaleInfReps(); geneCol must match rowData column name.
se_dtu <- isoformProportions(se_tx_dtu, geneCol = "gene")

rm(se_tx_dtu)
gc()

se_dtu <- computeInfRV(se_dtu)

set.seed(2)
y_dtu <- swish(se_dtu, x = "condition",
               pair  = pair_col,
               batch = batch_col,
               nperms = nperms)

dtu_df <- data.frame(
  transcript_id = rownames(y_dtu),
  gene_id       = rowData(y_dtu)$gene,
  log10mean     = mcols(y_dtu)$log10mean,
  log2FC        = mcols(y_dtu)$log2FC,
  stat          = mcols(y_dtu)$stat,
  pvalue        = mcols(y_dtu)$pvalue,
  qvalue        = mcols(y_dtu)$qvalue,
  stringsAsFactors = FALSE
)

write.csv(dtu_df[order(dtu_df$qvalue, -abs(dtu_df$log2FC)), ],
          file.path(results_dir, "DTU_full_results.csv"), row.names = FALSE)
cat("DTU full results written:", nrow(dtu_df), "transcripts tested\n")
sig_dtu <- dtu_df[!is.na(dtu_df$qvalue) & dtu_df$qvalue < alpha, ]
sig_dtu <- sig_dtu[order(sig_dtu$qvalue, -abs(sig_dtu$log2FC)), ]

# "increased_usage"/"decreased_usage" rather than "up"/"down": the quantity
# tested is an isoform proportion [0,1], not raw expression.
dtu_out <- write_sig_csvs(sig_dtu, "DTU", "increased_usage", "decreased_usage")

# =============================================================================
# 12. DTU PLOTS — plots/DTU_plots.pdf  (4 pages)
# =============================================================================

sig_dtu_idx <- !is.na(mcols(y_dtu)$qvalue) & mcols(y_dtu)$qvalue < alpha
hi_dtu <- which(sig_dtu_idx)[order(-abs(mcols(y_dtu)$log2FC[sig_dtu_idx]))]

pdf(file.path(plots_dir, "DTU_plots.pdf"), width = 8, height = 6)

hist(mcols(y_dtu)$pvalue,
     breaks = 40, col = "grey75", border = "white",
     main = paste0("DTU: P-value distribution\n(", cond_b, " vs ", cond_a, ")"),
     xlab = "P-value", ylab = "Frequency")

plotMASwish(y_dtu, alpha = alpha,
            main = paste0("DTU: MA plot — isoform proportion shifts\n(FDR < ", alpha, " coloured)"))

# After isoformProportions() the infRep assays hold proportions [0,1].
# plotInfReps() uses them automatically.
if (any(sig_dtu_idx)) {
  plotInfReps(y_dtu, idx = hi_dtu[1], x = "condition",
              main = paste0("DTU top hit (isoform proportions)\n",
                            rownames(y_dtu)[hi_dtu[1]]))
}

plot(mcols(y_dtu)$log2FC,
     -log10(pmax(mcols(y_dtu)$qvalue, 1e-300)),
     col = ifelse(sig_dtu_idx, "darkorchid", "grey70"),
     pch = 20, cex = 0.4,
     main = "DTU: Proportion shift vs significance",
     xlab = paste0("log2 Fold Change of isoform proportion (", cond_b, " / ", cond_a, ")"),
     ylab = "-log10(q-value)")
abline(v = 0, col = "navy", lty = 2)
legend("topright",
       legend = c(paste0("Sig (n=", sum(sig_dtu_idx), ")"), "Not sig"),
       col = c("darkorchid", "grey70"), pch = 20, bty = "n")

dev.off()
cat("DTU plots written to", file.path(plots_dir, "DTU_plots.pdf"), "\n\n")

rm(y_dtu, se_dtu, dtu_df, hi_dtu, sig_dtu_idx)
gc()

# =============================================================================
# 13. DGE — DIFFERENTIAL GENE EXPRESSION
# =============================================================================

# Subset to comparison samples AND valid genes BEFORE aggregation so that
# gene-level effective lengths are weighted only by the two groups compared.
valid <- !is.na(rowData(se)$gene)
se2   <- se[valid, se$condition %in% c(cond_a, cond_b)]
se2$condition <- droplevels(se2$condition)
se2$condition <- relevel(se2$condition, ref = cond_a)

# Carry pair / batch into the gene-level SE (fix #3)
gene_condition <- se2$condition
gene_pair  <- if (!is.null(pair_col))  se2[[pair_col]]  else NULL
gene_batch <- if (!is.null(batch_col)) se2[[batch_col]] else NULL

rm(se)
gc()

gene <- rowData(se2)$gene

# Expression-weighted harmonic mean of transcript lengths per gene:
#   gene_length = sum(counts) / sum(counts / tx_length)
# Matches tximport/tximeta summarizeToGene formula exactly.
tx_counts  <- assay(se2, "counts")
tx_lengths <- assay(se2, "length")
tx_lengths[tx_lengths == 0] <- 1L

tx_rate     <- tx_counts / tx_lengths
counts_gene <- rowsum(tx_counts, group = gene)
length_gene <- counts_gene / rowsum(tx_rate, group = gene)

n_nonfinite <- sum(!is.finite(length_gene))
if (n_nonfinite > 0) {
  med_len <- median(length_gene[is.finite(length_gene)], na.rm = TRUE)
  warning(n_nonfinite, " gene-sample entries have non-finite effective length",
          " — replacing with median (", round(med_len, 1), " bp)")
  length_gene[!is.finite(length_gene)] <- med_len
}

rm(tx_counts, tx_lengths, tx_rate)
gc()

inf_assay_names <- grep("^infRep", assayNames(se2), value = TRUE)

inf_gene_list <- lapply(inf_assay_names, function(a) {
  rowsum(assay(se2, a), group = gene)
})
names(inf_gene_list) <- inf_assay_names

rm(se2)
gc()

se_gene <- SummarizedExperiment(
  assays = c(
    list(counts = counts_gene,
         length = length_gene),
    inf_gene_list
  )
)

se_gene$condition <- gene_condition
if (!is.null(pair_col))  se_gene[[pair_col]]  <- gene_pair
if (!is.null(batch_col)) se_gene[[batch_col]] <- gene_batch

rm(counts_gene, length_gene, inf_gene_list, gene_condition, gene_pair, gene_batch)
gc()

se_gene <- scaleInfReps(se_gene)
se_gene <- labelKeep(se_gene, minCount = min_count, minN = min_n)

cat("DGE — genes before labelKeep:", nrow(se_gene), "\n")
cat("DGE — genes retained        :", sum(mcols(se_gene)$keep),
    sprintf("(%.1f%%)\n", sum(mcols(se_gene)$keep) / nrow(se_gene) * 100))

se_gene <- se_gene[mcols(se_gene)$keep, ]

set.seed(3)
y_gene <- swish(se_gene, x = "condition",
                pair  = pair_col,
                batch = batch_col,
                nperms = nperms)

dge_df <- data.frame(
  gene_id   = rownames(y_gene),
  log10mean = mcols(y_gene)$log10mean,
  log2FC    = mcols(y_gene)$log2FC,
  stat      = mcols(y_gene)$stat,
  pvalue    = mcols(y_gene)$pvalue,
  qvalue    = mcols(y_gene)$qvalue,
  stringsAsFactors = FALSE
)

write.csv(dge_df[order(dge_df$qvalue, -abs(dge_df$log2FC)), ],
          file.path(results_dir, "DGE_full_results.csv"), row.names = FALSE)
cat("DGE full results written:", nrow(dge_df), "genes tested\n")
sig_dge <- dge_df[!is.na(dge_df$qvalue) & dge_df$qvalue < alpha, ]
sig_dge <- sig_dge[order(sig_dge$qvalue, -abs(sig_dge$log2FC)), ]

dge_out <- write_sig_csvs(sig_dge, "DGE", "upregulated", "downregulated")

# =============================================================================
# 14. DGE PLOTS — plots/DGE_plots.pdf  (5 pages)
# =============================================================================

sig_dge_idx <- !is.na(mcols(y_gene)$qvalue) & mcols(y_gene)$qvalue < alpha
hi_gene <- which(sig_dge_idx)[order(-mcols(y_gene)$log2FC[sig_dge_idx])]
lo_gene <- which(sig_dge_idx)[order( mcols(y_gene)$log2FC[sig_dge_idx])]

pdf(file.path(plots_dir, "DGE_plots.pdf"), width = 8, height = 6)

hist(mcols(y_gene)$pvalue,
     breaks = 40, col = "grey75", border = "white",
     main = paste0("DGE: P-value distribution\n(", cond_b, " vs ", cond_a, ")"),
     xlab = "P-value", ylab = "Frequency")

plotMASwish(y_gene, alpha = alpha,
            main = paste0("DGE: MA plot (FDR < ", alpha, " coloured)"))

if (any(sig_dge_idx)) {
  plotInfReps(y_gene, idx = hi_gene[1], x = "condition",
              main = paste0("DGE top up-regulated\n", rownames(y_gene)[hi_gene[1]]))
  plotInfReps(y_gene, idx = lo_gene[1], x = "condition",
              main = paste0("DGE top down-regulated\n", rownames(y_gene)[lo_gene[1]]))
}

plot(mcols(y_gene)$log2FC,
     -log10(pmax(mcols(y_gene)$qvalue, 1e-300)),
     col = ifelse(sig_dge_idx, "forestgreen", "grey70"),
     pch = 20, cex = 0.5,
     main = "DGE: Effect size vs significance",
     xlab = paste0("log2 Fold Change (", cond_b, " / ", cond_a, ")"),
     ylab = "-log10(q-value)")
abline(v = 0, col = "navy", lty = 2)
legend("topright",
       legend = c(paste0("Sig (n=", sum(sig_dge_idx), ")"), "Not sig"),
       col = c("forestgreen", "grey70"), pch = 20, bty = "n")

dev.off()
cat("DGE plots written to", file.path(plots_dir, "DGE_plots.pdf"), "\n\n")

rm(se_gene, y_gene, dge_df, hi_gene, lo_gene, sig_dge_idx)
gc()

# =============================================================================
# 15. SUMMARY
# =============================================================================

cat("=================================================\n")
cat("SWISH RESULTS\n")
cat("Comparison:", cond_b, "vs", cond_a, "\n")
if (!is.null(pair_col))  cat("Paired design  : TRUE (auto-detected from 'pair' column)\n")
if (!is.null(batch_col)) cat("Batch variable : TRUE (auto-detected from 'batch' column)\n")
if (is.null(pair_col) && is.null(batch_col)) cat("Design         : unpaired, no batch\n")
cat("FDR threshold:", alpha, "\n")
cat("-------------------------------------------------\n")
cat("DTE transcripts :", nrow(sig_dte), "\n")
cat("DTU transcripts :", nrow(sig_dtu), "\n")
cat("DGE genes       :", nrow(sig_dge), "\n")
cat("=================================================\n\n")

cat("Results CSVs in:", results_dir, "\n")
cat("  DTE_full_results.csv  |  DTE_all_significant.csv  |  DTE_upregulated.csv  |  DTE_downregulated.csv\n")
cat("  DTU_full_results.csv  |  DTU_all_significant.csv  |  DTU_increased_usage.csv  |  DTU_decreased_usage.csv\n")
cat("  DGE_full_results.csv  |  DGE_all_significant.csv  |  DGE_upregulated.csv  |  DGE_downregulated.csv\n\n")

cat("Plots in:", plots_dir, "\n")
cat("  DTE_plots.pdf — 6 pages\n")
cat("  DTU_plots.pdf — 4 pages\n")
cat("  DGE_plots.pdf — 5 pages\n\n")

cat("Run finished:", format(Sys.time()), "\n")
cat("Log written to:", log_file, "\n")

sink(type = "message")
sink(type = "output")
close(log_con)

message("swish_analysis.R complete. Log: ", log_file)
