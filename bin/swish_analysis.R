#!/usr/bin/env Rscript
############################################################
# swish_analysis.R — CLI wrapper for LoReRNA pipeline
# Based on proven swish_pipeline.R (bash version)
############################################################

rm(list = ls())
gc()

# ── CLI args ─────────────────────────────────────────────
get_arg <- function(args, key, default = NULL) {
  idx <- which(args == paste0("--", key))
  if (!length(idx)) return(default)
  args[idx[1L] + 1L]
}

args         <- commandArgs(trailingOnly = TRUE)
quant_dir    <- get_arg(args, "quant_dir")
tx2gene_f    <- get_arg(args, "tx2gene")
cond_csv     <- get_arg(args, "conditions")
cond_a       <- get_arg(args, "condition_a", "Control")
cond_b       <- get_arg(args, "condition_b", "KD")
results_dir  <- get_arg(args, "results_dir", "results")
plots_dir    <- get_arg(args, "plots_dir",   "plots")
logs_dir     <- get_arg(args, "logs_dir",    "logs")
min_count    <- as.integer(get_arg(args, "min_count", "10"))
min_n        <- as.integer(get_arg(args, "min_n",     "3"))
nperms       <- as.integer(get_arg(args, "nperms",    "100"))
alpha        <- as.numeric(get_arg(args, "alpha",     "0.05"))

stopifnot(!is.null(quant_dir), !is.null(tx2gene_f), !is.null(cond_csv))

for (d in c(results_dir, plots_dir, logs_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, paste0("swish_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("SWISH analysis log\n")
cat("Run started :", format(Sys.time()), "\n")
cat("quant_dir   :", quant_dir, "\n")
cat("tx2gene     :", tx2gene_f, "\n")
cat("conditions  :", cond_csv, "\n")
cat("condition_a :", cond_a, "\n")
cat("condition_b :", cond_b, "\n")
cat("alpha       :", alpha, "\n\n")

############################################################
# 1. LIBRARIES
############################################################
library(data.table)
library(arrow)
library(SummarizedExperiment)
library(fishpond)

############################################################
# 2. CONDITIONS
############################################################
cond_tbl <- fread(cond_csv)
stopifnot(all(c("sample_id", "condition") %in% colnames(cond_tbl)))

samples    <- cond_tbl$sample_id
conditions <- cond_tbl$condition

# Detect paired/batch design — NA-safe checks
has_pair  <- "pair"  %in% colnames(cond_tbl) &&
             any(!is.na(cond_tbl$pair) & nchar(as.character(cond_tbl$pair)) > 0)
has_batch <- "batch" %in% colnames(cond_tbl) &&
             any(!is.na(cond_tbl$batch) & nchar(as.character(cond_tbl$batch)) > 0)

if (has_pair)  cat("Paired design detected from 'pair' column\n")
if (has_batch) cat("Batch design detected from 'batch' column\n")
if (has_pair && has_batch) {
  cat("WARNING: both pair and batch present - pair takes precedence\n")
  has_batch <- FALSE
}

############################################################
# 3. LOAD COUNTS + LENGTHS
############################################################
counts_list <- list()
length_list <- list()

for (s in samples) {
  qf <- file.path(quant_dir, s, paste0(s, ".quant"))
  if (!file.exists(qf)) stop("quant file not found: ", qf)
  q <- fread(qf)
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

rm(counts_list, length_list, q); gc()

############################################################
# 4. LOAD INFERENTIAL REPLICATES
############################################################
inf_list <- lapply(samples, function(s) {
  as.matrix(read_parquet(file.path(quant_dir, s, paste0(s, ".infreps.pq"))))
})
names(inf_list) <- samples

stopifnot(
  "infRep row count mismatch" = all(sapply(inf_list, nrow) == length(tx_names))
)
nboot_per_sample <- sapply(inf_list, ncol)
stopifnot(
  "differing bootstrap counts" = all(nboot_per_sample == nboot_per_sample[1])
)
nboot <- nboot_per_sample[1]
cat("Inferential replicates per sample:", nboot, "\n")

inf_assays <- lapply(seq_len(nboot), function(b) {
  mat <- do.call(cbind, lapply(inf_list, function(m) m[, b]))
  colnames(mat) <- samples
  mat
})
names(inf_assays) <- paste0("infRep", seq_len(nboot))

rm(inf_list, nboot_per_sample); gc()

############################################################
# 5. TX2GENE
############################################################
tx2gene <- fread(tx2gene_f, header = FALSE, col.names = c("tx", "gene"))

m           <- match(tx_names, tx2gene$tx)
n_unmatched <- sum(is.na(m))
cat("Unmatched transcripts:", n_unmatched,
    sprintf("(%.1f%%)\n", n_unmatched / length(m) * 100))
if (n_unmatched > 0)
  warning(n_unmatched, " transcripts have no tx2gene entry")

tx2gene_matched <- tx2gene[m, ]

############################################################
# 6. BUILD SE
############################################################
se <- SummarizedExperiment(
  assays = c(list(counts = counts, length = eff_length), inf_assays)
)
rownames(se)     <- tx_names
rowData(se)$gene <- tx2gene_matched$gene
se$condition     <- factor(conditions, levels = c(cond_a, cond_b))

if (has_pair)  se$pair  <- factor(cond_tbl$pair)
if (has_batch) se$batch <- factor(cond_tbl$batch)

rm(counts, eff_length, inf_assays, tx_names, tx2gene, tx2gene_matched, m); gc()

############################################################
# HELPER: write CSVs
############################################################
write_sig_csvs <- function(df, prefix, up_label = "upregulated", dn_label = "downregulated") {
  df$sig       <- !is.na(df$qvalue) & df$qvalue < alpha
  df$direction <- ifelse(!df$sig, "ns", ifelse(df$log2FC > 0, "up", "down"))
  df <- df[order(!df$sig, df$qvalue, -abs(df$log2FC)), ]

  sig_df <- df[df$sig, ]
  sig_up <- df[df$direction == "up",   ]
  sig_dn <- df[df$direction == "down", ]

  write.csv(df,     file.path(results_dir, paste0(prefix, "_full_results.csv")),    row.names = FALSE)
  write.csv(sig_df, file.path(results_dir, paste0(prefix, "_all_significant.csv")), row.names = FALSE)
  write.csv(sig_up, file.path(results_dir, paste0(prefix, "_", up_label, ".csv")),  row.names = FALSE)
  write.csv(sig_dn, file.path(results_dir, paste0(prefix, "_", dn_label, ".csv")),  row.names = FALSE)

  cat(prefix, "tested:", nrow(df), " sig:", nrow(sig_df),
      sprintf("(%.1f%%)", nrow(sig_df)/nrow(df)*100),
      " up:", nrow(sig_up), " down:", nrow(sig_dn), "\n")
  invisible(sig_df)
}

############################################################
# 7. DTE — matches working script sections 8-10
############################################################
cat("\n=== DTE ===\n")
se_tx <- se
se_tx <- scaleInfReps(se_tx)
se_tx <- labelKeep(se_tx, minCount = min_count, minN = min_n)

cat("DTE transcripts before filter:", nrow(se_tx),
    " retained:", sum(mcols(se_tx)$keep), "\n")

se_tx <- se_tx[mcols(se_tx)$keep, ]
se_tx <- computeInfRV(se_tx)

set.seed(1)
if (has_pair) {
  y_tx <- swish(se_tx, x = "condition", pair = "pair", nperms = nperms)
} else if (has_batch) {
  y_tx <- swish(se_tx, x = "condition", cov = "batch", nperms = nperms)
} else {
  y_tx <- swish(se_tx, x = "condition", nperms = nperms)
}

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
sig_dte <- write_sig_csvs(dte_df, "DTE")

# DTE plots
sig_dte_idx <- !is.na(mcols(y_tx)$qvalue) & mcols(y_tx)$qvalue < alpha
hi_tx <- which(sig_dte_idx)[order(-mcols(y_tx)$log2FC[sig_dte_idx])]
lo_tx <- which(sig_dte_idx)[order( mcols(y_tx)$log2FC[sig_dte_idx])]

pdf(file.path(plots_dir, "DTE_plots.pdf"), width = 8, height = 6)
hist(mcols(y_tx)$pvalue, breaks = 40, col = "grey75", border = "white",
     main = "DTE: P-value distribution", xlab = "P-value", ylab = "Frequency")
plotMASwish(y_tx, alpha = alpha, main = "DTE: MA plot (FDR < 5% coloured)")
if (any(sig_dte_idx)) {
  plotInfReps(y_tx, idx = hi_tx[1], x = "condition",
              main = paste0("DTE top up-regulated\n", rownames(y_tx)[hi_tx[1]]))
  plotInfReps(y_tx, idx = lo_tx[1], x = "condition",
              main = paste0("DTE top down-regulated\n", rownames(y_tx)[lo_tx[1]]))
}
plot(mcols(y_tx)$log2FC, -log10(pmax(mcols(y_tx)$qvalue, 1e-300)),
     col = ifelse(sig_dte_idx, "firebrick", "grey70"), pch = 20, cex = 0.4,
     main = "DTE: Effect size vs significance", xlab = "log2FC", ylab = "-log10(q)")
abline(v = 0, col = "navy", lty = 2)
legend("topright", legend = c(paste0("Sig (n=", sum(sig_dte_idx), ")"), "Not sig"),
       col = c("firebrick", "grey70"), pch = 20, bty = "n")
hist(log10(mcols(y_tx)$meanInfRV + 1e-6), breaks = 30, col = "steelblue", border = "white",
     main = "DTE: Inferential uncertainty (InfRV)", xlab = "log10(mean InfRV)")
dev.off()
cat("DTE plots written\n")

# Keep se_tx for DTU — do NOT rm it here
rm(y_tx, dte_df, hi_tx, lo_tx, sig_dte_idx); gc()

############################################################
# 8. DTU — matches working script sections 11-13
#    CRITICAL: reuses se_tx from DTE (already filtered)
#    isoformProportions requires genes with >1 transcript
############################################################
cat("\n=== DTU ===\n")

# Reuse se_tx from DTE — already scaleInfReps + labelKeep + filtered
se_tx_dtu <- se_tx[!is.na(rowData(se_tx)$gene), ]
se_dtu    <- isoformProportions(se_tx_dtu, geneCol = "gene")

rm(se_tx_dtu, se_tx); gc()

se_dtu <- computeInfRV(se_dtu)

set.seed(2)
if (has_pair) {
  y_dtu <- swish(se_dtu, x = "condition", pair = "pair", nperms = nperms)
} else if (has_batch) {
  y_dtu <- swish(se_dtu, x = "condition", cov = "batch", nperms = nperms)
} else {
  y_dtu <- swish(se_dtu, x = "condition", nperms = nperms)
}

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
sig_dtu <- write_sig_csvs(dtu_df, "DTU", "increased_usage", "decreased_usage")

# DTU plots
sig_dtu_idx <- !is.na(mcols(y_dtu)$qvalue) & mcols(y_dtu)$qvalue < alpha
hi_dtu <- which(sig_dtu_idx)[order(-abs(mcols(y_dtu)$log2FC[sig_dtu_idx]))]

pdf(file.path(plots_dir, "DTU_plots.pdf"), width = 8, height = 6)
hist(mcols(y_dtu)$pvalue, breaks = 40, col = "grey75", border = "white",
     main = "DTU: P-value distribution", xlab = "P-value", ylab = "Frequency")
plotMASwish(y_dtu, alpha = alpha, main = "DTU: MA plot - isoform proportion shifts (FDR < 5% coloured)")
if (any(sig_dtu_idx)) {
  plotInfReps(y_dtu, idx = hi_dtu[1], x = "condition",
              main = paste0("DTU top hit\n", rownames(y_dtu)[hi_dtu[1]]))
}
plot(mcols(y_dtu)$log2FC, -log10(pmax(mcols(y_dtu)$qvalue, 1e-300)),
     col = ifelse(sig_dtu_idx, "darkorchid", "grey70"), pch = 20, cex = 0.4,
     main = "DTU: Proportion shift vs significance",
     xlab = "log2FC proportion", ylab = "-log10(q)")
abline(v = 0, col = "navy", lty = 2)
legend("topright", legend = c(paste0("Sig (n=", sum(sig_dtu_idx), ")"), "Not sig"),
       col = c("darkorchid", "grey70"), pch = 20, bty = "n")
dev.off()
cat("DTU plots written\n")

rm(y_dtu, se_dtu, dtu_df, hi_dtu, sig_dtu_idx); gc()

############################################################
# 9. DGE — matches working script sections 14-16
#    Uses original 'se' (unfiltered) for gene-level aggregation
############################################################
cat("\n=== DGE ===\n")

valid <- !is.na(rowData(se)$gene)
se2   <- se[valid, ]
se2$condition <- droplevels(se2$condition)

gene_condition <- se2$condition
gene_pair      <- if (has_pair)  se2$pair  else NULL
gene_batch     <- if (has_batch) se2$batch else NULL

rm(se); gc()

gene <- rowData(se2)$gene

tx_counts  <- assay(se2, "counts")
tx_lengths <- assay(se2, "length")
tx_lengths[tx_lengths == 0] <- 1

tx_rate     <- tx_counts / tx_lengths
counts_gene <- rowsum(tx_counts, group = gene)
length_gene <- counts_gene / rowsum(tx_rate, group = gene)

n_nonfinite <- sum(!is.finite(length_gene))
if (n_nonfinite > 0) {
  med_len <- median(length_gene[is.finite(length_gene)], na.rm = TRUE)
  warning(n_nonfinite, " gene-sample entries non-finite length - using median (", round(med_len, 1), ")")
  length_gene[!is.finite(length_gene)] <- med_len
}

rm(tx_counts, tx_lengths, tx_rate); gc()

inf_assay_names <- grep("^infRep", assayNames(se2), value = TRUE)
inf_gene_list <- lapply(inf_assay_names, function(a) {
  rowsum(assay(se2, a), group = gene)
})
names(inf_gene_list) <- inf_assay_names

rm(se2); gc()

se_gene <- SummarizedExperiment(
  assays = c(list(counts = counts_gene, length = length_gene), inf_gene_list)
)
se_gene$condition <- gene_condition
if (has_pair)  se_gene$pair  <- gene_pair
if (has_batch) se_gene$batch <- gene_batch

rm(counts_gene, length_gene, inf_gene_list, gene_condition); gc()

se_gene <- scaleInfReps(se_gene)
se_gene <- labelKeep(se_gene, minCount = min_count, minN = min_n)

cat("DGE genes before filter:", nrow(se_gene),
    " retained:", sum(mcols(se_gene)$keep), "\n")

se_gene <- se_gene[mcols(se_gene)$keep, ]

set.seed(3)
if (has_pair) {
  y_gene <- swish(se_gene, x = "condition", pair = "pair", nperms = nperms)
} else if (has_batch) {
  y_gene <- swish(se_gene, x = "condition", cov = "batch", nperms = nperms)
} else {
  y_gene <- swish(se_gene, x = "condition", nperms = nperms)
}

dge_df <- data.frame(
  gene_id   = rownames(y_gene),
  log10mean = mcols(y_gene)$log10mean,
  log2FC    = mcols(y_gene)$log2FC,
  stat      = mcols(y_gene)$stat,
  pvalue    = mcols(y_gene)$pvalue,
  qvalue    = mcols(y_gene)$qvalue,
  stringsAsFactors = FALSE
)
sig_dge <- write_sig_csvs(dge_df, "DGE")

# DGE plots
sig_dge_idx <- !is.na(mcols(y_gene)$qvalue) & mcols(y_gene)$qvalue < alpha
hi_gene <- which(sig_dge_idx)[order(-mcols(y_gene)$log2FC[sig_dge_idx])]
lo_gene <- which(sig_dge_idx)[order( mcols(y_gene)$log2FC[sig_dge_idx])]

pdf(file.path(plots_dir, "DGE_plots.pdf"), width = 8, height = 6)
hist(mcols(y_gene)$pvalue, breaks = 40, col = "grey75", border = "white",
     main = "DGE: P-value distribution", xlab = "P-value", ylab = "Frequency")
plotMASwish(y_gene, alpha = alpha, main = "DGE: MA plot (FDR < 5% coloured)")
if (any(sig_dge_idx)) {
  plotInfReps(y_gene, idx = hi_gene[1], x = "condition",
              main = paste0("DGE top up-regulated\n", rownames(y_gene)[hi_gene[1]]))
  plotInfReps(y_gene, idx = lo_gene[1], x = "condition",
              main = paste0("DGE top down-regulated\n", rownames(y_gene)[lo_gene[1]]))
}
plot(mcols(y_gene)$log2FC, -log10(pmax(mcols(y_gene)$qvalue, 1e-300)),
     col = ifelse(sig_dge_idx, "forestgreen", "grey70"), pch = 20, cex = 0.5,
     main = "DGE: Effect size vs significance", xlab = "log2FC", ylab = "-log10(q)")
abline(v = 0, col = "navy", lty = 2)
legend("topright", legend = c(paste0("Sig (n=", sum(sig_dge_idx), ")"), "Not sig"),
       col = c("forestgreen", "grey70"), pch = 20, bty = "n")
dev.off()
cat("DGE plots written\n")

rm(se_gene, y_gene, dge_df, hi_gene, lo_gene, sig_dge_idx); gc()

############################################################
# 10. SUMMARY
############################################################
cat("\n============================\n")
cat("SWISH RESULTS\n")
cat("============================\n")
cat("DTE significant:", nrow(sig_dte), "\n")
cat("DTU significant:", nrow(sig_dtu), "\n")
cat("DGE significant:", nrow(sig_dge), "\n")
cat("============================\n\n")

cat("results/\n")
cat("  DTE_full_results.csv  |  DTE_all_significant.csv  |  DTE_upregulated.csv  |  DTE_downregulated.csv\n")
cat("  DTU_full_results.csv  |  DTU_all_significant.csv  |  DTU_increased_usage.csv  |  DTU_decreased_usage.csv\n")
cat("  DGE_full_results.csv  |  DGE_all_significant.csv  |  DGE_upregulated.csv  |  DGE_downregulated.csv\n\n")

cat("Run finished:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(log_con)

message("Analysis complete. Log: ", log_file)
