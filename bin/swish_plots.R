#!/usr/bin/env Rscript
# =============================================================================
# swish_plots.R — publication-quality figures for LoReRNA swish results
#
# Called by modules/local/swish_plots.nf after SWISH completes.
# Reads the full-results CSVs written by swish_analysis.R and the
# per-sample oarfish quant files for expression values.
#
# Produces one multi-page PDF per analysis type:
#   DTE_publication_plots.pdf  — 5 pages
#   DTU_publication_plots.pdf  — 5 pages
#   DGE_publication_plots.pdf  — 5 pages
#   summary_panel.pdf          — 1 page: overview across all three analyses
#
# Usage:
#   Rscript swish_plots.R \
#     --results_dir  results/       \  # from SWISH (contains *_full_results.csv)
#     --quant_dir    quant_input/   \  # oarfish per-sample quant files
#     --tx2gene      tx2gene.tsv    \
#     --conditions   conditions.csv \
#     --condition_a  WT             \
#     --condition_b  KD             \
#     --plots_dir    pub_plots/     \
#     --alpha        0.05           \
#     --top_n        30
# =============================================================================

# ── 0. CLI ────────────────────────────────────────────────────────────────────

get_arg <- function(args, key, default = NULL) {
  idx <- which(args == paste0("--", key))
  if (!length(idx)) return(default)
  args[idx[1L] + 1L]
}

args         <- commandArgs(trailingOnly = TRUE)
results_dir  <- get_arg(args, "results_dir",  "results")
quant_dir    <- get_arg(args, "quant_dir")
tx2gene_f    <- get_arg(args, "tx2gene")
cond_csv     <- get_arg(args, "conditions")
cond_a       <- get_arg(args, "condition_a",  "Control")
cond_b       <- get_arg(args, "condition_b",  "KD")
plots_dir    <- get_arg(args, "plots_dir",    "pub_plots")
alpha        <- as.numeric(get_arg(args, "alpha",   "0.05"))
top_n        <- as.integer(get_arg(args, "top_n",   "30"))

stopifnot(
  "--results_dir is required" = !is.null(results_dir),
  "--quant_dir is required"   = !is.null(quant_dir),
  "--tx2gene is required"     = !is.null(tx2gene_f),
  "--conditions is required"  = !is.null(cond_csv)
)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
message("swish_plots.R started — plots_dir: ", plots_dir)

# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(patchwork)
  library(scales)
  library(viridis)
  library(RColorBrewer)
})

THEME <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 9),
    axis.title    = element_text(size = 10),
    legend.position = "bottom"
  )

COND_COLOURS <- c(cond_a, cond_b)
names(COND_COLOURS) <- c(cond_a, cond_b)

# ── 2. Load support data ──────────────────────────────────────────────────────

tx2gene   <- fread(tx2gene_f, header = FALSE, col.names = c("tx_id", "gene"))
cond_tbl  <- fread(cond_csv)
samples   <- cond_tbl$sample_id

# Load oarfish counts: normalise to TPM per sample for heatmaps
load_tpm <- function(quant_dir, samples) {
  dt_list <- lapply(samples, function(s) {
    qf <- file.path(quant_dir, s, paste0(s, ".quant"))
    stopifnot(paste0("quant file not found: ", qf) = file.exists(qf))
    q <- fread(qf, select = c("tname", "num_reads", "len"))
    q[, tpm := {
      rate <- num_reads / pmax(len, 1)
      rate / sum(rate) * 1e6
    }]
    q[, .(tname, tpm)]
  })
  names(dt_list) <- samples

  tpm_wide <- Reduce(function(a, b) merge(a, b, by = "tname", all = TRUE), dt_list)
  setnames(tpm_wide, c("tname", samples))
  mat <- as.matrix(tpm_wide[, -1, with = FALSE])
  rownames(mat) <- tpm_wide$tname
  mat[is.na(mat)] <- 0
  mat
}

message("Loading expression data...")
tpm_mat <- load_tpm(quant_dir, samples)
message("  Loaded TPM matrix: ", nrow(tpm_mat), " transcripts × ", ncol(tpm_mat), " samples")

# Add gene names to tpm_mat row annotations
rowdata <- data.frame(
  tx_id = rownames(tpm_mat),
  stringsAsFactors = FALSE
)
rowdata <- merge(rowdata, tx2gene, by = "tx_id", all.x = TRUE)
rowdata$gene[is.na(rowdata$gene)] <- rowdata$tx_id[is.na(rowdata$gene)]

# ── 3. Helper functions ───────────────────────────────────────────────────────

# Load results CSV — uses *_full_results.csv (all tested features)
# Falls back to *_all_significant.csv if full results not available
load_results <- function(prefix) {
  full_path <- file.path(results_dir, paste0(prefix, "_full_results.csv"))
  sig_path  <- file.path(results_dir, paste0(prefix, "_all_significant.csv"))

  if (file.exists(full_path)) {
    dt <- fread(full_path)
    message("  Loaded full results: ", nrow(dt), " features (", prefix, ")")
  } else if (file.exists(sig_path)) {
    dt <- fread(sig_path)
    message("  NOTE: only significant results found for ", prefix,
            " — MA/volcano plots show sig features only")
  } else {
    stop("No results CSV found for prefix: ", prefix,
         "\nLooked for:\n  ", full_path, "\n  ", sig_path)
  }
  dt
}

# Add gene name column from tx2gene lookup
annotate_gene <- function(dt) {
  if (!"gene" %in% names(dt)) {
    id_col <- if ("transcript_id" %in% names(dt)) "transcript_id" else names(dt)[1]
    dt <- merge(dt, tx2gene, by.x = id_col, by.y = "tx_id", all.x = TRUE)
    dt$gene[is.na(dt$gene)] <- dt[[id_col]][is.na(dt$gene)]
  }
  dt
}

# Label top hits by |log2FC| and lowest qvalue for ggrepel
top_labels <- function(dt, n = 15, sig_col = "qvalue") {
  sig <- dt[!is.na(dt[[sig_col]]) & dt[[sig_col]] < alpha, ]
  if (nrow(sig) == 0) return(character(0))
  sig <- sig[order(sig[[sig_col]], -abs(sig$log2FC)), ]
  head(sig$label, n)
}

# ── 4. MA plot (ggplot2) ──────────────────────────────────────────────────────

plot_ma <- function(dt, title, subtitle = "", sig_col = "qvalue",
                    x_col = "log10mean", y_col = "log2FC", colour = "#E74C3C") {
  dt <- copy(dt)
  dt$sig   <- !is.na(dt[[sig_col]]) & dt[[sig_col]] < alpha
  dt$label <- ifelse(dt$sig & dt$gene %in% top_labels(dt), dt$gene, NA_character_)

  ggplot(dt, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(colour = sig), size = 0.6, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_text_repel(
      aes(label = label), size = 2.5, colour = "black",
      max.overlaps = 20, segment.size = 0.3, segment.colour = "grey60",
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey75", "TRUE" = colour),
      labels = c("FALSE" = "Not sig", "TRUE" = paste0("FDR < ", alpha)),
      name = NULL
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Mean expression (log10)",
      y        = paste0("log2FC (", cond_b, " / ", cond_a, ")")
    ) +
    THEME
}

# ── 5. Volcano plot ───────────────────────────────────────────────────────────

plot_volcano <- function(dt, title, subtitle = "", sig_col = "qvalue",
                         colour = "#E74C3C") {
  dt <- copy(dt)
  dt$sig   <- !is.na(dt[[sig_col]]) & dt[[sig_col]] < alpha
  dt$neglog10q <- -log10(pmax(dt[[sig_col]], 1e-300))
  dt$label <- ifelse(dt$sig & dt$gene %in% top_labels(dt), dt$gene, NA_character_)

  ggplot(dt, aes(x = log2FC, y = neglog10q)) +
    geom_point(aes(colour = sig), size = 0.6, alpha = 0.5) +
    geom_vline(xintercept = 0,    linetype = "dashed", colour = "grey30") +
    geom_hline(yintercept = -log10(alpha), linetype = "dotted", colour = colour) +
    geom_text_repel(
      aes(label = label), size = 2.5, colour = "black",
      max.overlaps = 20, segment.size = 0.3, segment.colour = "grey60",
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey75", "TRUE" = colour),
      labels = c("FALSE" = "Not sig", "TRUE" = paste0("FDR < ", alpha)),
      name = NULL
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = paste0("log2FC (", cond_b, " / ", cond_a, ")"),
      y        = expression(-log[10](q-value))
    ) +
    THEME
}

# ── 6. Heatmap of top-N hits (pheatmap) ──────────────────────────────────────

plot_heatmap <- function(sig_df, tpm_mat, cond_tbl, title, n = top_n,
                         id_col = "transcript_id") {

  if (nrow(sig_df) == 0) {
    message("  Heatmap: no significant features — skipping")
    return(invisible(NULL))
  }

  top_ids <- head(sig_df[[id_col]], n)
  top_ids <- top_ids[top_ids %in% rownames(tpm_mat)]

  if (length(top_ids) == 0) {
    message("  Heatmap: no features found in TPM matrix — skipping")
    return(invisible(NULL))
  }

  mat <- tpm_mat[top_ids, cond_tbl$sample_id, drop = FALSE]

  # Z-score per row
  mat <- t(scale(t(log1p(mat))))
  mat[!is.finite(mat)] <- 0

  # Row labels: gene name if available
  gene_labels <- rowdata$gene[match(rownames(mat), rowdata$tx_id)]
  gene_labels[is.na(gene_labels)] <- rownames(mat)[is.na(gene_labels)]
  rownames(mat) <- gene_labels

  # Column annotation
  ann_col <- data.frame(
    Condition = cond_tbl$condition,
    row.names = cond_tbl$sample_id
  )
  ann_colours <- list(
    Condition = setNames(
      c("#3498DB", "#E74C3C"),
      c(cond_a, cond_b)
    )
  )

  pheatmap(mat,
    main            = title,
    annotation_col  = ann_col,
    annotation_colors = ann_colours,
    color           = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    show_rownames   = nrow(mat) <= 50,
    fontsize_row    = max(4, 10 - nrow(mat) / 5),
    fontsize_col    = 9,
    silent          = TRUE
  )
}

# ── 7. p-value histogram ──────────────────────────────────────────────────────

plot_pval_hist <- function(dt, title, sig_col = "qvalue") {
  dt <- dt[!is.na(dt[["pvalue"]]), ]
  ggplot(dt, aes(x = pvalue)) +
    geom_histogram(bins = 40, fill = "grey65", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = alpha, linetype = "dashed", colour = "#E74C3C") +
    labs(title = title, x = "P-value", y = "Frequency") +
    THEME
}

# ── 8. Summary bar chart ──────────────────────────────────────────────────────

plot_summary <- function(dte_sig, dtu_sig, dge_sig) {
  dt <- data.frame(
    analysis  = rep(c("DTE", "DTU", "DGE"), each = 2),
    direction = rep(c("Up", "Down"), 3),
    n = c(
      sum(dte_sig$log2FC > 0, na.rm = TRUE),
      sum(dte_sig$log2FC < 0, na.rm = TRUE),
      sum(dtu_sig$log2FC > 0, na.rm = TRUE),
      sum(dtu_sig$log2FC < 0, na.rm = TRUE),
      sum(dge_sig$log2FC > 0, na.rm = TRUE),
      sum(dge_sig$log2FC < 0, na.rm = TRUE)
    )
  )
  dt$direction <- factor(dt$direction, levels = c("Up", "Down"))
  dt$n_signed  <- ifelse(dt$direction == "Down", -dt$n, dt$n)

  ggplot(dt, aes(x = analysis, y = n_signed, fill = direction)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, colour = "grey30") +
    scale_fill_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB")) +
    scale_y_continuous(labels = function(x) abs(x)) +
    labs(
      title    = paste0("Significant features (FDR < ", alpha, ")"),
      subtitle = paste0(cond_b, " vs ", cond_a),
      x        = NULL, y = "Count", fill = NULL
    ) +
    THEME +
    theme(legend.position = "right")
}

# ── 9. DTE plots ──────────────────────────────────────────────────────────────

message("Plotting DTE...")
dte <- annotate_gene(load_results("DTE"))
dte_sig_csv <- file.path(results_dir, "DTE_all_significant.csv")
dte_sig <- if (file.exists(dte_sig_csv)) annotate_gene(fread(dte_sig_csv)) else dte[dte$qvalue < alpha & !is.na(dte$qvalue), ]

dte$label <- dte$gene

pdf(file.path(plots_dir, "DTE_publication_plots.pdf"), width = 8, height = 6)

# Page 1: p-value histogram
print(plot_pval_hist(dte,
  title = paste0("DTE: P-value distribution\n(", cond_b, " vs ", cond_a, ")")))

# Page 2: MA plot
print(plot_ma(dte,
  title    = paste0("DTE: MA plot (n=", nrow(dte_sig), " sig. transcripts)"),
  subtitle = paste0(cond_b, " vs ", cond_a, "  |  FDR < ", alpha, " coloured"),
  colour   = "#E74C3C"))

# Page 3: Volcano
print(plot_volcano(dte,
  title    = paste0("DTE: Volcano plot"),
  subtitle = paste0(cond_b, " vs ", cond_a, "  |  top hits labelled"),
  colour   = "#E74C3C"))

# Page 4: Heatmap
if (nrow(dte_sig) > 0)
  plot_heatmap(dte_sig, tpm_mat, cond_tbl,
    title  = paste0("DTE: Top ", min(top_n, nrow(dte_sig)), " transcripts (Z-score log1p TPM)"),
    id_col = names(dte_sig)[1])

# Page 5: top-up vs top-down bar chart
n_dte_up   <- sum(dte_sig$log2FC > 0, na.rm = TRUE)
n_dte_down <- sum(dte_sig$log2FC < 0, na.rm = TRUE)
dte_dir_df <- data.frame(
  direction = c("Up", "Down"),
  n         = c(n_dte_up, n_dte_down)
)
print(
  ggplot(dte_dir_df, aes(x = direction, y = n, fill = direction)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB")) +
    labs(title = "DTE: Sig. transcripts by direction",
         x = NULL, y = "Count", fill = NULL) +
    THEME + theme(legend.position = "none")
)

dev.off()
message("  DTE plots written to ", file.path(plots_dir, "DTE_publication_plots.pdf"))

# ── 10. DTU plots ─────────────────────────────────────────────────────────────

message("Plotting DTU...")
dtu <- annotate_gene(load_results("DTU"))
dtu_sig_csv <- file.path(results_dir, "DTU_all_significant.csv")
dtu_sig <- if (file.exists(dtu_sig_csv)) annotate_gene(fread(dtu_sig_csv)) else dtu[dtu$qvalue < alpha & !is.na(dtu$qvalue), ]

dtu$label <- dtu$gene

pdf(file.path(plots_dir, "DTU_publication_plots.pdf"), width = 8, height = 6)

print(plot_pval_hist(dtu,
  title = paste0("DTU: P-value distribution (proportion shifts)\n(", cond_b, " vs ", cond_a, ")")))

print(plot_ma(dtu,
  title    = paste0("DTU: MA plot — isoform proportion shifts (n=", nrow(dtu_sig), " sig.)"),
  subtitle = paste0("log2FC of isoform proportions  |  FDR < ", alpha, " coloured"),
  colour   = "#8E44AD"))

print(plot_volcano(dtu,
  title  = "DTU: Proportion shift vs significance",
  colour = "#8E44AD"))

# Top-N DTU heatmap using isoform-level TPM
if (nrow(dtu_sig) > 0)
  plot_heatmap(dtu_sig, tpm_mat, cond_tbl,
    title  = paste0("DTU: Top ", min(top_n, nrow(dtu_sig)), " isoforms (Z-score log1p TPM)"),
    id_col = names(dtu_sig)[1])

# Scatter: proportion shift vs mean expression
print(
  ggplot(dtu, aes(x = log10mean, y = log2FC,
                  colour = !is.na(qvalue) & qvalue < alpha)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    scale_colour_manual(values = c("FALSE" = "grey75", "TRUE" = "#8E44AD"),
                        labels = c("Not sig", paste0("FDR < ", alpha)), name = NULL) +
    labs(title = "DTU: Mean expression vs proportion shift",
         x = "Mean expression (log10)", y = "log2 Proportion shift") +
    THEME
)

dev.off()
message("  DTU plots written to ", file.path(plots_dir, "DTU_publication_plots.pdf"))

# ── 11. DGE plots ─────────────────────────────────────────────────────────────

message("Plotting DGE...")
dge <- annotate_gene(load_results("DGE"))
setnames(dge, old = names(dge)[1], new = "gene_id", skip_absent = TRUE)
dge$gene <- if ("gene" %in% names(dge)) dge$gene else dge$gene_id

dge_sig_csv <- file.path(results_dir, "DGE_all_significant.csv")
dge_sig <- if (file.exists(dge_sig_csv)) annotate_gene(fread(dge_sig_csv)) else dge[dge$qvalue < alpha & !is.na(dge$qvalue), ]
dge$label <- dge$gene

# Aggregate TPM to gene level for DGE heatmap
gene_tpm <- rowsum(tpm_mat, group = rowdata$gene[match(rownames(tpm_mat), rowdata$tx_id)])

pdf(file.path(plots_dir, "DGE_publication_plots.pdf"), width = 8, height = 6)

print(plot_pval_hist(dge,
  title = paste0("DGE: P-value distribution\n(", cond_b, " vs ", cond_a, ")")))

print(plot_ma(dge,
  title    = paste0("DGE: MA plot (n=", nrow(dge_sig), " sig. genes)"),
  subtitle = paste0(cond_b, " vs ", cond_a, "  |  FDR < ", alpha, " coloured"),
  colour   = "#27AE60"))

print(plot_volcano(dge,
  title  = "DGE: Gene expression volcano",
  colour = "#27AE60"))

# Gene-level heatmap
if (nrow(dge_sig) > 0) {
  top_genes <- head(dge_sig$gene, top_n)
  top_genes <- top_genes[top_genes %in% rownames(gene_tpm)]
  if (length(top_genes) > 0) {
    gmat <- gene_tpm[top_genes, cond_tbl$sample_id, drop = FALSE]
    gmat <- t(scale(t(log1p(gmat))))
    gmat[!is.finite(gmat)] <- 0
    ann_col <- data.frame(Condition = cond_tbl$condition, row.names = cond_tbl$sample_id)
    pheatmap(gmat,
      main           = paste0("DGE: Top ", length(top_genes), " genes (Z-score log1p TPM)"),
      annotation_col = ann_col,
      annotation_colors = list(Condition = setNames(c("#3498DB", "#E74C3C"), c(cond_a, cond_b))),
      color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      cluster_cols   = FALSE,
      fontsize_row   = max(4, 10 - length(top_genes) / 5),
      silent         = TRUE
    )
  }
}

# Top-up and top-down genes, labelled bar
n_dge_up   <- sum(dge_sig$log2FC > 0, na.rm = TRUE)
n_dge_down <- sum(dge_sig$log2FC < 0, na.rm = TRUE)
dge_dir_df <- data.frame(direction = c("Up", "Down"), n = c(n_dge_up, n_dge_down))
print(
  ggplot(dge_dir_df, aes(x = direction, y = n, fill = direction)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = c("Up" = "#27AE60", "Down" = "#2980B9")) +
    labs(title = "DGE: Sig. genes by direction", x = NULL, y = "Count") +
    THEME + theme(legend.position = "none")
)

dev.off()
message("  DGE plots written to ", file.path(plots_dir, "DGE_publication_plots.pdf"))

# ── 12. Summary panel ─────────────────────────────────────────────────────────

message("Plotting summary panel...")

pdf(file.path(plots_dir, "summary_panel.pdf"), width = 10, height = 6)
print(plot_summary(dte_sig, dtu_sig, dge_sig))
dev.off()
message("  Summary panel written to ", file.path(plots_dir, "summary_panel.pdf"))

message("")
message("=================================================")
message("swish_plots.R complete")
message("  DTE_publication_plots.pdf")
message("  DTU_publication_plots.pdf")
message("  DGE_publication_plots.pdf")
message("  summary_panel.pdf")
message("All written to: ", plots_dir)
message("=================================================")
