#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(SummarizedExperiment)
  library(fishpond)
  library(arrow)
})

option_list <- list(
  make_option("--counts", type = "character"),
  make_option("--manifest", type = "character"),
  make_option("--tx2gene", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--condition-col", type = "character", default = "condition"),
  make_option("--pair-col", type = "character", default = NA_character_),
  make_option("--batch-col", type = "character", default = NA_character_),
  make_option("--min-count", type = "integer", default = 10),
  make_option("--min-n", type = "integer", default = 3),
  make_option("--nperms", type = "integer", default = 100),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--outdir", type = "character", default = ".")
)

opts <- parse_args(OptionParser(option_list = option_list))

counts_tbl <- read_tsv(opts$counts, show_col_types = FALSE)
manifest_tbl <- read_tsv(opts$manifest, show_col_types = FALSE)
metadata_tbl <- read_csv(opts$metadata, show_col_types = FALSE)
tx2gene_tbl <- read_tsv(opts$tx2gene, show_col_types = FALSE)

required_metadata_cols <- c("sample", opts$`condition-col`)
missing_cols <- setdiff(required_metadata_cols, colnames(metadata_tbl))
if (length(missing_cols) > 0) {
  stop(sprintf("Metadata is missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

sample_names <- metadata_tbl$sample
count_cols <- intersect(sample_names, colnames(counts_tbl))
if (length(count_cols) != length(sample_names)) {
  stop("Counts matrix does not contain all metadata sample columns.")
}

counts_mat <- counts_tbl %>%
  select(all_of(sample_names)) %>%
  as.matrix()
rownames(counts_mat) <- counts_tbl$transcript_id
mode(counts_mat) <- "numeric"

row_data <- counts_tbl %>%
  select(transcript_id, gene_id, gene_name) %>%
  left_join(tx2gene_tbl, by = "transcript_id", suffix = c("", ".tx2gene")) %>%
  mutate(
    gene_id = coalesce(gene_id, gene_id.tx2gene),
    gene_name = coalesce(gene_name, gene_name.tx2gene)
  ) %>%
  select(transcript_id, gene_id, gene_name)

build_lengths <- function() {
  length_mat <- matrix(1, nrow = nrow(counts_mat), ncol = ncol(counts_mat))
  rownames(length_mat) <- rownames(counts_mat)
  colnames(length_mat) <- colnames(counts_mat)

  for (sample in sample_names) {
    sample_manifest <- manifest_tbl %>% filter(sample == !!sample)
    if (nrow(sample_manifest) != 1) {
      stop(sprintf("Manifest requires exactly one row for sample '%s'", sample))
    }
    quant_tbl <- read_tsv(sample_manifest$quant[[1]], show_col_types = FALSE)
    id_col <- intersect(c("transcript_id", "target_id", "target_name", "name", "Name"), colnames(quant_tbl))
    id_col <- if (length(id_col) > 0) id_col[[1]] else colnames(quant_tbl)[[1]]
    len_col <- intersect(c("length", "eff_length", "effective_length", "len"), colnames(quant_tbl))
    if (length(len_col) == 0) {
      next
    }
    aligned <- quant_tbl %>%
      transmute(transcript_id = .data[[id_col]], length = .data[[len_col[[1]]]]) %>%
      right_join(tibble(transcript_id = rownames(counts_mat)), by = "transcript_id") %>%
      mutate(length = coalesce(as.numeric(length), 1))
    length_mat[, sample] <- aligned$length
  }

  length_mat
}

load_infreps <- function() {
  paths <- manifest_tbl$infreps
  if (any(is.na(paths) | paths == "")) {
    stop("Swish requires inferential replicates for all samples, but one or more infreps paths are missing.")
  }

  infrep_tables <- vector("list", length(sample_names))
  names(infrep_tables) <- sample_names
  n_boot <- NULL

  for (sample in sample_names) {
    sample_manifest <- manifest_tbl %>% filter(sample == !!sample)
    quant_tbl <- read_tsv(sample_manifest$quant[[1]], show_col_types = FALSE)
    inf_tbl <- read_parquet(sample_manifest$infreps[[1]]) %>% as_tibble()
    if (nrow(inf_tbl) != nrow(quant_tbl)) {
      stop(sprintf("Inferential replicate rows do not match quant rows for sample '%s'", sample))
    }
    id_col <- intersect(c("transcript_id", "target_id", "target_name", "name", "Name"), colnames(quant_tbl))
    id_col <- if (length(id_col) > 0) id_col[[1]] else colnames(quant_tbl)[[1]]
    inf_tbl <- bind_cols(tibble(transcript_id = quant_tbl[[id_col]]), inf_tbl)
    infrep_tables[[sample]] <- inf_tbl
    sample_boot <- ncol(inf_tbl) - 1
    if (is.null(n_boot)) {
      n_boot <- sample_boot
    } else if (sample_boot != n_boot) {
      stop("All samples must have the same number of inferential replicates.")
    }
  }

  assays <- list()
  for (idx in seq_len(n_boot)) {
    assay_name <- sprintf("infRep%d", idx)
    assay_mat <- matrix(0, nrow = nrow(counts_mat), ncol = ncol(counts_mat))
    rownames(assay_mat) <- rownames(counts_mat)
    colnames(assay_mat) <- colnames(counts_mat)
    for (sample in sample_names) {
      tbl <- infrep_tables[[sample]]
      values <- tbl[[idx + 1]]
      names(values) <- tbl$transcript_id
      assay_mat[, sample] <- as.numeric(values[rownames(counts_mat)])
    }
    assays[[assay_name]] <- assay_mat
  }
  assays
}

length_mat <- build_lengths()
inf_assays <- load_infreps()

metadata_tbl <- metadata_tbl %>%
  filter(sample %in% sample_names) %>%
  slice(match(sample_names, sample))

col_data <- DataFrame(metadata_tbl)
row_data_df <- DataFrame(row_data)
rownames(row_data_df) <- row_data$transcript_id

se <- SummarizedExperiment(
  assays = c(list(counts = counts_mat, length = length_mat), inf_assays),
  colData = col_data,
  rowData = row_data_df
)

se <- scaleInfReps(se, lengthCorrect = FALSE)
se <- labelKeep(se, minCount = opts$`min-count`, minN = opts$`min-n`, x = opts$`condition-col`)
se <- se[mcols(se)$keep, ]

design_args <- list(y = se, x = opts$`condition-col`, nperms = opts$`nperms`)
if (!is.na(opts$`pair-col`) && nzchar(opts$`pair-col`)) {
  design_args$pair <- opts$`pair-col`
}
if (!is.na(opts$`batch-col`) && nzchar(opts$`batch-col`)) {
  design_args$cov <- opts$`batch-col`
}

set.seed(1)
dte <- do.call(swish, design_args)

iso <- isoformProportions(se, geneCol = "gene_id", quiet = TRUE)
set.seed(1)
dtu <- do.call(swish, c(list(y = iso), design_args[names(design_args) != "y"]))

dte_res <- as_tibble(mcols(dte)) %>%
  mutate(transcript_id = rownames(dte), gene_id = rowData(dte)$gene_id) %>%
  arrange(qvalue, desc(abs(log2FC)))

dtu_res <- as_tibble(mcols(dtu)) %>%
  mutate(transcript_id = rownames(dtu), gene_id = rowData(dtu)$gene_id) %>%
  arrange(qvalue, desc(abs(log2FC)))

write_tsv(dte_res, file.path(opts$outdir, "swish.dte.tsv"))
write_tsv(dtu_res, file.path(opts$outdir, "swish.dtu.tsv"))
