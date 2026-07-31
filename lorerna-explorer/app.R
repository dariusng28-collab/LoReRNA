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
                   DT = "0.20", readr = "2.0.0")
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


options(shiny.maxRequestSize = 500 * 1024^2)

# Above this many non-significant points, thin them for display.
NS_CAP <- 15000

COL_UP   <- "#B03A48"
COL_DOWN <- "#35618F"
COL_NS   <- "#C4C9CE"
STATUS_COLS <- setNames(c(COL_UP, COL_DOWN, COL_NS), STATUS_LEVELS)

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
    nav_panel("Genes",    DTOutput("gene_table"))
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
