#########################################
# LUAD age-bin transcriptomic analysis
# Volcano plots and paired heatmaps
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(pheatmap)
  library(glue)
  library(tibble)
  library(stringr)
})

# Working directory
# Run this script from the project folder containing: dea/, prepro/, and meta/.
ROOT <- getwd()

DEA_DIR <- file.path(ROOT, "dea")
PREPRO  <- file.path(ROOT, "prepro")
META    <- file.path(ROOT, "meta")
FIGS    <- file.path(ROOT, "figs")
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
GSE_IDS <- c("GSE115002", "GSE75037")
BINS    <- c("young", "old")

BIN_LABELS <- c(
  young = "<60",
  old   = "≥60"
)

FDR_CUT <- 0.05
LFC_MIN <- 1
TOP_N   <- 50

`%||%` <- function(a, b) if (!is.null(a)) a else b

as_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x
}

# Load paired differential expression results
load_dea <- function(gse, bin) {
  candidate_files <- c(
    file.path(DEA_DIR, sprintf("%s___TvN_%s.tsv", gse, bin)),
    file.path(DEA_DIR, sprintf("%s__TvN_%s.tsv",  gse, bin))
  )
  f <- candidate_files[file.exists(candidate_files)][1]
  if (is.na(f)) stop("DEA file not found for ", gse, " / ", bin, " in: ", DEA_DIR)

  tt <- readr::read_tsv(f, show_col_types = FALSE)

  if (!"GENE" %in% names(tt)) {
    gcol <- c("Gene", "SYMBOL", "symbol", "gene", "ID")
    gcol <- gcol[gcol %in% names(tt)][1] %||% stop("No gene column found in: ", f)
    names(tt)[names(tt) == gcol] <- "GENE"
  }

  if (!"FDR" %in% names(tt)) {
    if ("adj.P.Val" %in% names(tt)) {
      tt$FDR <- tt$adj.P.Val
    } else if ("P.Value" %in% names(tt)) {
      tt$FDR <- p.adjust(tt$P.Value, method = "BH")
    } else {
      stop("No FDR/adj.P.Val/P.Value column found in: ", f)
    }
  }

  if (!"logFC" %in% names(tt)) stop("No logFC column found in: ", f)

  tt %>%
    dplyr::mutate(GENE = as_chr(GENE))
}

# Load SYMBOL-level expression matrix
load_expr <- function(gse) {
  f <- file.path(PREPRO, paste0(gse, "_expr_SYMBOL.tsv"))
  if (!file.exists(f)) stop("Expression file not found: ", f)

  df <- readr::read_tsv(f, show_col_types = FALSE)
  rn <- df$GENE
  as.matrix(data.frame(row.names = rn, df[, -1, drop = FALSE], check.names = FALSE))
}

# Load metadata
load_meta <- function(gse) {
  f <- file.path(META, paste0(gse, "_meta.csv"))
  if (!file.exists(f)) stop("Metadata file not found: ", f)

  readr::read_csv(f, show_col_types = FALSE) %>%
    dplyr::mutate(
      sample_id  = as_chr(sample_id),
      patient_id = as_chr(patient_id),
      tissue     = dplyr::case_when(
        grepl("^normal$", tissue, ignore.case = TRUE) ~ "Normal",
        grepl("^tumor$",  tissue, ignore.case = TRUE) ~ "Tumor",
        TRUE ~ as_chr(tissue)
      ),
      age_bin = as_chr(age_bin)
    )
}

# Keep strict 1 Normal + 1 Tumor pairs within each age bin
paired_subset <- function(E, M, bin_label) {
  mm0 <- M %>%
    dplyr::filter(
      !is.na(sample_id),
      !is.na(patient_id),
      tissue %in% c("Normal", "Tumor"),
      age_bin == bin_label
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)

  ids <- intersect(colnames(E), mm0$sample_id)
  if (length(ids) < 4) return(list(Emm = NULL, mm = NULL))

  Emm <- E[, ids, drop = FALSE]
  mm  <- mm0[match(ids, mm0$sample_id), , drop = FALSE]

  ct <- mm %>%
    dplyr::count(patient_id, tissue) %>%
    tidyr::pivot_wider(names_from = tissue, values_from = n, values_fill = 0)

  if (!"Normal" %in% names(ct)) ct$Normal <- 0L
  if (!"Tumor"  %in% names(ct)) ct$Tumor  <- 0L

  keep_p <- ct$patient_id[ct$Normal == 1L & ct$Tumor == 1L]
  mm <- mm %>% dplyr::filter(patient_id %in% keep_p)
  if (!nrow(mm)) return(list(Emm = NULL, mm = NULL))

  Emm <- Emm[, mm$sample_id, drop = FALSE]
  stopifnot(all(colnames(Emm) == mm$sample_id))

  list(Emm = Emm, mm = mm)
}

# Save paired heatmap and the plotted gene list
save_heatmap <- function(E, M, tt, gse, bin, out_dir = FIGS, top_n = TOP_N) {
  tt2 <- tt %>% dplyr::arrange(FDR, dplyr::desc(abs(logFC)))

  up <- tt2 %>%
    dplyr::filter(FDR <= FDR_CUT, logFC >= LFC_MIN) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(GENE)

  dn <- tt2 %>%
    dplyr::filter(FDR <= FDR_CUT, logFC <= -LFC_MIN) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(GENE)

  genes <- unique(c(up, dn))

  if (!length(genes)) {
    message("No genes passed FDR/logFC thresholds for heatmap: ", gse, " / ", bin)
    return(invisible(NULL))
  }

  sub <- paired_subset(E, M, bin)
  Emm <- sub$Emm
  mm  <- sub$mm

  if (is.null(Emm)) {
    message("No strict 1N-1T pairs available for heatmap: ", gse, " / ", bin)
    return(invisible(NULL))
  }

  genes <- intersect(genes, rownames(Emm))
  if (!length(genes)) {
    message("Top genes were not found in expression matrix: ", gse, " / ", bin)
    return(invisible(NULL))
  }

  X <- Emm[genes, , drop = FALSE]

  # Order columns by patient, with Normal followed by Tumor.
  mm <- mm %>%
    dplyr::mutate(Tissue = factor(tissue, levels = c("Normal", "Tumor"))) %>%
    dplyr::arrange(patient_id, Tissue)

  X <- X[, mm$sample_id, drop = FALSE]

  # Row-wise z-score
  z <- t(scale(t(X)))
  z[!is.finite(z)] <- 0

  ann_col <- data.frame(
    patient_id = mm$patient_id,
    tissue     = mm$tissue,
    age_bin    = mm$age_bin,
    row.names  = mm$sample_id
  )

  ann_colors <- list(
    tissue  = c(Normal = "#4daf4a", Tumor = "#e41a1c"),
    age_bin = c(young  = "#377eb8", old   = "#984ea3")
  )

  fn_png   <- file.path(out_dir, glue::glue("Heatmap_{gse}_{bin}_top{top_n}.png"))
  fn_pdf   <- file.path(out_dir, glue::glue("Heatmap_{gse}_{bin}_top{top_n}.pdf"))
  fn_genes <- file.path(out_dir, glue::glue("Heatmap_{gse}_{bin}_top{top_n}_genes.tsv"))

  readr::write_tsv(tibble::tibble(GENE = genes), fn_genes)

  grDevices::png(fn_png, width = 7.5, height = 7.5, units = "in", res = 300)
  pheatmap::pheatmap(
    z,
    color                    = colorRampPalette(c("navy", "white", "firebrick3"))(256),
    annotation_col           = ann_col,
    annotation_colors        = ann_colors,
    show_colnames            = FALSE,
    show_rownames            = TRUE,
    fontsize_row             = 6.5,
    clustering_distance_cols = "correlation",
    clustering_method        = "complete",
    cluster_cols             = FALSE,
    border_color             = NA,
    filename                 = NA
  )
  grDevices::dev.off()

  grDevices::pdf(fn_pdf, width = 7.5, height = 7.5)
  pheatmap::pheatmap(
    z,
    color                    = colorRampPalette(c("navy", "white", "firebrick3"))(256),
    annotation_col           = ann_col,
    annotation_colors        = ann_colors,
    show_colnames            = FALSE,
    show_rownames            = TRUE,
    fontsize_row             = 6.5,
    clustering_distance_cols = "correlation",
    clustering_method        = "complete",
    cluster_cols             = FALSE,
    border_color             = NA,
    filename                 = NA
  )
  grDevices::dev.off()

  message("Heatmap saved: ", fn_png)
}

# Optional labels for volcano plots
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (!have_ggrepel) message("Package ggrepel is not installed; volcano plots will be generated without labels.")

for (gse in GSE_IDS) {
  message("\n== ", gse, " ==")

  E <- load_expr(gse)
  M <- load_meta(gse)

  for (bin in BINS) {
    tt <- load_dea(gse, bin)

    df <- tt %>%
      dplyr::mutate(
        sig = dplyr::case_when(
          FDR <= FDR_CUT & logFC >=  LFC_MIN ~ "UP (Tumor)",
          FDR <= FDR_CUT & logFC <= -LFC_MIN ~ "DOWN (Tumor)",
          TRUE ~ "NS"
        ),
        sig = factor(sig, levels = c("UP (Tumor)", "DOWN (Tumor)", "NS")),
        negLog10FDR = -log10(pmax(FDR, .Machine$double.xmin))
      )

    lab_up <- df %>%
      dplyr::filter(sig == "UP (Tumor)") %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 15)

    lab_dn <- df %>%
      dplyr::filter(sig == "DOWN (Tumor)") %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 15)

    labs <- dplyr::bind_rows(lab_up, lab_dn)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = negLog10FDR)) +
      ggplot2::geom_point(ggplot2::aes(color = sig), alpha = 0.7, size = 0.8) +
      ggplot2::scale_color_manual(
        values = c(
          "UP (Tumor)"   = "#d62728",
          "DOWN (Tumor)" = "#08519c",
          "NS"           = "grey75"
        )
      ) +
      ggplot2::geom_vline(xintercept = c(-LFC_MIN, LFC_MIN),
                          linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = -log10(FDR_CUT),
                          linetype = "dashed", linewidth = 0.3) +
      ggplot2::labs(
        title = glue::glue("{gse} — {BIN_LABELS[[bin]]}"),
        x     = "log2 fold-change (Tumor vs Normal)",
        y     = "-log10(FDR)",
        color = NULL
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        plot.title      = ggplot2::element_text(face = "bold")
      )

    if (have_ggrepel && nrow(labs)) {
      p <- p + ggrepel::geom_text_repel(
        data = labs,
        ggplot2::aes(label = GENE),
        size = 2.8,
        max.overlaps       = Inf,
        box.padding        = 0.25,
        point.padding      = 0.15,
        min.segment.length = 0
      )
    }

    fn_png <- file.path(FIGS, glue::glue("Volcano_{gse}_{bin}.png"))
    fn_pdf <- file.path(FIGS, glue::glue("Volcano_{gse}_{bin}.pdf"))

    ggplot2::ggsave(fn_png, p, width = 5.2, height = 4.8, dpi = 300)
    ggplot2::ggsave(fn_pdf, p, width = 5.2, height = 4.8, dpi = 300)

    message("Volcano saved: ", fn_png)

    save_heatmap(E, M, tt, gse, bin, out_dir = FIGS, top_n = TOP_N)
  }
}

message("\nDone. Figures saved in: ", FIGS)
