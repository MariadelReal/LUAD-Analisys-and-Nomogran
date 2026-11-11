#########################################
# ssGSEA consensus summaries and plots
# LUAD paired tumor-normal analysis by age bin
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(purrr)
})

# Avoid namespace conflicts with Bioconductor packages
select     <- dplyr::select
filter     <- dplyr::filter
arrange    <- dplyr::arrange
mutate     <- dplyr::mutate
rename     <- dplyr::rename
group_by   <- dplyr::group_by
summarise  <- dplyr::summarise
slice_head <- dplyr::slice_head

# -----------------------------
# Configuration
# -----------------------------

# Run this script from the project directory that contains the ssgsea/ folder.
ROOT    <- getwd()
OUT_DIR <- file.path(ROOT, "ssgsea")
CONS    <- file.path(OUT_DIR, "consensus")
dir.create(CONS, recursive = TRUE, showWarnings = FALSE)

GSE_IDS <- c("GSE115002", "GSE75037")
BINS    <- c("young", "old")
collection_file <- file.path(OUT_DIR,"collections_used.txt")
if (!file.exists(collection_file)) stop("Ejecuta primero el paso GEO 11 de esta carpeta.")
COLL <- readLines(collection_file)
stopifnot(length(COLL)==2L, COLL[1]=="H", COLL[2] %in% c("Reactome","KEGG"))

FDR_CUT <- 0.05
TOP_N   <- 15
MIN_LFC <- 0.10

# -----------------------------
# Helper functions
# -----------------------------

dea_path <- function(gse, coll, bin) {
  file.path(OUT_DIR, sprintf("%s___ssgsea_DEA_%s_%s.tsv", gse, coll, bin))
}

read_dea <- function(gse, coll, bin, fdr = FDR_CUT, min_lfc = MIN_LFC) {
  f <- dea_path(gse, coll, bin)
  if (!file.exists(f)) stop("Falta la tabla del paso 11: ",f)

  df <- readr::read_tsv(f, show_col_types = FALSE)

  if (!"Pathway" %in% names(df)) {
    names(df)[1] <- "Pathway"
  }

  if (!"FDR" %in% names(df)) {
    if ("adj.P.Val" %in% names(df)) {
      df$FDR <- df$adj.P.Val
    } else if ("P.Value" %in% names(df)) {
      df$FDR <- p.adjust(df$P.Value, method = "BH")
    } else {
      stop("No FDR/adj.P.Val/P.Value column found in: ", f)
    }
  }

  if (!"logFC" %in% names(df)) {
    stop("No logFC column found in: ", f)
  }

  df %>%
    mutate(
      logFC = as.numeric(logFC),
      FDR   = as.numeric(FDR),
      GSE   = gse,
      Coll  = coll,
      Bin   = bin,
      sig   = FDR <= fdr & abs(logFC) >= min_lfc
    )
}

pretty_path <- function(x, width = 40) {
  x %>%
    str_remove("^(HALLMARK|KEGG|REACTOME|WP|PID)_") %>%
    str_replace_all("_", " ") %>%
    str_trim() %>%
    str_wrap(width)
}

make_signed_score <- function(df) {
  df %>%
    mutate(score = -log10(pmax(FDR, 1e-300)) * sign(logFC))
}

# -----------------------------
# 1. Load ssGSEA differential pathway results
# -----------------------------

combos <- expand.grid(
  GSE  = GSE_IDS,
  Coll = COLL,
  Bin  = BINS,
  stringsAsFactors = FALSE
)

dea_list <- purrr::pmap(combos, function(GSE, Coll, Bin) {
  read_dea(GSE, Coll, Bin)
})

dea_list <- purrr::compact(dea_list)

if (!length(dea_list)) {
  stop("No ssGSEA DEA files were loaded. Check the ssgsea/ directory and COLL setting.")
}

all_dea <- dplyr::bind_rows(dea_list)
stopifnot(nrow(all_dea) > 0)

# -----------------------------
# 2. Count significant pathways
# -----------------------------

counts <- all_dea %>%
  group_by(GSE, Coll, Bin) %>%
  summarise(
    n_sig = sum(sig, na.rm = TRUE),
    n_tot = n(),
    .groups = "drop"
  )

p_counts <- ggplot(counts, aes(x = Bin, y = n_sig, fill = GSE)) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(~Coll, scales = "free_y") +
  labs(
    title = "Significant pathways by age bin and dataset",
    y = sprintf("FDR <= %.02f and |logFC| >= %.2f", FDR_CUT, MIN_LFC),
    x = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(CONS, "counts_sig_paths.png"), p_counts, width = 10, height = 5, dpi = 200)
readr::write_tsv(counts, file.path(CONS, "counts_sig_paths.tsv"))

# -----------------------------
# 3. Consensus significant pathways between datasets
# -----------------------------

for (coll in COLL) {
  for (bin in BINS) {
    sub <- all_dea %>% filter(Coll == coll, Bin == bin, sig)

    x1 <- sub %>% filter(GSE == GSE_IDS[1]) %>% pull(Pathway) %>% unique()
    x2 <- sub %>% filter(GSE == GSE_IDS[2]) %>% pull(Pathway) %>% unique()

    consensus_paths <- intersect(x1, x2)

    readr::write_tsv(
      tibble(Pathway = consensus_paths),
      file.path(CONS, sprintf("consensus_%s_%s.tsv", coll, bin))
    )
  }
}

# -----------------------------
# 4. Heatmaps of top directed pathway scores
# -----------------------------

for (coll in COLL) {
  for (bin in BINS) {
    sub <- all_dea %>%
      filter(Coll == coll, Bin == bin) %>%
      make_signed_score()

    if (!nrow(sub)) next

    top1 <- sub %>%
      filter(GSE == GSE_IDS[1]) %>%
      arrange(desc(abs(score))) %>%
      slice_head(n = TOP_N) %>%
      pull(Pathway)

    top2 <- sub %>%
      filter(GSE == GSE_IDS[2]) %>%
      arrange(desc(abs(score))) %>%
      slice_head(n = TOP_N) %>%
      pull(Pathway)

    top_union <- unique(c(top1, top2))

    mat <- sub %>%
      filter(Pathway %in% top_union) %>%
      select(Pathway, GSE, score) %>%
      pivot_wider(names_from = GSE, values_from = score, values_fill = 0) %>%
      column_to_rownames("Pathway") %>%
      as.matrix()

    if (nrow(mat) <= 1) next

    mat <- mat[, GSE_IDS, drop = FALSE]
    mat <- mat[order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]

    row_labs <- pretty_path(rownames(mat), width = 42)
    lim <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) lim <- 1

    breaks <- seq(-lim, lim, length.out = 101)
    cols   <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)

    n_rows <- nrow(mat)
    png_h  <- max(4, min(18, 3 + 0.22 * n_rows))

    pheatmap(
      mat,
      color        = cols,
      breaks       = breaks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      labels_row   = row_labs,
      labels_col   = GSE_IDS,
      fontsize_row = if (n_rows > 45) 6 else if (n_rows > 30) 7 else 8,
      fontsize_col = 10,
      border_color = NA,
      main         = sprintf("%s | %s - directed score (-log10FDR x signFC)", coll, bin),
      filename     = file.path(CONS, sprintf("heatmap_%s_%s_top%d.png", coll, bin, TOP_N)),
      width        = 9,
      height       = png_h
    )
  }
}

# -----------------------------
# 5. Volcano plots for pathway-level DEA
# -----------------------------

plot_volcano <- function(df, title, outfile,
                         fdr_cut = FDR_CUT,
                         min_lfc = MIN_LFC,
                         n_labels_up = 5,
                         n_labels_down = 5,
                         label_lfc_min = 0.05,
                         wrap_width = 28) {

  df <- df %>%
    mutate(
      sig        = FDR <= fdr_cut & abs(logFC) >= min_lfc,
      mlog10     = -log10(pmax(FDR, 1e-300)),
      rank_score = abs(logFC) * mlog10,
      Path_clean = pretty_path(Pathway, width = wrap_width)
    )

  to_lab_up <- df %>%
    filter(sig, logFC > 0, abs(logFC) >= label_lfc_min) %>%
    arrange(desc(rank_score)) %>%
    slice_head(n = n_labels_up)

  to_lab_down <- df %>%
    filter(sig, logFC < 0, abs(logFC) >= label_lfc_min) %>%
    arrange(desc(rank_score)) %>%
    slice_head(n = n_labels_down)

  p <- ggplot(df, aes(x = logFC, y = mlog10)) +
    geom_point(data = subset(df, !sig), color = "grey75", size = 1, alpha = 0.7) +
    geom_point(data = subset(df, sig & logFC > 0), color = "#d73027", size = 1.6, alpha = 0.9) +
    geom_point(data = subset(df, sig & logFC < 0), color = "#1f78b4", size = 1.6, alpha = 0.9) +
    geom_vline(xintercept = c(-min_lfc, min_lfc), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dotted", color = "grey60") +
    labs(
      title = title,
      x = "logFC (Tumor vs Normal)",
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12) +
    ggrepel::geom_text_repel(
      data = to_lab_up,
      aes(label = Path_clean),
      seed = 1,
      size = 3.2,
      color = "#d73027",
      box.padding = 0.25,
      point.padding = 0.2
    ) +
    ggrepel::geom_text_repel(
      data = to_lab_down,
      aes(label = Path_clean),
      seed = 1,
      size = 3.2,
      color = "#1f78b4",
      box.padding = 0.25,
      point.padding = 0.2
    )

  ggsave(outfile, p, width = 9, height = 6, dpi = 220)

  readr::write_tsv(
    bind_rows(to_lab_up, to_lab_down) %>%
      select(Pathway, Path_clean, logFC, FDR, mlog10, rank_score),
    sub("\\.png$", "_labeled.tsv", outfile)
  )
}

for (coll in COLL) {
  for (bin in BINS) {
    for (gse in GSE_IDS) {
      sub <- all_dea %>% filter(Coll == coll, Bin == bin, GSE == gse)
      if (!nrow(sub)) next

      outfile <- file.path(CONS, sprintf("volcano_%s_%s_%s.png", coll, bin, gse))
      plot_volcano(sub, sprintf("%s | %s | %s", coll, bin, gse), outfile)
    }
  }
}

message("Contrast used: TissueTumor (Tumor - Normal). logFC > 0 indicates higher pathway score in tumor.")
message("Done. Results written to: ", CONS)
