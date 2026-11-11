#########################################
# LUAD age-bin analysis
# Consensus DEG lists and Venn plots
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggvenn)
  library(ggrepel)
})

# This script reads paired limma differential expression results for two GEO
# datasets and identifies concordant DEGs by age bin (<60 and >=60 years).
# Positive logFC corresponds to higher expression in Tumor vs Normal.

# Set the working directory to the folder containing the dea/ directory before
# running this script, or edit ROOT accordingly.
ROOT    <- getwd()
DEA_DIR <- file.path(ROOT, "dea")
FIG_DIR <- file.path(ROOT, "figs")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

GSE_A <- "GSE115002"
GSE_B <- "GSE75037"

BINS       <- c("young", "old")
BIN_LABELS <- c(young = "<60 years", old = ">=60 years")

FDR_CUT <- 0.05
LFC_MIN <- 1
N_LABEL <- 20

as_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x
}

load_dea <- function(gse, bin) {
  cand <- c(
    file.path(DEA_DIR, sprintf("%s___TvN_%s.tsv", gse, bin)),
    file.path(DEA_DIR, sprintf("%s__TvN_%s.tsv",  gse, bin))
  )
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) stop("DEA file not found for ", gse, " / ", bin)

  tt <- readr::read_tsv(f, show_col_types = FALSE)

  if (!"GENE" %in% names(tt)) {
    gcol <- c("Gene", "SYMBOL", "symbol", "gene", "ID")
    gcol <- gcol[gcol %in% names(tt)][1]
    if (is.na(gcol)) stop("No recognizable gene column in: ", f)
    names(tt)[names(tt) == gcol] <- "GENE"
  }

  if (!"logFC" %in% names(tt)) {
    stop("No logFC column in: ", f)
  }

  if (!"FDR" %in% names(tt)) {
    if ("adj.P.Val" %in% names(tt)) {
      tt$FDR <- tt$adj.P.Val
    } else if ("P.Value" %in% names(tt)) {
      tt$FDR <- p.adjust(tt$P.Value, method = "BH")
    } else {
      stop("No FDR, adj.P.Val, or P.Value column in: ", f)
    }
  }

  tt %>%
    mutate(
      GENE = as_chr(GENE),
      FDR = as.numeric(FDR),
      logFC = as.numeric(logFC)
    ) %>%
    filter(!is.na(GENE), is.finite(FDR), is.finite(logFC))
}

make_volcano <- function(tt, gse, bin, bin_label,
                         fdr_cut = FDR_CUT,
                         lfc_min = LFC_MIN,
                         n_label = N_LABEL,
                         out_dir = FIG_DIR) {
  df <- tt %>%
    mutate(
      neg_logFDR = -log10(pmax(FDR, 1e-300)),
      direction = case_when(
        FDR <= fdr_cut & logFC >=  lfc_min ~ "UP",
        FDR <= fdr_cut & logFC <= -lfc_min ~ "DOWN",
        TRUE                               ~ "NS"
      )
    )

  to_lab <- df %>%
    filter(direction != "NS") %>%
    arrange(desc(neg_logFDR)) %>%
    slice_head(n = n_label)

  p <- ggplot(df, aes(x = logFC, y = neg_logFDR)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1.3) +
    geom_vline(xintercept = c(-lfc_min, lfc_min),
               linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_cut),
               linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = to_lab,
      aes(label = GENE),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.2,
      point.padding = 0.15,
      seed = 1
    ) +
    scale_color_manual(
      values = c("UP" = "#d62728", "DOWN" = "#08519c", "NS" = "grey75"),
      name = NULL,
      breaks = c("UP", "DOWN", "NS"),
      labels = c(
        "UP (FDR <= 0.05, logFC >= 1)",
        "DOWN (FDR <= 0.05, logFC <= -1)",
        "Not significant"
      )
    ) +
    labs(
      title = sprintf("%s - Tumor vs Normal (%s)", gse, bin_label),
      subtitle = sprintf("Volcano plot (FDR <= %.2g, |logFC| >= %g)", fdr_cut, lfc_min),
      x = "log2 fold change (Tumor vs Normal)",
      y = expression(-log[10]("FDR"))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      legend.position = "right"
    )

  out_png <- file.path(out_dir, sprintf("Volcano_%s_%s.png", gse, bin))
  out_pdf <- sub("\\.png$", ".pdf", out_png)

  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
  ggsave(out_pdf, p, width = 6, height = 5, device = cairo_pdf)

  message("Volcano saved: ", out_png)
  invisible(p)
}

get_deg_tables <- function(tt, suffix) {
  up <- tt %>%
    filter(FDR <= FDR_CUT, logFC >= LFC_MIN) %>%
    select(gene = GENE, logFC = logFC, FDR = FDR)

  down <- tt %>%
    filter(FDR <= FDR_CUT, logFC <= -LFC_MIN) %>%
    select(gene = GENE, logFC = logFC, FDR = FDR)

  names(up)[names(up) == "logFC"] <- paste0("logFC_", suffix)
  names(up)[names(up) == "FDR"]   <- paste0("FDR_", suffix)
  names(down)[names(down) == "logFC"] <- paste0("logFC_", suffix)
  names(down)[names(down) == "FDR"]   <- paste0("FDR_", suffix)

  list(up = up, down = down)
}

save_venn_plot <- function(sets, bin, direction, bin_label, out_dir = FIG_DIR) {
  if (sum(lengths(sets)) == 0) {
    message("  [", bin, "] No ", direction, " genes detected; Venn plot skipped.")
    return(invisible(NULL))
  }

  fill_col <- if (direction == "UP") "#d62728" else "#08519c"

  p <- ggvenn(
    sets,
    fill_color = c(fill_col, fill_col),
    stroke_size = 0.5,
    set_name_size = 4
  ) +
    ggtitle(sprintf("%sregulated DEGs (%s)",
                    ifelse(direction == "UP", "Up", "Down"), bin_label))

  out_png <- file.path(out_dir, sprintf("Venn_%s_%s.png", bin, direction))
  ggsave(out_png, p, width = 4.5, height = 4.5, dpi = 300)

  message("Venn saved: ", out_png)
  invisible(p)
}

summary_rows <- list()
enriched_out <- list()

for (bin in BINS) {
  bin_label <- BIN_LABELS[[bin]]
  message("\n== Bin: ", bin, " (", bin_label, ") ==")

  tA <- load_dea(GSE_A, bin)
  tB <- load_dea(GSE_B, bin)

  if (bin == BINS[1]) {
    message("DEA columns in ", GSE_A, ":")
    print(colnames(tA))
  }

  make_volcano(tA, GSE_A, bin, bin_label)
  make_volcano(tB, GSE_B, bin, bin_label)

  A <- get_deg_tables(tA, "A")
  B <- get_deg_tables(tB, "B")

  common_up_tab <- inner_join(A$up, B$up, by = "gene") %>%
    mutate(group = "UP", logFC_mean = (logFC_A + logFC_B) / 2)

  common_dn_tab <- inner_join(A$down, B$down, by = "gene") %>%
    mutate(group = "DOWN", logFC_mean = (logFC_A + logFC_B) / 2)

  common_all_tab <- bind_rows(common_up_tab, common_dn_tab)

  message("  Common UP = ", nrow(common_up_tab),
          " | Common DOWN = ", nrow(common_dn_tab))

  venn_up <- setNames(list(A$up$gene, B$up$gene), c(GSE_A, GSE_B))
  venn_dn <- setNames(list(A$down$gene, B$down$gene), c(GSE_A, GSE_B))

  save_venn_plot(venn_up, bin, "UP", bin_label)
  save_venn_plot(venn_dn, bin, "DOWN", bin_label)

  writeLines(common_up_tab$gene,  file.path(FIG_DIR, sprintf("common_%s_UP.txt", bin)))
  writeLines(common_dn_tab$gene,  file.path(FIG_DIR, sprintf("common_%s_DOWN.txt", bin)))
  writeLines(common_all_tab$gene, file.path(FIG_DIR, sprintf("common_%s_ALL.txt", bin)))

  readr::write_tsv(common_up_tab,  file.path(FIG_DIR, sprintf("common_%s_UP.tsv", bin)))
  readr::write_tsv(common_dn_tab,  file.path(FIG_DIR, sprintf("common_%s_DOWN.tsv", bin)))
  readr::write_tsv(common_all_tab, file.path(FIG_DIR, sprintf("common_%s_ALL.tsv", bin)))

  enriched_out[[bin]] <- common_all_tab

  summary_rows[[bin]] <- data.frame(
    bin        = bin,
    up_115002  = nrow(A$up),
    up_75037   = nrow(B$up),
    dn_115002  = nrow(A$down),
    dn_75037   = nrow(B$down),
    common_up  = nrow(common_up_tab),
    common_dn  = nrow(common_dn_tab),
    common_all = nrow(common_all_tab),
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)
readr::write_tsv(summary_df, file.path(FIG_DIR, "venn_summary.tsv"))

if (!is.null(enriched_out[["young"]]) && nrow(enriched_out[["young"]]) > 0) {
  groups_young <- enriched_out[["young"]] %>%
    select(gene, group, logFC_A, logFC_B, logFC_mean, FDR_A, FDR_B)
  readr::write_tsv(groups_young, file.path(FIG_DIR, "groups_young_enriched.tsv"))
}

if (!is.null(enriched_out[["old"]]) && nrow(enriched_out[["old"]]) > 0) {
  groups_old <- enriched_out[["old"]] %>%
    select(gene, group, logFC_A, logFC_B, logFC_mean, FDR_A, FDR_B)
  readr::write_tsv(groups_old, file.path(FIG_DIR, "groups_old_enriched.tsv"))
}

message("\nDone. Output files written to: ", FIG_DIR)
