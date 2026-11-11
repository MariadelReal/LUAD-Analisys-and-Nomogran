#########################################
# LUAD age-bin transcriptomic workflow
# 02_initial_expression_qc.R
#
# Purpose:
#   Load GEO ExpressionSet objects for the selected GSE datasets and generate
#   initial QC outputs: density plots, sample-median histograms, boxplots, PCA,
#   data-processing summaries, and GSM audit tables.
#
# Expected local input:
#   GSE115002_series.tsv
#   GSE75037_series.tsv
#
# Main outputs:
#   qc0/*_precheck.txt
#   qc0/*_densities.png
#   qc0/*_medians_hist.png
#   qc0/*_boxplot_colored.png
#   qc0/*_PCA.png / .pdf
#   qc0/AUD_*
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase)
  library(ggplot2);  library(dplyr);  library(tidyr); library(readr)
  library(stringr);  library(limma)
})

# ---- Analysis directory ----
# Run this script from the folder containing the GEO Series Matrix TSV files
# produced in step 01_download_geo_series_matrix.R, or replace getwd() with that path.
ROOT <- getwd()
QC_DIR <- file.path(ROOT, "qc0")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------- GEO datasets ----------
# Local *_series.tsv files are used only to detect which GSE IDs to process.
# If they are not present yet, the default study list is used.
series_files <- list.files(ROOT, pattern = "^GSE\\d+_series\\.tsv$", full.names = FALSE)
GSE_IDS_LOCAL <- unique(sub("_series\\.tsv$", "", series_files))
GSE_IDS <- if (length(GSE_IDS_LOCAL)) GSE_IDS_LOCAL else c("GSE115002", "GSE75037")
rm(GSE_IDS_LOCAL, series_files)

# ---------- Plot colors ----------
BASE_PALETTE <- c("#1f78b4","#e31a1c","#33a02c","#ff7f00","#6a3d9a","#b15928","#a6cee3","#fb9a99")
GSE_BASE_COL <- setNames(BASE_PALETTE[seq_along(GSE_IDS)], GSE_IDS)

# ---------- Helpers ----------
is_log2_like <- function(M){
  rng <- range(M, na.rm = TRUE)
  (rng[2] <= 20 && rng[1] >= -5)
}

ks_summary <- function(M, max_samples = 60, per_pair_n = 1000, max_pairs = 300){
  M0 <- sweep(M, 2, apply(M, 2, median, na.rm = TRUE), FUN = "-")
  smp <- colnames(M0)
  if (length(smp) > max_samples) smp <- sample(smp, max_samples)
  if (length(smp) < 2) return(c(median = NA_real_, q90 = NA_real_))
  pairs <- t(combn(smp, 2))
  if (nrow(pairs) > max_pairs) pairs <- pairs[sample(seq_len(nrow(pairs)), max_pairs), , drop = FALSE]
  ks_vals <- numeric(nrow(pairs))
  for (k in seq_len(nrow(pairs))){
    a <- M0[, pairs[k,1]]; b <- M0[, pairs[k,2]]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) > per_pair_n) a <- sample(a, per_pair_n)
    if (length(b) > per_pair_n) b <- sample(b, per_pair_n)
    ks_vals[k] <- suppressWarnings(ks.test(a, b)$statistic)
  }
  c(median = stats::median(ks_vals, na.rm = TRUE),
    q90    = stats::quantile(ks_vals, 0.90, na.rm = TRUE, names = FALSE))
}

save_density_plot <- function(M, tag, out){
  keep <- colnames(M); if (length(keep) > 40) keep <- sample(keep, 40)
  df <- as.data.frame(M[, keep, drop = FALSE]) |>
    pivot_longer(everything(), names_to = "sample", values_to = "expr")
  p <- ggplot(df, aes(x = expr, group = sample)) +
    geom_density(alpha = 0.25) +
    labs(title = paste(tag, "— Densities (sample subset)"), x = "Expression", y = "Density") +
    theme_bw(base_size = 10)
  ggsave(out, p, width = 7.5, height = 4.8, dpi = 150)
}

save_medians_hist <- function(M, tag, out){
  meds <- apply(M, 2, median, na.rm = TRUE)
  df <- data.frame(median = meds)
  p <- ggplot(df, aes(x = median)) +
    geom_histogram(bins = 40) +
    labs(title = paste(tag, "— Histogram of sample medians"), x = "Median", y = "Frequency") +
    theme_bw(base_size = 10)
  ggsave(out, p, width = 6.2, height = 4.2, dpi = 150)
}

# ---- Per-block boxplots (GSE_GPL) ----
save_boxplot_colored <- function(M, tag, out, fill_color, ylim = NULL){
  df <- as.data.frame(M)
  df$FEATURE <- rownames(df)
  long <- df |>
    pivot_longer(-FEATURE, names_to = "GSM", values_to = "expr") |>
    mutate(GSM = factor(GSM, levels = colnames(M)))
  p <- ggplot(long, aes(x = GSM, y = expr)) +
    geom_boxplot(fill = fill_color, colour = "black", outlier.size = 0.35) +
    labs(title = paste0(tag, " | boxplot log2"), x = NULL, y = "Expression") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          panel.grid.minor = element_blank())
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  w <- max(1800, 14 * ncol(M)); h <- 900
  ggsave(out, p, width = w/150, height = h/150, dpi = 150, limitsize = FALSE)
}

# ---- Load all ExpressionSets for each GSE ----
load_esets_all <- function(gse){
  esets <- getGEO(gse, GSEMatrix = TRUE)
  if (!is.list(esets)) esets <- list(esets)
  names(esets) <- vapply(esets, function(e){
    gpl <- tryCatch(annotation(e), error = function(...) NA_character_)
    paste0(gse, "_", ifelse(is.na(gpl), "NA", gpl))
  }, character(1))
  esets
}

# ---- SYMBOL helpers for optional joint PCA ----
guess_symbol_column <- function(fdat){
  pri <- c("Gene.symbol","Gene Symbol","SYMBOL","Symbol","GENE_SYMBOL","gene_symbol",
           "GeneName","GENE","GENE_SYMBOLS")
  for (nm in pri) if (nm %in% colnames(fdat)) return(nm)
  alt <- grep("symbol|gene", colnames(fdat), ignore.case = TRUE, value = TRUE)
  if (length(alt)) return(alt[1])
  NA_character_
}
clean_symbols <- function(x){
  x <- as.character(x); x[is.na(x)] <- NA_character_
  x <- sapply(x, function(s){
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    s <- strsplit(s, "\\s*[/|;,+]+\\s*")[[1]][1]
    s <- toupper(trimws(s))
    if (s %in% c("", "NA", "---", "?", "N/A")) return(NA_character_)
    s
  }, USE.NAMES = FALSE)
  unname(x)
}
symbol_matrix <- function(eset){
  M  <- exprs(eset)
  if (!is_log2_like(M)) M <- log2(M + 1)
  fd <- Biobase::fData(eset)
  sym_col <- guess_symbol_column(fd)
  if (is.na(sym_col)) stop("No symbol column found in fData().")
  syms <- clean_symbols(fd[[sym_col]])
  keep <- !is.na(syms) & nzchar(syms)
  limma::avereps(M[keep, , drop = FALSE], ID = syms[keep])
}

# ---- One-panel boxplot across all blocks ----
save_box_onepanel <- function(M_list, out_png, out_pdf, fill_map){
  all_stats <- list()
  sample_order <- character(0)
  
  for (tag in names(M_list)){
    M <- M_list[[tag]]
    sample_order <- c(sample_order, colnames(M))
    qs <- apply(M, 2, function(x){
      x <- x[is.finite(x)]
      quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE, type = 7)
    })
    stats <- data.frame(
      sample = colnames(M),
      block  = tag,
      ymin   = qs[1,],
      lower  = qs[2,],
      middle = qs[3,],
      upper  = qs[4,],
      ymax   = qs[5,],
      stringsAsFactors = FALSE
    )
    all_stats[[tag]] <- stats
  }
  
  ST <- do.call(rbind, all_stats)
  ST$sample <- factor(ST$sample, levels = sample_order)
  ST$block  <- factor(ST$block,  levels = names(M_list))
  
  total_n  <- length(sample_order)
  width_in <- max(12, 0.09 * total_n)
  height_in <- 8
  
  p <- ggplot(ST, aes(x = sample, fill = block)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity", colour = "black", width = 0.7, size = 0.3
    ) +
    scale_fill_manual(values = fill_map, guide = guide_legend(title = "Block")) +
    labs(title = "Sample boxplots (one panel; color = GSE_GPL block)", x = NULL, y = "Expression") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
          panel.grid.minor = element_blank())
  
  ggsave(out_png, p, width = width_in, height = height_in, units = "in",
         dpi = 300, limitsize = FALSE)
  ggsave(out_pdf, p, width = width_in, height = height_in, units = "in", device = cairo_pdf)
}

# ---- PCA helpers ----
pca_of_matrix <- function(M){
  pc <- prcomp(t(M), scale. = TRUE)
  var_expl <- (pc$sdev^2) / sum(pc$sdev^2)
  list(scores = pc$x[,1:2, drop = FALSE], var = var_expl[1:2])
}
save_pca_points <- function(scores, var, labels, colors, title, out_png, out_pdf){
  df <- data.frame(PC1 = scores[,1], PC2 = scores[,2], group = labels)
  width_in <- 7.5; height_in <- 6.5
  p <- ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 1.8, alpha = 0.9) +
    scale_color_manual(values = colors) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", 100*var[1]),
         y = sprintf("PC2 (%.1f%%)", 100*var[2]),
         color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "right")
  ggsave(out_png, p, width = width_in, height = height_in, units = "in", dpi = 300)
  ggsave(out_pdf, p, width = width_in, height = height_in, units = "in", device = cairo_pdf)
}

# ---------- MAIN ----------
mat_list    <- list()
global_rng  <- c(Inf, -Inf)
dp_summary  <- list()
block_color <- c()

for (gse in GSE_IDS){
  message("\n== QC0 ", gse, " ==")
  
  # Carga TODOS los ExpressionSet (multi-GPL)
  esets <- tryCatch(load_esets_all(gse), error = function(e) NULL)
  if (is.null(esets)) { message("  Could not load ", gse); next }
  
  # GSE-level data_processing metadata, when available
  gse_obj <- tryCatch(GEOquery::getGEO(gse, GSEMatrix = FALSE), error = function(e) NULL)
  dp_lines <- character(0)
  if (!is.null(gse_obj)) {
    m <- Meta(gse_obj)
    if (!is.null(m$data_processing)) dp_lines <- c(dp_lines, m$data_processing)
  }
  txt_path <- file.path(QC_DIR, paste0(gse, "_data_processing.txt"))
  .dp_keywords <- c("rma","quantile","neqc","loess","vsn","lumi","background","log2","gc-rma","gcrma","robust multiarray")
  .detect_methods <- function(txt){
    if (is.na(txt) || !nzchar(txt)) return(character(0))
    found <- unique(.dp_keywords[grepl(paste0("\\b", .dp_keywords, "\\b"), txt, ignore.case = TRUE)])
    pretty <- c("rma"="RMA","gcrma"="GCRMA","gc-rma"="GCRMA",
                "quantile"="Quantile","neqc"="NEQC","loess"="LOESS",
                "vsn"="VSN","lumi"="lumi","background"="Background correction","log2"="log2")
    unname(ifelse(found %in% names(pretty), pretty[found], tools::toTitleCase(found)))
  }
  if (length(dp_lines)) {
    writeLines(dp_lines, txt_path)
    dp_txt_lower  <- tolower(paste(dp_lines, collapse = " "))
    methods_found <- .detect_methods(dp_txt_lower)
    dp_summary[[length(dp_summary)+1]] <- data.frame(
      GSE        = gse, has_report = TRUE,
      methods    = if (length(methods_found)) paste(methods_found, collapse = "; ") else "",
      n_lines    = length(dp_lines),
      snippet    = substr(gsub("\\s+", " ", dp_txt_lower), 1, 240),
      stringsAsFactors = FALSE
    )
  } else {
    writeLines("(no data_processing field found for this GSE)", txt_path)
    dp_summary[[length(dp_summary)+1]] <- data.frame(
      GSE = gse, has_report = FALSE, methods = "", n_lines = 0L, snippet = "",
      stringsAsFactors = FALSE
    )
  }
  
  # Process each GSE_GPL block
  for (es_nm in names(esets)){
    eset <- esets[[es_nm]]
    
    # Save sample-level data_processing columns when present
    dp_cols <- grep("data\\s*processing", names(pData(eset)), ignore.case = TRUE, value = TRUE)
    if (length(dp_cols)) {
      readr::write_tsv(
        data.frame(GSM = rownames(pData(eset)),
                   pData(eset)[, dp_cols, drop = FALSE],
                   check.names = FALSE),
        file.path(QC_DIR, paste0(es_nm, "_data_processing_pheno.tsv"))
      )
    }
    
    M <- exprs(eset)
    if (anyNA(M)) {
      keep <- rowSums(is.finite(M)) > 0
      M <- M[keep, , drop = FALSE]
    }
    mat_list[[es_nm]] <- M
    
    global_rng[1] <- min(global_rng[1], min(M, na.rm = TRUE))
    global_rng[2] <- max(global_rng[2], max(M, na.rm = TRUE))
    
    # Block color follows the GSE-level color
    gse_base <- sub("_.*$", "", es_nm)
    col_blk  <- GSE_BASE_COL[[gse_base]]
    if (is.null(col_blk)) col_blk <- "#999999"
    block_color[es_nm] <- col_blk
    
    # Precheck text report
    gpl <- tryCatch(annotation(eset), error = function(e) NA_character_)
    rng <- range(M, na.rm = TRUE)
    log2_flag <- is_log2_like(M)
    ks <- ks_summary(M)
    meds <- apply(M, 2, median, na.rm = TRUE)
    sd_meds <- stats::sd(meds); rng_meds <- diff(range(meds))
    txt <- file.path(QC_DIR, paste0(es_nm, "_precheck.txt"))
    cat(
      sprintf("BLOCK: %s\n", es_nm),
      sprintf("GSE: %s   GPL: %s\n", gse, ifelse(is.na(gpl), "(NA)", gpl)),
      sprintf("Dimensions: %d features × %d samples\n", nrow(M), ncol(M)),
      sprintf("Range: [%.3f, %.3f]\n", rng[1], rng[2]),
      sprintf("Log2-like scale?: %s\n", ifelse(log2_flag, "Yes", "No/undetermined")),
      sprintf("KS median: %s\n", ifelse(is.na(ks["median"]), "NA", sprintf("%.3f", ks["median"]))),
      sprintf("KS p90   : %s\n", ifelse(is.na(ks["q90"]), "NA", sprintf("%.3f", ks["q90"]))),
      sprintf("sd(sample medians): %.3f\n", sd_meds),
      sprintf("range(sample medians): %.3f\n", rng_meds),
      file = txt, sep = ""
    )
    
    # Per-block figures
    save_density_plot(M, es_nm, file.path(QC_DIR, paste0(es_nm, "_densities.png")))
    save_medians_hist(M, es_nm, file.path(QC_DIR, paste0(es_nm, "_medians_hist.png")))
    save_boxplot_colored(M, es_nm,
                         out = file.path(QC_DIR, paste0(es_nm, "_boxplot_colored.png")),
                         fill_color = block_color[es_nm], ylim = NULL)
    
    sc <- try(pca_of_matrix(if (!is_log2_like(M)) log2(M + 1) else M), silent = TRUE)
    if (!inherits(sc, "try-error")){
      save_pca_points(
        scores = sc$scores, var = sc$var,
        labels = rep(es_nm, nrow(sc$scores)),
        colors = setNames(block_color[es_nm], es_nm),
        title  = paste0(es_nm, " | PCA"),
        out_png = file.path(QC_DIR, paste0(es_nm, "_PCA.png")),
        out_pdf = file.path(QC_DIR, paste0(es_nm, "_PCA.pdf"))
      )
    }
    
    message("  -> ", es_nm, " precheck and figures generated")
  }
}

# ---------- Joint PCA using intersected probe IDs ----------
if (length(mat_list) >= 2){
  common <- Reduce(intersect, lapply(mat_list, rownames))
  if (length(common) >= 1000){
    Mj <- do.call(cbind, lapply(names(mat_list), function(nm) mat_list[[nm]][common, , drop = FALSE]))
    keep_var <- apply(Mj, 1, function(v){ v <- v[is.finite(v)]; if (!length(v)) return(FALSE); stats::var(v) > 0 })
    Mj <- Mj[keep_var, , drop = FALSE]
    if (nrow(Mj) > 30000) { set.seed(1L); Mj <- Mj[sample(seq_len(nrow(Mj)), 30000), , drop = FALSE] }
    
    pc <- try(prcomp(t(Mj), scale. = TRUE), silent = TRUE)
    if (!inherits(pc, "try-error")){
      var_expl <- (pc$sdev^2) / sum(pc$sdev^2)
      labs <- unlist(lapply(names(mat_list), function(nm) rep(nm, ncol(mat_list[[nm]][common, , drop = FALSE]))))
      dfp  <- data.frame(sample = colnames(Mj),
                         block  = factor(labs, levels = names(mat_list)),
                         PC1 = pc$x[,1], PC2 = pc$x[,2])
      
      png(file.path(QC_DIR, paste0("PCA_joint_", paste(names(mat_list), collapse = "_"), ".png")),
          width = 1300, height = 1000, res = 140)
      op <- par(mar = c(4,4,3,1))
      cols_vec <- block_color[levels(dfp$block)]
      plot(dfp$PC1, dfp$PC2, pch = 19, col = cols_vec[dfp$block],
           xlab = sprintf("PC1 (%.1f%%)", 100*var_expl[1]),
           ylab = sprintf("PC2 (%.1f%%)", 100*var_expl[2]),
           main = "Joint PCA (intersected probe IDs)")
      grid(col = "grey85"); box()
      legend("topright", inset = 0.01, legend = levels(dfp$block),
             col = cols_vec, pch = 19, bty = "n", cex = 0.9)
      par(op); dev.off()
      
      write_tsv(dfp, file.path(QC_DIR, "PCA_joint_coords.tsv"))
      write_tsv(
        tibble(PC = paste0("PC", seq_along(var_expl)), var_explained = var_expl),
        file.path(QC_DIR, "PCA_joint_var_explained.tsv")
      )
      message("  -> joint PCA generated")
    }
  } else {
    message("  (Skipped joint PCA: probe intersection < 1000)")
  }
}

# ---- Reported normalization summary at GSE level ----
if (length(dp_summary)) {
  DP_SUM <- dplyr::bind_rows(dp_summary)
  out_sum <- file.path(QC_DIR, "data_processing_summary.tsv")
  readr::write_tsv(DP_SUM, out_sum)
  print(DP_SUM[, c("GSE","has_report","methods","n_lines")], row.names = FALSE)
  message("data_processing summary written to ", out_sum)
} else {
  message("No data_processing summary was generated because no valid GSE was processed.")
}

# =============================================================
# Sample audit
# =============================================================
aud <- lapply(names(mat_list), function(nm){
  data.frame(block = nm, GSM = colnames(mat_list[[nm]]), stringsAsFactors = FALSE)
})
AUD <- dplyr::bind_rows(aud)
if (nrow(AUD)){
  readr::write_tsv(AUD, file.path(QC_DIR, "AUD_samples_by_block.tsv"))
  AUD$GSE <- sub("(_GPL.*)$", "", AUD$block)
  counts <- AUD |> dplyr::count(GSE, name = "n_in_qc")
  readr::write_tsv(counts, file.path(QC_DIR, "AUD_counts_by_GSE.tsv"))
}

all_pd <- list()
for (gse in GSE_IDS){
  esets0 <- tryCatch(getGEO(gse, GSEMatrix = TRUE), error = function(e) NULL)
  if (is.null(esets0)) next
  if (!is.list(esets0)) esets0 <- list(esets0)
  gsm_all <- unique(unlist(lapply(esets0, function(e) rownames(pData(e)))))
  all_pd[[gse]] <- data.frame(GSE = gse, GSM = gsm_all, stringsAsFactors = FALSE)
}
PD <- if (length(all_pd)) dplyr::bind_rows(all_pd) else data.frame(GSE=character(0), GSM=character(0))

missing <- if (nrow(PD) && nrow(AUD)) dplyr::anti_join(PD, AUD, by = c("GSE","GSM")) else PD
readr::write_tsv(missing, file.path(QC_DIR, "AUD_missing_samples.tsv"))
if (nrow(missing)){
  message("WARNING: some GSM accessions are present in GEO but absent from the QC matrices. Check AUD_missing_samples.tsv")
} else {
  message("Sample audit passed: all GEO GSM accessions are represented in at least one GSE_GPL block.")
}

# ---------- Global one-panel boxplot ----------
if (length(mat_list) >= 1){
  save_box_onepanel(
    M_list = mat_list,
    out_png = file.path(QC_DIR, "BOX_onepanel.png"),
    out_pdf = file.path(QC_DIR, "BOX_onepanel.pdf"),
    fill_map = block_color
  )
  message("  -> BOX_onepanel.png / .pdf")
}

message("\nDone. QC PNG/TSV outputs written to: ", QC_DIR)
# =============================================================
