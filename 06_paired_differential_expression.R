#########################################
# LUAD age-bin transcriptomic analysis
# Paired tumor-normal differential expression by age bin
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(limma)
  library(glue)
})

# Set working directory to the project folder containing prepro/ and meta/.
# Example: setwd("/path/to/Two_Agebins/Series_Matrix")
ROOT <- getwd()
PREPRO_DIR <- file.path(ROOT, "prepro")
META_DIR   <- file.path(ROOT, "meta")
OUT_DIR    <- file.path(ROOT, "dea")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Project directory: ", ROOT)
message("Output directory : ", OUT_DIR)

# Analysis parameters
BINS <- c("young", "old")
FDR_CUT <- 0.05
LFC_MIN <- 0

# -----------------------------
# Helper functions
# -----------------------------

as_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x
}

read_expr <- function(gse) {
  f <- file.path(PREPRO_DIR, paste0(gse, "_expr_SYMBOL.tsv"))
  if (!file.exists(f)) stop("Expression file not found: ", f)

  df <- readr::read_tsv(f, show_col_types = FALSE)
  if (!"GENE" %in% colnames(df)) stop("Missing GENE column in: ", f)

  rn <- df$GENE
  M <- as.matrix(data.frame(row.names = rn, df[, -1, drop = FALSE], check.names = FALSE))
  storage.mode(M) <- "double"
  M
}

assert_pairs_exact <- function(mm_bin, gse, bin_label) {
  chk <- mm_bin %>%
    count(patient_id, tissue) %>%
    tidyr::pivot_wider(names_from = tissue, values_from = n, values_fill = 0)

  if (!"Normal" %in% names(chk)) chk$Normal <- 0L
  if (!"Tumor"  %in% names(chk)) chk$Tumor  <- 0L

  chk <- chk %>% mutate(ok = (Normal == 1L & Tumor == 1L))

  if (!all(chk$ok)) {
    bad <- chk %>%
      filter(!ok) %>%
      transmute(
        patient_id,
        Normal = as.integer(Normal),
        Tumor = as.integer(Tumor)
      )

    msg <- paste0(
      "\n[", gse, " / ", bin_label, "] Non-strict pairing detected. ",
      "Each patient must have exactly 1 Normal and 1 Tumor sample:\n",
      paste(apply(bad, 1, function(r) {
        sprintf("  - %s (N=%s, T=%s)", r[["patient_id"]], r[["Normal"]], r[["Tumor"]])
      }), collapse = "\n"),
      "\nFix the metadata or exclude these patients before rerunning."
    )
    stop(msg, call. = FALSE)
  }

  invisible(TRUE)
}

# -----------------------------
# Paired DEA within one age bin
# -----------------------------

run_dea_bin <- function(E, meta_gse, gse, bin_label, fdr_cut = FDR_CUT, lfc_min = LFC_MIN) {
  # Keep only samples from the selected age bin with valid patient/tissue labels.
  mm0 <- meta_gse %>%
    transmute(
      sample_id  = as_chr(sample_id),
      patient_id = as_chr(patient_id),
      tissue = case_when(
        grepl("^normal$", tissue, ignore.case = TRUE) ~ "Normal",
        grepl("^tumor$",  tissue, ignore.case = TRUE) ~ "Tumor",
        TRUE ~ as_chr(tissue)
      ),
      age_bin = as_chr(age_bin)
    ) %>%
    filter(
      !is.na(sample_id),
      !is.na(patient_id),
      tissue %in% c("Normal", "Tumor"),
      age_bin == bin_label
    ) %>%
    distinct(sample_id, .keep_all = TRUE)

  if (nrow(mm0) == 0) {
    message(glue("  [{bin_label}] no samples found in metadata; skipped."))
    stop("El contraste está incompleto. Comprueba los identificadores antes de continuar.")
  }

  # Intersect expression and metadata sample IDs, then force the same order.
  if (any(duplicated(colnames(E)))) {
    E <- E[, !duplicated(colnames(E)), drop = FALSE]
  }

  ids <- intersect(colnames(E), mm0$sample_id)
  if (length(ids) < 4) {
    message(glue("  [{bin_label}] fewer than 4 matched expression columns; skipped."))
    stop("El contraste está incompleto. Comprueba los identificadores antes de continuar.")
  }

  Emm <- E[, ids, drop = FALSE]
  mm  <- mm0[match(ids, mm0$sample_id), , drop = FALSE]
  stopifnot(all(mm$sample_id == colnames(Emm)))

  # Strict 1N + 1T pairing is required within each age bin.
  assert_pairs_exact(mm, gse, bin_label)

  # Remove genes with zero variance.
  keep_var <- apply(Emm, 1, function(v) {
    v <- v[is.finite(v)]
    length(v) > 1 && stats::var(v) > 0
  })
  Emm <- Emm[keep_var, , drop = FALSE]

  if (nrow(Emm) < 50) {
    message(glue("  [{bin_label}] fewer than 50 genes remained after variance filtering; skipped."))
    stop("El contraste está incompleto. Comprueba los identificadores antes de continuar.")
  }

  # Paired design with intercept.
  # Coefficient TissueTumor corresponds to Tumor - Normal.
  mm$pair_id <- droplevels(factor(mm$patient_id))
  mm$Tissue  <- droplevels(factor(mm$tissue, levels = c("Normal", "Tumor")))
  design <- model.matrix(~ pair_id + Tissue, data = mm)
  stopifnot(nrow(design) == ncol(Emm))

  if (!"TissueTumor" %in% colnames(design)) {
    stop("Coefficient TissueTumor was not found in the design matrix.")
  }

  fit  <- limma::lmFit(Emm, design)
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  message("  Contrast: TissueTumor (Tumor - Normal); logFC > 0 indicates higher expression in Tumor")
  tt <- limma::topTable(fit2, coef = "TissueTumor", number = Inf, sort.by = "P")

  if (!"AveExpr" %in% names(tt)) {
    tt$AveExpr <- rowMeans(Emm, na.rm = TRUE)
  }

  if (!"FDR" %in% names(tt)) {
    if ("adj.P.Val" %in% names(tt)) {
      tt$FDR <- tt$adj.P.Val
    } else if ("P.Value" %in% names(tt)) {
      tt$FDR <- p.adjust(tt$P.Value, method = "BH")
    } else {
      stop("P.Value/adj.P.Val columns not found; cannot compute FDR.")
    }
  }

  tt$SE  <- with(tt, ifelse(is.finite(t) & t != 0, abs(logFC / t), NA_real_))
  tt$dir <- ifelse(tt$logFC > 0, "up_in_Tumor", "down_in_Tumor")

  out_file <- file.path(OUT_DIR, sprintf("%s__TvN_%s.tsv", gse, bin_label))
  readr::write_tsv(tibble::rownames_to_column(tt, "GENE"), out_file)

  n_pairs <- length(unique(mm$patient_id))
  n_sig   <- sum(tt$FDR <= fdr_cut & abs(tt$logFC) >= lfc_min, na.rm = TRUE)
  n_up    <- sum(tt$FDR <= fdr_cut &  tt$logFC >=  lfc_min, na.rm = TRUE)
  n_down  <- sum(tt$FDR <= fdr_cut &  tt$logFC <= -lfc_min, na.rm = TRUE)

  message(glue("  [{bin_label}] pairs = {n_pairs}, DEGs = {n_sig}; saved: {out_file}"))

  tibble(
    gse = gse,
    bin = bin_label,
    n_pairs = n_pairs,
    n_sig = n_sig,
    n_up = n_up,
    n_down = n_down,
    file = out_file
  )
}

# -----------------------------
# Optional sanity check
# -----------------------------

sanity_check_dea_direction <- function(gse, bin_label) {
  message("\nRunning optional sanity check for ", gse, " / ", bin_label)

  tt <- readr::read_tsv(
    file.path(OUT_DIR, sprintf("%s__TvN_%s.tsv", gse, bin_label)),
    show_col_types = FALSE
  )
  E <- read_expr(gse)
  M <- readr::read_csv(file.path(META_DIR, paste0(gse, "_meta.csv")), show_col_types = FALSE) %>%
    mutate(
      sample_id  = as.character(sample_id),
      patient_id = as.character(patient_id),
      tissue = case_when(
        grepl("^normal$", tissue, ignore.case = TRUE) ~ "Normal",
        grepl("^tumor$",  tissue, ignore.case = TRUE) ~ "Tumor",
        TRUE ~ as.character(tissue)
      ),
      age_bin = as.character(age_bin)
    ) %>%
    filter(age_bin == bin_label, tissue %in% c("Normal", "Tumor")) %>%
    distinct(sample_id, .keep_all = TRUE)

  ids <- intersect(colnames(E), M$sample_id)
  mm <- M[match(ids, M$sample_id), , drop = FALSE]

  ct <- mm %>%
    count(patient_id, tissue) %>%
    tidyr::pivot_wider(names_from = tissue, values_from = n, values_fill = 0)

  keep_pat <- with(ct, patient_id[Normal == 1 & Tumor == 1])
  mm2 <- mm %>%
    filter(patient_id %in% keep_pat) %>%
    mutate(Tissue = factor(tissue, levels = c("Normal", "Tumor"))) %>%
    arrange(patient_id, Tissue)

  X <- E[, mm2$sample_id, drop = FALSE]
  stopifnot(ncol(X) %% 2 == 0)

  D <- X[, seq(2, ncol(X), by = 2), drop = FALSE] -
       X[, seq(1, ncol(X), by = 2), drop = FALSE]
  mean_diff <- rowMeans(D, na.rm = TRUE)

  common <- intersect(names(mean_diff), tt$GENE)
  r <- cor(mean_diff[common], tt$logFC[match(common, tt$GENE)], use = "complete.obs")

  message("Correlation between paired mean Tumor-Normal difference and limma logFC: ", round(r, 4))
  invisible(r)
}

# -----------------------------
# Main loop
# -----------------------------

gse_expr <- sub("_expr_SYMBOL\\.tsv$", "", list.files(PREPRO_DIR, pattern = "_expr_SYMBOL\\.tsv$"))
gse_meta <- sub("_meta\\.csv$", "", list.files(META_DIR, pattern = "_meta\\.csv$"))
GSE_IDS  <- intersect(gse_expr, gse_meta)

if (!length(GSE_IDS)) {
  stop("No matching expression and metadata files were found in prepro/ and meta/.")
}

SUMMARY <- list()
for (gse in GSE_IDS) {
  message("\n== ", gse, " ==")
  E <- read_expr(gse)
  M <- readr::read_csv(file.path(META_DIR, paste0(gse, "_meta.csv")), show_col_types = FALSE)

  bin_results <- lapply(BINS, function(bin_label) {
    run_dea_bin(E, M, gse, bin_label)
  })

  SUMMARY[[gse]] <- dplyr::bind_rows(bin_results)
}

RES <- dplyr::bind_rows(SUMMARY)
out_summary <- file.path(OUT_DIR, "dea_pairs_summary.tsv")
readr::write_tsv(RES, out_summary)

message("\nSummary saved: ", out_summary)
print(RES)
message("\nPaired differential expression analysis complete.")

# To inspect directionality for one result, run for example:
# sanity_check_dea_direction("GSE115002", "young")
