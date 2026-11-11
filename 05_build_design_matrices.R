if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
## 05_build_design_matrices.R
## Purpose: build gene-signature and clinical design matrices for survival modeling.
## Inputs:
##   data/clin_survival_cohort.rds
##   data/expr_tumor_SYMBOL_filtered.rds
##   sig/CellCycle_LUAD_v1.csv
## Outputs:
##   data/05_preliminar/X_genes_signature_raw.rds
##   data/05_preliminar/X_genes_signature_scaled.rds
##   data/05_preliminar/X_clinic_noIntercept.rds
##   data/05_preliminar/clinical_for_modeling.rds
##   tables/05_clinic_design_columns.tsv
## Notes:
##   Age is kept in years and is not scaled.
##   Gene expression values are z-score scaled by gene for penalized Cox modeling.


set.seed(123L)

# ---------- Project directory ----------
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# ---------- Packages ----------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

# ---------- Paths ----------
data_dir <- file.path(base_dir, "data")
sig_dir  <- file.path(base_dir, "sig")
tab_dir  <- file.path(base_dir, "tables")
log_dir  <- file.path(base_dir, "logs")

dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "05_build_design_matrices.txt")
cat("# 05_build_design_matrices\n", file = log_file)

# ---------- Utilities ----------
norm_bar <- function(x) {
  x <- toupper(as.character(x))
  gsub("[.]", "-", x)
}

# ---------- Load final survival cohort and tumor expression matrix ----------
clin_path <- file.path(data_dir, "clin_survival_cohort.rds")
expr_path <- file.path(data_dir, "expr_tumor_SYMBOL_filtered.rds")

if (!file.exists(clin_path)) stop("Missing file: ", clin_path)
if (!file.exists(expr_path)) stop("Missing file: ", expr_path)

clin <- readRDS(clin_path)
expr <- readRDS(expr_path)  # first column = gene; remaining columns = barcodes

stopifnot("gene" %in% names(expr))
stopifnot("barcode" %in% names(clin))

# ---------- Align expression and clinical data ----------
all_bars_expr <- setdiff(colnames(expr), "gene")
colnames(expr)[match(all_bars_expr, colnames(expr))] <- norm_bar(all_bars_expr)

clin$barcode    <- norm_bar(clin$barcode)
clin$patient_id <- toupper(substr(clin$barcode, 1, 12))

cat("=== Expression-clinical alignment diagnostic ===\n", file = log_file, append = TRUE)
cat(sprintf("Expression barcodes: %d | clinical rows: %d\n",
            length(setdiff(colnames(expr), "gene")), nrow(clin)),
    file = log_file, append = TRUE)

expr_barcodes <- setdiff(colnames(expr), "gene")
n_match_bar <- sum(expr_barcodes %in% clin$barcode)
n_match_pid <- sum(substr(expr_barcodes, 1, 12) %in% clin$patient_id)

cat(sprintf("Match by barcode: %d | match by patient_id(12): %d\n",
            n_match_bar, n_match_pid),
    file = log_file, append = TRUE)

if (n_match_bar == 0 && n_match_pid > 0) {
  pid_from_expr <- substr(expr_barcodes, 1, 12)
  clin_rebuilt <- clin[match(pid_from_expr, clin$patient_id), , drop = FALSE]
  clin_rebuilt$barcode <- expr_barcodes
  keep <- !is.na(clin_rebuilt$patient_id)
  clin <- clin_rebuilt[keep, , drop = FALSE]
  expr <- expr[, c("gene", clin$barcode), drop = FALSE]
} else {
  bars_common <- intersect(expr_barcodes, clin$barcode)
  if (length(bars_common) == 0) {
    stop("No overlap between expression barcodes and clinical barcodes/patient IDs.")
  }
  expr <- expr[, c("gene", bars_common), drop = FALSE]
  clin <- clin[match(bars_common, clin$barcode), , drop = FALSE]
}

stopifnot(all(clin$barcode == setdiff(colnames(expr), "gene")))
cat(sprintf("OK: expression and clinical data aligned 1:1 by barcode. n = %d\n", nrow(clin)),
    file = log_file, append = TRUE)

# ---------- Load gene signature ----------
sig_file <- file.path(sig_dir, "CellCycle_LUAD_v1.csv")
if (!file.exists(sig_file)) stop("Missing signature file: ", sig_file)

sig <- readr::read_csv(sig_file, show_col_types = FALSE)
colnames(sig) <- tolower(colnames(sig))
stopifnot("gene_symbol" %in% names(sig))
genes_signature <- unique(sig$gene_symbol)

# ---------- Subset expression matrix to signature genes ----------
bars <- setdiff(colnames(expr), "gene")
present <- intersect(genes_signature, expr$gene)
missing <- setdiff(genes_signature, present)

cat(sprintf("Signature genes: total=%d | present=%d | missing=%d\n",
            length(genes_signature), length(present), length(missing)),
    file = log_file, append = TRUE)

if (length(missing)) {
  cat("Missing signature genes in expression matrix: ",
      paste(missing, collapse = ", "), "\n",
      file = log_file, append = TRUE)
}

if (length(present) == 0) {
  stop("No signature genes were found in the expression matrix.")
}

Xg <- t(as.matrix(expr[match(present, expr$gene), bars, drop = FALSE]))
colnames(Xg) <- present
stopifnot(!is.null(rownames(Xg)))

# ---------- Clinical modeling variables ----------
clin_mod <- clin %>%
  mutate(
    age     = as.numeric(age_years),
    gender  = factor(gender, levels = c("Female", "Male")),
    stage   = factor(stage_use, levels = c("I", "II", "III", "IV")),
    smoking = factor(smoking_status,
                     levels = c("Never", "Former", "Current", "Unknown"))
  )
clin_mod <- as.data.frame(clin_mod)
rownames(clin_mod) <- clin_mod$barcode

# ---------- Final row intersection ----------
common <- intersect(rownames(Xg), rownames(clin_mod))
Xg       <- Xg[common, , drop = FALSE]
clin_mod <- clin_mod[common, , drop = FALSE]

stopifnot(nrow(Xg) == nrow(clin_mod))
stopifnot(all(rownames(Xg) == rownames(clin_mod)))

# ---------- Clinical design matrix without intercept ----------
Xd <- model.matrix(~ age + gender + stage + smoking, data = clin_mod)[, -1, drop = FALSE]

# ---------- Scale genes for penalized modeling ----------
Xg_z <- scale(Xg)

# ---------- Save outputs ----------
preview_dir <- file.path(base_dir,"data","05_preliminar")
dir.create(preview_dir,showWarnings=FALSE,recursive=TRUE)
saveRDS(Xg,       file.path(preview_dir, "X_genes_signature_raw.rds"))
saveRDS(Xg_z,     file.path(preview_dir, "X_genes_signature_scaled.rds"))
saveRDS(Xd,       file.path(preview_dir, "X_clinic_noIntercept.rds"))
saveRDS(clin_mod, file.path(preview_dir, "clinical_for_modeling.rds"))

readr::write_tsv(tibble(variable = colnames(Xd)),
                 file.path(tab_dir, "05_clinic_design_columns.tsv"))

cat(sprintf("[OK] Design matrices ready. Rows=%d | signature genes=%d | clinical variables=%d\n",
            nrow(Xg_z), ncol(Xg_z), ncol(Xd)),
    file = log_file, append = TRUE)

message("[OK] Design matrices created.")
message("- Gene matrix (raw):    data/05_preliminar/X_genes_signature_raw.rds")
message("- Gene matrix (scaled): data/05_preliminar/X_genes_signature_scaled.rds")
message("- Clinical design:      data/05_preliminar/X_clinic_noIntercept.rds")
message("- Modeling clinical df: data/05_preliminar/clinical_for_modeling.rds")
