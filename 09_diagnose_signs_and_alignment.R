if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
## 09_diagnose_signs_and_alignment.R
## Purpose: diagnose unexpected Cox coefficient signs in the TCGA-LUAD
##          signature model by checking alignment, univariable effects,
##          ridge coefficients, collinearity, PCA, VIF, and cluster-based refit.
## Inputs:
##   data/X_genes_signature_scaled.rds
##   data/clinical_for_modeling.rds
## Outputs:
##   tables/Diag_univariable_by_gene.tsv
##   tables/Diag_ridge_multivariable_betas.tsv
##   tables/Diag_spearman_cor_matrix.tsv
##   tables/Diag_PCA_loadings_PC1.tsv
##   tables/Diag_gene_clusters.tsv
##   tables/Diag_multivariable_classic_clusterRepresentatives.tsv
##   logs/09_diagnose_signs_and_alignment.txt

set.seed(2025L)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(survival)
  library(glmnet)
  library(stats)
  library(Matrix)
})

# ---------- Paths ----------
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
tab_dir  <- file.path(base_dir, "tables")
log_dir  <- file.path(base_dir, "logs")
fig_dir  <- file.path(base_dir, "figures")

invisible(lapply(c(tab_dir, log_dir, fig_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

log_file <- file.path(log_dir, "09_diagnose_signs_and_alignment.txt")
cat("# 09_diagnose_signs_and_alignment\n", file = log_file)

# ---------- Helpers ----------
norm_bar <- function(x) {
  x <- toupper(as.character(x))
  gsub("[.]", "-", x)
}

safe_cox_summary <- function(x, y) {
  fit <- try(survival::coxph(y ~ x, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(c(beta = NA_real_, HR = NA_real_, lo95 = NA_real_,
             hi95 = NA_real_, pval = NA_real_))
  }
  s <- summary(fit)
  c(
    beta = unname(s$coef[1, "coef"]),
    HR   = unname(s$coef[1, "exp(coef)"]),
    lo95 = unname(s$conf.int[1, "lower .95"]),
    hi95 = unname(s$conf.int[1, "upper .95"]),
    pval = unname(s$coef[1, "Pr(>|z|)"])
  )
}

# ---------- Load modeling objects ----------
x_path    <- file.path(data_dir, "X_genes_signature_scaled.rds")
clin_path <- file.path(data_dir, "clinical_for_modeling.rds")

if (!file.exists(x_path)) stop("Missing gene matrix: ", x_path)
if (!file.exists(clin_path)) stop("Missing clinical modeling file: ", clin_path)

Xg_z   <- readRDS(x_path)      # rows = barcodes, columns = genes
clin_m <- readRDS(clin_path)   # includes OS.time, OS.event, barcode

stopifnot(all(c("OS.time", "OS.event", "barcode") %in% names(clin_m)))

Xg_z <- as.matrix(Xg_z)
storage.mode(Xg_z) <- "double"

# ---------- Strict expression-clinical alignment ----------
rownames(Xg_z) <- norm_bar(rownames(Xg_z))
clin_m$barcode <- norm_bar(clin_m$barcode)

ids <- rownames(Xg_z)
stopifnot(!anyDuplicated(ids), !anyDuplicated(clin_m$barcode), setequal(ids,clin_m$barcode))
if (!all(ids %in% clin_m$barcode)) {
  missing_ids <- setdiff(ids, clin_m$barcode)
  stop("Some expression barcodes are absent from clinical_for_modeling.rds. First missing IDs: ",
       paste(head(missing_ids, 10), collapse = ", "))
}

clin_use <- clin_m[match(ids, clin_m$barcode), , drop = FALSE]
stopifnot(all(rownames(Xg_z) == clin_use$barcode))

keep <- complete.cases(clin_use[, c("OS.time", "OS.event")]) &
  apply(Xg_z, 1, function(r) all(is.finite(r)))

Xg_z     <- Xg_z[keep, , drop = FALSE]
clin_use <- clin_use[keep, , drop = FALSE]

if (nrow(Xg_z) < 30) stop("Too few complete cases after filtering: n = ", nrow(Xg_z))
if (ncol(Xg_z) < 2)  stop("This diagnostic script expects at least two genes.")

surv_y <- survival::Surv(
  time  = as.numeric(clin_use$OS.time),
  event = as.numeric(clin_use$OS.event)
)

n_subjects <- nrow(Xg_z)
n_genes    <- ncol(Xg_z)
n_events   <- sum(as.numeric(clin_use$OS.event), na.rm = TRUE)

cat(sprintf("N=%d | events=%d | genes=%d\n", n_subjects, n_events, n_genes),
    file = log_file, append = TRUE)

# ============================================================
# 1) Univariable Cox model by gene
# ============================================================
univ_mat <- t(vapply(colnames(Xg_z), function(g) {
  safe_cox_summary(Xg_z[, g], surv_y)
}, numeric(5)))

univ_df <- tibble(
  gene = rownames(univ_mat),
  beta = as.numeric(univ_mat[, "beta"]),
  HR   = as.numeric(univ_mat[, "HR"]),
  lo95 = as.numeric(univ_mat[, "lo95"]),
  hi95 = as.numeric(univ_mat[, "hi95"]),
  pval = as.numeric(univ_mat[, "pval"])
) %>%
  arrange(pval)

readr::write_tsv(univ_df, file.path(tab_dir, "Diag_univariable_by_gene.tsv"))

# Interpretation note:
# If a gene already shows HR < 1 in the univariable model, the direction is not
# caused by multivariable collinearity. If it is HR > 1 univariably but flips in
# multivariable models, collinearity or suppressor effects are likely.

# ============================================================
# 2) Penalized multivariable Cox model: ridge
# ============================================================
cv_ridge <- glmnet::cv.glmnet(
  x = Xg_z,
  y = surv_y,
  family = "cox",
  alpha = 0,
  nfolds = 10,
  lambda.min.ratio = 1e-4
)

lambda_ridge <- cv_ridge$lambda.min
ridge_coef <- as.matrix(coef(cv_ridge, s = lambda_ridge))

ridge_df <- tibble(
  gene = rownames(ridge_coef),
  beta = as.numeric(ridge_coef[, 1])
) %>%
  arrange(desc(abs(beta)))

readr::write_tsv(ridge_df, file.path(tab_dir, "Diag_ridge_multivariable_betas.tsv"))
saveRDS(cv_ridge, file.path(data_dir, "ridge_diagnostic_cv_fit.rds"))

# ============================================================
# 3) Collinearity diagnostics: Spearman correlation, PCA, VIF
# ============================================================
cor_mat <- suppressWarnings(cor(Xg_z, method = "spearman"))
readr::write_tsv(
  tibble::as_tibble(cor_mat, rownames = "gene"),
  file.path(tab_dir, "Diag_spearman_cor_matrix.tsv")
)

pca_fit <- prcomp(Xg_z, center = TRUE, scale. = TRUE)
pca_var <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)

cat(sprintf("PC1 variance=%.1f%% | PC2 variance=%.1f%%\n",
            100 * pca_var[1], 100 * pca_var[2]),
    file = log_file, append = TRUE)

pca_loadings <- tibble(
  gene = colnames(Xg_z),
  loading_PC1 = as.numeric(pca_fit$rotation[, 1])
) %>%
  arrange(desc(abs(loading_PC1)))

readr::write_tsv(pca_loadings, file.path(tab_dir, "Diag_PCA_loadings_PC1.tsv"))

# VIF is optional because package 'car' may not be installed.
tryCatch({
  genes_for_vif <- if (ncol(Xg_z) > 20) colnames(Xg_z)[seq_len(20)] else colnames(Xg_z)
  vif_df <- as.data.frame(Xg_z[, genes_for_vif, drop = FALSE])
  colnames(vif_df) <- make.names(colnames(vif_df), unique = TRUE)

  fit_vif <- survival::coxph(surv_y ~ ., data = vif_df, ties = "efron")

  if ("car" %in% rownames(installed.packages())) {
    suppressPackageStartupMessages(library(car))
    vif_values <- car::vif(fit_vif)
    readr::write_tsv(
      tibble(variable = names(vif_values), VIF = as.numeric(vif_values)),
      file.path(tab_dir, "Diag_VIF_subset.tsv")
    )
  } else {
    cat("Package 'car' not installed; VIF was skipped.\n", file = log_file, append = TRUE)
  }
}, error = function(e) {
  cat("VIF step skipped due to error: ", conditionMessage(e), "\n",
      file = log_file, append = TRUE, sep = "")
})

# ============================================================
# 4) Classic multivariable Cox model with one representative per cluster
#    Clusters are defined from 1 - |Spearman rho|.
# ============================================================
dist_mat <- as.dist(1 - abs(cor_mat))
hc <- hclust(dist_mat, method = "average")

# h = 1 - 0.85 groups genes with approximately |rho| >= 0.85.
clusters <- cutree(hc, h = 1 - 0.85)
cluster_df <- tibble(
  gene = names(clusters),
  cluster = as.integer(clusters)
)

readr::write_tsv(cluster_df, file.path(tab_dir, "Diag_gene_clusters.tsv"))

rep_by_cluster <- cluster_df %>%
  left_join(univ_df %>% select(gene, beta, pval), by = "gene") %>%
  group_by(cluster) %>%
  slice_max(order_by = abs(beta), n = 1, with_ties = FALSE) %>%
  ungroup()

genes_rep <- rep_by_cluster$gene
X_rep <- as.data.frame(Xg_z[, genes_rep, drop = FALSE])
safe_names <- make.names(colnames(X_rep), unique = TRUE)
name_map <- tibble(safe_name = safe_names, gene = colnames(X_rep))
colnames(X_rep) <- safe_names

fit_rep <- survival::coxph(surv_y ~ ., data = X_rep, ties = "efron")
rep_summary <- summary(fit_rep)

coef_names <- gsub("^`|`$", "", rownames(rep_summary$coef))
rep_genes <- name_map$gene[match(coef_names, name_map$safe_name)]

rep_df <- tibble(
  gene = rep_genes,
  beta = as.numeric(rep_summary$coef[, "coef"]),
  HR   = as.numeric(rep_summary$coef[, "exp(coef)"]),
  lo95 = as.numeric(rep_summary$conf.int[, "lower .95"]),
  hi95 = as.numeric(rep_summary$conf.int[, "upper .95"]),
  pval = as.numeric(rep_summary$coef[, "Pr(>|z|)"])
) %>%
  arrange(pval)

readr::write_tsv(rep_df, file.path(tab_dir, "Diag_multivariable_classic_clusterRepresentatives.tsv"))

# ---------- Console report ----------
message("---- Diagnostic summary ----")
message(sprintf("N=%d, events=%d, genes=%d", n_subjects, n_events, n_genes))
message("Generated diagnostic tables in 'tables/':")
message("  - Diag_univariable_by_gene.tsv")
message("  - Diag_ridge_multivariable_betas.tsv")
message("  - Diag_spearman_cor_matrix.tsv")
message("  - Diag_PCA_loadings_PC1.tsv")
message("  - Diag_VIF_subset.tsv (if package 'car' is installed)")
message("  - Diag_gene_clusters.tsv")
message("  - Diag_multivariable_classic_clusterRepresentatives.tsv")
message("Log file: logs/09_diagnose_signs_and_alignment.txt")
