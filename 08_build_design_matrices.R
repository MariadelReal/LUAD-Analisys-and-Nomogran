if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
## 08_build_design_matrices.R (v2)
## Construye:
##  - Xg_z: matriz de expresión z-score (pacientes × genes candidatos)
##  - Xd: diseño clínico sin intercepto
##  - clin_mod: clínica lista para modelar (incluye OS.time / OS.event)
##
## Usa:
##  - data/clin_survival_cohort.rds
##  - data/expr_SIGNATURE_SYMBOL_z.rds   (salida de 07_signature_ingest...)

base_dir <- getwd()

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(stringr)
})

data_dir <- file.path(base_dir, "data")
tab_dir  <- file.path(base_dir, "tables"); dir.create(tab_dir, TRUE, TRUE)
log_dir  <- file.path(base_dir, "logs");   dir.create(log_dir, TRUE, TRUE)
log_file <- file.path(log_dir, "08_build_design.txt")
cat("# 08_build_design_matrices (v2)\n", file = log_file)

## --- Cargar cohorte clínica de supervivencia (de 04) ---
clin <- readRDS(file.path(data_dir, "clin_survival_cohort.rds"))

## --- Cargar matriz de firma z-score (de 07) ---
expr_sig_z <- readRDS(file.path(data_dir, "expr_SIGNATURE_SYMBOL_z.rds"))
expr_sig_z <- as.data.frame(expr_sig_z, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("gene" %in% colnames(expr_sig_z))

## Genes y barcodes de expresión
genes_expr <- expr_sig_z$gene
bars_expr  <- setdiff(colnames(expr_sig_z), "gene")

## Comprobar que la clínica trae 'barcode'
if (!"barcode" %in% colnames(clin)) {
  stop("clin_survival_cohort.rds does not contain 'barcode'. Check script 04.")
}

bars_clin <- as.character(clin$barcode)

## Intersección de barcodes
stopifnot(!anyDuplicated(names(expr_sig_z)), !anyDuplicated(bars_clin),
          !anyDuplicated(genes_expr), length(genes_expr)==30L, setequal(bars_expr,bars_clin))
common_bc <- bars_expr
cat(sprintf("Expr barcodes: %d | Clin barcodes: %d | Common: %d\n",
            length(bars_expr), length(bars_clin), length(common_bc)),
    file = log_file, append = TRUE)

if (length(common_bc) == 0) {
  stop("No common barcodes between expr_SIGNATURE_SYMBOL_z and clin_survival_cohort.")
}

## Reordenar expresión y clínica según common_bc
expr_sig_z_sub <- expr_sig_z[, c("gene", common_bc), drop = FALSE]
clin_sub       <- clin[match(common_bc, clin$barcode), , drop = FALSE]
stopifnot(all(colnames(expr_sig_z_sub)[-1] == clin_sub$barcode))

## --- Construir Xg_z (muestras × genes) ---
Xg_z <- t(as.matrix(expr_sig_z_sub[, -1, drop = FALSE]))
colnames(Xg_z) <- genes_expr
rownames(Xg_z) <- common_bc

cat(sprintf("Xg_z: %d samples × %d genes\n", nrow(Xg_z), ncol(Xg_z)),
    file = log_file, append = TRUE)

## --- Clínica para modelado ---
clin_mod <- clin_sub %>%
  mutate(
    age     = as.numeric(age_years),
    gender  = factor(gender, levels = c("Female", "Male")),                  # reference: Female
    stage   = factor(stage_use, levels = c("I", "II", "III", "IV")),         # reference: I
    smoking = factor(smoking_status, levels = c("Never","Former","Current","Unknown")) # ref: Never
  )

clin_mod <- as.data.frame(clin_mod)
rownames(clin_mod) <- clin_mod$barcode

## --- Matriz de diseño clínica SIN intercepto ---
clin_cc <- clin_mod[complete.cases(clin_mod[,c("age","gender","stage","smoking")]), ,drop=FALSE]
Xd <- model.matrix(~ age + gender + stage + smoking, data = clin_cc)[, -1, drop = FALSE]
stopifnot(identical(rownames(Xd),clin_cc$barcode), identical(rownames(Xg_z),clin_mod$barcode),
          all(is.finite(Xg_z)), all(is.finite(Xd)))

cat(sprintf("Xd: %d samples × %d clinical covariates\n",
            nrow(Xd), ncol(Xd)),
    file = log_file, append = TRUE)

## Guardar nombres de columnas clínicas
readr::write_tsv(
  tibble(variable = colnames(Xd)),
  file.path(tab_dir, "08_clinic_design_columns.tsv")
)

## --- Guardados para modelado ---
saveRDS(Xg_z,     file.path(data_dir, "X_genes_signature_scaled.rds"))  # z-score already
saveRDS(Xd,       file.path(data_dir, "X_clinic_noIntercept.rds"))
saveRDS(clin_mod, file.path(data_dir, "clinical_for_modeling.rds"))

cat(sprintf("[OK] 08 done. Samples = %d | Genes = %d | Clinical vars = %d\n",
            nrow(Xg_z), ncol(Xg_z), ncol(Xd)),
    file = log_file, append = TRUE)
cat(sprintf("[OK] Paso 08 | Matriz genética: %d x %d | Matriz clínica: %d x %d | Muertes en casos completos: %d\n",nrow(Xg_z),ncol(Xg_z),nrow(Xd),ncol(Xd),sum(clin_cc$OS.event)))
