if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
## 06_cox_univariate_30genes.R
## Cox univariada para 30 genes mitóticos en TCGA-LUAD
## Usa:
##   - data/clin_survival_cohort.rds   (salida de 04)
##   - data/expr_tumor_SYMBOL_filtered.rds (salida de 04)
##
## Salidas:
##   - tables/Cox_univariate_30genes.tsv
##   - tables/Cox_univariate_30genes_with_FDR.tsv
##   - figures/Cox_univariate_30genes_forest.png / .pdf

base_dir <- getwd()

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
  library(survival); library(ggplot2)
})

data_dir <- file.path(base_dir, "data")
tab_dir  <- file.path(base_dir, "tables");  dir.create(tab_dir, TRUE, TRUE)
fig_dir  <- file.path(base_dir, "figures"); dir.create(fig_dir, TRUE, TRUE)
log_dir  <- file.path(base_dir, "logs");    dir.create(log_dir, TRUE, TRUE)
log_file <- file.path(log_dir, "06_cox_univariate_30genes.txt")
cat("# 06_cox_univariate_30genes\n", file = log_file)

## --------- Lista de genes de interés ---------
genes_interest <- c(
  "CDCA3","TPX2","AURKB","FOXM1","PRC1","CCNB2","BIRC5","KIFC1",
  "CENPA","KIF20A","BUB1","PBK","KIF4A","KIF2C","CHEK1","UBE2C",
  "NUSAP1","TOP2A","MKI67","TTK","CDCA8","MELK","ASPM","AURKA",
  "CENPF","RRM2","CCNB1","KIF11","CCNA2","MAD2L1"
)

## --------- Cargar clínica y expresión (outputs de 04) ---------
clin <- readRDS(file.path(data_dir, "clin_survival_cohort.rds"))
expr <- readRDS(file.path(data_dir, "expr_tumor_SYMBOL_filtered.rds"))

## me aseguro de que expr sea data.frame y tenga 'gene'
expr <- as.data.frame(expr, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene" %in% colnames(expr)) {
  stop("expr_tumor_SYMBOL_filtered.rds NO tiene columna 'gene'. Revisa el script 04.")
}

cat(sprintf("N clínica (cohorte supervivencia): %d\n", nrow(clin)),
    file = log_file, append = TRUE)

## expr: fila = gen, columnas = barcodes
genes_present <- intersect(genes_interest, expr$gene)
genes_missing <- setdiff(genes_interest, expr$gene)

cat("Genes presentes en expresión: ",
    paste(genes_present, collapse = ", "), "\n",
    file = log_file, append = TRUE)
if (length(genes_missing)) {
  cat("Genes AUSENTES en la matriz: ",
      paste(genes_missing, collapse = ", "), "\n",
      file = log_file, append = TRUE)
}

if (length(genes_missing)) stop("Faltan genes: ",paste(genes_missing,collapse=", "))
stopifnot(!anyDuplicated(expr$gene))

## Subset solo a esos genes
expr_sub <- expr[expr$gene %in% genes_present, , drop = FALSE]

## Asegurar orden consistente de genes
expr_sub$gene <- factor(expr_sub$gene, levels = genes_present)
expr_sub <- expr_sub[order(expr_sub$gene), , drop = FALSE]
expr_sub$gene <- as.character(expr_sub$gene)

## Extraer matriz de expresión (sin la columna gene)
mat <- as.matrix(expr_sub[, setdiff(colnames(expr_sub), "gene"), drop = FALSE])

## Transponer: filas = muestras (barcodes), columnas = genes
expr_mat <- t(mat)  # dim: N_muestras × G_genes
colnames(expr_mat) <- expr_sub$gene   # AQUÍ fijamos los nombres de columnas

## Comprobaciones
cat(sprintf("Dim expr_mat: %d muestras × %d genes\n",
            nrow(expr_mat), ncol(expr_mat)),
    file = log_file, append = TRUE)
cat("Colnames(expr_mat)[1:5]: ",
    paste(colnames(expr_mat)[1:min(5, ncol(expr_mat))], collapse = ", "), "\n",
    file = log_file, append = TRUE)

## alineamos con clínica vía barcode
bar_expr <- rownames(expr_mat)

if (!"barcode" %in% colnames(clin)) {
  stop("clin_survival_cohort.rds NO tiene columna 'barcode'. Revisa el script 04.")
}
bar_clin <- clin$barcode
common_bc <- intersect(bar_expr, bar_clin)

cat(sprintf("Barcodes en expresión: %d | en clínica: %d | comunes: %d\n",
            length(bar_expr), length(bar_clin), length(common_bc)),
    file = log_file, append = TRUE)

stopifnot(!anyDuplicated(bar_expr), !anyDuplicated(bar_clin), setequal(bar_expr,bar_clin))

## reordenar expr_mat y clin_sub según common_bc
expr_mat <- expr_mat[common_bc, , drop = FALSE]
clin_sub <- clin[match(common_bc, clin$barcode), ]
stopifnot(all(rownames(expr_mat) == clin_sub$barcode))

# La matriz del 02 ya está en log2(TPM+1).
stopifnot(is.numeric(expr_mat), all(is.finite(expr_mat)), all(clin_sub$OS.event %in% 0:1),
          all(is.finite(clin_sub$OS.time)), all(clin_sub$OS.time > 0))

## --------- Cox univariada ---------
res_list <- vector("list", length = ncol(expr_mat))
names(res_list) <- colnames(expr_mat)

for (g in colnames(expr_mat)) {
  x <- as.numeric(expr_mat[, g])
  df <- data.frame(
    OS.time  = as.numeric(clin_sub$OS.time),
    OS.event = as.numeric(clin_sub$OS.event),
    gene_exp = x
  )
  
  keep <- complete.cases(df)
  df   <- df[keep, , drop = FALSE]
  
  if (nrow(df) < 20 || length(unique(df$gene_exp)) < 3) {
    res_list[[g]] <- tibble(
      gene    = g,
      n       = nrow(df),
      events  = sum(df$OS.event, na.rm = TRUE),
      HR      = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      z       = NA_real_,
      pval    = NA_real_
    )
    next
  }
  
  fit <- coxph(Surv(OS.time, OS.event) ~ gene_exp, data = df)
  s   <- summary(fit)
  
  HR    <- unname(s$coef[,"exp(coef)"][1])
  z     <- unname(s$coef[,"z"][1])
  p     <- unname(s$coef[,"Pr(>|z|)"][1])
  lower <- unname(s$conf.int[,"lower .95"][1])
  upper <- unname(s$conf.int[,"upper .95"][1])
  
  res_list[[g]] <- tibble(
    gene    = g,
    n       = s$n,
    events  = s$nevent,
    HR      = HR,
    lower95 = lower,
    upper95 = upper,
    z       = z,
    pval    = p
  )
}

cox_tbl <- bind_rows(res_list)

if (!"pval" %in% names(cox_tbl)) {
  stop("La tabla de resultados no tiene columna 'pval'. Algo fue mal en res_list.")
}

cox_tbl <- cox_tbl %>%
  mutate(FDR = p.adjust(pval, method = "BH")) %>%
  arrange(pval)

write_tsv(cox_tbl, file.path(tab_dir, "Cox_univariate_30genes.tsv"))
write_tsv(cox_tbl, file.path(tab_dir, "Cox_univariate_30genes_with_FDR.tsv"))

cat("\n[OK] Cox univariada completada para 30 genes.\n",
    file = log_file, append = TRUE)

## --------- Forest plot ---------
cox_plot <- cox_tbl %>%
  filter(!is.na(HR)) %>%
  mutate(
    gene = factor(gene, levels = rev(genes_interest[genes_interest %in% gene]))
  )

p_forest <- ggplot(cox_plot,
                   aes(x = HR, y = gene)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point() +
  geom_segment(aes(x = lower95, xend = upper95, yend = gene)) +
  scale_x_log10() +
  labs(
    x = "Hazard ratio per unit of log2(TPM + 1) (95% CI)",
    y = "Gene",
    title = "Univariate Cox analysis of 30 mitotic genes (TCGA-LUAD)"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(fig_dir, "Cox_univariate_30genes_forest.png"),
       p_forest, width = 7, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "Cox_univariate_30genes_forest.pdf"),
       p_forest, width = 7, height = 6)
cat(sprintf("[OK] Paso 06 | Genes evaluados: %d | Pacientes: %d | Muertes: %d | FDR < 0.05: %d\n",sum(is.finite(cox_tbl$pval)),nrow(clin_sub),sum(clin_sub$OS.event),sum(cox_tbl$FDR < 0.05,na.rm=TRUE)))
print(cox_tbl[cox_tbl$gene == "CCNB1", ])
