if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
# Paso 07. Panel preespecificado CellCycle_LUAD_v1 del repositorio TCGA_Nomogram.
base_dir <- getwd()
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(stringr); library(ggplot2)
})

data_dir <- file.path(base_dir, "data")
tab_dir  <- file.path(base_dir, "tables");  dir.create(tab_dir, TRUE, TRUE)
fig_dir  <- file.path(base_dir, "figures"); dir.create(fig_dir, TRUE, TRUE)
log_dir  <- file.path(base_dir, "logs");    dir.create(log_dir, TRUE, TRUE)
log_file <- file.path(log_dir, "07_signature.txt")
cat("# 07_signature_ingest_and_validation\n", file = log_file)

# ----- Cargar firma -----
genes_sig <- c(
  "KIF20A", "RRM2", "BUB1", "CCNB1", "CENPA", "AURKB",
  "AURKA", "PBK", "TOP2A", "CENPF", "CDCA3", "NUSAP1",
  "UBE2C", "CHEK1", "MELK", "ASPM", "KIF11", "MAD2L1",
  "MKI67", "KIFC1", "FOXM1", "CCNA2", "KIF4A", "TTK",
  "CCNB2", "TPX2", "PRC1", "CDCA8", "KIF2C", "BIRC5"
)
dir.create(file.path(base_dir,"sig"),showWarnings=FALSE)
sig_csv <- file.path(base_dir,"sig","CellCycle_LUAD_v1.csv")
if (!file.exists(sig_csv)) write_csv(tibble(gene_symbol=genes_sig),sig_csv)
existing <- read_csv(sig_csv,show_col_types=FALSE)
stopifnot("gene_symbol" %in% names(existing), length(existing$gene_symbol)==30,
          setequal(existing$gene_symbol,genes_sig))
cat(sprintf("Genes en firma: %d\n", length(genes_sig)), file=log_file, append=TRUE)

# ----- Cargar expresión (cohorte supervivencia) -----
expr <- readRDS(file.path(data_dir,"expr_tumor_SYMBOL_filtered.rds"))  # col1=gene
stopifnot("gene" %in% names(expr))
expr$gene <- toupper(expr$gene)

# ----- Subset por firma y reporte de cobertura -----
present_mask  <- expr$gene %in% genes_sig
expr_sig      <- expr[present_mask, , drop=FALSE]
missing_genes <- setdiff(genes_sig, expr$gene)

cat(sprintf("Presentes: %d | Faltantes: %d\n", nrow(expr_sig), length(missing_genes)),
    file=log_file, append=TRUE)

# Guardar lista de faltantes y presencia
write_tsv(tibble(gene_symbol = sort(missing_genes)), file.path(base_dir,"sig","Signature_missing.tsv"))
write_tsv(tibble(gene_symbol = genes_sig,
                 present     = genes_sig %in% expr$gene),
          file.path(base_dir,"sig","Signature_resolved.tsv"))

clin_full <- readRDS(file.path(data_dir,"clin_survival_cohort.rds"))
stopifnot(length(missing_genes)==0L, nrow(expr_sig)==30L, !anyDuplicated(expr_sig$gene),
          !anyDuplicated(clin_full$barcode), identical(names(expr_sig)[-1],clin_full$barcode))

# ----- Z-score por gen (cohorte supervivencia) -----
mat <- as.matrix(expr_sig[,-1, drop=FALSE])
row_means <- rowMeans(mat, na.rm=TRUE)
row_sds   <- apply(mat, 1, sd, na.rm=TRUE)
stopifnot(all(is.finite(mat)), all(row_sds > 0), all(is.finite(row_sds)))
write_tsv(tibble(gene=expr_sig$gene, center=row_means, sd=row_sds),file.path(tab_dir,"07_standardization.tsv"))
mat_z <- (mat - row_means) / row_sds
expr_sig_z <- cbind(gene=expr_sig$gene, as.data.frame(mat_z, check.names=FALSE))

saveRDS(expr_sig_z, file.path(data_dir,"expr_SIGNATURE_SYMBOL_z.rds"))
write_tsv(expr_sig_z, file.path(data_dir,"expr_SIGNATURE_SYMBOL_z.tsv"))

# ----- También matrices z-score para las subcohortes complete-case (si existen) -----
build_z_for_subset <- function(clin_file, out_stub){
  f <- file.path(data_dir, clin_file)
  if (!file.exists(f)) return(invisible(NULL))
  clin <- readRDS(f)
  bars <- clin$barcode
  subm <- expr_sig_z[, c("gene", bars), drop=FALSE]
  saveRDS(subm, file.path(data_dir, paste0(out_stub,".rds")))
  write_tsv(subm, file.path(data_dir, paste0(out_stub,".tsv")))
}

# (A) CC clásico (edad+sexo+estadio)
build_z_for_subset("clin_complete_case.rds", "expr_SIGNATURE_SYMBOL_z_complete_case")
# (B) CC + smoking (si lo generaste en 04)
if (file.exists(file.path(data_dir,"clin_complete_case_w_smoking.rds"))) {
  build_z_for_subset("clin_complete_case_w_smoking.rds", "expr_SIGNATURE_SYMBOL_z_complete_case_w_smoking")
}

# ----- Figura simple de cobertura de firma -----
cover_tbl <- tibble(
  category = c("Genes presentes", "Genes faltantes"),
  n = c(nrow(expr_sig), length(missing_genes))
)
write_tsv(cover_tbl, file.path(tab_dir,"07_signature_coverage.tsv"))

p_cov <- ggplot(cover_tbl, aes(x=category, y=n, label=n)) +
  geom_col(width=.6) + geom_text(vjust=-.2, size=4) +
  labs(title="Cobertura de la firma en la matriz de expresión", x=NULL, y="Genes") +
  theme_bw(12)
ggsave(file.path(fig_dir,"07_signature_coverage.png"), p_cov, width=6, height=4, dpi=300)

cat("[OK] 07 listo. Matrices z-score y cobertura exportadas.\n", file=log_file, append=TRUE)

cat(sprintf("[OK] Paso 07 | Genes presentes: %d/30 | Pacientes: %d | Faltantes: %d\n", nrow(mat_z), ncol(mat_z), length(missing_genes)))
