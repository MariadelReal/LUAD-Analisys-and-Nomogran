if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
## 01_download_gdc_expression_clinical.R
## Propósito: descargar TCGA-LUAD (RNA-seq counts + clínica) desde GDC y registrar conteos.

# ---------- Parámetros mínimos ----------
base_dir <- getwd()
project  <- "TCGA-LUAD"
workflow <- "STAR - Counts"

# ---------- Cargar paquetes ----------
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---------- Rutas de salida ----------
raw_dir  <- file.path(base_dir, "data_raw")
log_dir  <- file.path(base_dir, "logs")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_counts <- file.path(log_dir, "01_counts.txt")
log_dl     <- file.path(log_dir, "01_download.txt")

cat("# Log de descarga TCGA-LUAD (GDC)\n", file = log_dl)

message("Carpeta de la réplica: ", base_dir)

# ---------- Descarga EXPRESIÓN (SummarizedExperiment) ----------
cat(sprintf("%s | Iniciando consulta de expresión (%s, %s)\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), project, workflow),
    file = log_dl, append = TRUE)

query_exp <- GDCquery(project = project,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = workflow)

# Punto de revisión antes de descargar los archivos
print(nrow(getResults(query_exp)))

# ---------- Ejecutar la descarga ----------
GDCdownload(query_exp, directory = raw_dir)

se_exp <- GDCprepare(query_exp,
                     directory = raw_dir,
                     save      = FALSE)

saveRDS(se_exp, file.path(raw_dir, "expression_SE_LUAD.rds"))

# ---------- Conteos por tipo de muestra (robusto, sin asumir columnas) ----------
cd <- as.data.frame(SummarizedExperiment::colData(se_exp))

# Asegurar 'barcode'
if (!"barcode" %in% names(cd)) {
  cd$barcode <- rownames(cd)
}

# Derivar short_letter_code desde el barcode (posiciones 14-15)
# 01 = Primary Tumor, 11 = Solid Tissue Normal, etc.
cd$short_letter_code <- substr(cd$barcode, 14, 15)

# Obtener sample_type; si no existe en colData, derivarlo por código
if (!"sample_type" %in% names(cd)) {
  code_map <- c(
    "01" = "Primary Tumor",
    "02" = "Recurrent Tumor",
    "03" = "Primary Blood-Derived Cancer - Peripheral Blood",
    "06" = "Metastatic",
    "11" = "Solid Tissue Normal",
    "12" = "Blood Derived Normal"
  )
  cd$sample_type <- ifelse(cd$short_letter_code %in% names(code_map),
                           code_map[cd$short_letter_code],
                           "Unknown")
}

st <- cd[, c("barcode","sample_type","short_letter_code")]
st[] <- lapply(st, as.character)

# Conteos y exportes
counts <- st %>%
  dplyr::count(sample_type, name = "n") %>%
  dplyr::arrange(dplyr::desc(n))

readr::write_tsv(counts, file.path(log_dir, "01_counts_by_sample_type.tsv"))

# Guardar resumen legible
cat("# Conteos por tipo de muestra (EXPRESIÓN)\n", file = log_counts)
apply(counts, 1, function(r) {
  cat(sprintf("- %s: %s\n", r[1], r[2]), file = log_counts, append = TRUE)
})

# ---------- Descarga CLÍNICA cruda ----------
cat(sprintf("%s | Descargando clínica cruda\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    file = log_dl, append = TRUE)

clin_raw <- GDCquery_clinic(project = project, type = "clinical")
readr::write_tsv(as_tibble(clin_raw), file.path(raw_dir, "clinical_LUAD.tsv"))

# ---------- Cierre ----------
cat(sprintf("%s | Descarga completa. Archivos en: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), raw_dir),
    file = log_dl, append = TRUE)

cat("\n[OK] Descarga completada.\n",
    "- SummarizedExperiment: data_raw/expression_SE_LUAD.rds\n",
    "- Clínica cruda:        data_raw/clinical_LUAD.tsv\n",
    "- Conteos por tipo:     logs/01_counts_by_sample_type.tsv + logs/01_counts.txt\n",
    "- Log de descarga:      logs/01_download.txt\n", sep = "")

# Verificación interactiva opcional
if (interactive()) {
  clin_raw_chk <- readr::read_tsv(file.path(raw_dir, "clinical_LUAD.tsv"),
                                  show_col_types = FALSE)
  print(names(clin_raw_chk)[1:25])
  glimpse(clin_raw_chk)
}

file.exists(file.path(raw_dir, c("expression_SE_LUAD.rds", "clinical_LUAD.tsv")))
