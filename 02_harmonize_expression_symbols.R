if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
# Paso 02 | Matriz de expresión TCGA-LUAD por símbolo génico
# Ejecutar en RStudio después del 01. También se puede iniciar en sesión nueva.
base_dir <- getwd()
se_path <- file.path(base_dir, "data_raw", "expression_SE_LUAD.rds")
out_dir <- file.path(base_dir, "data")
log_dir <- file.path(base_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Usar la copia guardada en el paso 01 (también funciona tras reiniciar R).
se_02 <- readRDS(se_path)
stopifnot("tpm_unstrand" %in% SummarizedExperiment::assayNames(se_02))
tpm <- SummarizedExperiment::assay(se_02, "tpm_unstrand")
annotation <- SummarizedExperiment::rowData(se_02)
symbol_field <- intersect(c("gene_name", "external_gene_name"), names(annotation))
if (!length(symbol_field)) stop("La anotación no contiene una columna de símbolos génicos.")
symbols <- as.character(annotation[[symbol_field[1]]])
rm(se_02, annotation)
keep <- !is.na(symbols) & nzchar(trimws(symbols))
genes <- sort(unique(symbols[keep]), method = "radix")
samples <- colnames(tpm)
stopifnot(length(genes) > 0, length(samples) == ncol(tpm), !anyNA(samples),
          all(nzchar(samples)), !anyDuplicated(samples))
group <- match(symbols[keep], genes)

# Misma regla del 02 original: primero log2(TPM+1), después media por símbolo.
# Procesar por bloques de muestras reduce las matrices temporales en memoria.
expr_sym <- data.frame(gene = genes, check.names = FALSE)
for (start in seq.int(1L, ncol(tpm), by = 20L)) {
  cols <- start:min(start + 19L, ncol(tpm))
  block <- as.matrix(tpm[keep, cols, drop = FALSE])
  if (any(!is.finite(block)) || any(block < 0)) stop("Hay TPM no finitos o negativos.")
  sums <- rowsum(log2(block + 1), group, reorder = TRUE)
  averages <- sums / tabulate(group, nbins = length(genes))
  expr_sym[samples[cols]] <- as.data.frame(averages)
}
stopifnot(identical(names(expr_sym)[-1], samples), !anyDuplicated(expr_sym$gene))

saveRDS(expr_sym, file.path(out_dir, "expr_SYMBOL_matrix.rds"))
readr::write_tsv(expr_sym, file.path(out_dir, "expr_SYMBOL_matrix.tsv"))
resumen <- c("Ensayo: tpm_unstrand | Escala: log2(TPM+1)",
             sprintf("Filas iniciales: %d | Sin símbolo: %d", length(symbols), sum(!keep)),
             sprintf("Símbolos únicos: %d | Muestras: %d", length(genes), length(samples)))
writeLines(resumen, file.path(log_dir, "02_gene_mapping.txt"))
cat(resumen, "[OK] Paso 02 completado.", sep = "\n")
