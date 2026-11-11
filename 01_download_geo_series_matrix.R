#########################################
# LUAD age-bin transcriptomic workflow
# 01_download_geo_series_matrix.R
#
# Purpose:
#   Download GEO Series Matrix ExpressionSet objects for GSE75037 and
#   GSE115002, then export expression and phenotype tables as TSV files.
#
# Outputs:
#   <ROOT>/GSE75037_series.tsv
#   <ROOT>/GSE75037_pheno.tsv
#   <ROOT>/GSE115002_series.tsv
#   <ROOT>/GSE115002_pheno.tsv
#
# Notes:
#   Both datasets were verified as single-platform datasets in this workflow:
#     GSE75037  : GPL6884
#     GSE115002 : GPL13497
#   If GEOquery returns a list of ExpressionSet objects, the object with the
#   largest number of features is selected.
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(readr)
})

# Set analysis directory.
# Run this script from the folder where the GEO TSV files should be written,
# or replace getwd() with your local analysis directory.
ROOT <- getwd()
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)

GSE_IDS <- c("GSE75037", "GSE115002")

pick_main_eset <- function(gse_object) {
  if (is.list(gse_object)) {
    n_features <- sapply(gse_object, function(x) nrow(Biobase::exprs(x)))
    gse_object[[which.max(n_features)]]
  } else {
    gse_object
  }
}

download_and_save_series <- function(gse_id, outdir = ROOT) {
  message("\n==> Downloading ", gse_id, " ...")

  gse_object <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  eset <- pick_main_eset(gse_object)

  expr_matrix <- Biobase::exprs(eset)
  pheno_table <- Biobase::pData(eset)

  expr_out <- file.path(outdir, paste0(gse_id, "_series.tsv"))
  pheno_out <- file.path(outdir, paste0(gse_id, "_pheno.tsv"))

  readr::write_tsv(
    data.frame(FEATURE = rownames(expr_matrix), expr_matrix, check.names = FALSE),
    expr_out
  )
  readr::write_tsv(as.data.frame(pheno_table), pheno_out)

  message("Saved: ", basename(expr_out),
          " (", nrow(expr_matrix), " features x ", ncol(expr_matrix), " samples)")
  message("Saved: ", basename(pheno_out),
          " (", nrow(pheno_table), " rows)")

  invisible(list(expression = expr_out, phenotype = pheno_out))
}

invisible(lapply(GSE_IDS, download_and_save_series))
message("\nDone. Files written to:\n", ROOT)
