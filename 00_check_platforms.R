#########################################
# LUAD age-bin transcriptomic workflow
# 00_check_platforms.R
#
# Purpose:
# Check whether each GEO series contains one or more GPL platforms and verify
# that the locally saved *_series.tsv files contain the same GSM samples as GEO.
#
# Expected local input:
#   GSE115002_series.tsv
#   GSE75037_series.tsv
#
# Output, only if samples are missing:
#   qc0/<GSE>_GSM_missing_from_series.tsv
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(readr)
})

# Set the directory containing the *_series.tsv files.
# If running from the original analysis folder, keep ROOT <- getwd().
# Otherwise, replace this with the local path to the Series Matrix folder.
ROOT <- getwd()

GSE_IDS <- c("GSE115002", "GSE75037")

check_gse_platforms <- function(gse_id, root_dir = ROOT) {
  cat("\n== ", gse_id, " ==\n", sep = "")

  es <- try(GEOquery::getGEO(gse_id, GSEMatrix = TRUE), silent = TRUE)
  if (inherits(es, "try-error") || is.null(es)) {
    cat("  Could not load data from GEO.\n")
    return(invisible(NULL))
  }

  es_list <- if (is.list(es)) es else list(es)

  gpls <- vapply(
    es_list,
    function(x) tryCatch(Biobase::annotation(x), error = function(e) NA_character_),
    character(1)
  )

  cat("  GPL platforms detected: ", paste(unique(gpls), collapse = ", "), "\n", sep = "")

  # GSM accessions reported by GEO.
  gsm_geo <- unique(unlist(lapply(es_list, function(x) rownames(Biobase::pData(x)))))
  n_geo <- length(gsm_geo)

  # GSM accessions present in local *_series.tsv files.
  # This accepts files named either GSE_series.tsv or GSE_GPLxxx_series.tsv.
  pattern <- paste0("^", gse_id, "(_GPL[^_]*)?_series\\.tsv$")
  local_files <- list.files(root_dir, pattern = pattern, full.names = TRUE)

  if (!length(local_files)) {
    cat("  No local *_series.tsv files were found in: ", root_dir, "\n", sep = "")
    return(invisible(NULL))
  }

  gsm_local <- unique(unlist(lapply(local_files, function(f) {
    df <- suppressMessages(readr::read_tsv(f, show_col_types = FALSE, n_max = 2))
    setdiff(colnames(df), "FEATURE")
  })))
  n_local <- length(gsm_local)

  cat("  GSM in GEO       : ", n_geo, "\n", sep = "")
  cat("  GSM in local TSV : ", n_local, "  (files: ", length(local_files), ")\n", sep = "")

  if (setequal(gsm_local, gsm_geo)) {
    cat("  Complete local series files; no multi-GPL preprocessing was required.\n")
  } else {
    missing_gsm <- setdiff(gsm_geo, gsm_local)
    cat("  Missing GSM samples in local TSV files: ", length(missing_gsm), "\n", sep = "")

    out_dir <- file.path(root_dir, "qc0")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    out_file <- file.path(out_dir, paste0(gse_id, "_GSM_missing_from_series.tsv"))
    readr::write_tsv(tibble::tibble(GSM = missing_gsm), out_file)
    cat("  Missing-sample report written to: ", out_file, "\n", sep = "")
  }

  invisible(NULL)
}

invisible(lapply(GSE_IDS, check_gse_platforms))
cat("\nDone.\n")
