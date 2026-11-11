if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
#########################################
# TCGA-LUAD Cox/LASSO workflow
# 00_env_check_and_sessioninfo.R
#
# Purpose:
#   Prepare the project folder structure, set a fixed random seed,
#   check/install required packages, and export session information.
#
# Notes:
#   Run this script from the root directory of the project/repository.
#   Outputs are written to logs/.
#########################################


set.seed(123L)

# ---------- Project directory ----------
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
seed_val <- 123L

# ---------- Folder structure ----------
subdirs <- c(
  "config", "data_raw", "data", "sig", "logs", "models",
  "tables", "figures", "reports", "scripts"
)

invisible(sapply(
  file.path(base_dir, subdirs),
  dir.create,
  showWarnings = FALSE,
  recursive = TRUE
))

# ---------- Required packages ----------
cran_pkgs <- c(
  "dplyr", "readr", "tibble", "sessioninfo"
)

bioc_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment"
)

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran)) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, dependencies = TRUE)
}

missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc)) {
  if (!"BiocManager" %in% rownames(installed.packages())) {
    install.packages("BiocManager")
  }
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# ---------- Load packages ----------
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
  library(sessioninfo)
})

# ---------- Export session information ----------
log_dir   <- file.path(base_dir, "logs")
sess_file <- file.path(log_dir, "sessionInfo.txt")
pkg_file  <- file.path(log_dir, "package_versions.tsv")

sink(sess_file)
print(sessionInfo())
sink()

pkg_df <- as.data.frame(
  installed.packages()[, c("Package", "Version")],
  stringsAsFactors = FALSE
)
pkg_df <- tibble::as_tibble(pkg_df[order(pkg_df$Package), ])
readr::write_tsv(pkg_df, pkg_file)

cat(
  "\n[OK] Environment prepared.\n",
  "- Base dir: ", base_dir, "\n",
  "- Seed: ", seed_val, "\n",
  "- sessionInfo: ", sess_file, "\n",
  "- package versions: ", pkg_file, "\n",
  sep = ""
)
