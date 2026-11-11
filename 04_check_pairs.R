#########################################
# LUAD age-bin workflow
# Check and export strict tumor-normal pairs
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

# Set working directory to the folder that contains the meta/ directory.
# Expected inputs:
#   meta/GSE115002_meta.csv
#   meta/GSE75037_meta.csv
ROOT     <- getwd()
META_DIR <- file.path(ROOT, "meta")
OUT_DIR  <- file.path(META_DIR, "pairs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GSE_IDS <- c("GSE115002", "GSE75037")

# Age bins used in the paired analysis.
is_young <- function(a) !is.na(a) & a < 60
is_old   <- function(a) !is.na(a) & a >= 60

as_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x
}

# Patient-level age consensus between matched normal and tumor samples.
# If ages differ by <=1 year, the mean is used; otherwise, the modal rounded age is used.
age_consensus <- function(a, b) {
  aa <- suppressWarnings(as.numeric(c(a, b)))
  aa <- aa[!is.na(aa)]

  if (!length(aa)) return(NA_real_)
  if (length(aa) == 1) return(aa)

  if (abs(diff(range(aa))) <= 1) return(round(mean(aa), 1))

  tbl <- sort(table(round(aa, 1)), decreasing = TRUE)
  as.numeric(names(tbl)[1])
}

bin_from_age <- function(a) {
  if (is.na(a)) return(NA_character_)
  if (is_young(a)) return("young")
  if (is_old(a))   return("old")
  NA_character_
}

make_pairs_for_gse <- function(gse) {
  f_meta <- file.path(META_DIR, paste0(gse, "_meta.csv"))

  if (!file.exists(f_meta)) {
    message("Metadata file not found: ", f_meta)
    return(invisible(NULL))
  }

  meta <- readr::read_csv(f_meta, show_col_types = FALSE) %>%
    mutate(
      sample_id  = as_chr(sample_id),
      patient_id = as_chr(patient_id),
      tissue = case_when(
        grepl("^normal$", tissue, ignore.case = TRUE) ~ "Normal",
        grepl("^tumor$",  tissue, ignore.case = TRUE) ~ "Tumor",
        TRUE ~ as_chr(tissue)
      ),
      age     = suppressWarnings(as.numeric(age)),
      age_bin = as_chr(age_bin),
      title   = as_chr(title)
    ) %>%
    filter(!is.na(sample_id), !is.na(patient_id), tissue %in% c("Normal", "Tumor"))

  if (nrow(meta) == 0) {
    message("[", gse, "] No valid normal/tumor samples after filtering metadata.")
    return(invisible(NULL))
  }

  # Count normal and tumor samples per patient.
  pair_counts <- meta %>%
    count(patient_id, tissue) %>%
    pivot_wider(names_from = tissue, values_from = n, values_fill = 0)

  if (!"Normal" %in% names(pair_counts)) pair_counts$Normal <- 0L
  if (!"Tumor"  %in% names(pair_counts)) pair_counts$Tumor  <- 0L

  pair_issues <- pair_counts %>%
    filter(Normal != 1 | Tumor != 1)

  # Keep only strict 1 normal + 1 tumor pairs.
  keep_ids <- pair_counts %>%
    filter(Normal == 1, Tumor == 1) %>%
    pull(patient_id)

  normal_samples <- meta %>%
    filter(tissue == "Normal", patient_id %in% keep_ids) %>%
    select(patient_id, N_sample_id = sample_id, N_age = age, N_age_bin = age_bin)

  tumor_samples <- meta %>%
    filter(tissue == "Tumor", patient_id %in% keep_ids) %>%
    select(patient_id, T_sample_id = sample_id, T_age = age, T_age_bin = age_bin)

  pairs <- inner_join(normal_samples, tumor_samples, by = "patient_id") %>%
    rowwise() %>%
    mutate(
      age_consensus = age_consensus(N_age, T_age),
      bin_consensus = bin_from_age(age_consensus),
      agree_N = (!is.na(N_age) & bin_from_age(N_age) == bin_consensus),
      agree_T = (!is.na(T_age) & bin_from_age(T_age) == bin_consensus)
    ) %>%
    ungroup() %>%
    select(
      patient_id, N_sample_id, T_sample_id,
      N_age, T_age, age_consensus, bin_consensus,
      N_age_bin, T_age_bin, agree_N, agree_T
    )

  out_pairs <- file.path(OUT_DIR, paste0(gse, "_pairs.tsv"))
  readr::write_tsv(pairs, out_pairs)

  if (nrow(pair_issues)) {
    out_issues <- file.path(OUT_DIR, paste0(gse, "_pair_issues.tsv"))
    readr::write_tsv(pair_issues, out_issues)
  }

  summary_by_bin <- pairs %>%
    count(bin_consensus, name = "n_pairs") %>%
    arrange(factor(bin_consensus, levels = c("young", "old", NA)))

  message("\n== ", gse, " ==")
  message("Strict 1N+1T pairs: ", nrow(pairs))

  if (nrow(pair_issues)) {
    message("Pairing issues: ", file.path(OUT_DIR, paste0(gse, "_pair_issues.tsv")))
  } else {
    message("No pairing issues detected.")
  }

  message("Pairs by age bin:")
  print(summary_by_bin, n = Inf)
  message("Pair table: ", out_pairs)

  invisible(list(pairs = pairs, issues = pair_issues, summary_by_bin = summary_by_bin))
}

res <- lapply(GSE_IDS, make_pairs_for_gse)

message("\nDone. Output directory: ", OUT_DIR)
