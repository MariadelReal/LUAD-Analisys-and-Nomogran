#########################################
# LUAD age-bin transcriptomic workflow
# 03_build_metadata_and_pairs.R
#
# Purpose:
#   Build harmonized sample metadata for GSE75037 and GSE115002 from GEO
#   phenotype TSV files, assign age bins (<60 and >=60 years), validate
#   patient-level tumor-normal pairing, and export strict 1N-1T pair tables.
#
# Expected local input:
#   GSE75037_pheno.tsv
#   GSE115002_pheno.tsv
#
# Main outputs:
#   meta/GSE75037_meta.csv
#   meta/GSE115002_meta.csv
#   meta/summary.tsv
#   meta/pairs/GSE75037_pairs.tsv
#   meta/pairs/GSE115002_pairs.tsv
#   meta/*_pair_issues.tsv, if pairing issues are detected
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(glue); library(tibble)
})

# ---- Analysis directory ----
# Run this script from the folder containing the GEO phenotype TSV files
# produced in step 01_download_geo_series_matrix.R, or replace getwd() with that path.
ROOT <- getwd()
META_DIR <- file.path(ROOT, "meta")
OUT_DIR  <- file.path(META_DIR, "pairs")
dir.create(META_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)

# ---- Age-bin rules (<60 vs >=60 years) ----
is_young <- function(a) !is.na(a) & a < 60
is_old   <- function(a) !is.na(a) & a >= 60

# ---- Patient-level age/bin consensus for 1N-1T pairs ----
age_consensus <- function(a, b){
  aa <- suppressWarnings(as.numeric(c(a, b)))
  aa <- aa[!is.na(aa)]
  if (!length(aa)) return(NA_real_)
  if (length(aa) == 1) return(aa)
  # If the paired ages differ by <=1 year, use the rounded mean; otherwise use the mode.
  if (abs(diff(range(aa))) <= 1) return(round(mean(aa), 1))
  tbl <- sort(table(round(aa, 1)), decreasing = TRUE)
  as.numeric(names(tbl)[1])
}

bin_from_age <- function(a){
  if (is.na(a)) return(NA_character_)
  if (is_young(a)) return("young")
  if (is_old(a))   return("old")
  NA_character_
}

# ---- General utilities ----
`%||%` <- function(a,b) if(!is.null(a)) a else b
as_chr <- function(x){ x <- as.character(x); x[is.na(x)] <- NA_character_; x }
norm_id <- function(x){
  x <- tolower(as_chr(x))
  x <- str_replace_all(x, "[^a-z0-9]+","_")
  x <- str_replace_all(x, "_+","_")
  x <- str_replace(x, "^_|_$","")
  ifelse(x=="", NA_character_, x)
}

# Generic age parser kept for compatibility; not used directly below.
parse_age_strict <- function(...){
  s <- tolower(paste(..., sep=" | "))
  s <- gsub("\\b(gsm|gse)[0-9]+"," ", s)
  yr <- stringr::str_match(s, "\\b([0-9]{1,3})\\s*(years?|yrs?|yr)\\b")[,2]
  ag <- stringr::str_match(s, "age[^0-9]*([0-9]{1,3})")[,2]
  mo <- stringr::str_match(s, "\\b([0-9]{1,3})\\s*(months?|mo)\\b")[,2]
  out <- rep(NA_real_, length(s))
  if (!all(is.na(yr))) {
    yy <- suppressWarnings(as.numeric(yr)); out[!is.na(yy)] <- yy[!is.na(yy)]
  }
  fill <- is.na(out) & !is.na(ag); if (any(fill)) out[fill] <- suppressWarnings(as.numeric(ag[fill]))
  fill <- is.na(out) & !is.na(mo); if (any(fill)) out[fill] <- suppressWarnings(as.numeric(mo[fill]))/12
  out[out < 0 | out > 120] <- NA_real_
  out
}

# ---- Check pairing structure by patient ----
assert_pairs_strict <- function(meta, gse){
  ct <- meta %>%
    dplyr::filter(tissue %in% c("Tumor","Normal")) %>%
    dplyr::count(patient_id, tissue) %>%
    tidyr::pivot_wider(names_from=tissue, values_from=n, values_fill=0)
  prob <- ct %>% dplyr::filter(Normal>1 | Tumor>1 | (Normal==0 & Tumor==1) | (Normal==1 & Tumor==0))
  if (nrow(prob)==0){
    message(glue("[{gse}] Pareo estricto OK (máx 1N/1T por paciente)."))
  } else {
    message(glue("[{gse}] Pareo NO estricto para estos patient_id:"))
    apply(prob, 1, function(r) message(glue("  - {r[['patient_id']]} (N={r[['Normal']]}, T={r[['Tumor']]})")))
  }
}

# ---- Dataset-specific metadata parsers ----

# ---- GSE115002 parser ----
build_meta_115002 <- function(P){
  norm_simple <- function(x){
    x <- tolower(as_chr(x))
    x <- stringr::str_replace_all(x, "[^a-z0-9]+","_")
    x <- stringr::str_replace_all(x, "_+","_")
    x <- stringr::str_replace(x, "^_|_$","")
    ifelse(x=="", NA_character_, x)
  }
  # Remove final T/N suffixes, with or without separator, and tumor/normal suffixes.
  strip_nt_root <- function(tt){
    tt <- tolower(as_chr(tt))
    tt <- sub("\\s*\\((?:tumou?r|normal)\\)\\s*$", "", tt, perl=TRUE)
    tt <- sub("\\s*[-_ ](?:t|n)\\s*$", "", tt, perl=TRUE)  # "-T", "_N", " T"
    tt <- sub("(?<=\\w)(?:t|n)\\s*$", "", tt, perl=TRUE)   # "22T"/"22N" pegado
    norm_simple(tt)
  }
  src   <- if ("source_name_ch1" %in% names(P)) as_chr(P$source_name_ch1) else rep(NA_character_, nrow(P))
  title <- if ("title"          %in% names(P)) as_chr(P$title)            else as_chr(P$geo_accession)
  tlow  <- tolower(title)
  
  tissue <- dplyr::case_when(
    grepl("[-_ ]n\\s*$", tlow) | grepl("(?<=\\w)n\\s*$", tlow, perl=TRUE) ~ "Normal",
    grepl("[-_ ]t\\s*$", tlow) | grepl("(?<=\\w)t\\s*$", tlow, perl=TRUE) ~ "Tumor",
    grepl("\\(normal\\)\\s*$",  tlow) ~ "Normal",
    grepl("\\(tumou?r\\)\\s*$", tlow) ~ "Tumor",
    grepl("normal", src, ignore.case=TRUE) ~ "Normal",
    grepl("tumou?r", src, ignore.case=TRUE) ~ "Tumor",
    TRUE ~ NA_character_
  )
  
  age <- if ("age:ch1" %in% names(P)) suppressWarnings(readr::parse_number(as_chr(P$`age:ch1`))) else rep(NA_real_, nrow(P))
  age[age < 0 | age > 120] <- NA_real_
  
  pid_digits <- stringr::str_match(tlow, "(\\d{8,12})")[,2]
  pid_root   <- strip_nt_root(title)
  pid        <- ifelse(!is.na(pid_digits), pid_digits, pid_root)
  
  tibble::tibble(
    sample_id  = as_chr(P$geo_accession),
    patient_id = as_chr(pid),
    tissue     = as_chr(tissue),
    age        = as.numeric(age),
    age_bin    = dplyr::case_when(is_young(age)~"young", is_old(age)~"old", TRUE~NA_character_),
    title      = as_chr(title)
  )
}

# ---- GSE75037 parser ----
build_meta_75037 <- function(P){
  sn    <- if ("source_name_ch1" %in% names(P)) as_chr(P$source_name_ch1) else rep(NA_character_, nrow(P))
  title <- if ("title" %in% names(P)) as_chr(P$title) else as_chr(P$geo_accession)
  tissue <- dplyr::case_when(
    grepl("non[- ]?malignant\\s+lung", sn, ignore.case=TRUE) ~ "Normal",
    grepl("lung\\s+cancer", sn, ignore.case=TRUE) ~ "Tumor",
    TRUE ~ NA_character_
  )
  age <- if ("age (yrs):ch1" %in% names(P)) suppressWarnings(as.numeric(P$`age (yrs):ch1`)) else rep(NA_real_, nrow(P))
  age[age < 0 | age > 120] <- NA_real_
  pid <- sub("_(N|T)\\s*$","", title, ignore.case=TRUE)
  pid <- norm_id(pid)
  tibble::tibble(
    sample_id  = as_chr(P$geo_accession),
    patient_id = as_chr(pid),
    tissue     = as_chr(tissue),
    age        = as.numeric(age),
    age_bin    = dplyr::case_when(is_young(age)~"young", is_old(age)~"old", TRUE~NA_character_),
    title      = as_chr(title)
  )
}

# ---- Main: build metadata and summary tables ----
GSE_IDS <- c("GSE75037","GSE115002")
summary_rows <- list()
safe_pick <- function(df, col) if (col %in% names(df)) df[[col]] else 0L

for (gse in GSE_IDS){
  pheno_path <- file.path(ROOT, paste0(gse, "_pheno.tsv"))
  if (!file.exists(pheno_path)) { message("File not found: ", pheno_path, ". Did you run script 01?"); next }
  P <- readr::read_tsv(pheno_path, show_col_types = FALSE)
  stopifnot("geo_accession" %in% names(P))
  
  meta <- switch(gse,
                 "GSE115002" = build_meta_115002(P),
                 "GSE75037"  = build_meta_75037(P))
  
  # Standardize tissue labels
  meta <- meta %>% dplyr::mutate(
    tissue = dplyr::case_when(
      grepl("^normal$", tissue, ignore.case=TRUE) ~ "Normal",
      grepl("^tumor$",  tissue, ignore.case=TRUE) ~ "Tumor",
      TRUE ~ tissue
    ),
    dplyr::across(c(sample_id, patient_id, tissue, title, age_bin), as_chr)
  )
  
  # Check pairing structure and report issues
  assert_pairs_strict(meta, gse)
  
  # Pair counts and pairing issues
  ct0 <- meta %>%
    dplyr::filter(tissue %in% c("Tumor","Normal")) %>%
    dplyr::count(patient_id, tissue) %>%
    tidyr::pivot_wider(names_from = tissue, values_from = n, values_fill=0)
  n_pairs <- sum(ct0$Tumor>0 & ct0$Normal>0)
  issues <- ct0 %>% dplyr::filter(!(Tumor==1 & Normal==1) & !(Tumor==0 & Normal==0))
  if (nrow(issues)){
    out_issues <- file.path(META_DIR, paste0(gse,"_pair_issues.tsv"))
    write.table(issues, file=out_issues, sep="\t", row.names=FALSE, quote=FALSE)
    message(glue("{gse}: pares T/N={n_pairs}  -> Detalle: {out_issues}"))
  } else {
    message(glue("{gse}: pares T/N={n_pairs}"))
  }
  
  # Write metadata table
  out_meta <- file.path(META_DIR, paste0(gse,"_meta.csv"))
  readr::write_csv(meta, out_meta)
  message("  -> meta: ", out_meta)
  
  # Summary by age bin and tissue
  cbt <- meta %>%
    dplyr::filter(tissue %in% c("Normal","Tumor")) %>%
    dplyr::mutate(age_bin = dplyr::case_when(
      is_young(as.numeric(age))~"young",
      is_old(as.numeric(age))~"old",
      TRUE~NA_character_)
    ) %>%
    dplyr::count(age_bin, tissue) %>%
    tidyr::pivot_wider(names_from=c(age_bin, tissue), values_from=n)
  
  young_nrm <- safe_pick(cbt, "young_Normal"); young_tum <- safe_pick(cbt, "young_Tumor")
  old_nrm   <- safe_pick(cbt, "old_Normal");   old_tum  <- safe_pick(cbt, "old_Tumor")
  
  tb_tissue <- table(meta$tissue, useNA="no")
  n_tumor   <- if ("Tumor"  %in% names(tb_tissue)) as.integer(tb_tissue[["Tumor"]])  else 0L
  n_normal  <- if ("Normal" %in% names(tb_tissue)) as.integer(tb_tissue[["Normal"]]) else 0L
  
  summary_rows[[length(summary_rows)+1]] <- tibble::tibble(
    gse=gse, n_total=nrow(meta), n_tumor=n_tumor, n_normal=n_normal, n_pairs=n_pairs,
    young_normal=young_nrm, young_tumor=young_tum, old_normal=old_nrm, old_tumor=old_tum
  )
}

if (length(summary_rows)){
  SUM <- dplyr::bind_rows(summary_rows) %>% dplyr::arrange(gse)
  readr::write_tsv(SUM, file.path(META_DIR, "summary.tsv"))
  message("Summary -> ", file.path(META_DIR, "summary.tsv"))
}

## ---- GSE115002 ID correction and final report ----
META_115 <- file.path(META_DIR, "GSE115002_meta.csv")
if (file.exists(META_115)) {
  m115 <- readr::read_csv(META_115, show_col_types = FALSE) %>%
    dplyr::mutate(
      patient_id = as_chr(patient_id),
      # Remove final t/n attached to a numeric patient ID (e.g., 22t -> 22).
      patient_id = ifelse(grepl("(?<=\\d)[tn]$", patient_id, perl=TRUE),
                          sub("(?<=\\d)[tn]$", "", patient_id, perl=TRUE),
                          patient_id)
    )
  readr::write_csv(m115, META_115)
  
  # Report
  m115$age <- suppressWarnings(as.numeric(m115$age))
  ct <- m115 %>% dplyr::filter(tissue %in% c("Tumor","Normal")) %>%
    dplyr::count(patient_id, tissue) %>% tidyr::pivot_wider(names_from=tissue, values_from=n, values_fill=0)
  n_pairs <- sum(ct$Tumor>0 & ct$Normal>0)
  only_t  <- ct %>% dplyr::filter(Tumor>0, Normal==0) %>% dplyr::pull(patient_id)
  only_n  <- ct %>% dplyr::filter(Normal>0, Tumor==0) %>% dplyr::pull(patient_id)
  
  pick_age <- function(v){
    v <- na.omit(as.numeric(v)); if (!length(v)) return(NA_real_)
    if (length(v) == 1) return(v)
    as.numeric(names(sort(table(round(v,3)), decreasing=TRUE))[1])
  }
  
  pairs_by_bin <- m115 %>%
    dplyr::filter(tissue %in% c("Tumor","Normal")) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(
      has_t = any(tissue=="Tumor"),
      has_n = any(tissue=="Normal"),
      age_patient = pick_age(age),
      .groups="drop"
    ) %>%
    dplyr::mutate(
      age_bin = dplyr::case_when(
        !is.na(age_patient) & age_patient < 60 ~ "young",
        !is.na(age_patient) & age_patient >= 60 ~ "old",
        TRUE ~ NA_character_
      ),
      is_pair = has_t & has_n
    ) %>%
    dplyr::filter(is_pair) %>%
    dplyr::count(age_bin, name="n_pairs") %>%
    dplyr::arrange(factor(age_bin, levels=c("young","old","NA")))
  
  message("\n===== FINAL REPORT GSE115002 =====")
  message("Total paired patients: ", n_pairs)
  message("Pairs by age bin:"); print(pairs_by_bin)
  message("Tumor-only patients: ", ifelse(length(only_t), paste(only_t, collapse=", "), "(none)"))
  message("Normal-only patients: ", ifelse(length(only_n), paste(only_n, collapse=", "), "(none)"))
  message("================================================\n")
}

# ---- GSE75037 final report ----
META_750 <- file.path(META_DIR, "GSE75037_meta.csv")
if (file.exists(META_750)) {
  m750 <- readr::read_csv(META_750, show_col_types = FALSE)
  
  m750$age <- suppressWarnings(as.numeric(m750$age))
  ct <- m750 %>% dplyr::filter(tissue %in% c("Tumor","Normal")) %>%
    dplyr::count(patient_id, tissue) %>% tidyr::pivot_wider(names_from=tissue, values_from=n, values_fill=0)
  n_pairs <- sum(ct$Tumor>0 & ct$Normal>0)
  only_t  <- ct %>% dplyr::filter(Tumor>0, Normal==0) %>% dplyr::pull(patient_id)
  only_n  <- ct %>% dplyr::filter(Normal>0, Tumor==0) %>% dplyr::pull(patient_id)
  
  pick_age <- function(v){
    v <- na.omit(as.numeric(v)); if (!length(v)) return(NA_real_)
    if (length(v) == 1) return(v)
    as.numeric(names(sort(table(round(v,3)), decreasing=TRUE))[1])
  }
  
  pairs_by_bin <- m750 %>%
    dplyr::filter(tissue %in% c("Tumor","Normal")) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(
      has_t = any(tissue=="Tumor"),
      has_n = any(tissue=="Normal"),
      age_patient = pick_age(age),
      .groups="drop"
    ) %>%
    dplyr::mutate(
      age_bin = dplyr::case_when(
        !is.na(age_patient) & age_patient < 60 ~ "young",
        !is.na(age_patient) & age_patient >= 60 ~ "old",
        TRUE ~ NA_character_
      ),
      is_pair = has_t & has_n
    ) %>%
    dplyr::filter(is_pair) %>%
    dplyr::count(age_bin, name="n_pairs") %>%
    dplyr::arrange(factor(age_bin, levels=c("young","old","NA")))
  
  message("\n===== FINAL REPORT GSE75037 =====")
  message("Total paired patients: ", n_pairs)
  message("Pairs by age bin:"); print(pairs_by_bin)
  message("Tumor-only patients: ", ifelse(length(only_t), paste(only_t, collapse=", "), "(none)"))
  message("Normal-only patients: ", ifelse(length(only_n), paste(only_n, collapse=", "), "(none)"))
  message("================================================\n")
}

# ---- Strict 1N-1T pair tables using consensus age/bin ----

make_pairs_for_gse <- function(gse){
  f_meta <- file.path(META_DIR, paste0(gse, "_meta.csv"))
  if (!file.exists(f_meta)) {
    message("File not found: ", f_meta, " (¿falló el paso de meta?).")
    return(invisible(NULL))
  }
  
  M <- readr::read_csv(f_meta, show_col_types = FALSE) %>%
    dplyr::mutate(
      sample_id  = as_chr(sample_id),
      patient_id = as_chr(patient_id),
      tissue     = dplyr::case_when(
        grepl("^normal$", tissue, ignore.case=TRUE) ~ "Normal",
        grepl("^tumor$",  tissue, ignore.case=TRUE) ~ "Tumor",
        TRUE ~ as_chr(tissue)
      ),
      age     = suppressWarnings(as.numeric(age)),
      age_bin = as_chr(age_bin),
      title   = as_chr(title)
    ) %>%
    dplyr::filter(
      !is.na(sample_id),
      !is.na(patient_id),
      tissue %in% c("Normal","Tumor")
    )
  
  if (nrow(M) == 0){
    message("[", gse, "] metadata is empty after filtering.")
    return(invisible(NULL))
  }
  
  # Count by patient_id and tissue
  ct <- M %>%
    dplyr::count(patient_id, tissue) %>%
    tidyr::pivot_wider(names_from = tissue, values_from = n, values_fill = 0)
  
  issues <- ct %>% dplyr::filter(Normal != 1 | Tumor != 1)
  
  # Strict 1N-1T pairs
  keep_ids <- ct %>%
    dplyr::filter(Normal == 1, Tumor == 1) %>%
    dplyr::pull(patient_id)
  
  N <- M %>%
    dplyr::filter(tissue == "Normal", patient_id %in% keep_ids) %>%
    dplyr::select(patient_id, N_sample_id = sample_id, N_age = age, N_age_bin = age_bin)
  
  Tt <- M %>%
    dplyr::filter(tissue == "Tumor",  patient_id %in% keep_ids) %>%
    dplyr::select(patient_id, T_sample_id = sample_id, T_age = age, T_age_bin = age_bin)
  
  pairs <- dplyr::inner_join(N, Tt, by = "patient_id") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      age_consensus = age_consensus(N_age, T_age),
      bin_consensus = bin_from_age(age_consensus),
      agree_N = (!is.na(N_age) & bin_from_age(N_age) == bin_consensus),
      agree_T = (!is.na(T_age) & bin_from_age(T_age) == bin_consensus)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      patient_id, N_sample_id, T_sample_id,
      N_age, T_age, age_consensus, bin_consensus,
      N_age_bin, T_age_bin, agree_N, agree_T
    )
  
  # Write pair table
  out_pairs <- file.path(OUT_DIR, paste0(gse, "_pairs.tsv"))
  readr::write_tsv(pairs, out_pairs)
  
  # Write pairing issues if present
  if (nrow(issues)){
    out_issues <- file.path(OUT_DIR, paste0(gse, "_pair_issues.tsv"))
    readr::write_tsv(issues, out_issues)
  }
  
  # Summary by consensus bin
  sum_bin <- pairs %>%
    dplyr::count(bin_consensus, name = "n_pairs") %>%
    dplyr::arrange(factor(bin_consensus, levels = c("young","old", NA)))
  
  message("\n== ", gse, " ==")
  message("Strict 1N+1T pairs: ", nrow(pairs))
  if (nrow(issues)){
    message("Pairing issues: ", file.path(OUT_DIR, paste0(gse, "_pair_issues.tsv")))
  } else {
    message("No pairing issues.")
  }
  message("By consensus age bin:"); print(sum_bin, n = Inf)
  message("Pair table -> ", out_pairs)
  
  invisible(list(pairs = pairs, issues = issues, sum_bin = sum_bin))
}

GSE_IDS_PAIRS <- c("GSE115002","GSE75037")
res_pairs <- lapply(GSE_IDS_PAIRS, make_pairs_for_gse)

message("\nDone. Pair tables are in: ", OUT_DIR)
