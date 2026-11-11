#########################################
# LUAD age-bin transcriptomic analysis
# Collapse probe-level expression matrices to gene symbols
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(readr)
  library(dplyr)
  library(stringr)
  library(limma)
  library(tidyr)
  library(tibble)
})

# Set working directory to the folder containing GSE*_series.tsv files.
# Example: setwd("/path/to/Two_Agebins/Series_Matrix")
ROOT <- getwd()
SERIES_DIR <- ROOT
OUT_DIR <- file.path(ROOT, "prepro")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Project directory: ", ROOT)
message("Output directory : ", OUT_DIR)

# -----------------------------
# Helper functions
# -----------------------------

guess_symbol_column <- function(fdat) {
  candidates <- colnames(fdat)
  priority <- c(
    "Gene.symbol", "Gene Symbol", "SYMBOL", "Symbol", "GENE_SYMBOL",
    "gene_symbol", "GeneName", "GENE", "GENE_SYMBOLS"
  )

  for (nm in priority) {
    if (nm %in% candidates) return(nm)
  }

  alt <- grep("symbol", candidates, ignore.case = TRUE, value = TRUE)
  if (length(alt)) return(alt[1])

  NA_character_
}

clean_symbols <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- as.character(x)
  x[is.na(x)] <- NA_character_

  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)

    # If multiple symbols are listed for one probe, keep the first one.
    if (toupper(trimws(s)) %in% c("N/A","NA","---","?","")) return(NA_character_)
    s <- strsplit(s, "\\s*[/|;,+]+\\s*")[[1]][1]
    s <- toupper(trimws(s))

    if (s %in% c("", "NA", "---", "?", "N/A")) return(NA_character_)
    s
  }, character(1), USE.NAMES = FALSE)
}

needs_log2 <- function(M) {
  rng <- range(M, na.rm = TRUE)
  !(rng[2] <= 20 && rng[1] >= -5)
}

build_symbol_map_all_gpls <- function(gse_id) {
  esets <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  if (!is.list(esets)) esets <- list(esets)

  maps <- lapply(esets, function(es) {
    fd <- Biobase::fData(es)
    if (is.null(fd) || nrow(fd) == 0) return(NULL)

    sym_col <- guess_symbol_column(fd)
    if (is.na(sym_col)) return(NULL)

    data.frame(
      PROBE = rownames(fd),
      SYMBOL = clean_symbols(fd[[sym_col]]),
      GPL = tryCatch(annotation(es), error = function(e) NA_character_),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  maps <- Filter(Negate(is.null), maps)
  if (!length(maps)) {
    stop("No gene symbol column was found in fData() for ", gse_id)
  }

  MAP <- dplyr::bind_rows(maps)

  # If a probe appears more than once, keep the first non-missing symbol.
  MAP <- MAP %>%
    arrange(factor(GPL)) %>%
    group_by(PROBE) %>%
    summarise(
      SYMBOL = {
        s <- SYMBOL[!is.na(SYMBOL) & nzchar(SYMBOL)]
        if (length(s)) s[1] else NA_character_
      },
      GPL = GPL[1],
      .groups = "drop"
    )

  MAP
}

# -----------------------------
# Process one dataset
# -----------------------------

process_one_gse <- function(gse_id) {
  message("\n== ", gse_id, " ==")
  log_applied <- FALSE

  series_file <- file.path(SERIES_DIR, paste0(gse_id, "_series.tsv"))
  if (!file.exists(series_file)) stop("File not found: ", series_file)

  df <- readr::read_tsv(series_file, show_col_types = FALSE)
  if (!"FEATURE" %in% colnames(df)) {
    stop("The expression file must contain a FEATURE column: ", series_file)
  }

  feature_ids <- df$FEATURE
  M <- as.matrix(data.frame(row.names = feature_ids, df[, -1, drop = FALSE], check.names = FALSE))

  n_features_in <- nrow(M)
  n_samples_in <- ncol(M)

  # Build probe-to-symbol map from GEO platform annotations.
  MAP <- build_symbol_map_all_gpls(gse_id)

  common <- intersect(rownames(M), MAP$PROBE)
  if (length(common) < 1000) {
    warning("Few probes matched fData() annotation for ", gse_id, ": ", length(common),
            " of ", n_features_in)
  }

  M <- M[common, , drop = FALSE]
  MAPc <- MAP[match(common, MAP$PROBE), , drop = FALSE]

  if (needs_log2(M)) {
    message("Applying log2(x + 1) transformation")
    M <- log2(M + 1)
    log_applied <- TRUE
  }

  symbols <- MAPc$SYMBOL
  keep <- !is.na(symbols) & nzchar(symbols)
  M_keep <- M[keep, , drop = FALSE]
  symbols_keep <- symbols[keep]

  if (!nrow(M_keep)) {
    stop("No probes with valid gene symbols remained after filtering for ", gse_id)
  }

  # Collapse multiple probes mapping to the same symbol by average expression.
  M_gene <- limma::avereps(M_keep, ID = symbols_keep)
  stopifnot(nrow(M_gene) > 0, ncol(M_gene) == ncol(M))

  expr_out <- file.path(OUT_DIR, paste0(gse_id, "_expr_SYMBOL.tsv"))
  map_out <- file.path(OUT_DIR, paste0(gse_id, "_symbol_map.tsv"))
  audit_out <- file.path(OUT_DIR, paste0(gse_id, "_prepro_audit.txt"))

  readr::write_tsv(
    data.frame(GENE = rownames(M_gene), M_gene, check.names = FALSE),
    expr_out
  )
  readr::write_tsv(MAP, map_out)

  dropped_no_map <- setdiff(feature_ids, MAP$PROBE)
  n_no_map <- length(dropped_no_map)
  n_common <- length(common)
  n_with_symbol <- sum(keep)

  cat(
    sprintf("GSE                       : %s\n", gse_id),
    sprintf("Features in TSV           : %d\n", n_features_in),
    sprintf("Samples in TSV            : %d\n", n_samples_in),
    sprintf("Probes with fData() map   : %d\n", n_common),
    sprintf("Probes with valid SYMBOL  : %d\n", n_with_symbol),
    sprintf("Genes after collapse      : %d\n", nrow(M_gene)),
    sprintf("Probes without map        : %d\n", n_no_map),
    sprintf("Log2 transformation       : %s\n", if (log_applied) "yes" else "no"),
    file = audit_out,
    sep = ""
  )

  message("Saved expression matrix: ", expr_out)
  message("Saved probe-symbol map : ", map_out)
  message("Saved audit report     : ", audit_out)

  invisible(list(
    n_genes = nrow(M_gene),
    n_samples = ncol(M_gene),
    out = expr_out,
    n_features_in = n_features_in,
    n_common = n_common,
    n_with_symbol = n_with_symbol,
    n_no_map = n_no_map
  ))
}

# -----------------------------
# Run all available GSE datasets
# -----------------------------

series_files <- list.files(SERIES_DIR, pattern = "^GSE\\d+_series\\.tsv$", full.names = FALSE)
GSE_IDS <- sub("_series\\.tsv$", "", series_files)

if (!length(GSE_IDS)) {
  stop("No GSE*_series.tsv files were found in: ", SERIES_DIR)
}

res <- lapply(GSE_IDS, process_one_gse)

summary_table <- tibble::tibble(
  gse = GSE_IDS,
  n_genes = sapply(res, `[[`, "n_genes"),
  n_samples = sapply(res, `[[`, "n_samples"),
  features_in = sapply(res, `[[`, "n_features_in"),
  probes_mapped = sapply(res, `[[`, "n_common"),
  with_symbol = sapply(res, `[[`, "n_with_symbol"),
  probes_no_map = sapply(res, `[[`, "n_no_map"),
  file = sapply(res, `[[`, "out")
)

summary_file <- file.path(OUT_DIR, "prepro_summary.tsv")
readr::write_tsv(summary_table, summary_file)

message("\nPreprocessing summary saved to: ", summary_file)
print(summary_table, n = Inf)
