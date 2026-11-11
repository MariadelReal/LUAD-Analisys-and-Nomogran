#########################################
# LUAD age-bin workflow
# ssGSEA in ENTREZ space followed by paired pathway-level limma
#########################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(limma)
  library(GSVA)
  library(Biobase)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Avoid common conflicts between dplyr and Bioconductor packages.
select     <- dplyr::select
mutate     <- dplyr::mutate
distinct   <- dplyr::distinct
filter     <- dplyr::filter
count      <- dplyr::count
arrange    <- dplyr::arrange
inner_join <- dplyr::inner_join
transmute  <- dplyr::transmute

# Set working directory.
# The expected structure is:
#   prepro/GSE*_expr_SYMBOL.tsv
#   meta/GSE*_meta.csv
# Outputs are written to:
#   ssgsea/
ROOT     <- getwd()
PREPRO   <- file.path(ROOT, "prepro")
META_DIR <- file.path(ROOT, "meta")
OUT_DIR  <- file.path(ROOT, "ssgsea")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GSE_IDS <- c("GSE115002", "GSE75037")
BINS    <- c("young", "old")
FDR_CUT <- 0.05
MIN_SZ  <- 10
MAX_SZ  <- 2000

as_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x
}

# ---------- Read expression matrices ----------
read_expr_symbol <- function(gse) {
  f <- file.path(PREPRO, paste0(gse, "_expr_SYMBOL.tsv"))
  if (!file.exists(f)) stop("Expression file not found: ", f)

  df <- readr::read_tsv(f, show_col_types = FALSE)
  if (!"GENE" %in% colnames(df)) stop("File does not contain a 'GENE' column: ", f)

  rn <- df$GENE
  M  <- as.matrix(data.frame(row.names = rn, df[, -1, drop = FALSE], check.names = FALSE))
  storage.mode(M) <- "double"
  M
}

# ---------- Map SYMBOL expression matrix to ENTREZ identifiers ----------
expr_symbol_to_entrez <- function(M) {
  syms <- rownames(M)

  map <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(syms),
      keytype = "SYMBOL",
      columns = c("ENTREZID")
    )
  )

  map <- map %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  keep <- syms %in% map$SYMBOL
  M2   <- M[keep, , drop = FALSE]
  ent  <- setNames(map$ENTREZID, map$SYMBOL)[rownames(M2)]

  M3 <- limma::avereps(M2, ID = ent)
  rownames(M3) <- as.character(rownames(M3))
  M3
}

# ---------- Strict paired-sample check ----------
assert_pairs_exact <- function(mm_bin, gse, bin_label) {
  chk <- mm_bin %>%
    count(patient_id, tissue) %>%
    tidyr::pivot_wider(names_from = tissue, values_from = n, values_fill = 0)

  if (!"Normal" %in% names(chk)) chk$Normal <- 0L
  if (!"Tumor"  %in% names(chk)) chk$Tumor  <- 0L

  chk <- chk %>% mutate(ok = (Normal == 1L & Tumor == 1L))

  if (!all(chk$ok)) {
    bad <- chk %>%
      filter(!ok) %>%
      transmute(patient_id, Normal = as.integer(Normal), Tumor = as.integer(Tumor))

    msg <- paste0(
      "\n[", gse, " / ", bin_label, "] Non-strict pairing detected. Required: 1 Normal + 1 Tumor per patient.\n",
      paste(
        apply(
          bad,
          1,
          function(r) sprintf("  - %s (N=%s, T=%s)", r[["patient_id"]], r[["Normal"]], r[["Tumor"]])
        ),
        collapse = "\n"
      ),
      "\nFix meta/*.csv or exclude problematic patients before rerunning."
    )
    stop(msg, call. = FALSE)
  }

  invisible(TRUE)
}

# ---------- Robust msigdbr wrapper ----------
.msigdbr_fetch <- function(species = "Homo sapiens", collection, sub = NULL) {
  res <- try(msigdbr(species = species, collection = collection, subcollection = sub), silent = TRUE)
  if (!inherits(res, "try-error") && is.data.frame(res) && nrow(res)) return(res)

  res <- try(msigdbr(species = species, collection = collection, subcategory = sub), silent = TRUE)
  if (!inherits(res, "try-error") && is.data.frame(res) && nrow(res)) return(res)

  res <- try(msigdbr(species = species, category = collection, subcategory = sub), silent = TRUE)
  if (!inherits(res, "try-error") && is.data.frame(res) && nrow(res)) return(res)

  res2 <- try(msigdbr(species = species, collection = collection), silent = TRUE)
  if (!inherits(res2, "try-error") && is.data.frame(res2) && nrow(res2)) {
    if (!is.null(sub)) {
      key <- intersect(c("gs_subcollection","gs_subcat","gs_subcategory"),names(res2))
      if (!length(key)) return(tibble())
      res2 <- res2[res2[[key[1]]] == sub,,drop=FALSE]
    }
    if (nrow(res2)) return(res2)
  }

  tibble()
}

# ---------- Build MSigDB gene sets in ENTREZ space ----------
build_msig_entrez_sets <- function(collection, subcollection = NULL, species = "Homo sapiens") {
  tbl <- .msigdbr_fetch(species = species, collection = collection, sub = subcollection)
  if (!nrow(tbl)) return(list())
  tbl <- tbl %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  if (!nrow(tbl)) return(list())

  map <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(tbl$gene_symbol),
      keytype = "SYMBOL",
      columns = c("ENTREZID")
    )
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  tblE <- inner_join(tbl, map, by = c("gene_symbol" = "SYMBOL")) %>%
    select(gs_name, ENTREZID) %>%
    distinct()

  sets <- split(tblE$ENTREZID, tblE$gs_name)
  lapply(sets, function(v) unique(as.character(v)))
}

# ---------- Gene sets ----------
hallmark_sets <- build_msig_entrez_sets(collection = "H")
kegg_sets     <- build_msig_entrez_sets(collection = "C2", subcollection = "CP:KEGG")

if (!length(kegg_sets)) {
  message("KEGG was not available through msigdbr. Using Reactome fallback (C2:CP:REACTOME).")
  path_sets  <- build_msig_entrez_sets(collection = "C2", subcollection = "CP:REACTOME")
  PATH_LABEL <- "Reactome"
} else {
  path_sets  <- kegg_sets
  PATH_LABEL <- "KEGG"
}

hallmark_sets <- hallmark_sets[vapply(hallmark_sets, length, integer(1)) >= MIN_SZ &
                                 vapply(hallmark_sets, length, integer(1)) <= MAX_SZ]
path_sets     <- path_sets[vapply(path_sets, length, integer(1)) >= MIN_SZ &
                             vapply(path_sets, length, integer(1)) <= MAX_SZ]

# ---------- Mapping audit ----------
report_mapping <- function(syms, map_df, tag, out_dir) {
  map_dir <- file.path(out_dir, "mapping")
  dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)

  mapped_syms <- unique(map_df$SYMBOL)
  unmapped    <- setdiff(unique(syms), mapped_syms)

  readr::write_tsv(
    dplyr::distinct(map_df, SYMBOL, ENTREZID, .keep_all = FALSE),
    file.path(map_dir, paste0(tag, "_SYMBOL_to_ENTREZ_mapped.tsv"))
  )

  readr::write_lines(
    unmapped,
    file.path(map_dir, paste0(tag, "_SYMBOL_unmapped.txt"))
  )

  tot <- length(unique(syms))
  ok  <- length(mapped_syms)
  message(sprintf(
    "[%s] SYMBOL to ENTREZ mapping: %d/%d (%.1f%%). Details: %s",
    tag, ok, tot, 100 * ok / tot, map_dir
  ))

  invisible(list(mapped = mapped_syms, unmapped = unmapped))
}

# ---------- Run ssGSEA and paired pathway-level differential analysis ----------
run_ssgsea_one <- function(E_entrez, meta_gse, bin_label, setlist, setname, out_prefix, gse_id) {
  mm0 <- meta_gse %>%
    transmute(
      sample_id  = as_chr(sample_id),
      patient_id = as_chr(patient_id),
      tissue     = as_chr(tissue),
      age_bin    = as_chr(age_bin)
    ) %>%
    filter(tissue %in% c("Normal", "Tumor"), age_bin == bin_label) %>%
    distinct(sample_id, .keep_all = TRUE)

  if (nrow(mm0) == 0) {
    message("  [", bin_label, " / ", setname, "] 0 samples; skipped.")
    stop("El contraste está incompleto. Comprueba los identificadores antes de continuar.")
  }

  ids <- base::intersect(colnames(E_entrez), mm0$sample_id)
  if (length(ids) < 4) {
    message("  [", bin_label, " / ", setname, "] <4 columns after expression/metadata intersection; skipped.")
    stop("El contraste está incompleto. Comprueba los identificadores antes de continuar.")
  }

  Emm <- E_entrez[, ids, drop = FALSE]
  mm  <- mm0[match(ids, mm0$sample_id), , drop = FALSE]
  stopifnot(identical(as.character(mm$sample_id), colnames(Emm)))

  assert_pairs_exact(mm, gse = gse_id, bin_label = bin_label)

  sizes <- vapply(setlist, function(g) length(base::intersect(rownames(Emm), g)), integer(1))
  med <- unname(stats::median(sizes))
  q25 <- unname(stats::quantile(sizes, 0.25, names = FALSE))
  q75 <- unname(stats::quantile(sizes, 0.75, names = FALSE))
  mx  <- max(sizes)

  cat(sprintf(
    "  [%s / %s] gene set coverage — median: %.0f  P25: %.0f  P75: %.0f  max: %d\n",
    bin_label, setname, med, q25, q75, mx
  ))

  keep <- which(sizes >= MIN_SZ & sizes <= MAX_SZ)
  if (!length(keep)) {
    message("  [", bin_label, " / ", setname, "] 0 sets after size filtering; skipped.")
    stop("No quedan conjuntos para ",gse_id," / ",bin_label," / ",setname)
  }

  set_f <- setlist[keep]

  Emm <- as.matrix(Emm)
  storage.mode(Emm) <- "double"
  ssg <- NULL

  # Prefer current GSVA API when available, while keeping compatibility with older releases.
  if (exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    sp_fun <- get("ssgseaParam", asNamespace("GSVA"))
    sp <- sp_fun(
      exprData = Emm,
      geneSets = set_f,
      normalize = TRUE,
      minSize = MIN_SZ,
      maxSize = MAX_SZ
    )
    ssg <- GSVA::gsva(sp)

    if (inherits(ssg, "ExpressionSet")) {
      ssg <- Biobase::exprs(ssg)
    } else if (!is.matrix(ssg) && requireNamespace("SummarizedExperiment", quietly = TRUE) &&
               inherits(ssg, "SummarizedExperiment")) {
      ssg <- as.matrix(SummarizedExperiment::assay(ssg))
    }
  } else {
    ssg <- tryCatch(
      GSVA::gsva(Emm, set_f, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE),
      error = function(e) NULL
    )

    if (is.null(ssg)) {
      eset <- Biobase::ExpressionSet(assayData = Emm)
      ssg <- tryCatch(
        GSVA::gsva(eset, set_f, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE),
        error = function(e) NULL
      )

      if (!is.null(ssg)) {
        if (inherits(ssg, "ExpressionSet")) {
          ssg <- Biobase::exprs(ssg)
        } else if (!is.matrix(ssg) && requireNamespace("SummarizedExperiment", quietly = TRUE) &&
                   inherits(ssg, "SummarizedExperiment")) {
          ssg <- as.matrix(SummarizedExperiment::assay(ssg))
        }
      }
    }
  }

  if (!is.matrix(ssg)) stop("Could not obtain a GSVA score matrix; check GSVA version/API.")

  # Save ssGSEA scores.
  out_mat <- paste0(out_prefix, "__ssgsea_", setname, "_", bin_label, ".tsv")
  readr::write_tsv(data.frame(Pathway = rownames(ssg), ssg, check.names = FALSE), out_mat)

  # Paired limma on pathway scores. With the intercept design, TissueTumor = Tumor - Normal.
  mm$pair_id <- droplevels(factor(mm$patient_id))
  mm$Tissue  <- droplevels(factor(mm$tissue, levels = c("Normal", "Tumor")))
  design <- model.matrix(~ pair_id + Tissue, mm)

  if (nrow(design) != ncol(ssg)) {
    stop("Design matrix and ssGSEA score matrix are not aligned: nrow(design) != ncol(ssg).")
  }

  if (!identical(colnames(ssg), mm$sample_id)) {
    ssg <- ssg[, mm$sample_id, drop = FALSE]
  }

  if (!"TissueTumor" %in% colnames(design)) {
    stop("Coefficient 'TissueTumor' was not found in the design matrix. Check Tissue levels.")
  }

  fit  <- limma::lmFit(ssg, design)
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  tt <- limma::topTable(fit2, coef = "TissueTumor", number = Inf, sort.by = "P")
  tt$FDR <- if ("adj.P.Val" %in% names(tt)) tt$adj.P.Val else p.adjust(tt$P.Value, "BH")
  tt$dir <- ifelse(tt$logFC > 0, "UP_in_Tumor", "DOWN_in_Tumor")

  out_de <- paste0(out_prefix, "__ssgsea_DEA_", setname, "_", bin_label, ".tsv")
  readr::write_tsv(tibble::rownames_to_column(tt, "Pathway"), out_de)

  message(sprintf(
    "  [%s / %s] %d pathways with FDR <= %.2f -> %s",
    bin_label, setname, sum(tt$FDR <= FDR_CUT, na.rm = TRUE), FDR_CUT, out_de
  ))

  invisible(list(scores = ssg, dea = tt))
}

# ---------- Main loop ----------
manifest <- file.path(OUT_DIR,"collections_used.txt")
if (file.exists(manifest)) unlink(manifest)
stopifnot(length(hallmark_sets)>0L,length(path_sets)>0L)
for (gse in GSE_IDS) {
  message("\n== ", gse, " ==")

  E_sym <- read_expr_symbol(gse)

  syms <- rownames(E_sym)
  raw_map <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(syms),
      keytype = "SYMBOL",
      columns = c("ENTREZID")
    )
  ) %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)

  report_mapping(syms, raw_map, tag = gse, out_dir = OUT_DIR)

  M <- readr::read_csv(file.path(META_DIR, paste0(gse, "_meta.csv")), show_col_types = FALSE) %>%
    mutate(
      sample_id  = as_chr(sample_id),
      patient_id = as_chr(patient_id),
      tissue     = as_chr(tissue),
      age_bin    = as_chr(age_bin)
    )

  cat(sprintf("Genes in SYMBOL matrix: %d | Samples: %d\n", nrow(E_sym), ncol(E_sym)))
  cat("Example genes: ", paste(head(rownames(E_sym), 8), collapse = ", "), "\n")

  E_entrez <- expr_symbol_to_entrez(E_sym)
  cat(sprintf("Genes mapped to ENTREZ after collapsing: %d\n", nrow(E_entrez)))

  out_prefix <- file.path(OUT_DIR, paste0(gse, "_"))

  for (bin in BINS) {
    invisible(run_ssgsea_one(E_entrez, M, bin, hallmark_sets, "H", out_prefix, gse_id = gse))
    invisible(run_ssgsea_one(E_entrez, M, bin, path_sets, PATH_LABEL, out_prefix, gse_id = gse))
  }
}

message("Contrast used: TissueTumor (Tumor - Normal). logFC > 0 means higher pathway score in tumor.")
cat("\nDone. Results written to: ", OUT_DIR, "\n")

writeLines(c("H",PATH_LABEL),file.path(OUT_DIR,"collections_used.txt"))
