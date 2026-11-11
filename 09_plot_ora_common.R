#########################################
# LUAD age-bin analysis
# ORA for concordant DEG lists
#########################################

#Nota: esta es la corrección más actual, incluye KEGG legacy
set.seed(123L)

required_pkgs <- c("readr", "dplyr", "tibble", "ggplot2", "clusterProfiler",
                   "org.Hs.eg.db", "msigdbr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Paquetes faltantes: ", paste(missing_pkgs, collapse = ", "),
       ". Instalalos antes de continuar.", call. = FALSE)
}

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
DEA_DIR <- file.path(ROOT, "dea")
FIG_DIR <- file.path(ROOT, "figs")
OUT_DIR <- file.path(ROOT, "enrich_common")
REF_DIR <- file.path(OUT_DIR, "reference_sets")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

GSES <- c("GSE115002", "GSE75037")
BINS <- c("young", "old")
DIRECTIONS <- c("UP", "DOWN")
FDR_CUT <- 0.05
LFC_MIN <- 1
TOP_N <- 5L
EXPECTED_MSIGDB <- "2026.1.Hs"

# Un solo lugar para elegir colecciones; no usar CP:KEGG (obsoleta).
COLLECTIONS <- c(Reactome = "CP:REACTOME", KEGG = "CP:KEGG_LEGACY")

# Registro de ejecucion y entorno.
writeLines(c(
  sprintf("ROOT: %s", ROOT),
  sprintf("Run at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Thresholds: FDR <= %.2f and abs(logFC) >= %.1f", FDR_CUT, LFC_MIN),
  sprintf("Expected MSigDB: %s", EXPECTED_MSIGDB),
  sprintf("Collections: %s", paste(names(COLLECTIONS), COLLECTIONS, collapse = "; ")),
  "",
  capture.output(sessionInfo())
), file.path(OUT_DIR, "09_session_and_parameters.txt"))

dea_path <- function(gse, bin) {
  candidates <- file.path(DEA_DIR,
                          c(sprintf("%s__TvN_%s.tsv", gse, bin),
                            sprintf("%s___TvN_%s.tsv", gse, bin)))
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) stop("Falta DEA para ", gse, " / ", bin,
                              ". Busque en: ", DEA_DIR, call. = FALSE)
  if (length(existing) > 1L && length(unique(unname(tools::md5sum(existing)))) > 1L) {
    stop("Hay dos tablas DEA distintas para ", gse, "/", bin, ": ", paste(existing, collapse=", "))
  }
  existing[1]
}

load_dea <- function(gse, bin) {
  path <- dea_path(gse, bin)
  tab <- readr::read_tsv(path, show_col_types = FALSE)
  if (!"GENE" %in% names(tab)) {
    alt <- intersect(c("Gene", "SYMBOL", "symbol", "gene", "ID"), names(tab))
    if (!length(alt)) stop("No hay columna de genes en ", path, call. = FALSE)
    names(tab)[names(tab) == alt[1]] <- "GENE"
  }
  if (!"FDR" %in% names(tab)) {
    if ("adj.P.Val" %in% names(tab)) {
      tab$FDR <- tab$adj.P.Val
    } else if ("P.Value" %in% names(tab)) {
      tab$FDR <- p.adjust(tab$P.Value, method = "BH")
    } else stop("No hay FDR / adj.P.Val / P.Value en ", path, call. = FALSE)
  }
  if (!"logFC" %in% names(tab)) stop("No hay logFC en ", path, call. = FALSE)
  tab <- dplyr::mutate(tab,
                       GENE = trimws(as.character(GENE)),
                       FDR = as.numeric(FDR),
                       logFC = as.numeric(logFC))
  tab <- dplyr::filter(tab, !is.na(GENE), nzchar(GENE),
                       is.finite(FDR), is.finite(logFC))
  if (anyDuplicated(tab$GENE)) {
    stop("Genes duplicados en DEA: ", path,
         ". Resolver antes de ORA.", call. = FALSE)
  }
  tab
}

# No confiar ciegamente en listas antiguas de figs/. Recalcular siempre
# directamente a partir de los cuatro archivos del paso 06.
get_common_list <- function(bin, direction, tabA, tabB) {
  if (direction == "UP") {
    A <- tabA$GENE[tabA$FDR <= FDR_CUT & tabA$logFC >= LFC_MIN]
    B <- tabB$GENE[tabB$FDR <= FDR_CUT & tabB$logFC >= LFC_MIN]
  } else {
    A <- tabA$GENE[tabA$FDR <= FDR_CUT & tabA$logFC <= -LFC_MIN]
    B <- tabB$GENE[tabB$FDR <= FDR_CUT & tabB$logFC <= -LFC_MIN]
  }
  genes <- sort(intersect(unique(A), unique(B)), method = "radix")
  listfile <- file.path(FIG_DIR, sprintf("common_%s_%s.txt", bin, direction))
  if (file.exists(listfile)) {
    old <- unique(readLines(listfile, warn = FALSE))
    old <- old[nzchar(old)]
    if (!setequal(genes, old)) {
      stop("La lista en figs/ NO coincide con DEA para ", bin, "/", direction,
           ". Revisa el paso 08 y evita usar archivos mezclados: ", listfile,
           call. = FALSE)
    }
  }
  writeLines(genes, listfile)
  writeLines(genes, file.path(OUT_DIR,
                              sprintf("common_%s_%s.txt", bin, direction)))
  genes
}

# Colecciones congeladas en TSV. La primera ejecucion consulta msigdbr;
# ejecuciones posteriores usan la misma instantanea, incluso si el paquete
# es actualizado. Conservar reference_sets/ al compartir el repositorio.
get_msig_t2g <- function(subcollection) {
  snapshot <- file.path(REF_DIR, paste0("MSigDB_", EXPECTED_MSIGDB, "_",
                                        gsub(":", "_", subcollection), ".tsv"))
  if (file.exists(snapshot)) {
    tab <- readr::read_tsv(snapshot, show_col_types = FALSE,
                           col_types = readr::cols(.default = "c"))
    message("[MSigDB] Usando instantanea: ", basename(snapshot))
  } else {
    arg_names <- names(formals(msigdbr::msigdbr))
    args <- if ("collection" %in% arg_names) {
      list(collection = "C2", subcollection = subcollection)
    } else {
      list(category = "C2", subcategory = subcollection)
    }
    db <- do.call(msigdbr::msigdbr,
                  c(list(species = "Homo sapiens"), args))
    if (!nrow(db)) stop("MSigDB devolvio cero filas: ", subcollection)
    if (!"db_version" %in% names(db) ||
        !identical(unique(as.character(db$db_version)), EXPECTED_MSIGDB)) {
      stop("La version de MSigDB no coincide para ", subcollection)
    }
    tab <- unique(data.frame(gs_name = as.character(db$gs_name),
                             gene_symbol = as.character(db$gene_symbol),
                             stringsAsFactors = FALSE))
    tab <- tab[order(tab$gs_name, tab$gene_symbol), , drop = FALSE]
    
  }
  if (!all(c("gs_name", "gene_symbol") %in% names(tab)) || !nrow(tab)) {
    stop("Instantanea invalida: ", snapshot, call. = FALSE)
  }
  prefix <- if (subcollection == "CP:REACTOME") "REACTOME_" else "KEGG_"
  if (!all(startsWith(tab$gs_name, prefix))) {
    stop("Nombres de rutas inesperados en: ", snapshot, call. = FALSE)
  }
  stopifnot(!anyNA(tab$gs_name), !anyNA(tab$gene_symbol))
  tab <- dplyr::distinct(tab, gs_name, gene_symbol)
  if (!file.exists(snapshot)) readr::write_tsv(tab,snapshot)
  message("[MSigDB] ", subcollection, ": ", nrow(tab),
          " asociaciones / ", length(unique(tab$gs_name)), " rutas")
  tab
}

# Cargar ambas colecciones y detenerse si alguna falla.
reactome <- get_msig_t2g(COLLECTIONS[["Reactome"]])
kegg <- get_msig_t2g(COLLECTIONS[["KEGG"]])

# Misma logica original de GO: SYMBOL -> ENTREZ, IDs unicos.
to_entrez <- function(symbols) {
  symbols <- unique(toupper(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (!length(symbols)) return(character())
  mapped <- suppressMessages(clusterProfiler::bitr(
    symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
  mapped <- dplyr::distinct(mapped, SYMBOL, .keep_all = TRUE)
  unique(as.character(mapped$ENTREZID[!is.na(mapped$ENTREZID)]))
}

run_go <- function(genes, universe, ontology) {
  if (length(genes) < 10L || length(universe) < 100L) return(NULL)
  tryCatch(clusterProfiler::enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID", ont = ontology,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    readable = TRUE
  ), error = function(e) {
    stop("Error en GO-", ontology, ": ", conditionMessage(e), call. = FALSE)
  })
}

run_pathway <- function(genes, universe, term2gene, label) {
  # Importante: identificadores SYMBOL en gene, universe y TERM2GENE.
  genes <- intersect(unique(genes), unique(universe))
  mapped <- intersect(genes, unique(term2gene$gene_symbol))
  message("  ", label, ": ", length(mapped), "/", length(genes),
          " genes asociados a alguna ruta")
  if (length(mapped) < 10L) {
    warning(label, ": menos de 10 genes asociados; tabla vacia.")
    return(NULL)
  }
  tryCatch(clusterProfiler::enricher(
    gene = genes,
    universe = universe,
    TERM2GENE = term2gene[, c("gs_name", "gene_symbol")],
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500
  ), error = function(e) {
    stop("Error en ", label, ": ", conditionMessage(e), call. = FALSE)
  })
}

empty_result <- function() {
  data.frame(ID = character(), Description = character(),
             GeneRatio = character(), BgRatio = character(),
             pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
             geneID = character(), Count = integer())
}

save_table <- function(result, filepath) {
  dat <- if (is.null(result)) empty_result() else as.data.frame(result)
  if (!nrow(dat)) dat <- empty_result()
  if (nrow(dat) && "p.adjust" %in% names(dat)) {
    dat <- dat[order(dat$p.adjust, dat$ID), , drop = FALSE]
  }
  readr::write_tsv(dat, filepath)
  nrow(dat)
}

quick_dot <- function(result, out_png, title, topn = TOP_N) {
  dat <- if (is.null(result)) empty_result() else as.data.frame(result)
  dat <- dat[is.finite(dat$p.adjust) & dat$p.adjust <= FDR_CUT,,drop=FALSE]
  if (!nrow(dat)) {
    if (file.exists(out_png)) unlink(out_png)
    return(invisible(NULL))
  }
  dat <- head(dat[order(dat$p.adjust), , drop = FALSE], topn)
  p <- ggplot2::ggplot(dat, ggplot2::aes(
    x = reorder(Description, -p.adjust),
    y = -log10(pmax(p.adjust, 1e-300)), size = Count,
    color = -log10(pmax(p.adjust, 1e-300))
  )) + ggplot2::geom_point() + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "-log10(FDR)", title = title,
                  color = "-log10(FDR)") +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(out_png, p, width = 6, height = 4, dpi = 150)
}

# Un resumen siempre presente, incluso cuando no hay rutas devueltas.
run_log <- list()
for (bin in BINS) {
  t115 <- load_dea("GSE115002", bin)
  t750 <- load_dea("GSE75037", bin)
  universe <- sort(intersect(t115$GENE, t750$GENE), method = "radix")
  universe_entrez <- to_entrez(universe)
  
  for (direction in DIRECTIONS) {
    message("\n== ", bin, " / ", direction, " ==")
    genes <- get_common_list(bin, direction, t115, t750)
    message("  Genes concordantes: ", length(genes),
            "; universo: ", length(universe))
    genes_entrez <- to_entrez(genes)
    
    go_bp <- run_go(genes_entrez, universe_entrez, "BP")
    go_cc <- run_go(genes_entrez, universe_entrez, "CC")
    go_mf <- run_go(genes_entrez, universe_entrez, "MF")
    pw_reactome <- run_pathway(genes, universe, reactome, "Reactome")
    pw_kegg <- run_pathway(genes, universe, kegg, "KEGG_LEGACY")
    
    pref <- file.path(OUT_DIR, sprintf("%s_%s_ORA_common", bin, direction))
    outputs <- list(GO_BP = go_bp, GO_CC = go_cc, GO_MF = go_mf,
                    PATHWAY = pw_reactome, KEGG = pw_kegg)
    
    for (category in names(outputs)) {
      res <- outputs[[category]]
      out_path <- paste0(pref, "_", category, ".tsv")
      n_total <- save_table(res, out_path)
      if (is.null(res)) {
        n_sig <- 0L
      } else {
        dd <- as.data.frame(res)
        n_sig <- if (!nrow(dd)) 0L else sum(dd$p.adjust <= 0.05, na.rm = TRUE)
      }
      message("  ", category, ": ", n_total,
              " rutas en tabla; ", n_sig, " con FDR <= 0.05 -> ",
              basename(out_path))
      run_log[[length(run_log) + 1L]] <- data.frame(
        age_bin = bin, direction = direction, collection = category,
        n_input_genes = length(genes), n_universe = length(universe),
        n_terms_reported = n_total, n_fdr_005 = n_sig,
        output_file = basename(out_path), stringsAsFactors = FALSE
      )
    }
    
    quick_dot(go_bp, paste0(pref, "_GO_BP.png"), paste(bin, direction, "GO:BP"))
    quick_dot(go_cc, paste0(pref, "_GO_CC.png"), paste(bin, direction, "GO:CC"))
    quick_dot(go_mf, paste0(pref, "_GO_MF.png"), paste(bin, direction, "GO:MF"))
    quick_dot(pw_reactome, paste0(pref, "_PATHWAY.png"),
              paste(bin, direction, "Reactome"))
    quick_dot(pw_kegg, paste0(pref, "_KEGG.png"),
              paste(bin, direction, "KEGG_LEGACY"))
  }
}

run_log_df <- dplyr::bind_rows(run_log)
readr::write_tsv(run_log_df, file.path(OUT_DIR, "09_ORA_run_summary.tsv"))

input_paths <- unlist(lapply(BINS, function(bin) vapply(GSES, dea_path, character(1), bin=bin)))
snapshot_names <- paste0("MSigDB_",EXPECTED_MSIGDB,"_",gsub(":", "_",COLLECTIONS),".tsv")
input_paths <- c(input_paths,file.path(REF_DIR,snapshot_names))
stopifnot(all(file.exists(input_paths)))
md5_table <- data.frame(file = basename(input_paths),
                        md5 = unname(tools::md5sum(input_paths)),
                        stringsAsFactors = FALSE)
readr::write_tsv(md5_table, file.path(OUT_DIR, "09_input_checksums.tsv"))
message("\n[OK] Paso 09 completado. Resultados en: ", OUT_DIR)
message("Resumen: enrich_common/09_ORA_run_summary.tsv")
