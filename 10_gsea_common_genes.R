############################################################
# GSEA of common genes across LUAD age bins
#
# This script combines gene-level paired Tumor-vs-Normal
# statistics from GSE115002 and GSE75037 using a Stouffer
# z-score ranking, then runs fgsea on MSigDB Hallmark pathways.
#
# Inputs:
#   dea/<GSE>__TvN_<bin>.tsv
#
# Outputs:
#   gsea_common_genes/fgsea_H_common_<bin>.tsv
#   gsea_common_genes/*ES*.png/pdf
#   gsea_common_genes/panel*.png/pdf
#   gsea_common_genes/leading_edge_genes_H_panel2x2.tsv
############################################################

if (!dir.exists("dea")) stop("Abre GEO.Rproj antes de ejecutar este script.")
set.seed(123)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(fgsea); library(ggplot2);   library(patchwork)
  library(glue);  library(grid);   library(BiocParallel)
})

## --- ejecución en serie (determinista) ---
BiocParallel::register(BiocParallel::SerialParam())
set.seed(123)
options(bitmapType = "cairo")

## ------------------ Configuración básica ------------------

ROOT    <- getwd()
DEA_DIR <- file.path(ROOT, "dea")
OUTDIR  <- file.path(ROOT, "gsea_common_genes")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## Datasets
GSE_A <- "GSE115002"
GSE_B <- "GSE75037"

## Pesos relativos (por nº de pares T/N si quieres; ahora 1:1)
weights <- c(A = 1, B = 1)

## Cómo tratar genes con signo opuesto entre GSEs:
## "keep" = combina tal cual, "zero" = penaliza poniéndolos a 0, "drop" = excluye
CONCORDANCE <- "keep"

## Colección MSigDB para GSEA:
CATEGORY <- "H"   # Hallmark

## Parámetros fgsea / reporting
MIN_SZ  <- 15
MAX_SZ  <- 5000
TOPK    <- 10          # nº de vías con curvas ES por bin
FDR_CUT <- 0.05        # corte para paneles
NES_MIN <- 0           # |NES| mínimo para panel (0 = sin filtro)

## Etiquetas bonitas de cohorte para las figuras
BIN_LABELS <- c(
  Young = "<60 years",
  Old   = "≥60 years"
)

## ------------------ Utilidades ------------------

## Lee TvN por bin; tolera __TvN_ y ___TvN_ y columnas varias de gen/P
read_tt <- function(gse, bin){
  cand <- c(
    file.path(DEA_DIR, sprintf("%s___TvN_%s.tsv", gse, bin)),
    file.path(DEA_DIR, sprintf("%s__TvN_%s.tsv",  gse, bin))
  )
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) stop("No encuentro TvN para ", gse, " / ", bin)
  
  tt <- readr::read_tsv(f, show_col_types = FALSE)
  
  ## columna de gen
  gcol <- c("GENE","Gene","SYMBOL","symbol","gene","ID")
  gcol <- gcol[gcol %in% names(tt)][1]
  if (is.na(gcol)) stop("No encuentro columna de genes en: ", f)
  names(tt)[names(tt) == gcol] <- "GENE"
  
  ## columna de P.Value
  if (!"P.Value" %in% names(tt)) {
    pcol <- c("P.Value","PValue","pvalue","PVAL","P.value")
    pcol <- pcol[pcol %in% names(tt)][1]
    if (is.na(pcol)) stop("No encuentro P.Value en: ", f)
    names(tt)[names(tt) == pcol] <- "P.Value"
  }
  
  ## columna t (limma típico); si no hay, construimos un estadístico sustituto
  if (!"t" %in% names(tt)) {
    if (all(c("logFC", "P.Value") %in% names(tt))) {
      eps <- 1e-300
      tt <- tt |>
        mutate(t = sign(logFC) * (-log10(pmax(P.Value, eps))))
    } else {
      stop("No encuentro 't' ni (logFC + P.Value) en: ", f)
    }
  }
  
  tt |>
    mutate(GENE = as.character(GENE)) |>
    group_by(GENE) |>
    summarise(
      t       = mean(t, na.rm = TRUE),
      P.Value = pmax(1e-300, min(P.Value, na.rm = TRUE)),
      .groups = "drop"
    )
}

## Ranking combinado por Stouffer (z-score) limitado a genes comunes
stouffer_rank <- function(ttA, ttB, weights = c(A=1,B=1),
                          concordance = c("keep","zero","drop")){
  concordance <- match.arg(concordance)
  need <- c("GENE","t","P.Value")
  stopifnot(all(need %in% names(ttA)), all(need %in% names(ttB)))
  
  G <- intersect(ttA$GENE, ttB$GENE)
  if (!length(G)) stop("No hay genes comunes entre los dos GSE en este bin.")
  
  ttA <- ttA[match(G, ttA$GENE), ]
  ttB <- ttB[match(G, ttB$GENE), ]
  
  ## Z de dos colas con signo del t
  zA <- sign(ttA$t) * qnorm(ttA$P.Value/2, lower.tail = FALSE)
  zB <- sign(ttB$t) * qnorm(ttB$P.Value/2, lower.tail = FALSE)
  
  samesign <- sign(zA) == sign(zB)
  if (concordance == "zero") {
    zA[!samesign] <- 0; zB[!samesign] <- 0
  }
  if (concordance == "drop") {
    keep <- which(samesign)
    zA <- zA[keep]; zB <- zB[keep]; G <- G[keep]
  }
  
  wA <- weights["A"]; wB <- weights["B"]
  zc <- (wA*zA + wB*zB) / sqrt(wA^2 + wB^2)
  names(zc) <- G
  
  zc <- zc[is.finite(zc)]
  sort(zc, decreasing = TRUE)
}

## Obtiene pathways MSigDB restringidos al universo del ranking
gmt <- "h.all.v2025.1.Hs.symbols.gmt"
if (!file.exists(gmt)) {
  url <- paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/",gmt)
  pending <- tempfile(tmpdir=".")
  download.file(url,pending,mode="wb")
  stopifnot(length(fgsea::gmtPathways(pending))==50L,file.rename(pending,gmt))
}
get_pathways <- function(category, genes_universe, min_sz=MIN_SZ, max_sz=MAX_SZ){
  pw <- fgsea::gmtPathways(gmt)
  stopifnot(length(pw)==50L,all(startsWith(names(pw),"HALLMARK_")))
  pwU  <- lapply(pw, function(g) intersect(g, genes_universe))
  pwU  <- pwU[lengths(pwU) >= min_sz & lengths(pwU) <= max_sz]
  pwU
}

## fgsea ordenado por FDR y |NES|
run_fgsea <- function(pathways, stats){
  stats <- stats[is.finite(stats)]
  fg <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats    = stats,
    BPPARAM  = BiocParallel::SerialParam()
  )
  fg[order(fg$padj, -abs(fg$NES)), ]
}

## Etiqueta bonita a partir de ID Hallmark/Reactome/etc.
pretty_label <- function(id){
  x <- gsub("^HALLMARK_", "", id)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

nice_title <- function(fg_tbl, pathway_id){
  row <- fg_tbl[fg_tbl$pathway == pathway_id, ][1, ]
  sprintf("%s\nNES=%.2f   FDR=%s",
          pretty_label(pathway_id),
          row$NES,
          format(row$padj, digits = 2, scientific = TRUE))
}

## Curvas ES para TOPK vías
save_es_curves <- function(fg_tbl, stats, pathways, bin_label = "young"){
  if (!nrow(fg_tbl)) return(invisible(NULL))
  top_paths <- unique(c(head(fg_tbl$pathway, TOPK), intersect(c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1"), fg_tbl$pathway)))
  
  for (nm in top_paths){
    row <- fg_tbl[fg_tbl$pathway == nm, ]
    p_pretty <- fgsea::plotEnrichment(pathways[[nm]], stats) +
      ggplot2::labs(title = nice_title(fg_tbl, nm),
                    subtitle = sprintf("p = %.3g", row$pval))
    safe <- gsub("[^A-Za-z0-9._-]+", "_", nm)
    outfile_png <- file.path(
      OUTDIR,
      sprintf("%s_ES_%s_common_%s.png", bin_label, safe, CATEGORY)
    )
    ggsave(outfile_png, p_pretty, width = 9, height = 7, dpi = 300)
    ggsave(sub("\\.png$", ".pdf", outfile_png), p_pretty,
           width = 9, height = 7, device = cairo_pdf)
  }
}

## ------------------ Panel 2x2 |NES| (ambos bins juntos) ------------------

barpanel_2x2 <- function(tbl_y, tbl_o, out_png, out_pdf,
                         topk=TOPK, fdr=FDR_CUT, nes_min=NES_MIN){
  make_df <- function(tbl, bin_code){
    bin_label <- BIN_LABELS[[bin_code]]
    tbl |>
      filter(padj <= fdr, abs(NES) >= nes_min) |>
      mutate(
        Direction = ifelse(NES >= 0, "UP", "DOWN"),
        score     = abs(NES),
        Label     = pretty_label(pathway),
        Bin       = bin_label
      ) |>
      arrange(desc(score)) |>
      group_by(Direction) |>
      slice_head(n = topk) |>
      ungroup()
  }
  
  dy <- make_df(tbl_y, "Young")
  do <- make_df(tbl_o, "Old")
  if (!nrow(dy) && !nrow(do)) {
    warning("Nada que graficar en el panel 2x2 (revisa FDR_CUT / NES_MIN).")
    return(invisible(NULL))
  }
  
  pal_dir <- c("UP"="#E41A1C","DOWN"="#08519c")
  
  make_p <- function(subdf, title){
    if (!nrow(subdf)) return(ggplot() + theme_void() + ggtitle(title))
    subdf <- subdf |>
      mutate(Label = factor(Label, levels = rev(unique(Label))))
    ggplot(subdf, aes(x = Label, y = score, fill = Direction)) +
      geom_col(width = 0.78, alpha = 0.95) +
      geom_text(aes(label = sprintf("NES=%.2f", ifelse(Direction=="UP", score, -score))),
                hjust = -0.05, size = 3.2, colour = "grey20") +
      coord_flip(clip = "off") +
      scale_fill_manual(values = pal_dir, name = "Direction") +
      scale_y_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0, .08))) +
      labs(x = NULL, y = "|NES|", title = title) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position    = "right",
        plot.title.position = "plot",
        plot.margin        = margin(8, 36, 8, 14)
      )
  }
  
  p_y_up   <- make_p(dy |> filter(Direction=="UP"),
                     sprintf("%s — Top UP", BIN_LABELS[["Young"]]))
  p_y_down <- make_p(dy |> filter(Direction=="DOWN"),
                     sprintf("%s — Top DOWN", BIN_LABELS[["Young"]]))
  p_o_up   <- make_p(do |> filter(Direction=="UP"),
                     sprintf("%s — Top UP", BIN_LABELS[["Old"]]))
  p_o_down <- make_p(do |> filter(Direction=="DOWN"),
                     sprintf("%s — Top DOWN", BIN_LABELS[["Old"]]))
  
  panel <- (p_y_up | p_y_down) / (p_o_up | p_o_down) +
    plot_layout(guides = "collect") &
    theme(
      legend.position    = "right",
      legend.box         = "vertical",
      legend.box.spacing = unit(18, "pt"),
      legend.margin      = margin(10,10,10,10),
      legend.box.margin  = margin(0,24,0,8),
      legend.key.height  = unit(14, "pt"),
      legend.key.width   = unit(18, "pt")
    )
  
  ggsave(out_png, panel, width = 15, height = 10, dpi = 300)
  ggsave(out_pdf, panel, width = 15, height = 10, device = cairo_pdf)
  message("Panel 2x2 → ", out_png, " ; ", out_pdf)
  invisible(panel)
}

## ------------------ Paneles separados por cohorte ------------------

barpanel_split_bins <- function(tbl_y, tbl_o, out_prefix,
                                topk=TOPK, fdr=FDR_CUT, nes_min=NES_MIN) {
  pal_dir <- c("UP"="#E41A1C","DOWN"="#08519c")
  
  make_df <- function(tbl, bin_code){
    bin_label <- BIN_LABELS[[bin_code]]
    tbl |>
      filter(padj <= fdr, abs(NES) >= nes_min) |>
      mutate(
        Direction = ifelse(NES >= 0, "UP", "DOWN"),
        score     = abs(NES),
        Label     = pretty_label(pathway),
        Bin       = bin_label
      ) |>
      arrange(desc(score)) |>
      group_by(Direction) |>
      slice_head(n = topk) |>
      ungroup()
  }
  
  dy <- make_df(tbl_y, "Young")
  do <- make_df(tbl_o, "Old")
  
  make_p <- function(subdf, bin_label, dir_label){
    title_txt <- sprintf("%s — Top %s", bin_label, dir_label)
    if (!nrow(subdf)) return(ggplot() + theme_void() + ggtitle(title_txt))
    subdf <- subdf |>
      mutate(Label = factor(Label, levels = rev(unique(Label))))
    
    ggplot(subdf, aes(x = Label, y = score, fill = Direction)) +
      geom_col(width = 0.78, alpha = 0.95) +
      geom_text(aes(label = sprintf("NES=%.2f", ifelse(Direction=="UP", score, -score))),
                hjust = -0.05, size = 3.0, colour = "grey20") +
      coord_flip(clip = "off") +
      scale_fill_manual(values = pal_dir, name = "Direction") +
      scale_y_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0, .08))) +
      labs(x = NULL, y = "|NES|", title = title_txt) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position    = "right",
        plot.title.position = "plot",
        plot.margin        = margin(8, 36, 8, 14)
      )
  }
  
  ## ---- Panel <60 years ----
  bin_y_lab <- BIN_LABELS[["Young"]]
  p_y_up    <- make_p(dy |> filter(Direction=="UP"),   bin_y_lab, "UP")
  p_y_down  <- make_p(dy |> filter(Direction=="DOWN"), bin_y_lab, "DOWN")
  
  panel_y <- (p_y_up | p_y_down) +
    plot_layout(guides = "collect") &
    theme(
      legend.position    = "right",
      legend.box         = "vertical",
      legend.box.spacing = unit(18, "pt"),
      legend.margin      = margin(10,10,10,10),
      legend.key.height  = unit(14, "pt"),
      legend.key.width   = unit(18, "pt")
    )
  
  out_png_y <- paste0(out_prefix, "_lt60.png")
  out_pdf_y <- paste0(out_prefix, "_lt60.pdf")
  ggsave(out_png_y, panel_y, width = 13, height = 6, dpi = 300)
  ggsave(out_pdf_y, panel_y, width = 13, height = 6, device = cairo_pdf)
  message("Panel <60 years → ", out_png_y)
  
  ## ---- Panel ≥60 years ----
  bin_o_lab <- BIN_LABELS[["Old"]]
  p_o_up    <- make_p(do |> filter(Direction=="UP"),   bin_o_lab, "UP")
  p_o_down  <- make_p(do |> filter(Direction=="DOWN"), bin_o_lab, "DOWN")
  
  panel_o <- (p_o_up | p_o_down) +
    plot_layout(guides = "collect") &
    theme(
      legend.position    = "right",
      legend.box         = "vertical",
      legend.box.spacing = unit(18, "pt"),
      legend.margin      = margin(10,10,10,10),
      legend.key.height  = unit(14, "pt"),
      legend.key.width   = unit(18, "pt")
    )
  
  out_png_o <- paste0(out_prefix, "_ge60.png")
  out_pdf_o <- paste0(out_prefix, "_ge60.pdf")
  ggsave(out_png_o, panel_o, width = 13, height = 6, dpi = 300)
  ggsave(out_pdf_o, panel_o, width = 13, height = 6, device = cairo_pdf)
  message("Panel ≥60 years → ", out_png_o)
  
  invisible(list(lt60 = panel_y, ge60 = panel_o))
}

## ------------------ Pipeline por bin ------------------

run_bin <- function(bin = c("young","old")){
  bin <- match.arg(bin)
  
  message("\n== Bin: ", bin, " ==")
  
  ## 1) Lee TvN de ambos GSE
  ttA <- read_tt(GSE_A, bin)
  ttB <- read_tt(GSE_B, bin)
  
  ## 2) Ranking combinado SOLO con genes comunes
  z_rank <- stouffer_rank(ttA, ttB, weights = weights, concordance = CONCORDANCE)
  if (!length(z_rank)) stop("Ranking vacío para bin=", bin)
  
  ## 3) Pathways (universo restringido al ranking)
  pathways <- get_pathways(
    CATEGORY,
    genes_universe = names(z_rank),
    min_sz = MIN_SZ,
    max_sz = MAX_SZ
  )
  if (!length(pathways)) {
    stop("Sin pathways tras filtrar por tamaño en bin=", bin,
         ". Revisa MIN_SZ / MAX_SZ / CATEGORY.")
  }
  
  ## 4) fgsea
  set.seed(123)
  fg <- run_fgsea(pathways, z_rank)
  
  ## 5) Exporta tabla
  out_tab <- file.path(OUTDIR,
                       sprintf("fgsea_%s_common_%s.tsv", CATEGORY, bin))
  readr::write_tsv(as.data.frame(fg) |> mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))), out_tab)
  message("fgsea ", bin, " → ", out_tab)
  
  ## 6) Curvas ES para TOPK vías
  bin_lab <- ifelse(bin == "young", "Young", "Old")
  save_es_curves(fg, z_rank, pathways, bin_label = bin_lab)
  
  list(fg = fg, pathways = pathways, stats = z_rank)
}

## ------------------ Ejecuta: young y old ------------------

message("\n== Running GSEA (young) ==")
res_y <- run_bin("young")

message("\n== Running GSEA (old) ==")
res_o <- run_bin("old")

## ------------------ Panel 2×2 + paneles separados ------------------

panel_png <- file.path(OUTDIR, sprintf("panel2x2_%s_both_bins.png", CATEGORY))
panel_pdf <- file.path(OUTDIR, sprintf("panel2x2_%s_both_bins.pdf", CATEGORY))

barpanel_2x2(res_y$fg, res_o$fg, panel_png, panel_pdf,
             topk = TOPK, fdr = FDR_CUT, nes_min = NES_MIN)

split_prefix <- file.path(OUTDIR, sprintf("panel_%s_bins", CATEGORY))
barpanel_split_bins(res_y$fg, res_o$fg,
                    out_prefix = split_prefix,
                    topk = TOPK, fdr = FDR_CUT, nes_min = NES_MIN)

message("\nListo. Salidas en: ", OUTDIR,
        "\n- fgsea_", CATEGORY, "_common_{young,old}.tsv",
        "\n- {Young,Old}_ES_<pathway>_common_", CATEGORY, ".(png/pdf)",
        "\n- panel2x2_", CATEGORY, "_both_bins.(png/pdf)",
        "\n- panel_", CATEGORY, "_bins_{lt60,ge60}.(png/pdf)\n")

## ============================================================
## Leading-edge genes SOLO para los pathways del panel 2×2
## ============================================================

## 1) Recrear la tabla de pathways que aparecen en el panel 2×2
panel_make_df <- function(tbl, bin_code){
  bin_label <- BIN_LABELS[[bin_code]]
  tbl %>%
    filter(padj <= FDR_CUT, abs(NES) >= NES_MIN) %>%
    mutate(
      Direction = ifelse(NES >= 0, "UP", "DOWN"),
      score     = abs(NES),
      Bin       = bin_label
    ) %>%
    arrange(desc(score)) %>%
    group_by(Direction) %>%
    slice_head(n = TOPK) %>%
    ungroup() %>%
    select(Bin, Direction, pathway, NES, padj)
}

panel_y <- panel_make_df(res_y$fg, "Young")
panel_o <- panel_make_df(res_o$fg, "Old")
panel_paths <- bind_rows(panel_y, panel_o)

## 2) Función robusta para extraer leading-edge
extract_leading_table <- function(fg_tbl, bin_label){
  df <- as.data.frame(fg_tbl)
  le_col <- grep("leading", names(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(le_col)) {
    stop("No encuentro columna tipo 'leadingEdge' en la tabla fgsea.")
  }
  df <- df %>%
    mutate(Bin = bin_label) %>%
    select(Bin, pathway, NES, padj, all_of(le_col)) %>%
    tidyr::unnest(cols = all_of(le_col))
  names(df)[names(df) == le_col] <- "GENE"
  df[order(df$Bin, df$pathway, df$GENE), ]
}

leading_young_all <- extract_leading_table(res_y$fg, BIN_LABELS[["Young"]])
leading_old_all   <- extract_leading_table(res_o$fg, BIN_LABELS[["Old"]])

leading_all <- bind_rows(leading_young_all, leading_old_all)

## 4) Filtrar SOLO pathways del panel 2×2
panel_leading <- leading_all %>%
  inner_join(panel_paths,
             by = c("Bin", "pathway", "NES", "padj")) %>%
  arrange(Bin, Direction, pathway, GENE)

## 5) Guardar tabla final
out_panel_leading <- file.path(
  OUTDIR,
  sprintf("leading_edge_genes_%s_panel2x2.tsv", CATEGORY)
)

write_tsv(panel_leading, out_panel_leading)
message("Leading-edge (solo pathways del panel 2×2) → ", out_panel_leading)

# Paneles de la figura del artículo con los objetos de ESTA ejecución.
bars <- function(result, label) {
  tab <- dplyr::filter(result$fg,NES>0,padj<=FDR_CUT) |> dplyr::arrange(dplyr::desc(NES)) |> head(10)
  ggplot(tab,aes(x=NES,y=reorder(pathway,NES))) + geom_col(fill="#205776") +
    scale_y_discrete(labels=function(x) gsub("_"," ",sub("^HALLMARK_","",x))) +
    labs(x="Normalized enrichment score", y=NULL,title=label) + theme_bw(base_size=10)
}
curves <- function(result,label) {
  ids <- paste0("HALLMARK_",c("E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1"))
  lapply(ids,function(id) {
    row <- result$fg[result$fg$pathway==id,]
    stopifnot(nrow(row)==1)
    fgsea::plotEnrichment(result$pathways[[id]],result$stats) +
      labs(title=pretty_label(id), subtitle=sprintf("%s   NES %.2f   FDR %.2g",label,row$NES,row$padj)) +
      theme_bw(base_size=10)
  })
}
c_panel <- bars(res_y,"<60 years") + bars(res_o,"≥60 years")
d_panel <- patchwork::wrap_plots(c(curves(res_y,"<60 years"),curves(res_o,"≥60 years")),ncol=3)
combined <- patchwork::wrap_plots(c_panel,d_panel,ncol=1,heights=c(1,1.4))
for (ext in c("pdf","tiff")) {
  ggsave(file.path(OUTDIR,paste0("Fig2cd_GSEA.",ext)), combined,width=13,height=10,dpi=300)
}
writeLines(capture.output(sessionInfo()),file.path(OUTDIR,"session.txt"))
