if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
## 04_map_samples_and_apply_survival_filters.R
## Mapea barcodes -> patient_id, filtra Primary Tumor, aplica filtros OS, deduplica pacientes.
## Añadidos:
## - Arrastra smoking_status y crea smoking_bin (Ever/Never, Unknown->NA)
## - Flow de Ns por paso (tabla + figura) + Badge
## - Distribuciones post-filtros (género, estadio, tabaquismo)
## - Dos cohorts complete-case: sin y con tabaquismo (excluye Unknown)

# ----- Parámetros -----
base_dir <- getwd()
min_os_days  <- 15
max_os_days  <- 3650   # ~10 años

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
  library(ggplot2)
})

out_dir <- file.path(base_dir, "data");    dir.create(out_dir, TRUE, TRUE)
tab_dir <- file.path(base_dir, "tables");  dir.create(tab_dir, TRUE, TRUE)
fig_dir <- file.path(base_dir, "figures"); dir.create(fig_dir, TRUE, TRUE)
log_dir <- file.path(base_dir, "logs");    dir.create(log_dir, TRUE, TRUE)
log_file <- file.path(log_dir, "04_filters.txt")
cat("# 04_map_samples_and_apply_survival_filters\n", file = log_file)

# ---------- Utils ----------
mk_dist <- function(x, nm) {
  tb <- as.data.frame(table(x, useNA = "ifany"))
  names(tb) <- c(nm, "n")
  tb[[nm]] <- as.character(tb[[nm]])
  tb$perc <- round(100 * tb$n / sum(tb$n), 1)
  tb
}

plot_counts <- function(df, var, title, fname, order_levels = NULL, width=7, height=5) {
  df$.x <- df[[var]]
  if (!is.null(order_levels)) df$.x <- factor(df$.x, levels = order_levels)
  p <- ggplot(df, aes(x = .x)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.25, size = 3) +
    labs(x = var, y = "Conteo", title = title) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))
  ggsave(file.path(fig_dir, fname), p, width=width, height=height, dpi=300)
}

# ---------- Cargar clínica (de 03) ----------
clin <- readRDS(file.path(out_dir, "clinical_OS.rds"))

# El paso 02 conservó los barcodes del SE como nombres de columnas.
expr_df <- readRDS(file.path(out_dir, "expr_SYMBOL_matrix.rds"))
stopifnot(names(expr_df)[1] == "gene", !anyDuplicated(expr_df$gene),
          !anyNA(clin$patient_id), !anyDuplicated(clin$patient_id))
barcodes <- names(expr_df)[-1]
stopifnot(!anyDuplicated(barcodes), all(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}",barcodes)))
sample_code <- substr(barcodes,14,15)
keep_barcodes <- barcodes[sample_code == "01"]
stopifnot(length(keep_barcodes) > 0)
cat(sprintf("Expresión: %d | Tumor primario: %d | Tejido normal sólido: %d\n",
            length(barcodes),length(keep_barcodes),sum(sample_code=="11")),file=log_file,append=TRUE)

# ---------- Mapear a clínica (por patient_id) ----------
map_df  <- tibble(barcode = keep_barcodes, patient_id = substr(keep_barcodes, 1, 12))
clin_map <- map_df %>% left_join(clin, by = "patient_id")

# ---------- Filtros de supervivencia (flow) ----------
flow <- tibble(step = character(), n = integer())
flow <- add_row(flow, step = "S0_Expr_PrimaryTumor_barcodes", n = nrow(clin_map))

clin1 <- clin_map %>% filter(!is.na(OS.time), !is.na(OS.event))
flow  <- add_row(flow, step = "S1_OS_defined", n = nrow(clin1))

clin2 <- clin1 %>% filter(OS.time >= min_os_days)
flow  <- add_row(flow, step = sprintf("S2_OS_ge_%dd", min_os_days), n = nrow(clin2))

clin3 <- clin2 %>% filter(OS.time <= max_os_days)
flow  <- add_row(flow, step = sprintf("S3_OS_le_%dd", max_os_days), n = nrow(clin3))

# Deduplicar por paciente (elige barcode alfabéticamente primero)
clin3 <- clin3 %>%
  arrange(patient_id, barcode) %>%
  group_by(patient_id) %>% slice(1) %>% ungroup()
flow <- add_row(flow, step = "S4_unique_patient", n = nrow(clin3))

# ---------- Badge de N clave ----------
n_badge <- tibble(
  item = c("Primary (barcodes)", "OS definido",
           sprintf("OS ≥ %dd", min_os_days),
           sprintf("OS ≤ %dd", max_os_days),
           "Pacientes únicos"),
  N = flow$n
)
write_tsv(flow,   file.path(tab_dir, "Flow_04_counts.tsv"))
write_tsv(n_badge,file.path(tab_dir, "Flow04_badge.tsv"))

# Figuras flow + badge
p_flow <- ggplot(flow, aes(x = factor(step, levels = step), y = n)) +
  geom_col() + geom_text(aes(label = n), vjust = -0.3, size = 3) +
  labs(title = "Flujo de filtros (Paso 04)", x = "Paso", y = "N muestras/pacientes") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 28, hjust = 1))
ggsave(file.path(fig_dir, "Flow04_bar.png"), p_flow, width = 10, height = 4.8, dpi = 300)
ggsave(file.path(fig_dir, "Flow04_bar.pdf"), p_flow, width = 10, height = 4.8)

p_badge <- ggplot(n_badge, aes(x = reorder(item, N), y = N, label = N)) +
  geom_col(width = 0.6) +
  geom_text(vjust = -0.25, size = 4) +
  coord_flip() +
  labs(title = "Badge de N por paso (04)", x = NULL, y = "N") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(fig_dir, "Flow04_badge.png"), p_badge, width = 6, height = 4, dpi = 300)

cat("\nResumen Flow 04:\n", file = log_file, append = TRUE)
capture.output(flow, file = log_file, append = TRUE)

# ---------- Cohorte final (agrega tabaquismo antes de exportar) ----------
final_barcodes   <- clin3$barcode
expr_tumor_final <- expr_df[, c("gene", final_barcodes), drop = FALSE]
stopifnot(identical(names(expr_tumor_final)[-1], clin3$barcode), !anyDuplicated(clin3$patient_id), nrow(clin3)>0)

clin3 <- clin3 %>%
  mutate(
    smoking_status = factor(smoking_status, levels = c("Never","Former","Current","Unknown")),
    smoking_bin = case_when(
      smoking_status == "Never" ~ "Never",
      smoking_status %in% c("Former","Current") ~ "Ever",
      TRUE ~ NA_character_
    ),
    smoking_bin = factor(smoking_bin, levels = c("Never","Ever"))
  )

# ---------- Exportes principales ----------
saveRDS(expr_tumor_final, file.path(out_dir, "expr_tumor_SYMBOL_filtered.rds"))
write_tsv(expr_tumor_final, file.path(out_dir, "expr_tumor_SYMBOL_filtered.tsv"))
saveRDS(clin3, file.path(out_dir, "clin_survival_cohort.rds"))
write_tsv(clin3, file.path(out_dir, "clin_survival_cohort.tsv"))

cat("\n[OK] Cohorte de supervivencia exportada (incluye smoking_status y smoking_bin).\n",
    file = log_file, append = TRUE)

# ---------- Distribuciones post-filtro ----------
write_tsv(mk_dist(clin3$gender, "gender"),                 file.path(tab_dir, "AfterFilters_Dist_gender.tsv"))
write_tsv(mk_dist(clin3$stage_use, "stage"),               file.path(tab_dir, "AfterFilters_Dist_stage.tsv"))
write_tsv(mk_dist(clin3$smoking_status, "smoking_status"), file.path(tab_dir, "AfterFilters_Dist_smoking.tsv"))
write_tsv(mk_dist(clin3$smoking_bin, "smoking_bin"),       file.path(tab_dir, "AfterFilters_Dist_smokingBIN.tsv"))

plot_counts(clin3, "gender",         "Después de filtros: Género",           "AfterFilters_counts_gender.png",
            c("Female","Male"))
plot_counts(clin3, "stage_use",      "Después de filtros: Estadio (I–IV)",   "AfterFilters_counts_stage.png",
            c("I","II","III","IV"))
plot_counts(clin3, "smoking_status", "Después de filtros: Tabaquismo (4 niv.)","AfterFilters_counts_smoking.png",
            c("Never","Former","Current","Unknown"), width=8)
plot_counts(clin3, "smoking_bin",    "Después de filtros: Tabaquismo (bin)", "AfterFilters_counts_smokingBIN.png",
            c("Never","Ever"), width=6)

# ---------- OS (años) ----------
p_os_hist <- ggplot(clin3, aes(x = OS.time/365.25)) +
  geom_histogram(bins = 30) +
  labs(title = "Distribución del tiempo de seguimiento (años)", x = "Años", y = "Conteo") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "AfterFilters_OS_hist_years.png"), p_os_hist, width = 7, height = 5, dpi = 300)

# ---------- Complete-case sets ----------
# (A) edad + género + estadio
clin_cc <- clin3 %>%
  mutate(age_years = as.numeric(age_years),
         gender    = as.factor(gender),
         stage_use = as.factor(stage_use)) %>%
  filter(!is.na(age_years), !is.na(gender), !is.na(stage_use))
expr_cc <- expr_tumor_final[, c("gene", clin_cc$barcode), drop = FALSE]
saveRDS(clin_cc, file.path(out_dir, "clin_complete_case.rds"))
write_tsv(clin_cc, file.path(out_dir, "clin_complete_case.tsv"))
saveRDS(expr_cc, file.path(out_dir, "expr_tumor_SYMBOL_complete_case.rds"))
write_tsv(expr_cc, file.path(out_dir, "expr_tumor_SYMBOL_complete_case.tsv"))
cat(sprintf("\nComplete-case (edad+género+estadio): n=%d\n", nrow(clin_cc)), file = log_file, append = TRUE)

# (B) + tabaquismo (excluye Unknown)
clin_cc_sm <- clin3 %>%
  mutate(age_years = as.numeric(age_years),
         gender    = as.factor(gender),
         stage_use = as.factor(stage_use),
         smoking_status = as.factor(smoking_status),
         smoking_bin    = as.factor(smoking_bin)) %>%
  filter(!is.na(age_years), !is.na(gender), !is.na(stage_use),
         !is.na(smoking_status), smoking_status != "Unknown")
expr_cc_sm <- expr_tumor_final[, c("gene", clin_cc_sm$barcode), drop = FALSE]
saveRDS(clin_cc_sm, file.path(out_dir, "clin_complete_case_w_smoking.rds"))
write_tsv(clin_cc_sm, file.path(out_dir, "clin_complete_case_w_smoking.tsv"))
saveRDS(expr_cc_sm, file.path(out_dir, "expr_tumor_SYMBOL_complete_case_w_smoking.rds"))
write_tsv(expr_cc_sm, file.path(out_dir, "expr_tumor_SYMBOL_complete_case_w_smoking.tsv"))
cat(sprintf("Complete-case + smoking: n=%d (excluye 'Unknown')\n", nrow(clin_cc_sm)), file = log_file, append = TRUE)

cat("\n[OK] 04_map_samples_and_apply_survival_filters COMPLETADO.\n", file = log_file, append = TRUE)

# ===== Donuts (composición final tras filtros) =====
plot_donut <- function(df, var, title, file_stub, order_levels = NULL, w=5, h=5, dpi=300) {
  v <- df[[var]]
  if (!is.null(order_levels)) v <- factor(v, levels = order_levels)
  tb <- as.data.frame(table(v, useNA = "ifany"), stringsAsFactors = FALSE)
  names(tb) <- c("cat","n")
  tb$cat <- as.character(tb$cat)
  # Etiqueta NA visible si existe
  tb$cat[is.na(tb$cat)] <- "NA"
  # Filtra categorías con n=0 (evita sectores vacíos)
  tb <- tb[tb$n > 0, , drop = FALSE]
  tb$perc <- tb$n / sum(tb$n)
  tb$lbl  <- sprintf("%s\n%d (%.1f%%)", tb$cat, tb$n, 100*tb$perc)
  
  p <- ggplot(tb, aes(x = 2, y = n, fill = cat)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    theme_void(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(title = title) +
    geom_text(aes(label = lbl), position = position_stack(vjust = 0.5), size = 3)
  
  ggsave(file.path(fig_dir, paste0(file_stub, ".png")), p, width = w, height = h, dpi = dpi)
  ggsave(file.path(fig_dir, paste0(file_stub, ".pdf")), p, width = w, height = h)
}

# Llamadas (ajusto el orden para que se vean bonitos)
plot_donut(clin3, "gender",
           "Composición final: Género",
           "AfterFilters_donut_gender",
           order_levels = c("Female","Male", NA))

plot_donut(clin3, "stage_use",
           "Composición final: Estadio (I–IV)",
           "AfterFilters_donut_stage",
           order_levels = c("I","II","III","IV", NA))

plot_donut(clin3, "smoking_status",
           "Composición final: Tabaquismo (Never/Former/Current/Unknown)",
           "AfterFilters_donut_smoking",
           order_levels = c("Never","Former","Current","Unknown", NA), w = 5.5, h = 5.5)

plot_donut(clin3, "smoking_bin",
           "Composición final: Tabaquismo (Never/Ever)",
           "AfterFilters_donut_smokingBIN",
           order_levels = c("Never","Ever", NA), w = 5, h = 5)
print(flow)
cat(sprintf("[OK] Paso 04 | Tumores primarios: %d | Pacientes únicos: %d | Muertes: %d | Censuras: %d\n",length(keep_barcodes),nrow(clin3),sum(clin3$OS.event == 1),sum(clin3$OS.event == 0)))
