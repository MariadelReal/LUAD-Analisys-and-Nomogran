if (!file.exists("TCGA.Rproj")) stop("Abre TCGA.Rproj antes de ejecutar el paso.")
## 03_clinical_qc_and_os_construction.R
## Objetivo:
## - Clínica por paciente (sin intersectar con expresión)
## - OS.time / OS.event numéricos
## - Edad continua (años) con rastro: age_source = diagnosis/index/missing
## - Stage I/II/III/IV (ajcc/tumor/pathologic/clinical)
## - Tabaquismo robusto (Never / Former / Current / Unknown) + binario (Never/Ever)
## - Tablas y figuras para inspección visual de N

# ===== Parámetros =====
base_dir <- getwd()
EXCLUDE_OBFUSCATED_AGE <- FALSE   # pon TRUE si quieres eliminar age_is_obfuscated==TRUE

suppressPackageStartupMessages({
  library(dplyr);   library(readr);  library(tibble)
  library(stringr); library(ggplot2)
})

`%||%` <- function(a,b) if(!is.null(a)) a else b

raw_dir <- file.path(base_dir, "data_raw")
out_dir <- file.path(base_dir, "data");    dir.create(out_dir, TRUE, TRUE)
tab_dir <- file.path(base_dir, "tables");  dir.create(tab_dir, TRUE, TRUE)
fig_dir <- file.path(base_dir, "figures"); dir.create(fig_dir, TRUE, TRUE)
log_dir <- file.path(base_dir, "logs");    dir.create(log_dir, TRUE, TRUE)
log_file <- file.path(log_dir, "03_clinical_qc.txt")
cat("# 03_clinical_qc_and_os_construction\n", file = log_file)

clin_path <- file.path(raw_dir, "clinical_LUAD.tsv")
stopifnot(file.exists(clin_path))
clin_raw <- readr::read_tsv(clin_path, show_col_types = FALSE)

# ---------- Utilidades ----------
norm_stage_I_IV <- function(x) {
  s <- tolower(as.character(x))
  s <- gsub("^stage\\s+", "", s)
  s <- gsub("\\s", "", s)
  s <- toupper(s)
  s[grepl("^I([AB])?$",    s)] <- "I"
  s[grepl("^II([AB])?$",   s)] <- "II"
  s[grepl("^III([ABC])?$", s)] <- "III"
  s[grepl("^IV([AB])?$",   s)] <- "IV"
  s[!s %in% c("I","II","III","IV")] <- NA
  factor(s, levels = c("I","II","III","IV"))
}

norm_gender <- function(x) {
  g <- tolower(as.character(x))
  g[g %in% c("female","f")] <- "Female"
  g[g %in% c("male","m")]   <- "Male"
  g[!g %in% c("Female","Male")] <- NA
  factor(g, levels = c("Female","Male"))
}

# Unificar nombres entre versiones de la tabla clínica.
campo <- function(candidatos) {
  out <- rep(NA_character_, nrow(clin_raw))
  for (nm in intersect(candidatos,names(clin_raw))) {
    v <- trimws(as.character(clin_raw[[nm]]))
    ok <- is.na(out) & !is.na(v) & nzchar(v) & !tolower(v) %in% c("unknown","not reported","not available","--")
    out[ok] <- v[ok]
  }
  out
}
pid <- campo(c("submitter_id","bcr_patient_barcode","case_submitter_id"))
stopifnot(!anyNA(pid), all(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$",pid)))
status <- tolower(campo("vital_status"))
evento <- match(status,c("alive","dead")) - 1L
stopifnot(any(!is.na(evento)))
for (nm in c("days_to_death","days_to_last_follow_up")) {
  if (!nm %in% names(clin_raw)) stop("Falta la columna ",nm,". Revisa names(clin_raw).")
}
tiempo <- ifelse(evento == 1L,as.numeric(campo("days_to_death")),as.numeric(campo("days_to_last_follow_up")))
sexo <- campo(c("gender","sex_at_birth"))
tabaco <- tolower(campo("tobacco_smoking_status"))
smk_status_vec <- factor(dplyr::case_when(
  tabaco == "lifelong non-smoker" ~ "Never", tabaco == "current smoker" ~ "Current",
  grepl("^current reformed smoker",tabaco) ~ "Former", TRUE ~ "Unknown"),
  levels=c("Never","Former","Current","Unknown"))
# ---------- EDAD: construir vectores primero (diagnosis -> index) ----------
age_diag_num  <- suppressWarnings(as.numeric(campo("age_at_diagnosis")))
age_index_num <- suppressWarnings(as.numeric(campo("age_at_index")))

age_diag_years  <- ifelse(!is.na(age_diag_num),  age_diag_num / 365.25, NA_real_)
age_index_years <- ifelse(!is.na(age_index_num), age_index_num,        NA_real_)

age_years_vec <- ifelse(!is.na(age_diag_years), age_diag_years, NA_real_)
use_index     <- is.na(age_years_vec) & !is.na(age_index_years)
age_years_vec[use_index] <- age_index_years[use_index]

age_source_vec <- ifelse(!is.na(age_diag_years), "diagnosis",
                         ifelse(!is.na(age_index_years), "index", "missing"))
age_source_vec <- factor(age_source_vec, levels = c("diagnosis","index","missing"))

age_obf <- if ("age_is_obfuscated" %in% names(clin_raw)) {
  tolower(as.character(clin_raw$age_is_obfuscated)) %in% c("true","1","yes","y")
} else {
  rep(NA, nrow(clin_raw))
}

# ---------- Stage: consolidar columna cruda antes del mutate ----------
stage_raw <- campo(c("ajcc_pathologic_stage","tumor_stage","pathologic_stage","ajcc_clinical_stage","clinical_stage"))
# ---------- Construcción clínica por paciente ----------
clin0 <- clin_raw %>%
  mutate(
    patient_id = pid, OS.time = tiempo, OS.event = evento,
    age_years  = age_years_vec,
    age_source = age_source_vec,
    age_is_obfuscated = age_obf,
    
    stage_use = norm_stage_I_IV(stage_raw),
    
    gender = norm_gender(sexo),
    smoking_status = smk_status_vec
  ) %>%
  mutate(
    smoking_bin = case_when(
      smoking_status == "Never" ~ "Never",
      smoking_status %in% c("Former","Current") ~ "Ever",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Never","Ever"))
  ) %>%
  select(
    patient_id, OS.time, OS.event,
    age_years, age_source, age_is_obfuscated,
    gender, stage_use, smoking_status, smoking_bin,
    any_of(c("ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m"))
  ) %>%
  distinct()

if (anyDuplicated(clin0$patient_id)) stop("Hay registros clínicos discordantes para un mismo paciente. Revisa clin0 antes de seguir.")

# (Opcional) excluir ofuscados
n_excl_obf <- 0L
if (EXCLUDE_OBFUSCATED_AGE && any(!is.na(clin0$age_is_obfuscated))) {
  n_excl_obf <- sum(clin0$age_is_obfuscated %in% TRUE, na.rm = TRUE)
  clin0 <- clin0 %>% filter(!(age_is_obfuscated %in% TRUE))
}
cat(sprintf("Excluidos por age_is_obfuscated==TRUE: %d\n",
            n_excl_obf),
    file = log_file, append = TRUE)

# ---------- QC y exportes base ----------
qc <- list(
  n_rows              = nrow(clin0),
  n_OS_nonmissing     = sum(!is.na(clin0$OS.time) & !is.na(clin0$OS.event)),
  n_age_nonmissing    = sum(!is.na(clin0$age_years)),
  n_age_from_diag     = sum(clin0$age_source == "diagnosis", na.rm = TRUE),
  n_age_from_index    = sum(clin0$age_source == "index",     na.rm = TRUE),
  n_gender_nonmissing = sum(!is.na(clin0$gender)),
  n_stage_nonmissing  = sum(!is.na(clin0$stage_use)),
  n_smok_nonmissing   = sum(!is.na(clin0$smoking_status) & clin0$smoking_status != "Unknown"),
  n_smok_bin_defined  = sum(!is.na(clin0$smoking_bin))
)
capture.output(str(qc), file = log_file, append = TRUE)

# Resumen QC en tabla
qc_tbl <- tibble(
  metric = names(qc),
  value  = unlist(qc, use.names = FALSE)
)
readr::write_tsv(qc_tbl, file.path(tab_dir, "Clinical_QC_summary.tsv"))

saveRDS(clin0, file.path(out_dir, "clinical_OS.rds"))
readr::write_tsv(clin0, file.path(out_dir, "clinical_OS.tsv"))

# ---------- Tablas de distribución ----------
mk_dist <- function(x, nm){
  tb <- as.data.frame(table(x, useNA = "ifany"), stringsAsFactors = FALSE)
  colnames(tb) <- c(nm, "n")
  tb[[nm]] <- as.character(tb[[nm]])
  tb$perc <- round(100 * tb$n / sum(tb$n), 1)
  tb
}

readr::write_tsv(
  mk_dist(clin0$gender, "gender"),
  file.path(tab_dir, "Clinical_Distribution_gender.tsv")
)
readr::write_tsv(
  mk_dist(clin0$stage_use, "stage"),
  file.path(tab_dir, "Clinical_Distribution_stage.tsv")
)
readr::write_tsv(
  mk_dist(clin0$smoking_status, "smoking4"),
  file.path(tab_dir, "Clinical_Distribution_smoking4.tsv")
)
readr::write_tsv(
  mk_dist(clin0$smoking_bin, "smoking2"),
  file.path(tab_dir, "Clinical_Distribution_smoking2.tsv")
)

# Fuente de la edad (auditoría para “a ojo”)
readr::write_tsv(
  clin0 %>% select(patient_id, age_years, age_source, age_is_obfuscated),
  file.path(tab_dir, "Clinical_Age_source.tsv")
)
readr::write_tsv(
  tibble(
    total_rows          = nrow(clin0),
    n_age_from_diag     = qc$n_age_from_diag,
    n_age_from_index    = qc$n_age_from_index,
    n_age_missing       = sum(is.na(clin0$age_years)),
    excl_age_obfuscated = n_excl_obf
  ),
  file.path(tab_dir, "Clinical_Age_source_counts.tsv")
)

# ---------- Figuras rápidas ----------
plot_counts <- function(df, var, title, fname){
  df$.x <- df[[var]]
  p <- ggplot(df, aes(x = .x)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)),
              vjust = -0.3, size = 3) +
    labs(x = var, y = "Conteo", title = title) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))
  ggsave(file.path(fig_dir, fname), p, width = 7, height = 5, dpi = 300)
}

plot_counts(clin0, "gender",
            "Distribución por género",
            "ClinicalQC_counts_gender.png")

plot_counts(clin0, "stage_use",
            "Distribución por estadio (I–IV)",
            "ClinicalQC_counts_stage.png")

plot_counts(clin0, "smoking_status",
            "Distribución por tabaquismo (4 niv.)",
            "ClinicalQC_counts_smoking4.png")

plot_counts(clin0, "smoking_bin",
            "Tabaquismo binario (Never/Ever)",
            "ClinicalQC_counts_smoking2.png")

# Edad: histograma con N y missing en el título
age_vec   <- suppressWarnings(as.numeric(clin0$age_years))
age_clean <- age_vec[!is.na(age_vec)]
ttl <- sprintf("Edad (años) | N=%d (Missing=%d)",
               length(age_clean), sum(is.na(age_vec)))

p_age_hist <- ggplot(data.frame(age = age_clean), aes(x = age)) +
  geom_histogram(bins = 30, alpha = .85) +
  geom_density(aes(y = after_stat(density))) +
  labs(title = ttl, x = "Edad (años)", y = "Frecuencia / Densidad") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "ClinicalQC_age_hist.png"),
       p_age_hist, width = 7, height = 5, dpi = 300)

# Caja/violín por estadio (edad continua)
p_age_stage <- ggplot(clin0, aes(x = stage_use, y = age_years)) +
  geom_violin(trim = TRUE, na.rm = TRUE) +
  geom_boxplot(width = .15, outlier.size = .8, na.rm = TRUE) +
  labs(title = "Edad por estadio", x = "Estadio (I–IV)", y = "Edad (años)") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "ClinicalQC_age_by_stage.png"),
       p_age_stage, width = 7, height = 5, dpi = 300)

cat(
  "\n[OK] 03 listo.\n",
  "- Pacientes totales: ", nrow(clin0), "\n",
  "- OS definidos:      ", qc$n_OS_nonmissing, "\n",
  "- Edad definida:     ", qc$n_age_nonmissing,
  " (diag=", qc$n_age_from_diag, ", index=", qc$n_age_from_index, ")\n",
  "- Género definido:   ", qc$n_gender_nonmissing, "\n",
  "- Stage definido:    ", qc$n_stage_nonmissing, "\n",
  "- Tabaquismo 4 niv.: ", qc$n_smok_nonmissing, "\n",
  "- Tabaquismo bin.:   ", qc$n_smok_bin_defined, "\n",
  "Tablas: Clinical_QC_summary.tsv, Clinical_Age_source*.tsv, Clinical_Distribution_*.tsv\n",
  "Figuras: ClinicalQC_counts_*.png, ClinicalQC_age_*.png\n",
  file = log_file, append = TRUE
)

cat(sprintf("[OK] Paso 03 | Pacientes: %d | OS definida: %d | Muertes con tiempo: %d | Censuras con tiempo: %d\n", nrow(clin0), qc$n_OS_nonmissing, sum(clin0$OS.event == 1 & !is.na(clin0$OS.time),na.rm=TRUE), sum(clin0$OS.event == 0 & !is.na(clin0$OS.time),na.rm=TRUE)))
