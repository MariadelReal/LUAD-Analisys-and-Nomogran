## 13_LASSO_coef_and_contribution_plots.R
## - Usa SIEMPRE los coeficientes LASSO de lambda.min
## - Genera:
##    1) Forest / coef-plot (beta por gen)
##    2) Barplot de contribución relativa por gen (|beta| / sum |beta|)
## - Exporta también tabla con contribuciones relativas

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

## ---------- paths ----------
base_dir <- getwd()
tab_dir  <- file.path(base_dir, "tables")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, TRUE, TRUE)

lambda_lab <- "lambda.min"
coef_file  <- file.path(base_dir,"results","lambda.min.tsv")

if (!file.exists(coef_file)) {
  stop("No encontré ", coef_file, ". Corre primero 10_lasso_cox.R")
}

## ---------- cargar coeficientes no cero ----------
coefs <- readr::read_tsv(coef_file, show_col_types = FALSE)

stopifnot(all(c("gene","beta") %in% names(coefs)))

coefs <- coefs %>%
  filter(is.finite(beta), beta != 0) %>%
  filter(!is.na(gene), gene != "(Intercept)") %>%
  arrange(beta) %>%
  mutate(gene = as.character(gene))

if (nrow(coefs) == 0L) {
  stop("lambda.min.tsv no tiene betas distintos de 0.")
}

## ---------- 1) Forest / coef-plot ----------
coefs <- coefs %>%
  mutate(gene = factor(gene, levels = gene))   # orden: de menor a mayor beta

p_coef <- ggplot(coefs, aes(x = beta, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = beta, y = gene, yend = gene),
               linewidth = 0.7) +
  geom_point(size = 2.5) +
  labs(
    title = "LASSO coefficients (lambda.min)",
    x     = "Coefficient (beta)",
    y     = "Gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.title.y = element_text(vjust = 1)
  )

fig_coef_png <- file.path(fig_dir, "FIG_LASSO_coefplot_lambda_min.png")
fig_coef_pdf <- file.path(fig_dir, "FIG_LASSO_coefplot_lambda_min.pdf")

ggsave(fig_coef_png, p_coef, width = 7, height = 5, dpi = 300)
ggsave(fig_coef_pdf, p_coef, width = 7, height = 5, device = cairo_pdf)

message("Coef-plot guardado → ", fig_coef_png, " & ", fig_coef_pdf)

## ---------- 2) Contribución relativa por gen ----------
coefs_contrib <- coefs %>%
  mutate(
    abs_beta    = abs(beta),
    rel_weight  = abs_beta / sum(abs_beta),
    rel_percent = 100 * rel_weight
  ) %>%
  arrange(desc(rel_weight)) %>%
  mutate(gene = factor(gene, levels = gene))   # el más importante arriba

out_tbl_path <- file.path(tab_dir, "LASSO_relative_contribution_lambda_min.tsv")
readr::write_tsv(coefs_contrib, out_tbl_path)
message("Tabla de contribuciones relativas → ", out_tbl_path)

p_contrib <- ggplot(coefs_contrib, aes(x = gene, y = rel_percent)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", rel_percent)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Relative magnitude of LASSO coefficients",
    x     = "Gene",
    y     = "Share of total absolute coefficient (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.title.y = element_text(vjust = 1)
  )

fig_contrib_png <- file.path(fig_dir, "FIG_LASSO_relative_contribution_lambda_min.png")
fig_contrib_pdf <- file.path(fig_dir, "FIG_LASSO_relative_contribution_lambda_min.pdf")

ggsave(fig_contrib_png, p_contrib, width = 7, height = 5, dpi = 300)
ggsave(fig_contrib_pdf, p_contrib, width = 7, height = 5, device = cairo_pdf)

message("Barplot de contribución relativa → ", fig_contrib_png, " & ", fig_contrib_pdf)
message("Listo: coef-plot + relative contribution usando lambda.min.")