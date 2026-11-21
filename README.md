# LUAD age-stratified analysis and prognostic nomogram

This repository contains paired GEO transcriptomic analyses (GSE75037 and GSE115002), age-stratified enrichment and network prioritisation, and TCGA-LUAD survival modelling using a 30-gene candidate panel. The original fitted LASSO signature contains 10 genes at lambda.min; lambda.1se selects none. The prognostic model requires external validation.

## Run one script at a time in RStudio

Keep this complete directory structure. For separate GEO and TCGA output directories, make two complete copies named `01_GEO` and `02_TCGA`. Open GEO.Rproj in the first and TCGA.Rproj in the second. The Rproj file sets the working directory; all `source()` calls are issued from that project root. Scripts with the same numeric prefix belong to different workflows, not interchangeable versions. `helpers.R` is loaded by scripts and is not a numbered analysis step.

Example (adapt the path only):
```r
setwd("C:/Users/delre/OneDrive/LUAD Royal Society/Scripts/02_TCGA")
stopifnot(file.exists("TCGA.Rproj"))
source("00_env_check_and_sessioninfo.R", encoding="UTF-8")
# After it finishes, run the next file from the TCGA list below.
```

GEO sequence:

1. `00_check_platforms.R`
2. `01_download_geo_series_matrix.R`
3. `02_initial_expression_qc.R`
4. `03_build_metadata_and_pairs.R`
5. `04_check_pairs.R`
6. `05_collapse_probes_to_symbols.R`
7. `06_paired_differential_expression.R`
8. `07_volcano_and_heatmap.R`
9. `08_consensus_deg_lists.R`
10. `09_plot_ora_common.R`
11. `10_gsea_common_genes.R`
12. `11_ssgsea_paired_entrez.R`
13. `12_ssgsea_consensus_plots.R`

TCGA sequence:

1. `00_env_check_and_sessioninfo.R`
2. `01_download_gdc_expression_clinical.R`
3. `02_harmonize_expression_symbols.R`
4. `03_clinical_qc_and_os_construction.R`
5. `04_map_samples_and_apply_survival_filters.R`
6. `05_build_design_matrices.R`
7. `06_cox_univariate_30genes.R`
8. `07_signature_ingest_and_validation.R`
9. `08_build_design_matrices.R`
10. `09_diagnose_signs_and_alignment.R`
11. `10_lasso_cox.R`
12. `11_cox_firma_clinica.R`
13. `12_nomogram.R`
14. `13_LASSO_coef_and_contribution_plots.R`
15. `14_ROC.R`
16. `15_Bootstrap_seleccion.R`

Step 05 TCGA is retained as a preliminary diagnostic builder and writes only to `data/05_preliminar/`. It does not replace step 08, whose matrices feed the final models. The pre-specified 30-gene CSV is included in `sig/`, so step 05 need not wait for step 07 to create that file. Step 09 is diagnostic. No numbered script has been removed or renumbered.

Do not advance beyond TCGA step 01 until `data_raw/expression_SE_LUAD.rds` and `data_raw/clinical_LUAD.tsv` exist and can be read. An empty directory after a GDC timeout is not a completed download. New live GDC releases may change sample counts; compare them with the archived cohort instead of forcing a particular N.

Dependencies used across the workflows include CRAN packages readr, dplyr, tibble, tidyr, stringr, ggplot2, glue, patchwork, pheatmap, ggrepel, ggvenn, purrr, scales, glmnet, survival, rms, timeROC, msigdbr and sessioninfo; Bioconductor packages GEOquery, Biobase, limma, GSVA, AnnotationDbi, org.Hs.eg.db, clusterProfiler, fgsea, BiocParallel, TCGAbiolinks and SummarizedExperiment. `car` is optional for the VIF diagnostic. Use the supplied session records to identify the versions from the reported runs; installing current packages is not a frozen historical environment.

## Reference inputs and figure reproduction

`referencia/clinical.tsv` is the supplied 491-patient clinical reference used by step 11 to compare a new cohort. It is not a substitute for downloading clinical data. `referencia/folds.tsv` preserves the assignment used by step 10. The script aligns folds by barcode; if a newly downloaded cohort differs, its fallback generates new folds and issues a warning, so results need not be identical. `sig/CellCycle_LUAD_v1.csv` contains the pre-specified 30 candidates, not the ten fitted coefficients.

The reported cohort contains 491 patients/180 deaths; clinical complete cases contain 474/175. Ten patients lack age and seven lack stage, without overlap; 13 complete cases retain Unknown smoking as a category. The supplied CYS sessions confirm MCODE modules of 67 nodes/2052 edges and 76 nodes/2578 edges, 35 MCC hubs per stratum and 30 shared canonical symbols. KIF4A corresponds to STRING 9606.ENSP00000363524; the legacy internal shared name KIF4 in the young session is documented in the audit.

To reproduce the submitted figures without re-running enrichment or models:
- `reproducir_figura2/`: run `01_regenerar_ORA.R`, then `02_montar_figura2_y_suplemento.R` in that directory. The latter preserves the supplied GSEA curve PDF; the full ranking is needed to recompute curves.
- `reproducir_figura3/`: run `python montar_figura3.py`; output TIFF/PNG is 600 ppi. Sources are original Cytoscape PNG exports.
- `reproducir_figura4/`: run `python montar_figura4.py`; it reassembles the supplied model-output panels at 300 dpi and separates overlapping labels.

The figure-assembly code is additional to the numbered analysis scripts. Bootstrap concordance validation repeats gene selection and tuning, whereas the displayed five-year calibration holds the fitted signature fixed. Seven of 200 successful selection-bootstrap replicates selected no genes; their integrated model reduced to the clinical model. The original fitted signature still has ten genes.

ORA Reactome/KEGG Legacy is configured for MSigDB 2026.1.Hs; GSEA Hallmark uses 2025.1.Hs. GO uses org.Hs.eg.db. Package versions and gene-set releases are different identifiers. This archive does not contain every raw download or every historical gene-set snapshot and does not establish a clean end-to-end rerun in the recorded Windows environment.

## Provenance and publication status

Numbered scripts were obtained from `MariadelReal/LUAD-Analisys-and-Nomogran`, commit `21f5598a473f01828d96768dacf6ba02a4d164c0`. The additions in this delivery are README/project setup, the candidate-gene CSV, clinical reference and figure reproduction folders. The numbered analysis scripts are preserved. These reproducibility additions are published alongside the original numbered scripts.

A software licence has not been selected on behalf of the authors. Authors should choose and publish an appropriate licence before release; this README does not grant one. The original public datasets retain their own access and reuse terms.
