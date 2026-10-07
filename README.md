# Air Pollution & Plasma Metabolomics in SALSA

A metabolome-wide association study (MWAS) of long-term air pollution exposure and the plasma metabolome among older Mexican Americans in the Sacramento Area Latino Study on Aging (SALSA).

## Overview

Air pollution is a widespread and regulable exposure tied to cardiometabolic and neurodegenerative disease, yet the circulating molecular signatures that connect it to health in older Latinos — an underrepresented group in omics research — remain poorly characterized. This project profiles the plasma metabolome of SALSA participants and tests, feature by feature, its association with long-term ambient and traffic-related air pollution. The aims are to:

1. identify circulating metabolites and metabolic pathways associated with individual pollutants and with the air-pollution mixture, and
2. establish the basis for downstream work relating air-pollution–associated metabolites to cardiometabolic and cognitive aging, including cross-cohort alignment with WHICAP.

## Methods in brief

- **Cohort.** SALSA, a prospective cohort of 1,789 community-dwelling Mexican Americans aged 60+ followed from 1998 to 2007; the analytic sample is the participants with untargeted plasma metabolomics.
- **Exposure.** Long-term ambient and traffic-related air pollution (e.g., NOx, ozone, traffic indicators) assigned to residence, modeled both as single pollutants and as a **mixture** via weighted quantile sum (WQS) regression and quantile g-computation (QGcomp).
- **Metabolome.** Untargeted high-resolution metabolomics on two complementary LC-HRMS platforms — C18 (negative ionization) and HILIC (positive ionization).
- **Association testing.** Per-feature empirical-Bayes linear models (`limma`) with within-participant correlation modeled via `duplicateCorrelation`, complemented by a PLS variable-importance ranking (`mixOmics`); significance at FDR < 0.05.
- **Interpretation.** Metabolite annotation (in-house reference-standard library + xMSannotator) and pathway enrichment (mummichog and metapone).
- **Replication / downstream.** Cross-cohort alignment with the WHICAP cohort and supplementary sensitivity analyses.

## Analysis pipeline (`scripts/`)

The pipeline runs end to end via **`0-run_all.R`**, which sources the numbered scripts in order.

| Script | Purpose |
|---|---|
| `0-run_all.R` | Entry point; runs the numbered scripts in sequence. |
| `1-functions.R` | Loads R packages (via `pacman`) and defines the shared helper functions used throughout. |
| `2-load_data.R` | Loads air-pollution exposures, SALSA cohort files, sample-to-participant linkage, and metabolomic feature tables. |
| `3-clean_data.R` | Cleans and imputes covariates, derives the exposures, and assembles the analysis data frames. |
| `4-mwas_analysis.R` | Runs the MWAS (`limma` + `duplicateCorrelation`) and PLS-VIP across exposures, platforms, population strata, and covariate sets. |
| `5-annotation.R` | Annotates MWAS-significant features (in-house library + xMSannotator). |
| `6-pathway_analysis.R` | Pathway enrichment with mummichog (MetaboAnalystR) and metapone. |
| `7-visualization.R` | Volcano, scatter, heatmap, and pathway-enrichment plots. |
| `8-create_table1.R` | Descriptive / Table 1 statistics (`gtsummary`). |
| `9-salsa_whicap_alignment.R` | Aligns SALSA findings with the WHICAP cohort. |
| `10-supplement_tables.R`, `10b-supplement_tables_concise.R` | Supplementary result tables. |
| `11-supplement_docx.R` | Builds the Supporting Information `.docx` from the supplement-table workbook. |

Supporting helpers: `wqs_modified.R`, `qgcomp_modified.R` (mixture-exposure models), `mummichog_pathway.R`, `metapone_pathway.R`, `pathway_viz_augment.R`.

## Revision analyses (`scripts/revision/`)

Scripts written in response to peer review (Reviewer 1 and Reviewer 2). They mirror the primary pipeline but write **only** under `revision_output/` (same sub-layout as `tables/`, `figures/`, `data/processed/`, `metaboAnalyst/`, `Metapone/`), so the submitted results stay untouched as the record of the original analysis. Each script sources `scripts/1-functions.R` and then `R1-revision_functions.R`.

The central change is an **outcome-free exposure definition**. The submitted WQS/QGcomp indices were estimated by regressing dementia/CIND on the pollutant mixture, so the same participants informed both the weights and the metabolite discovery. The revision replaces them with:

- **Unsupervised PCA index** (`comp_pca_all`, primary) — PC1 of the same quartile-scored pollutant matrix; never touches the outcome.
- **Feature-wise QGcomp** (`comp_qgcomp_fw_all`, secondary) — each metabolite is the outcome and the pollutants the exposures, so there is no outcome model at any stage.
- **Cross-fitted WQS / QGcomp** (`comp_*_cf_*`) — weights estimated in K−1 participant-level folds and applied to the held-out fold.

Single pollutants are refit under the same covariates, imputation, and `duplicateCorrelation` as the composites, at 5-, 3-, and 10-year exposure windows.

| Script | Mirrors | Purpose |
|---|---|---|
| `R0-run_revision.R` | `0-run_all.R` | Runs R2–R8, each in its own R process with a log; one failing step does not abort the run. |
| `R1-revision_functions.R` | `1-functions.R` | Shared paths (`rev_here()`), exposure/covariate/population configuration, `make_pca_index()`, `crossfit_composite()`. |
| `R2-exposure_windows_revision.R` | `3-clean_data.R` | Exposure-window sensitivity: 1-, 3-, and 10-year averages, their PCA indices, and between-window correlations (R1 comment 13). |
| `R3-composites_revision.R` | `3-clean_data.R` | Builds the PCA and cross-fitted indices at each window, plus agreement and window-concordance diagnostics (R1 comment 1, R2 comment 1). |
| `R4-mwas_revision.R` | `4-mwas_analysis.R` | limma MWAS (`duplicateCorrelation` + `lmFit` + `eBayes`) and PLS-VIP across exposures, strata, and covariate sets; restartable from per-cell checkpoints. |
| `R5-annotation_revision.R` | `5-annotation.R` | Attaches the existing xMSannotator / in-house annotations with traceable identification confidence, proposed ion, and adduct (R1 comment 8). The feature space is unchanged, so annotation is reused rather than re-run. |
| `R6-pathway_revision.R` | `6-pathway_analysis.R` | Mummichog (MetaboAnalystR) and metapone, metapone parallelised across cells. |
| `R7-visualization_revision.R` | `7-visualization.R` | Volcano, Manhattan, VIP, comparison scatter, heatmap, pathway plots, and the revised main-text combined panel (R1 comment 14). |
| `R8-dag_figure.R` | — | Covariate DAG in dagitty syntax; computes the minimal sufficient adjustment set and stops if it disagrees with the MWAS covariates (R1 comment 12, R2 comment 5). |
| `R9-whicap_concordance.R` | `9-salsa_whicap_alignment.R` | SALSA–WHICAP concordance restricted to PM2.5 vs PM2.5, so both cohorts estimate the same contrast (R1 comment 10, R2 comment 8). |
| `R10-imputation_revision.R` | — | Multiple imputation with Rubin's rules (m = 10) for the primary MWAS, replacing single imputation, plus a complete-case sensitivity analysis (R1 comment 6, R2 comment 5). |

Run the full revision (R2–R8), or a subset of steps:

```bash
Rscript scripts/revision/R0-run_revision.R
Rscript scripts/revision/R0-run_revision.R R6,R7
```

`R9` and `R10` run standalone (`Rscript scripts/revision/R9-whicap_concordance.R`) after `R4` has produced the revision MWAS results. For a quick smoke test (not reportable), set `SALSA_REVISION_QUICK=true` and point `SALSA_REVISION_ROOT` at a throwaway folder; `SALSA_METAPONE_WORKERS` controls metapone parallelism in `R6`.

## Repository layout

```
scripts/            Numbered R pipeline (0–11) plus mixture-model and pathway helpers
scripts/revision/   Peer-review revision analyses (R0–R10); outputs go to revision_output/
qmd/                Quarto report (AP_SALSA.qmd)
LICENSE
```

## Notes

- Raw data, intermediate `.RData`, and exported tables/figures (including everything under `revision_output/`) are produced locally and excluded from version control.
- Do not re-run the full pipeline to inspect results — read the saved outputs instead.
- The analytical infrastructure (covariates, metabolomics pipeline, MWAS code) is shared with the companion `sleep_metabolomics_salsa` project.
