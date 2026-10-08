# Air Toxicant Exposures & Plasma Metabolomics in SALSA

A metabolome-wide association study (MWAS) of long-term residential air toxicant exposure and the plasma metabolome among older Mexican Americans in the Sacramento Area Latino Study on Aging (SALSA).

## Overview

Air pollution is a widespread and regulable exposure linked to cardiometabolic and neurodegenerative disease, but the circulating molecular signatures associated with it in older Latinos, a group underrepresented in omics research, are poorly characterized. This project profiles the plasma metabolome of SALSA participants and tests each metabolic feature for association with long-term exposure to a mixture of eight air toxicants and with each pollutant individually. Associations are also examined by dementia/CIND status.

## Methods in brief

- **Cohort.** SALSA is a prospective cohort of 1,789 community-dwelling Mexican Americans aged 60 and older, enrolled in 1998–1999 and followed through 2007. The analytic sample is 952 participants contributing 1,546 EDTA plasma specimens (one to four per participant). Of these participants, 42 had prevalent and 99 incident dementia/CIND, and 811 remained cognitively intact.
- **Exposure.** Eight pollutants were assigned to residential address histories: five hazardous air pollutants (benzene, 1,3-butadiene, chromium, nickel, lead), traffic-derived NOx (CALINE4), and PM2.5 and NO2 (land-use regression). Each specimen's exposure is the average over the five calendar years before its draw, so exposure varies over time. The primary mixture exposure is an unsupervised index: the first principal component (PC1) of the quartile-scored pollutant matrix, which does not use the outcome. A cross-fitted WQS index serves as a sensitivity analysis. Single pollutants are IQR-scaled.
- **Metabolome.** Untargeted LC-HRMS on C18 (negative ESI) and HILIC (positive ESI) columns in technical triplicate. Features were retained if detected in ≥ 75% of samples (9,703 C18 and 10,422 HILIC), then log2-transformed, median-centered, and ComBat batch-corrected.
- **Association testing.** Pooled repeated-measures MWAS with `limma`, using `duplicateCorrelation` to block on participant. Benjamini–Hochberg FDR < 0.05 is applied within each platform. Covariates with missing values are multiply imputed (m = 10) and pooled with Rubin's rules. PLS-VIP (`mixOmics`) is reported only as an exploratory ranking, not as a selection criterion.
- **Interpretation.** Annotation uses the in-house reference-standard library and xMSannotator, with Schymanski confidence levels. Pathway enrichment uses mummichog (primary) and metapone.
- **Sensitivity analyses.** Cross-fitted WQS index, 3- and 10-year exposure windows, a pre-diagnosis-only population, a sensitivity covariate set, cognitive strata, and complete-case analysis.
- **Cross-cohort comparison.** A fixed-effect meta-analysis with the WHICAP PM2.5 MWAS on aligned features. It is reported as a descriptive comparison, not a replication.

## Analysis pipeline (`scripts/`)

The original pipeline runs end to end via **`0-run_all.R`**, which sources the numbered scripts in order.

| Script | Purpose |
|---|---|
| `0-run_all.R` | Entry point; runs the numbered scripts in sequence. |
| `1-functions.R` | Loads R packages (via `pacman`) and defines the shared helper functions used throughout. |
| `2-load_data.R` | Loads air toxicant exposures, SALSA cohort files, sample-to-participant linkage, and metabolomic feature tables. |
| `3-clean_data.R` | Cleans and imputes covariates, derives the exposures, and assembles the analysis data frames. |
| `4-mwas_analysis.R` | Runs the MWAS (`limma` + `duplicateCorrelation`) and PLS-VIP across exposures, platforms, population strata, and covariate sets. |
| `5-annotation.R` | Annotates MWAS-significant features (in-house library + xMSannotator). |
| `6-pathway_analysis.R` | Pathway enrichment with mummichog (MetaboAnalystR) and metapone. |
| `7-visualization.R` | Volcano, scatter, heatmap, and pathway-enrichment plots. |
| `8-create_table1.R` | Descriptive / Table 1 statistics (`gtsummary`). |
| `9-salsa_whicap_alignment.R` | Aligns SALSA and WHICAP features (m/z and retention time crosswalks). |
| `10-supplement_tables.R`, `10b-supplement_tables_concise.R` | Supplementary result tables. |
| `11-supplement_docx.R` | Builds the Supporting Information `.docx` from the supplement-table workbook. |

Supporting helpers: `wqs_modified.R`, `qgcomp_modified.R` (mixture-exposure models), `mummichog_pathway.R`, `metapone_pathway.R`, `pathway_viz_augment.R`.

## Revision analyses (`scripts/revision/`)

These scripts implement the analyses for the response to reviewers. They mirror the original pipeline but write only under `revision_output/`, which has the same sub-layout as `tables/`, `figures/`, `data/processed/`, `metaboAnalyst/`, and `Metapone/`. The submitted results therefore remain untouched. Each script sources `scripts/1-functions.R` and then `R1-revision_functions.R`.

The submitted WQS and QGcomp indices were estimated by regressing dementia/CIND on the pollutant mixture. The same participants therefore informed both the exposure weights and the metabolite discovery. The revision makes these changes:

- **Primary exposure.** The unsupervised PC1 index (`comp_pca_all`) replaces the outcome-informed indices. The primary analysis uses the full cohort, the primary covariate set, and the five-year window.
- **WQS and QGcomp.** WQS and QGcomp indices are re-estimated with five-fold cross-fitting at the participant level. The cross-fitted WQS index (r = 0.94 with PC1) is kept as a sensitivity analysis. The QGcomp indices did not hold up under cross-fitting, so they were withdrawn from the manuscript; the code still builds them as a diagnostic.
- **Single pollutants.** All eight are refit with the same covariates, imputation, and `duplicateCorrelation` as the mixture index.
- **Covariates and imputation.** `batch` is dropped from the covariates because ComBat already removes it. Single imputation is replaced by multiple imputation with Rubin's rules.

| Script | Mirrors | Purpose |
|---|---|---|
| `R0-run_revision.R` | `0-run_all.R` | Runs R2–R8, each in its own R process with a log. If one step fails, the others still run. |
| `R1-revision_functions.R` | `1-functions.R` | Shared paths (`rev_here()`); exposure, covariate, and population configuration; `make_pca_index()` and `crossfit_composite()`. |
| `R2-exposure_windows_revision.R` | `3-clean_data.R` | Exposure-window sensitivity: alternative averaging windows (3- and 10-year reported; NOx is unavailable for windows shorter than five years), their PC1 indices, and between-window correlations. |
| `R3-composites_revision.R` | `3-clean_data.R` | Builds the PC1 and cross-fitted indices at each window. Also computes index–pollutant correlations and fold-to-fold weight stability. |
| `R4-mwas_revision.R` | `4-mwas_analysis.R` | Pooled repeated-measures MWAS (`duplicateCorrelation` + `lmFit` + `eBayes`) across exposures, cognitive strata, the pre-diagnosis population, and both covariate sets. Significance is FDR only. PLS-VIP is fit on covariate-residualized data for exploratory ranking. The script can restart from per-cell checkpoints. |
| `R5-annotation_revision.R` | `5-annotation.R` | Attaches the existing xMSannotator and in-house annotations. Each annotation gets a Schymanski confidence level, its candidate count, adduct, and ppm error. Candidates implausible for ESI (poly-halogenated formulas, salts) are flagged. The feature table is unchanged, so annotation is reused rather than re-run. |
| `R6-pathway_revision.R` | `6-pathway_analysis.R` | Mummichog (MetaboAnalystR) and metapone at 1,000 permutations; metapone runs in parallel across cells. |
| `R7-visualization_revision.R` | `7-visualization.R` | Volcano, Manhattan, comparison scatter, heatmap, and pathway figures. Also builds the consolidated main-text MWAS figure: volcano, Manhattan, Level 1 effect estimates beside the single pollutants, and exposure–response. |
| `R8-dag_figure.R` | — | Covariate DAG written in dagitty syntax. Computes the minimal sufficient adjustment set and stops if it disagrees with the MWAS covariates. |
| `R9-whicap_concordance.R` | `9-salsa_whicap_alignment.R` | SALSA–WHICAP comparison on 3,325 aligned feature pairs, standardized per SD of exposure and of abundance. Runs a fixed-effect meta-analysis twice: SALSA PC1 with WHICAP PM2.5, and PM2.5 with PM2.5. Reports direction concordance across all aligned features. |
| `R10-imputation_revision.R` | — | Multiple imputation (m = 10) with Rubin's rules and Barnard–Rubin degrees of freedom for the primary MWAS, plus a complete-case sensitivity analysis. |

To run the revision (R2–R8), or only some steps:

```bash
Rscript scripts/revision/R0-run_revision.R
Rscript scripts/revision/R0-run_revision.R R6,R7
```

`R9` and `R10` are run separately, after `R4` has produced the revision MWAS results. For a quick smoke test (results not reportable), set `SALSA_REVISION_QUICK=true` and point `SALSA_REVISION_ROOT` at a throwaway folder. `SALSA_METAPONE_WORKERS` sets the number of parallel metapone workers in `R6`.

## Repository layout

```
scripts/            Numbered R pipeline (0–11) plus mixture-model and pathway helpers
scripts/revision/   Revision analyses (R0–R10); outputs go to revision_output/
qmd/                Quarto report (AP_SALSA.qmd)
LICENSE
```

## Notes

- Raw data, intermediate `.RData`, and exported tables/figures (including everything under `revision_output/`) are produced locally and excluded from version control.
- The metabolomics data are deposited in Metabolomics Workbench (study ST005208). SALSA cohort data are archived at NACDA (ICPSR 22760).
