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

Supporting helpers: `wqs_modified.R`, `qgcomp_modified.R` (mixture-exposure models), `mummichog_pathway.R`, `metapone_pathway.R`, `pathway_viz_augment.R`.

## Repository layout

```
scripts/   Numbered R pipeline (0–10) plus mixture-model and pathway helpers
qmd/       Quarto report (AP_SALSA.qmd)
LICENSE
```

## Notes

- Raw data, intermediate `.RData`, and exported tables/figures are produced locally and excluded from version control.
- Do not re-run the full pipeline to inspect results — read the saved outputs instead.
- The analytical infrastructure (covariates, metabolomics pipeline, MWAS code) is shared with the companion `sleep_metabolomics_salsa` project.
