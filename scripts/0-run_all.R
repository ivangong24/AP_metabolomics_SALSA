## ---------------------------
##
## Script name: 0-run_all.R
## Purpose of script: Master script to run the entire MWAS analysis pipeline
##                    for air toxicants and metabolomics in SALSA
##
## Author: Yufan Gong
##
## Date Created: 2026-01-29
##
## Date Modified: 2026-04-02
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: This master script runs all analysis scripts in sequence.
##        Run this script from the project root directory.
##
##        Pipeline Overview:
##        1-functions.R       - Load packages and define helper functions
##        2-load_data.R       - Load all raw data (exposure, covariates, metabolomics)
##        3-clean_data.R      - Clean and preprocess raw data
##        4-mwas_analysis.R   - MWAS analysis (limma + PLS with VIP scores)
##        5-annotation.R      - Annotate metabolic features (Emory inhouse, HMDB)
##        6-pathway_analysis.R - Pathway analysis (Mummichog via MetaboAnalystR)
##        7-visualization.R   - Create visualizations (volcano, heatmap, scatter)
##        8-create_table1.R   - Create Table 1 (demographic characteristics)
##
## ---------------------------

# Clear environment
rm(list = ls())
gc()

# Set working directory to project root
# setwd(here::here())

# Record start time
start_time <- Sys.time()

cat("\n")
cat("=============================================================\n")
cat("  AIR TOXICANTS & METABOLOMICS ANALYSIS PIPELINE - SALSA\n")
cat("=============================================================\n")
cat("\n")
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")


# =============================================================================
# STEP 1: Load Functions
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 1: Loading functions and packages...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "1-functions.R"))

cat("Functions loaded successfully.\n\n")


# =============================================================================
# STEP 2: Load Data
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 2: Loading all data...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "2-load_data.R"))

cat("Data loaded successfully.\n\n")


# =============================================================================
# STEP 3: Clean Data
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 3: Cleaning SALSA and exposure data...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "3-clean_data.R"))

cat("Data cleaning completed.\n\n")


# =============================================================================
# STEP 4: MWAS Analysis
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 4: Running MWAS analysis (limma + PLS)...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "4-mwas_analysis.R"))

cat("MWAS analysis completed.\n\n")


# =============================================================================
# STEP 5: Annotation
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 5: Annotating metabolic features...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "5-annotation.R"))

cat("Annotation completed.\n\n")


# =============================================================================
# STEP 6: Pathway Analysis
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 6: Running pathway analysis (Mummichog)...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "6-pathway_analysis.R"))

cat("Pathway analysis completed.\n\n")


# =============================================================================
# STEP 7: Create Visualizations
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 7: Creating visualizations...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "7-visualization.R"))

cat("Visualizations created.\n\n")


# =============================================================================
# STEP 8: Create Table 1
# =============================================================================

cat("-------------------------------------------------------------\n")
cat("STEP 8: Creating Table 1 (demographic characteristics)...\n")
cat("-------------------------------------------------------------\n")

source(here::here("scripts", "8-create_table1.R"))

cat("Table 1 created.\n\n")


# =============================================================================
# SUMMARY
# =============================================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("=============================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("=============================================================\n")
cat("\n")
cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total duration:", round(as.numeric(duration), 2), "minutes\n")
cat("\n")
cat("Output locations:\n")
cat("  - Tables: tables/mwas_results/\n")
cat("  - Figures: figures/mwas/\n")
cat("  - Pathway results: metaboAnalyst/\n")
cat("  - Processed data: data/metabolomics/processed/\n")
cat("\n")
cat("Key results files:\n")
cat("  - mwas_results_all.RData\n")
cat("  - mwas_annotation.RData\n")
cat("  - metapone_results_all.RData\n")
cat("  - Table 1 (demographic characteristics)\n")
cat("\n")

# Save session info for reproducibility
session_info <- sessionInfo()
save(session_info,
     file = here::here("data", "session_info.RData"))

cat("Session info saved for reproducibility.\n")
cat("=============================================================\n")

#--------------------------------End of the code--------------------------------
