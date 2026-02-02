## ---------------------------
##
## Script name: 5-get_residuals.R
## Purpose of script: To generate covariate-adjusted metabolomics residuals
##                    using empirical Bayes linear model (WGCNA)
##
## Author: Yufan Gong
##
## Date Created: 2026-01-29
##
## Date Modified: 2026-01-29
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: This script generates residualized metabolomics data after adjusting
##        for covariates using WGCNA's empiricalBayesLM function. The residuals
##        can then be used in downstream MWAS analysis.
##
##        Dependencies: Run scripts 1-4 before this script.
## ---------------------------

# Load required packages
library(WGCNA)
library(tidyverse)

# Define covariates for adjustment -------------------------------------------

# Covariates to adjust for in the metabolomics data
# These should NOT include the exposure of interest (NOx, air toxicants)
myvars_covar <- quote_all(blage, gender, edu_year, smoking_status)


# Prepare covariate matrices -------------------------------------------------

## Merge cleaned SALSA data with link files to match metabolomics samples

### For C18 negative mode
covar_c18 <- salsa_clean |>
  dplyr::mutate(rand_id = as.character(rand_id)) |>
  dplyr::inner_join(
    salsa_mapping_c18neg_surveylinked_links |>
      dplyr::rename_all(str_to_lower) |>
      dplyr::mutate(rand_id = as.character(rand_id)),
    by = "rand_id"
  ) |>
  dplyr::select(rand_id, fullrunname, all_of(myvars_covar)) |>
  # Convert categorical variables to numeric for empiricalBayesLM
  dplyr::mutate(
    gender = as.numeric(as.factor(gender)),
    smoking_status = as.numeric(as.factor(smoking_status))
  ) |>
  dplyr::select(-rand_id) |>
  tibble::column_to_rownames(var = "fullrunname")

### For HILIC positive mode
covar_hilic <- salsa_clean |>
  dplyr::mutate(rand_id = as.character(rand_id)) |>
  dplyr::inner_join(
    salsa_mapping_hilicpos_surveylinked_links |>
      dplyr::rename_all(str_to_lower) |>
      dplyr::mutate(rand_id = as.character(rand_id)),
    by = "rand_id"
  ) |>
  dplyr::select(rand_id, fullrunname, all_of(myvars_covar)) |>
  # Convert categorical variables to numeric for empiricalBayesLM
  dplyr::mutate(
    gender = as.numeric(as.factor(gender)),
    smoking_status = as.numeric(as.factor(smoking_status))
  ) |>
  dplyr::select(-rand_id) |>
  tibble::column_to_rownames(var = "fullrunname")


# Prepare metabolomics matrices ----------------------------------------------

## Load processed metabolomics data if not already loaded
## Assuming metabolomics matrices are stored with samples as columns and features as rows

### For C18 negative mode - subset to matching samples
metabolome_c18 <- matrix_c18_processed |>  # Replace with actual metabolomics matrix variable name

  as.data.frame() |>
  dplyr::select(all_of(rownames(covar_c18)))

### For HILIC positive mode - subset to matching samples
metabolome_hilic <- matrix_hilic_processed |>  # Replace with actual metabolomics matrix variable name
  as.data.frame() |>
  dplyr::select(all_of(rownames(covar_hilic)))


# Verify sample matching -----------------------------------------------------

## Check if rownames of covariate dataset match column names of metabolome matrix
message("C18 sample matching:")
print(table(rownames(covar_c18) == colnames(metabolome_c18)))

message("HILIC sample matching:")
print(table(rownames(covar_hilic) == colnames(metabolome_hilic)))


# Run empirical Bayes linear model -------------------------------------------

## This adjusts metabolomics data for covariates using WGCNA's empiricalBayesLM
## Parameters:
##   - data: transposed metabolomics matrix (samples as rows, features as columns)
##   - removedCovariates: covariate matrix (samples as rownames)
##   - automaticWeights = "bicov": use biweight midcorrelation for robust weighting
##   - aw.maxPOutliers = 0.01: maximum proportion of outliers

message("Running empirical Bayes linear model for C18...")
system.time({
  metabolome_residual_c18 <- WGCNA::empiricalBayesLM(
    data = t(metabolome_c18),
    removedCovariates = covar_c18,
    automaticWeights = "bicov",
    aw.maxPOutliers = 0.01
  )
})

message("Running empirical Bayes linear model for HILIC...")
system.time({
  metabolome_residual_hilic <- WGCNA::empiricalBayesLM(
    data = t(metabolome_hilic),
    removedCovariates = covar_hilic,
    automaticWeights = "bicov",
    aw.maxPOutliers = 0.01
  )
})


# Extract adjusted data ------------------------------------------------------

## The adjusted data is in the $adjustedData component of the result
combined_residual_c18 <- metabolome_residual_c18$adjustedData |>
  as.data.frame() |>
  t() |>
  as.data.frame()

combined_residual_hilic <- metabolome_residual_hilic$adjustedData |>
  as.data.frame() |>
  t() |>
  as.data.frame()


# Verify output dimensions ---------------------------------------------------

message("C18 residual matrix dimensions:")
print(dim(combined_residual_c18))

message("HILIC residual matrix dimensions:")
print(dim(combined_residual_hilic))


# Create ID link datasets ----------------------------------------------------

## These link sample IDs to SALSA participant IDs for downstream analysis
sample_link_c18 <- salsa_clean |>
  dplyr::mutate(rand_id = as.character(rand_id)) |>
  dplyr::inner_join(
    salsa_mapping_c18neg_surveylinked_links |>
      dplyr::rename_all(str_to_lower) |>
      dplyr::mutate(rand_id = as.character(rand_id)),
    by = "rand_id"
  ) |>
  dplyr::filter(fullrunname %in% colnames(combined_residual_c18)) |>
  dplyr::select(rand_id, fullrunname) |>
  dplyr::arrange(rand_id)

sample_link_hilic <- salsa_clean |>
  dplyr::mutate(rand_id = as.character(rand_id)) |>
  dplyr::inner_join(
    salsa_mapping_hilicpos_surveylinked_links |>
      dplyr::rename_all(str_to_lower) |>
      dplyr::mutate(rand_id = as.character(rand_id)),
    by = "rand_id"
  ) |>
  dplyr::filter(fullrunname %in% colnames(combined_residual_hilic)) |>
  dplyr::select(rand_id, fullrunname) |>
  dplyr::arrange(rand_id)


# Save residual data ---------------------------------------------------------

## Save residualized metabolomics data for downstream MWAS analysis
save(combined_residual_c18,
     file = here::here("data", "metabolomics", "processed",
                       "combined_residual_c18.RData"))

save(combined_residual_hilic,
     file = here::here("data", "metabolomics", "processed",
                       "combined_residual_hilic.RData"))

## Save sample link files
save(sample_link_c18, sample_link_hilic,
     file = here::here("data", "metabolomics", "processed",
                       "sample_links.RData"))

## Save covariate matrices for reference
save(covar_c18, covar_hilic,
     file = here::here("data", "metabolomics", "processed",
                       "covar_matrices.RData"))

message("Residual matrices saved successfully!")

#--------------------------------End of the code--------------------------------
