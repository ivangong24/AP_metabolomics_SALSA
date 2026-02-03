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
##        Dependencies: Run scripts 1-3 before this script.
## ---------------------------

# Load required packages
library(WGCNA)
library(tidyverse)

# Define covariates for adjustment -------------------------------------------

# Covariates to adjust for in the metabolomics data
# These should NOT include the exposure of interest (NOx, air toxicants)
myvars_covar <- quote_all(blage, gender, edu_year, mh62)


# Prepare covariate matrices -------------------------------------------------

## Merge cleaned SALSA data with link files to match metabolomics samples

list(
  list(salsa_mapping_c18neg_surveylinked_03jun2025,
       salsa_mapping_hilicpos_surveylinked_03jun2025),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |> 
  purrr::pmap(function(data1, data2){
    salsa_clean |>
      dplyr::mutate(rand_id = as.character(rand_id)) |>
      dplyr::left_join(salsa_id |>
                         dplyr::mutate(rand_id = as.character(rand_id),
                                       id = as.character(id)), 
                       by = "rand_id") |>
      dplyr::inner_join(
        data1 |>
          dplyr::rename_all(str_to_lower) |>
          dplyr::mutate(id = as.character(id)),
        by = "id"
      ) |>
      dplyr::select(file.name_new, all_of(myvars_covar)) |>
      dplyr::filter(file.name_new %in% colnames(data2)) |>
      # dplyr::select(-rand_id) |>
      tibble::column_to_rownames(var = "file.name_new")
  }) %>% 
  set_names("covar_c18", "covar_hilic") |>
  list2env(.GlobalEnv)


# Prepare metabolomics matrices ----------------------------------------------

## limit metabolite to SALSA samples (remove QCs)

list(
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed),
  list(covar_c18,
       covar_hilic)
) |> 
  pmap(function(data1, data2){
    data1 |>
      as.data.frame() |>
      dplyr::select(any_of(rownames(data2)))  # select only samples present in covariate matrix))
  }) |>
  set_names(c("metabolome_c18", "metabolome_hilic")) |>
  list2env(.GlobalEnv)


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


list(
  list("C18", "HILIC"),
  list(metabolome_c18, metabolome_hilic),
  list(covar_c18, covar_hilic)
) |>
  purrr::pmap(function(mode, metabo, covar){
    message(paste0("Running empirical Bayes linear model for ", mode, "..."))
    WGCNA::empiricalBayesLM(
      data = t(metabo),
      removedCovariates = covar,
      automaticWeights = "bicov",
      aw.maxPOutliers = 0.01
    )
  }) |>
  set_names("metabolome_residual_c18", "metabolome_residual_hilic") |>
  list2env(.,envir = .GlobalEnv)



# Extract adjusted data ------------------------------------------------------

## The adjusted data is in the $adjustedData component of the result
list(metabolome_residual_c18, metabolome_residual_hilic) |>
  map(function(data){
    data$adjustedData |> 
      as.data.frame() |>
      t() |>
      as.data.frame()
  }) %>% 
  set_names("combined_residual_c18", 
            "combined_residual_hilic") %>% 
  list2env(.,envir = .GlobalEnv)



# Verify output dimensions ---------------------------------------------------

list(
  list("C18", "HILIC"),
  list(combined_residual_c18, combined_residual_hilic)
) |>
  purrr::pmap(function(mode, data){
    message(paste0(mode, " residual matrix dimensions:"))
    print(dim(data))
  }) |>
  invisible()

# Create ID link datasets ----------------------------------------------------

## These link sample IDs to SALSA participant IDs for downstream analysis

list(
  list(salsa_mapping_c18neg_surveylinked_03jun2025,
       salsa_mapping_hilicpos_surveylinked_03jun2025),
  list(combined_residual_c18, combined_residual_hilic)
) |>
  pmap(function(mapping, metabo_residual){
    salsa_clean |>
      dplyr::mutate(rand_id = as.character(rand_id)) |>
      dplyr::left_join(salsa_id |>
                         dplyr::mutate(rand_id = as.character(rand_id),
                                       id = as.character(id)), 
                       by = "rand_id") |>
      dplyr::inner_join(
        mapping |>
          dplyr::rename_all(str_to_lower) |>
          dplyr::mutate(id = as.character(id)),
        by = "id"
      ) |>
      dplyr::filter(file.name_new %in% colnames(metabo_residual)) |>
      dplyr::select(rand_id, file.name_new) |>
      dplyr::arrange(rand_id)
  }) |>
  set_names("sample_link_c18", "sample_link_hilic") |>
  list2env(.GlobalEnv)

## check unique rand_id counts with metabolic samples

list(
  list("C18", "HILIC"),
  list(sample_link_c18, sample_link_hilic)
) |>
  purrr::pmap(function(mode, link_data){
    n_rand_id <- n_distinct(link_data$rand_id)
    n_samples <- nrow(link_data)
    message(paste0(mode, " unique rand_id count: ", n_rand_id))
    message(paste0(mode, " total sample count: ", n_samples))
  }) |>
  invisible()



# Save residual data ---------------------------------------------------------

## Save residualized metabolomics and link data for downstream MWAS analysis
list(
  list(combined_residual_c18,
       combined_residual_hilic),
  list(sample_link_c18,
       sample_link_hilic),
  list("c18", "hilic")
) |>
  pmap(function(data1, data2, mode){
    save(data1,
         file = here::here("data", "metabolomics", "processed",
                           paste0("combined_residual_", mode, ".Rdata")))
    message(paste0("Residual matrices for ", mode, " saved successfully!"))
    
    save(data2,
         file = here::here("data", "links", "processed",
                           paste0("sample_link_", mode, ".Rdata")))
    message(paste0("Sample link data for ", mode, " saved successfully!"))
  }) |>
  invisible()


#--------------------------------End of the code--------------------------------
