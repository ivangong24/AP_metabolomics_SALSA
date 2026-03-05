## ---------------------------
##
## Script name: 6-mwas_analysis.R
## Purpose of script: To perform Metabolome-Wide Association Study (MWAS)
##                    for air toxicants and metabolomic profiles in SALSA
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
## Notes: This script performs MWAS analysis using:
##        1. limma for linear model fitting with empirical Bayes
##        2. PLS with VIP scores for feature selection
##
##        Dependencies: Run scripts 1-4 before this script.
##        Key input: Residualized metabolomics data from script 4
## ---------------------------

# Load required packages -----------------------------------------------------

library(limma)
library(mixOmics)
library(tidyverse)
library(writexl)
library(future)
library(furrr)

# Load residualized metabolomics data ----------------------------------------

load(here::here("data", "metabolomics", "processed", "combined_residual_c18.RData"))
load(here::here("data", "metabolomics", "processed", "combined_residual_hilic.RData"))
# load(here::here("data", "links", "processed", "Sample_links.RData"))
load(here::here("data", "links", "processed", "salsa_clean.RData"))
load(here::here("data", "processed", "air_toxicants_exposure.RData"))

# Prepare exposure data for MWAS ---------------------------------------------

## finalize the metabolite matrix

list(
  list(salsa_blood_date_c18, salsa_blood_date_hilic),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |> 
  purrr::pmap(function(data1, data2){
    data1 |> 
      select(rand_id, blood_date, file.name_new) |> 
      filter(file.name_new %in% colnames(data2)) |> 
      filter(!is.na(rand_id)) |> 
      distinct()
  }) |>
  set_names("sample_link_c18", "sample_link_hilic") |>
  list2env(envir = .GlobalEnv)



list(
  list(sample_link_c18, sample_link_hilic),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |> 
  purrr::pmap(function(sample_link, metabo){
    combined_data_list_new[["total"]] |> 
      purrr::map(function(data){
        link <- data |> 
          dplyr::select(rand_id, blood_date) |>
          dplyr::inner_join(sample_link, by = c("rand_id", "blood_date")) |> 
          dplyr::arrange(match(file.name_new, colnames(metabo)))
        
        metabo |> 
          dplyr::select(all_of(link$file.name_new))
      })
  }) |> 
  purrr::set_names("metabo_list_c18_final", "metabo_list_hilic_final") |>
  list2env(.GlobalEnv)
  

## Merge air toxicants exposure with sample link files

list(
  list(sample_link_c18, sample_link_hilic),
  list(metabo_list_c18_final, metabo_list_hilic_final)
) |> 
  purrr::pmap(function(sample_link, metabo_list){
    list(combined_data_list_new[["total"]],
         metabo_list) |> 
      purrr::pmap(function(data, metabo){
        data |> 
          # dplyr::mutate(rand_id = as.character(rand_id)) |>
          dplyr::inner_join(sample_link, by = c("rand_id", "blood_date")) |> 
          dplyr::arrange(match(file.name_new, colnames(metabo)))
      })
  }) |>
  purrr::set_names("combined_data_list_c18", "combined_data_list_hilic") |>
  list2env(.GlobalEnv)

# list(
#   list(sample_link_c18, sample_link_hilic),
#   list(combined_residual_c18, combined_residual_hilic)
# ) |>
#   purrr::pmap(function(sample_link, metabo_residual){
#     sample_link |>
#       dplyr::left_join(air_toxicants_avg_ztrans |>
#                          dplyr::mutate(rand_id = as.character(rand_id)),
#                        by = "rand_id") |>
#       dplyr::arrange(match(file.name_new, colnames(metabo_residual)))
#   }) |>
#   purrr::set_names("exposure_c18", "exposure_hilic") |>
#   list2env(.GlobalEnv)



# Verify sample ordering -----------------------------------------------------

list(
  list("C18", "HILIC"),
  list(combined_data_list_c18, combined_data_list_hilic),
  list(metabo_list_c18_final, metabo_list_hilic_final)
) |>
  purrr::pmap(function(mode, combined_data_list, metabo_list) {
    message(paste0(mode, " sample ordering check:"))
    list(combined_data_list, metabo_list) |> 
      purrr::pmap(function(combined_data, metabo) {
        print(table(combined_data$file.name_new == colnames(metabo)))
      })
  }) |>
  invisible()

# Get exposure variable names ------------------------------------------------

## Extract all exposure variables (air toxicants)
exposure_vars <- combined_data_list_c18[["all"]] |>
  dplyr::select(starts_with("comp_")) |>
  colnames()

covars <- combined_data_list_c18[["all"]] |>
  dplyr::select(all_of(myvars_covar)) |>
  select(-demcind) |> 
  colnames()

message("Exposure variables for MWAS:")
print(exposure_vars)

message("Covariates for MWAS:")
print(covars)

# =============================================================================
# SECTION 1: LIMMA-BASED MWAS
# =============================================================================

# Create design matrices for each exposure -----------------------------------

## Function to create design matrix for a single exposure
create_design_matrix <- function(combined_data, exposure_var, covars) {
  formula_matrix <- as.formula(str_c("~ ", exposure_var, " + ", 
                                     paste(covars, collapse = " + ")))
  model.matrix(formula_matrix, data = combined_data)
}

formula_matrix <- as.formula(str_c("~ ", "comp_wqs_all", " + ", 
                                   paste(covars, collapse = " + ")))

test <- model.matrix(formula_matrix, data = combined_data_list_c18[["all"]])

## Create design matrices for all exposures
list(combined_data_list_c18, combined_data_list_hilic) |>
  purrr::map(function(combined_data_list){
    combined_data_list |> 
       purrr::map(function(combined_data){
         exposure_vars |>
           purrr::set_names() |>
           purrr::map(~ create_design_matrix(combined_data, .x, covars))
       })
  }) |>
  purrr::set_names("design_c18_list", "design_hilic_list") |>
  list2env(.GlobalEnv)


# Estimate within-subject correlation for duplicate measures -------------------

## Set up parallel backend
n_workers <- parallelly::availableCores() - 1
message(paste0("Setting up parallel plan with ", n_workers, " workers..."))
plan(multisession, workers = n_workers)

system.time({
  list(
    list("C18", "HILIC"),
    list(metabo_list_c18_final, metabo_list_hilic_final),
    list(design_c18_list, design_hilic_list),
    list(combined_data_list_c18, combined_data_list_hilic)
  ) |>
    purrr::pmap(function(mode, metabo_list, design_list, combined_data_list) {
      list(design_list, metabo_list, 
           combined_data_list, names(combined_data_list)) |> 
        purrr::pmap(function(designls, metabo, combined_data, population){
          message(paste0("Calculating correlations for ", mode, 
                         " in ", population, " ..."))
          block <- combined_data$rand_id
          stopifnot(
            length(block) == ncol(metabo),
            all(combined_data$file.name_new == colnames(metabo)))
          
          designls |>
            furrr::future_map(function(design) {
              limma::duplicateCorrelation(metabo, design, block = block)
            }, .options = furrr_options(seed = TRUE))
        })
    }) |>
    purrr::set_names("dupcor_c18_list", "dupcor_hilic_list") |>
    list2env(.GlobalEnv)
})

save(dupcor_c18_list, dupcor_hilic_list,
     file = here::here("data", "processed", 
                       "duplicate_correlation_results.RData"))

# Fit limma models -----------------------------------------------------------

## Function to fit limma model with duplicate correlation
fit_limma <- function(metabolome_matrix, design_matrix, block, correlation) {
  fit <- limma::lmFit(metabolome_matrix, design_matrix,
                      block = block, correlation = correlation)
  fit <- limma::eBayes(fit)
  return(fit)
}

list(
  list("C18", "HILIC"),
  list(design_c18_list, design_hilic_list),
  list(combined_residual_c18, combined_residual_hilic),
  list(dupcor_c18, dupcor_hilic),
  # list(dupcor_c18_list, dupcor_hilic_list),
  list(exposure_c18, exposure_hilic)
) |>
  purrr::pmap(function(mode, design_list, metabo_residual,
                        dupcor, exposure_data) {
    message(paste0("Fitting limma models for ", mode, "..."))
    block <- exposure_data$rand_id
    
    design_list |>
      purrr::map(function(design) {
        fit_limma(metabo_residual, design,
                  block = block,
                  correlation = dupcor$consensus.correlation)
      })
    # list(design_list, dupcor_list) |>
    #   furrr::future_pmap(function(design, dupcor) {
    #     fit_limma(metabo_residual, design,
    #               block = block,
    #               correlation = dupcor$consensus.correlation)
    #   }, .options = furrr_options(seed = TRUE))
  }) |>
  purrr::set_names("limma_fit_c18", "limma_fit_hilic") |>
  list2env(.GlobalEnv)

## Reset to sequential plan
plan(sequential)

save(limma_fit_c18, file = here::here("data", "metabolomics", 
                                      "results", "limma_fit_c18.RData"))

save(limma_fit_hilic, file = here::here("data", "metabolomics", 
                                       "results", "limma_fit_hilic.RData"))

# Extract MWAS results -------------------------------------------------------

## Function to extract topTable results
extract_toptable <- function(fit, design, metabolome_matrix) {
  limma::topTable(
    fit,
    coef = ncol(design),
    sort.by = "p",
    number = nrow(metabolome_matrix),
    adjust.method = "BH"  # Benjamini-Hochberg FDR adjustment
  )
}

list(
  list("C18", "HILIC"),
  list(limma_fit_c18, limma_fit_hilic),
  list(design_c18_list, design_hilic_list),
  list(combined_residual_c18, combined_residual_hilic)
) |>
  purrr::pmap(function(mode, fit_list, design_list, metabo_residual){
    message(paste0("Extracting MWAS results for ", mode, "..."))
    list(fit_list, design_list) |>
      purrr::pmap(function(fit, design) {
        extract_toptable(fit, design, metabo_residual)
      }) |>
      purrr::set_names(exposure_vars)
  }) |>
  purrr::set_names("mwas_results_c18", "mwas_results_hilic") |>
  list2env(.GlobalEnv) |>
  invisible()


# Summarize significant results ----------------------------------------------

## Function to summarize significant metabolites
summarize_significant <- function(results_list, fdr_threshold = 0.05) {
  results_list |>
    purrr::imap(function(tbl, exp_name) {
      sig_count <- sum(tbl$adj.P.Val < fdr_threshold, na.rm = TRUE)
      tibble(
        exposure = exp_name,
        n_significant = sig_count,
        n_total = nrow(tbl),
        pct_significant = round(sig_count / nrow(tbl) * 100, 2)
      )
    }) |>
    purrr::list_rbind()
}

message("\nC18 MWAS Summary (FDR < 0.05):")
print(summarize_significant(mwas_results_c18))

message("\nHILIC MWAS Summary (FDR < 0.05):")
print(summarize_significant(mwas_results_hilic))


# =============================================================================
# SECTION 2: PLS-DA WITH VIP SCORES
# =============================================================================

# Prepare metabolomics matrices for PLS --------------------------------------

## Transpose residual matrices for PLS (samples as rows)
list(
  list("C18", "HILIC"),
  list(combined_residual_c18, combined_residual_hilic),
  list(sample_link_c18, sample_link_hilic)
) |>
  purrr::pmap(function(mode, metabo_residual, sample_link){
    message(paste0("Preparing metabolomics matrix for ", mode, "..."))
    metabo_residual |>
      t() |>
      as.data.frame() |>
      tibble::rownames_to_column("file.name_new") |>
      dplyr::right_join(sample_link, by = "file.name_new") |>
      dplyr::select(-rand_id) |>
      tibble::column_to_rownames("file.name_new") |>
      as.matrix()
  }) |>
  purrr::set_names("metabo_matrix_c18", "metabo_matrix_hilic") |>
  list2env(.GlobalEnv)


# Fit PLS models for each exposure -------------------------------------------

## Function to fit PLS and extract VIP
fit_pls_vip <- function(X, Y, ncomp = 3) {
  # Remove NA values
  valid_idx <- !is.na(Y)
  X_valid <- X[valid_idx, ]
  Y_valid <- Y[valid_idx]

  # Fit PLS
  pls_fit <- mixOmics::pls(X_valid, Y_valid, ncomp = ncomp)

  # Extract VIP scores
  vip_scores <- mixOmics::vip(pls_fit) |>
    as.data.frame() |>
    arrange(desc(comp1))

  return(list(pls_fit = pls_fit, vip = vip_scores))
}

## Prepare exposure vectors for PLS

# check if the rownames of the metabolome matrix match the file.name_new of the exposure data
table(exposure_c18$file.name_new==rownames(metabo_matrix_c18))
table(exposure_hilic$file.name_new==rownames(metabo_matrix_hilic))

list(
  list("C18", "HILIC"),
  list(exposure_c18, exposure_hilic),
  list(metabo_matrix_c18, metabo_matrix_hilic)
) |>
  purrr::pmap(function(mode, exposure_data, metabo_matrix){
    message(paste0("Fitting PLS models for ", mode, "..."))
    exposure_vars |>
      purrr::set_names() |>
      purrr::map(function(exp_var) {
        fit_pls_vip(
          X = metabo_matrix,
          Y = exposure_data[[exp_var]],
          ncomp = 3
          )
        })
  }) |>
  purrr::set_names("pls_results_c18", "pls_results_hilic") |>
  list2env(.GlobalEnv)



# Extract VIP scores ---------------------------------------------------------

vip_c18 <- pls_results_c18 |>
  purrr::map(~ .x$vip)

vip_hilic <- pls_results_hilic |>
  purrr::map(~ .x$vip)


# Identify high-VIP metabolites ----------------------------------------------

## Function to count VIP > 2 metabolites
count_high_vip <- function(vip_list, threshold = 2) {
  vip_list |>
    purrr::imap(function(vip_df, exp_name) {
      high_vip_count <- sum(vip_df$comp1 > threshold, na.rm = TRUE)
      tibble(
        exposure = exp_name,
        n_vip_gt_2 = high_vip_count,
        n_total = nrow(vip_df)
      )
    }) |>
    purrr::list_rbind()
}

message("\nC18 VIP > 2 Summary:")
print(count_high_vip(vip_c18))

message("\nHILIC VIP > 2 Summary:")
print(count_high_vip(vip_hilic))


# =============================================================================
# SECTION 3: COMBINE LIMMA AND VIP RESULTS
# =============================================================================

# Function to combine MWAS and VIP results -----------------------------------

combine_mwas_vip <- function(mwas_results, vip_results, annotation_df = NULL) {
  list(mwas_results, vip_results) |>
    purrr::pmap(function(mwas, vip) {
      mwas |>
        tibble::rownames_to_column("met") |>
        dplyr::left_join(
          vip |>
            tibble::rownames_to_column("met") |>
            dplyr::select(met, VIP_comp1 = comp1, VIP_comp2 = comp2, VIP_comp3 = comp3),
          by = "met"
        ) |>
        dplyr::arrange(adj.P.Val)
    })
}

## Combine results for C18
combined_results_c18 <- combine_mwas_vip(mwas_results_c18, vip_c18) |>
  purrr::set_names(exposure_vars)

## Combine results for HILIC
combined_results_hilic <- combine_mwas_vip(mwas_results_hilic, vip_hilic) |>
  purrr::set_names(exposure_vars)


# Filter significant metabolites ---------------------------------------------

## Function to filter significant metabolites
filter_significant <- function(combined_results, fdr_thresh = 0.05, 
                               vip_thresh = 2) {
  combined_results |>
    purrr::map(function(df) {
      df |>
        dplyr::filter(adj.P.Val < fdr_thresh | VIP_comp1 > vip_thresh) |>
        dplyr::arrange(adj.P.Val)
    })
}

significant_c18 <- filter_significant(combined_results_c18)
significant_hilic <- filter_significant(combined_results_hilic)



# =============================================================================
# SECTION 4: SAVE RESULTS
# =============================================================================

# Create output directories --------------------------------------------------

dir.create(here::here("tables", "mwas_results"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("data", "metabolomics", "results"), showWarnings = FALSE, recursive = TRUE)

# Save MWAS results to Excel -------------------------------------------------

## Save combined MWAS + VIP results

list(
  list(combined_results_c18, combined_results_hilic),
  list("c18", "hilic")
) |>
  pmap(function(datalist, mode){
    message(paste0("Saving MWAS results for ", mode, " ..."))
    datalist |>
      purrr::imap(function(df, exp_name) {
        writexl::write_xlsx(
          df,
          path = here::here("tables", "mwas_results",
                            glue::glue("mwas_{mode}_{exp_name}.xlsx"))
        )
      })
  })


# Save significant metabolites -----------------------------------------------

list(
  list(significant_c18, significant_hilic),
  list("c18", "hilic")
) |>
  pmap(function(datalist, mode){
    message(paste0("Saving significant MWAS results for ", mode, " ..."))
    datalist |>
      purrr::imap(function(df, exp_name) {
        if (nrow(df) > 0) {
          writexl::write_xlsx(
            df,
            path = here::here("tables", "mwas_results",
                              glue::glue("mwas_{mode}_{exp_name}_sig.xlsx"))
          )
        }
      })
  })



# Save R objects for downstream analysis -------------------------------------

save(mwas_results_c18, mwas_results_hilic,
     vip_c18, vip_hilic,
     combined_results_c18, combined_results_hilic,
     significant_c18, significant_hilic,
     limma_fit_c18, limma_fit_hilic,
     file = here::here("data", "metabolomics", "results",
                       "mwas_results_all.RData"))

message("MWAS analysis completed! Results saved to tables/mwas_results/")


# =============================================================================
# SECTION 5: SUMMARY TABLE
# =============================================================================

# Create summary table for all exposures -------------------------------------

create_summary_table <- function(mwas_c18, mwas_hilic, vip_c18, vip_hilic) {
  exposure_vars <- names(mwas_c18)

  exposure_vars |>
    purrr::map(function(exp) {
      tibble(
        exposure = exp,
        c18_total = nrow(mwas_c18[[exp]]),
        c18_sig_fdr05 = sum(mwas_c18[[exp]]$adj.P.Val < 0.05, na.rm = TRUE),
        c18_sig_fdr10 = sum(mwas_c18[[exp]]$adj.P.Val < 0.10, na.rm = TRUE),
        c18_vip_gt2 = sum(vip_c18[[exp]]$comp1 > 2, na.rm = TRUE),
        hilic_total = nrow(mwas_hilic[[exp]]),
        hilic_sig_fdr05 = sum(mwas_hilic[[exp]]$adj.P.Val < 0.05, na.rm = TRUE),
        hilic_sig_fdr10 = sum(mwas_hilic[[exp]]$adj.P.Val < 0.10, na.rm = TRUE),
        hilic_vip_gt2 = sum(vip_hilic[[exp]]$comp1 > 2, na.rm = TRUE)
      )
    }) |>
    purrr::list_rbind()
}

summary_table <- create_summary_table(
  mwas_results_c18, mwas_results_hilic,
  vip_c18, vip_hilic
)

message("\nMWAS Summary Table:")
print(summary_table)

# writexl::write_xlsx(
#   summary_table,
#   path = here::here("tables", "mwas_results", "mwas_summary_table.xlsx")
# )

#--------------------------------End of the code--------------------------------
