## ---------------------------
##
## Script name: get_residuals.R
## Purpose of script: To get residuals from the linear regression model
##
## Author: Yufan Gong
##
## Date Created: 2025-12-08
##
## Date Modified: 2025-12-08
##
## Copyright (c) Yufan Gong, 2025
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: 


# Metabolome matrix residuals ---------------------------------------------

## Y ~ covariates

myvars_covar <- quote_all(blage, gender, bl_bmi)

#select covariates

covar_salsa <- salsa2_ap |> 
  dplyr::select(rand_id, all_of(myvars_covar))


#Generate the metabolome matrix for the study sample (n = 1564)

list(
  list(
    covar_c18_combined,
    covar_hilic_combined
  ),
  list(matrix_c18_f_lod_med_combat_peg_date, 
       matrix_hilic_f_lod_med_combat_peg_date_pc)
) %>% 
  pmap(function(covarlist, metabo_matrix){
    covarlist %>% 
      map(function(covar){
        metabo_matrix %>% 
          as.data.frame() %>%
          select(all_of(covar %>% rownames()))
      })
  }) %>% 
  set_names("metabolome_list_c18", "metabolome_list_hilic") %>%
  list2env(.,envir = .GlobalEnv)

#check if rownames of covar dataset match column names of metabolome matrix
table(rownames(covar_c18_combined[["case"]])==colnames(metabolome_list_c18[[1]]))
table(rownames(covar_c18_combined[["ctrl"]])==colnames(metabolome_list_c18[[2]]))
table(rownames(covar_c18_combined[["total"]])==colnames(metabolome_list_c18[[3]]))

table(rownames(covar_hilic_combined[["case"]])==colnames(metabolome_list_hilic[[1]]))
table(rownames(covar_hilic_combined[["ctrl"]])==colnames(metabolome_list_hilic[[2]]))
table(rownames(covar_hilic_combined[["total"]])==colnames(metabolome_list_hilic[[3]]))

####run empirical Bayes linear model####
library(WGCNA)


# run empirical Bayes linear model (parallel computation needed)
system.time({
  list(
    list(metabolome_list_c18, 
         covar_c18_combined
    ),
    list(metabolome_list_hilic, 
         covar_hilic_combined)
  ) %>% 
    map(function(datalist){
      datalist %>% 
        pmap(function(metabo_matrix, covar){
          empiricalBayesLM(
            data = t(metabo_matrix),
            removedCovariates = covar,
            automaticWeights = "bicov",
            aw.maxPOutliers = 0.01)
        })
    }) %>% 
    set_names("metabolome_residual_list_c18", "metabolome_residual_list_hilic") %>% 
    list2env(.,envir = .GlobalEnv)
})

list(metabolome_residual_list_c18, metabolome_residual_list_hilic) %>% 
  map(function(data){
    data %>% 
      map(function(df){
        as.data.frame(df$adjustedData)
      })
  }) %>% 
  set_names("adjusteddf_residual_list_c18", 
            "adjusteddf_residual_list_hilic") %>% 
  list2env(.,envir = .GlobalEnv)

list(adjusteddf_residual_list_c18, adjusteddf_residual_list_hilic) %>% 
  map(function(datalist){
    datalist %>% 
      map(function(data){
        data %>% 
          bind_cols() %>% 
          t() %>% 
          as.data.frame()
      })
  }) %>% 
  set_names("combined_residual_list_c18", "combined_residual_list_hilic") %>% 
  list2env(.,envir = .GlobalEnv)

table(rownames(covar_c18_combined[["case"]])==colnames(combined_residual_list_c18[["case"]]))

save(combined_residual_list_c18, file = "combined_residual_list_c18_more_nopd.RData")
save(combined_residual_list_hilic, file = "combined_residual_list_hilic_more_nopd.RData")