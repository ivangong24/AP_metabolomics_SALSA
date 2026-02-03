## ---------------------------
##
## Script name: 3-clean_data.R
## Purpose of script: To clean and preprocess raw data for analysis
##
## Author: Yufan Gong
##
## Date Created: 2025-11-12
##
## Date Modified: 2025-11-12
##
## Copyright (c) Yufan Gong, 2025
## Email: ivangong@ucla.edu
##
## ---------------------------

## clean the original salsa data
{
  salsa_clean <- salsa_data_04212016 |> 
    # removed people without baseline visit, n = 3
    dplyr::filter(!is.na(bl_date)) |> 
    # removed people with CIND at baseline, n = 115
    dplyr::filter(
      !(demcind == 1 & dcyear == 0) 
    # | (demcind == 1 & blage >= ageatcind) 
    # | (demcind == 1 & blage >= ageatdem)
    ) |> 
    # removed no-follow-ups & survival time = 0, n = 57
    dplyr::filter(!(dplyr::if_all(av1_date:fv6_date, is.na) & dcst == 0)) |> 
    # n = 1614 for now, need to further restrict to people who provided all necessary information (n= 53)
    # but not sure what variables are needed yet
    # this final number to this step may change after checking the variables needed
    ############################################################################################
    dplyr::select(rand_id, bl_date, enrollment, blage, birth_date, gender, ageatcind, ageatdem,
    ses3, cind, demcind, mh62) |> 
    # check all variables contains "smoke", if any of them is not NA, then classify as "ever smoker"
    dplyr::rename(
      edu_year = ses3
    ) |> 
    dplyr::mutate(
      bl_date = if_else(is.na(bl_date), enrollment, bl_date),
      bl_date = lubridate::ymd(bl_date)
    ) |>
    # impute missing data using mice, method = predictive mean matching (pmm)
    mice::mice(m = 5, maxit = 50, method = "pmm", seed = 42) |> 
    mice::complete(1) |> 
    # dplyr::mutate(timediff_cind = ageatcind - blage,
    #   timediff_demcind = ageatdem - blage
    #   # index_cind = 
    # ) |> 
    set_variable_labels(
      mh62 = "Smoking status",
      edu_year = "Years of education",
      gender = "Gender",
      blage = "Baseline age",
      ageatcind = "Age at CIND diagnosis",
      ageatdem = "Age at dementia diagnosis",
      cind = "Cognitive impairment no dementia (CIND)",
      demcind = "Dementia or CIND",
      enrollment = "Enrollment date"
    ) |>
    set_value_labels(
      gender = c("Male" = 1, "Female" = 2),
      mh62 = c("Never smoker" = 1, "Former smoker" = 2, "Current smoker" = 3),
      demcind = c("No Dementia/CIND" = 0, "Dementia/CIND" = 1),
      cind = c("No CIND" = 0, "CIND" = 1)
    ) |>
    modify_if(is.labelled, to_factor)
  
  skimr::skim(salsa_clean)
  
}


## clean air toxicants exposure data

### clean the caline nox data for 1998-2002

nox <- caline1789_nox_1998_2002 |> 
  select(-unique_id) |>
  tidyr::pivot_longer(
    cols = starts_with("nox_"),
    names_to = "date",
    names_prefix = "nox_",
    values_to = "nox"
  ) |>
  dplyr::mutate(
    .source_dir = "nox",
    date = lubridate::ymd(
      stringr::str_c(date, "-01-01")
    )
  )

mean(salsa_data_04212016$dcst, na.rm = TRUE)
### get the 5-yr average exposure prior to baseline for each air toxicant
{
  col_common <- quote_all(rand_id, date, .source_dir)
  
  exp_data_names_new <- c(exp_data_names, "nox")

  air_toxicants_avg_ztrans <- ls(envir = .GlobalEnv, all.names = TRUE) |> 
    purrr::keep(~ str_detect(.x, regex(str_c(exp_data_names_new, collapse = "|"), 
                                       ignore_case = TRUE))) |> 
    purrr::discard(~ str_detect(.x, regex(str_c("caline|final"),
                                 ignore_case = TRUE))) |>
    purrr::keep(~ exists(.x, envir = .GlobalEnv, inherits = TRUE)) |> 
    mget(envir = .GlobalEnv, inherits = TRUE) |> 
    (\(x) x[order(names(x))])() |> 
    purrr::map(function(data){
      data |> 
        dplyr::rename(value = all_of(setdiff(names(data), col_common)),
               toxicant = .source_dir
            ) |> 
        dplyr::mutate(value = extreme_remove_percentile_win(value)) # winsorize the extreme values
    }) |> 
    purrr::list_rbind() |> 
    dplyr::left_join(salsa_data_04212016 |> 
                select(rand_id, enrollment), by = "rand_id") |>
    dplyr::mutate(baseline_year = lubridate::year(enrollment),
                  year = lubridate::year(date)) |>
    dplyr::filter(year >= baseline_year - 5 | 
                    toxicant == "nox") |>
    dplyr::group_by(rand_id, toxicant) |> 
    dplyr::summarize(
      avg_exp = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    tidyr::pivot_wider(
      names_from = toxicant,
      values_from = avg_exp,
      names_prefix = "exp_"
    ) |> 
    dplyr::rename_all(str_to_lower) |> 
    dplyr::mutate(across(starts_with("exp_"), 
                         ~ scale(.x, center = TRUE) %>% as.vector()))

}

cor(air_toxicants_avg_ztrans %>% 
      select(-rand_id),
    use = "complete.obs")

### get the 1- 10-yr mean exposure prior to baseline for each air toxicant 
{
  lag_windows <- c(1, 3, 5, 10)

  air_toxicants_yearly <- lag_windows |> 
    purrr::map(function(lag_yrs){
      ls(envir = .GlobalEnv, all.names = TRUE) |> 
    purrr::keep(~ str_detect(.x, regex(str_c(exp_data_names, collapse = "|"), 
                                       ignore_case = TRUE))) |> 
    purrr::discard(~ str_detect(.x, regex(str_c("no2|o3|pm2_5"),
                                 ignore_case = TRUE))) |> 
    purrr::keep(~ exists(.x, envir = .GlobalEnv, inherits = TRUE)) |> 
    mget(envir = .GlobalEnv, inherits = TRUE) |> 
    (\(x) x[order(names(x))])() |> 
    purrr::map(function(data){
      data |> 
        dplyr::rename(value = all_of(setdiff(names(data), col_common)),
               toxicant = .source_dir
            ) |> 
        dplyr::mutate(value = extreme_remove_percentile_win(value) # winsorize the extreme values
      ) 
    }) |> 
    purrr::list_rbind() |> 
    dplyr::mutate(year = lubridate::year(date)) |>
    dplyr::select(-date) |> 
    dplyr::relocate(year, .after = rand_id) |> 
    dplyr::left_join(salsa_clean |> 
                dplyr::select(rand_id, bl_date, enrollment), by = "rand_id") |> 
    dplyr::mutate(baseline_year = lubridate::year(bl_date),
                  lag = lag_yrs) |>
    dplyr::filter(
      year <= baseline_year - lag
      )    
    }) |> 
    purrr::list_rbind() |> 
    dplyr::group_by(rand_id, toxicant, lag) |>
    dplyr::summarise(
      exp_mean = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    # dplyr::mutate(var_name = stringr::str_c(toxicant, "_", lag, "y")) |>
    # dplyr::select(rand_id, var_name, exp_mean) |>
    tidyr::pivot_wider(
      names_from  = c(toxicant, lag),
      values_from = exp_mean
    ) |> 
    dplyr::rename_all(str_to_lower) |> 
    dplyr::mutate(across(-rand_id,
                         ~ scale(.x, center = TRUE) %>% as.vector())) 
}

## get NOx IQR data

### first look at the cleaned caline nox data from Dr. Paul, i.e., salsa2_ap

# the nox_iqr variable is basically the rescaled version of nox (i.e., nox divided by IQR(nox)), with subtle difference
# I guess the we should use the IQR of nox in the total study population (n = 1789) to serve as the denominator for rescaling
# but what exactly is the nox variable here? 
# is it the average nox exposure during the enrollment year (https://pmc.ncbi.nlm.nih.gov/articles/PMC7591265/)? or some other period?
# we need to verify this with Dr. Paul

test_nox <- salsa2_ap |> 
  dplyr::select(rand_id, nox, nox_iqr) |> 
  dplyr::mutate(
    nox_iqr_check = nox / IQR(nox, na.rm = TRUE),
    nox_iqr_check_origin = nox / 2.31 # pre-calculated IQR value from Dr. Paul's paper
  )

test <- caline1789_nox_1998_2002 |> 
  mutate(nox_avg = rowMeans(across(nox_1998:nox_2002), na.rm = TRUE))
IQR(test$nox_avg, na.rm = TRUE)

### check the 2002 caline nox data

caline_nox_2002 <- caline_2002 |> 
  # merge with salsa_geocode_unique_ca to get rand_id
  dplyr::left_join(salsa_geocode_unique_ca |> 
    dplyr::select(rand_id, unique_id),
    by = "unique_id"
  ) |>
  # filter(rand_id %in% salsa2_ap$rand_id) |> 
  dplyr::group_by(rand_id) |>
  dplyr::summarise(
    nox_2002_avg = mean(nox, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  dplyr::mutate(
    nox_iqr_check_2002 = nox_2002_avg / IQR(nox_2002_avg, na.rm = TRUE)
  ) 

# obviously, the nox_iqr_check_2002 is different from the nox_iqr in salsa2_ap
# and it doesn't make sense to use only 2002 data to calculate the IQR for the entire study period

### calculate the yearly average nox exposure from caline data up to enrollment date
caline_nox_long <- salsa2_ap |>
  dplyr::select(rand_id, unique_id, starts_with("nox_19"),
                starts_with("nox_20")) |>
  # make the dataset long format
  tidyr::pivot_longer(
    cols = starts_with("nox_"),
    names_to = c("year", "month"),
    names_pattern = "nox_(\\d{4})_(\\d{1,2})",
    values_to = "value") |> 
  dplyr::left_join(
    salsa_data_04212016 |> 
      dplyr::select(rand_id, enrollment),
    by = "rand_id"
  ) |> 
  dplyr::mutate(
    enroll_year = lubridate::year(enrollment),
    enroll_month = lubridate::month(enrollment),
    year = as.numeric(year),
    month = as.numeric(month)
  ) |> 
  dplyr::group_by(rand_id, unique_id) |> 
  dplyr::filter(
    year == enroll_year,
    month <= enroll_month
  ) |>
  dplyr::ungroup()



caline_nox_yearly <- caline_nox_long %>%
  dplyr::group_by(rand_id, unique_id) %>%
  dplyr::summarise(
    yearly_avg = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

nox_iqr_value <- IQR(caline_nox_yearly$yearly_avg, na.rm = TRUE)
nox_median_value <- median(caline_nox_yearly$yearly_avg, na.rm = TRUE)

caline_nox_iqr <- caline_nox_yearly %>%
  dplyr::mutate(
    monthly_std = yearly_avg/nox_iqr_value
  ) |>
  # group_by(rand_id, unique_id) %>%
  # summarise(
  #   monthly_avg_std = mean(monthly_std, na.rm = TRUE),
  #   .groups = "drop"
  # ) |>
  dplyr::left_join(
    salsa2_ap |>
      dplyr::select(rand_id, unique_id, nox_iqr),
        by = c("rand_id", "unique_id")
  )

# Process metabolomic data ------------------------------------------------

## clean the link data

salsa_mapping_list <- ls(envir = .GlobalEnv, all.names = TRUE) |>
  keep(~ str_detect(.x, regex("salsa_mapping", ignore_case = TRUE))) |>
  discard(~ str_detect(.x, regex("list",ignore_case = TRUE))) |>
  keep(~ exists(.x, envir = .GlobalEnv, inherits = TRUE)) |>
  mget(envir = .GlobalEnv, inherits = TRUE) |>
  (\(x) x[order(names(x))])()

salsa_mapping_list |>
  purrr::map(function(data){
    data %>% 
      dplyr::mutate(file.name_new = str_c(file.name, "mzXML", sep = ".")) |> 
      dplyr::mutate(file.name_new = str_to_lower(file.name_new))
  }) %>% 
  purrr::set_names(names(salsa_mapping_list)) |>
  list2env(.,envir = .GlobalEnv)

list(raw_mzcalibrated_untargeted_mediansummarized_featuretable_c18neg,
     raw_mzcalibrated_untargeted_mediansummarized_featuretable_hilicpos) |>
  purrr::map(function(data){
    data %>%
      dplyr::mutate(met_new = str_c("mz_rt", round(mz, 4), 
                                 round(time, 4), sep = "_"),
                    met = str_c("mz_rt", round(mz, 4), 
                                     round(time, 1), sep = "_")) |>
      select(met, met_new)
  }) %>%
  purrr::set_names("met_c18", "met_hilic") %>%
  list2env(.,envir = .GlobalEnv)


# =============================================================================
# SECTION: PREPARE FINAL NOx EXPOSURE DATA
# =============================================================================

## Create final NOx exposure dataset with IQR scaling


nox_exposure_final <- salsa2_ap |>
  dplyr::select(rand_id, nox, nox_iqr) |>
  dplyr::rename(
    exp_nox = nox,
    exp_nox_iqr = nox_iqr
  ) |>
  dplyr::distinct(rand_id, .keep_all = TRUE)

# nox_exposure_final <- caline_nox_iqr |>
#   dplyr::select(rand_id, yearly_avg, monthly_std) |>
#   dplyr::rename(
#     exp_nox = yearly_avg,
#     exp_nox_iqr = monthly_std
#   ) |>
#   dplyr::distinct(rand_id, .keep_all = TRUE)


## Merge NOx with other air toxicants exposure data


air_toxicants_all <- air_toxicants_avg_ztrans |>
  dplyr::left_join(nox_exposure_final, by = "rand_id") |>
  dplyr::mutate(
    exp_nox_z = scale(exp_nox, center = TRUE)[,1]
  )


# =============================================================================
# SECTION: CREATE FINAL ANALYSIS DATASET
# =============================================================================

## Merge SALSA covariates with all exposure data

salsa_data_final <- salsa_clean |>
  dplyr::mutate(rand_id = as.character(rand_id)) |>
  dplyr::left_join(
    air_toxicants_all |>
      dplyr::mutate(rand_id = as.character(rand_id)),
    by = "rand_id"
  )

## Summary of final dataset
message("\nFinal analysis dataset summary:")
message(paste0("  Total participants: ", nrow(salsa_data_final)))
message(paste0("  Participants with air toxicants data: ",
               sum(!is.na(salsa_data_final$exp_benzene))))
message(paste0("  Participants with NOx data: ",
               sum(!is.na(salsa_data_final$exp_nox))))


## Create dataset with lagged exposures for Cox models

salsa_data_lagged <- salsa_clean |>
  dplyr::mutate(rand_id = as.character(rand_id)) |>
  dplyr::left_join(
    air_toxicants_yearly |>
      dplyr::mutate(rand_id = as.character(rand_id)),
    by = "rand_id"
  )


# =============================================================================
# SECTION: SAVE CLEANED DATA
# =============================================================================

## Create output directory
dir.create(here::here("data", "processed"), 
           showWarnings = FALSE, recursive = TRUE)

## Save cleaned datasets for downstream analysis

save(salsa_clean,
     file = here::here("data", "processed", "salsa_clean.RData"))

save(air_toxicants_avg_ztrans, 
     file = here::here("data", "processed", "air_toxicants_exposure.RData"))

save(met_c18, met_hilic,
     file = here::here("data", "processed", "metabolomics_met_link.RData"))

# save(nox_exposure_final, caline_nox_iqr,
#      file = here::here("data", "processed", "nox_exposure.RData"))
# 
# save(salsa_data_final, salsa_data_lagged,
#      file = here::here("data", "processed", "salsa_final_analysis.RData"))

message("\nCleaned data saved to data/processed/")
message("  - salsa_clean.RData")
message("  - air_toxicants_exposure.RData")
message("  - nox_exposure.RData")
message("  - salsa_final_analysis.RData")

#--------------------------------End of the code--------------------------------