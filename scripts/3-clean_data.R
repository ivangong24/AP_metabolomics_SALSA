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

# clean the original salsa data -------------------------------------------
{
  salsa_clean_total <- salsa_data_04212016 |> 
    # # removed people without baseline visit, n = 3
    # dplyr::filter(!is.na(bl_date)) |> 
    # # removed people with CIND at baseline, n = 115
    # dplyr::filter(
    #   !(demcind == 1 & dcyear == 0) 
    # # | (demcind == 1 & blage >= ageatcind) 
    # # | (demcind == 1 & blage >= ageatdem)
    # ) |> 
    # # removed no-follow-ups & survival time = 0, n = 57
    # dplyr::filter(!(dplyr::if_all(av1_date:fv6_date, is.na) & dcst == 0)) |> 
    # n = 1614 for now, need to further restrict to people who provided all necessary information (n= 53)
    # but not sure what variables are needed yet
    # join the wave variable to get the blood draw date
    # dplyr::left_join(salsa_mapping_c18neg_surveylinked_03jun2025 |> 
    #                    rename(
    #                      batch_c18 = batch,
    #                      wave_c18 = wave) |> 
    #                    select(id, batch_c18, wave_c18), by = "id") |>
    #  dplyr::left_join(salsa_mapping_hilicpos_surveylinked_03jun2025 |>
    #                    rename(
    #                      batch_hilic = batch,
    #                      wave_hilic = wave) |> 
    #                    select(id, batch_hilic, wave_hilic), by = "id")
    # this final number to this step may change after checking the variables needed

    dplyr::select(rand_id, bl_date, enrollment, birth_date, 
                  gender, ageatcind, ageatdem, ageatdc, 
                  ses3, dem, cind, demcind, mh62) |> 
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
    # join the linked id
    dplyr::left_join(salsa_id, by = "rand_id") |>
    # restrict to people who have metabolomics data
    dplyr::filter(id %in% salsa_mapping_c18neg_surveylinked_03jun2025$id |
                    id %in% salsa_mapping_hilicpos_surveylinked_03jun2025$id) |>
    # join the follow up dates back
    dplyr::left_join(salsa_data_04212016 |>
                       dplyr::select(rand_id, dyear, dcyear, dcst, demst,
                                     contains("mse3"),
                                     ends_with("date"), ends_with("diab"), 
                                     ends_with("bmi"), ends_with("age"), 
                                     ends_with("weight"), ends_with("height")
                                     ) |> 
                       dplyr::select(-matches("sa|bl_date|birth_date")),
                     by = "rand_id") |> 
    # transform weight from lbs to kg, height from inches to cm, and calculate BMI for each visit
    dplyr::mutate(
      across(ends_with("weight"), ~ .x * 0.45359237, .names = "{.col}_kg"),
      across(ends_with("height"), ~ .x * 2.54, .names = "{.col}_cm")) |>
    dplyr::mutate(
      bl_bmi_new = bl_weight_kg / (bl_height_cm/100)^2,
      av1_bmi_new = av1_weight_kg / (av1_height_cm/100)^2,
      # fv2_bmi_new = fv2_weight / (fv2_height/100)^2,
      fv3_bmi_new = fv3_weight_kg / (fv3_height_cm/100)^2,
      fv4_bmi_new = fv4_weight_kg / (fv4_height_cm/100)^2,
      fv5_bmi_new = fv5_weight_kg / (fv5_height_cm/100)^2,
      fv6_bmi_new = fv6_weight_kg / (fv6_height_cm/100)^2,
      bmi_ultimate = coalesce(fv6_weight_kg, fv5_weight_kg, fv4_weight_kg, 
                              fv3_weight_kg, av1_weight_kg, bl_weight_kg) /
        (coalesce(fv6_height_cm, fv5_height_cm, fv4_height_cm, 
                  fv3_height_cm, av1_height_cm, bl_height_cm)/100)^2
    ) |> 
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
      ageatdc = "Age at dementia/CIND",
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
  
  skimr::skim(salsa_clean_total)
  
  salsa_clean_cox <- salsa_clean_total |> 
    # removed people with demcind at baseline, n = 42
    dplyr::filter(!(demcind == "Dementia/CIND" & dcyear == 0)) |>
    # removed people without follow-up, n = 2
    dplyr::filter(!(dplyr::if_all(av1_date:fv6_date, is.na)))
    # # removed people without 3MSE data during follow-ups, n = 5
    # dplyr::filter(!(dplyr::if_all(av1_mse3_new:fv6_mse3_new, is.na)))
  
}
look_for(salsa_data_04212016, "mse")

skim(salsa_data_04212016$dcyear)
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
### get the 5-yr average exposure prior to blood draw for each air toxicant
{
  # link blood draw wave to visit date to determine the time window for calculating average exposure
  list(salsa_mapping_c18neg_surveylinked_03jun2025,
       salsa_mapping_hilicpos_surveylinked_03jun2025) |> 
    purrr::map(function(data){
      data |> 
        dplyr::left_join(salsa_clean_total |> 
                    select(rand_id, id, ends_with("date")), by = "id") |> 
        dplyr::mutate(blood_date = case_when(
          wave == 0 ~ bl_date,
          wave == 1 ~ av1_date,
          wave == 2 ~ fv2_date,
          wave == 3 ~ fv3_date,
          wave == 4 ~ fv4_date,
          wave == 5 ~ fv5_date,
          wave == 6 ~ fv6_date,
          TRUE ~ NA_Date_
        ))
    }) |> 
    purrr::set_names("salsa_blood_date_c18", "salsa_blood_date_hilic") |>
    list2env(.GlobalEnv)
  
  # check if the blood date for c18 and hilic are the same
  
   all(salsa_blood_date_c18$blood_date == salsa_blood_date_hilic$blood_date, 
       na.rm = TRUE)
   
   # join the blood_date to salsa_clean and make it a long format dataset for later use
   
   list(salsa_blood_date_c18, salsa_blood_date_hilic) |>
     purrr::map(function(data){
       list(salsa_clean_total, salsa_clean_cox) |> 
         purrr::map(function(df){
           df |>
             dplyr::left_join(data |>
                                dplyr::select(rand_id, file.name_new,
                                              blood_date), by = "rand_id") |>
             dplyr::mutate(
               final_diab_status = dplyr::coalesce(fv6_diab, fv5_diab, fv4_diab,
                                                   fv3_diab, fv2_diab, av1_diab, 
                                                   diab),
               bmi_at_blooddraw = dplyr::case_when(
                 blood_date == bl_date ~ bl_bmi,
                 blood_date == av1_date ~ dplyr::coalesce(av1_bmi, bl_bmi),
                 blood_date == fv2_date ~ dplyr::coalesce(av1_bmi, bl_bmi),
                 blood_date == fv3_date ~ dplyr::coalesce(fv3_bmi, av1_bmi, 
                                                          bl_bmi),
                 blood_date == fv4_date ~ dplyr::coalesce(fv4_bmi, fv3_bmi, 
                                                          av1_bmi, bl_bmi),
                 blood_date == fv5_date ~ dplyr::coalesce(fv5_bmi, fv4_bmi, 
                                                          fv3_bmi, av1_bmi,
                                                          bl_bmi),
                 blood_date == fv6_date ~ dplyr::coalesce(fv6_bmi, fv5_bmi, 
                                                          fv4_bmi, fv3_bmi,
                                                          av1_bmi, bl_bmi),
                 TRUE ~ NA
               ),
               bmi_at_blooddraw = dplyr::if_else(
                 is.na(bmi_at_blooddraw), 
                 bmi_ultimate, bmi_at_blooddraw),
               diab_at_blooddraw = dplyr::case_when(
                 blood_date == bl_date ~ diab,
                 blood_date == av1_date ~ av1_diab,
                 blood_date == fv2_date ~ fv2_diab,
                 blood_date == fv3_date ~ fv3_diab,
                 blood_date == fv4_date ~ fv4_diab,
                 blood_date == fv5_date ~ fv5_diab,
                 blood_date == fv6_date ~ fv6_diab,
                 TRUE ~ NA
               ),
               diab_at_blooddraw = dplyr::if_else(
                 is.na(diab_at_blooddraw), 
                 final_diab_status, diab_at_blooddraw),
               mse3_at_blooddraw = dplyr::case_when(
                 blood_date == bl_date ~ bl_mse3_new,
                 blood_date == av1_date ~ coalesce(av1_mse3_new, bl_mse3_new),
                 blood_date == fv2_date ~ coalesce(fv2_mse3_new, av1_mse3_new, 
                                                   bl_mse3_new),
                 blood_date == fv3_date ~ coalesce(fv3_mse3_new, fv2_mse3_new, 
                                                   av1_mse3_new, bl_mse3_new),
                 blood_date == fv4_date ~ coalesce(fv4_mse3_new, fv3_mse3_new, 
                                                   fv2_mse3_new, av1_mse3_new, 
                                                   bl_mse3_new),
                 blood_date == fv5_date ~ coalesce(fv5_mse3_new, fv4_mse3_new, 
                                                   fv3_mse3_new, fv2_mse3_new,
                                                   av1_mse3_new, bl_mse3_new),
                 blood_date == fv6_date ~ coalesce(fv6_mse3_new, fv5_mse3_new,
                                                   fv4_mse3_new, fv3_mse3_new, 
                                                   fv2_mse3_new, av1_mse3_new, 
                                                   bl_mse3_new),
                 TRUE ~ NA
               ),
               age_at_blooddraw = dplyr::case_when(
                 blood_date == bl_date ~ blage,
                 blood_date == av1_date ~ av1age,
                 blood_date == fv2_date ~ fv2age,
                 blood_date == fv3_date ~ fv3age,
                 blood_date == fv4_date ~ fv4age,
                 blood_date == fv5_date ~ fv5age,
                 blood_date == fv6_date ~ fv6age,
                 TRUE ~ NA
               ),
             ) |>
             dplyr::select(-starts_with("av"), -starts_with("fv"), 
                           -starts_with("bl_"), 
                           -c(diab, blage, bl_date, dyear, dcyear)) |> 
             dplyr::mutate(
               bmi_at_blooddraw = dplyr::if_else(
                 rand_id == "296898", median(bmi_at_blooddraw, na.rm = TRUE), 
                 bmi_at_blooddraw),
               diab_at_blooddraw = dplyr::if_else(
                 rand_id == "296898", 
                 median(diab_at_blooddraw, na.rm = TRUE), 
                 diab_at_blooddraw),
               bmi_ultimate = dplyr::if_else(
                 rand_id == "296898", 
                 median(bmi_ultimate, na.rm = TRUE), 
                 bmi_ultimate),
               final_diab_status = dplyr::if_else(
                 rand_id == "296898", 
                 median(final_diab_status, na.rm = TRUE), 
                 final_diab_status))
         })
       
     }) |>
     purrr::set_names("salsa_clean_long_list_c18", 
                      "salsa_clean_long_list_hilic") |>
     list2env(.GlobalEnv)

   all(salsa_clean_long_list_c18[[1]]$blood_date == 
                              salsa_clean_long_list_c18[[1]]$blood_date, 
       na.rm = TRUE)
   
   salsa_clean_new_list <- salsa_clean_long_list_c18 |> 
     purrr::map(function(data){
       total <- data |> 
         dplyr::select(-c(file.name_new)) |> 
         dplyr::relocate(blood_date, .after = rand_id) |>
         dplyr::relocate(lastdate, .after = enrollment) |>
         dplyr::distinct()
       
       stratified <- total |> 
         group_by(demcind) |> 
         group_split()
       
       combined <- c(list(total), stratified) |> 
         set_names("all", "no demcind", "demcind")
     }) |> 
     set_names("total", "cox")

   
   skimr::skim(salsa_clean_new_list[[1]][[1]])

  col_common <- quote_all(rand_id, date, .source_dir)
  
  exp_data_names_new <- c(exp_data_names, "nox")
  
  # system.time({
  #   salsa_clean_new_list |> 
  #     purrr::map(function(dflist){
  #       dflist |> 
  #         purrr::map(function(salsa_data){
  #           ls(envir = .GlobalEnv, all.names = TRUE) |> 
  #             purrr::keep(~ str_detect(.x, regex(str_c(exp_data_names_new, 
  #                                                      collapse = "|"), 
  #                                                ignore_case = TRUE))) |> 
  #             purrr::discard(~ str_detect(.x, regex(str_c("caline|final"),
  #                                                   ignore_case = TRUE))) |>
  #             purrr::keep(~ exists(.x, envir = .GlobalEnv, inherits = TRUE)) |> 
  #             mget(envir = .GlobalEnv, inherits = TRUE) |> 
  #             (\(x) x[order(names(x))])() |> 
  #             purrr::map(function(data){
  #               data |> 
  #                 dplyr::rename(value = all_of(setdiff(names(data), 
  #                                                      col_common)),
  #                               toxicant = .source_dir
  #                 ) |>
  #                 dplyr::mutate(value = extreme_remove_percentile_win(value)) # winsorize the extreme values
  #             }) |> 
  #             purrr::list_rbind() |> 
  #             dplyr::right_join(salsa_data |> 
  #                                 select(rand_id, id, blood_date), 
  #                               by = "rand_id") |>
  #             dplyr::mutate(blood_year = lubridate::year(blood_date),
  #                           year = lubridate::year(date)) |>
  #             dplyr::filter(year >= blood_year - 5 & year < blood_year) |>
  #             dplyr::group_by(rand_id, year, blood_date, toxicant) |> 
  #             dplyr::summarize(
  #               value = mean(value, na.rm = TRUE),
  #               .groups = "drop"
  #             ) |> 
  #             dplyr::group_by(rand_id, blood_date, toxicant) |> 
  #             dplyr::summarize(
  #               avg_exp = mean(value, na.rm = TRUE),
  #               .groups = "drop"
  #             ) |> 
  #             tidyr::pivot_wider(
  #               names_from = toxicant,
  #               values_from = avg_exp,
  #               names_prefix = "exp_"
  #             ) |> 
  #             dplyr::rename_all(str_to_lower)
  #         })
  #     }) -> air_toxicants_avg_list
  # })
  
  
  # ---- parallel plan + progress handlers ----
  future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")  # or "cli"
  # change the global option for future to allow larger objects 
  # (default is 500MB, set to 8GB here)
  options(future.globals.maxSize = 8000 * 1024^2)
  system.time({
    
    # 1) Precompute the exposure object names ONCE (instead of inside every loop)
    exp_pattern <- stringr::regex(stringr::str_c(exp_data_names_new, collapse = "|"),
                                  ignore_case = TRUE)
    
    drop_pattern <- stringr::regex("caline|final", ignore_case = TRUE)
    
    exp_obj_names <- ls(envir = .GlobalEnv, all.names = TRUE) |>
      purrr::keep(~ stringr::str_detect(.x, exp_pattern)) |>
      purrr::discard(~ stringr::str_detect(.x, drop_pattern)) |>
      purrr::keep(~ exists(.x, envir = .GlobalEnv, inherits = TRUE)) |>
      sort()
    
    # 2) Load those exposure objects ONCE (named list)
    exp_obj_list <- mget(exp_obj_names, envir = .GlobalEnv, inherits = TRUE)
    
    # 3) Pre-transform exposure datasets ONCE:
    #    - rename to (value, toxicant)
    #    - winsorize
    #    - bind into one long table
    exp_long <- exp_obj_list |>
      purrr::imap(function(data, nm) {
        data |>
          dplyr::rename(
            value = dplyr::all_of(setdiff(names(data), col_common)),
            toxicant = .source_dir
          ) |>
          dplyr::mutate(value = extreme_remove_percentile_win(value))
      }) |>
      purrr::list_rbind()
    
    # helper: process a single salsa_data using the precomputed exp_long
    process_one <- function(salsa_data) {
      exp_long |>
        dplyr::right_join(
          salsa_data |> dplyr::select(rand_id, id, blood_date),
          by = "rand_id"
        ) |>
        dplyr::mutate(
          blood_year = lubridate::year(blood_date),
          year = lubridate::year(date)
        ) |>
        dplyr::filter(year >= blood_year - 5, year < blood_year) |>
        dplyr::group_by(rand_id, year, blood_date, toxicant) |>
        dplyr::summarize(value = mean(value, na.rm = TRUE), .groups = "drop") |>
        dplyr::group_by(rand_id, blood_date, toxicant) |>
        dplyr::summarize(avg_exp = mean(value, na.rm = TRUE), 
                         .groups = "drop") |>
        tidyr::pivot_wider(
          names_from = toxicant,
          values_from = avg_exp,
          names_prefix = "exp_"
        ) |>
        # if NA after pivot_wider, it means no exposure data in the 5-year window, 
        # so we can impute it with 0 (assuming no exposure)
        dplyr::mutate(across(starts_with("exp_"), 
                             ~ if_else(is.na(.x), 0, .x))) |>
        dplyr::rename_all(stringr::str_to_lower)
    }
    
    # 4) Parallelize the inner map (over salsa_data) with progress
    air_toxicants_avg_list <-
      progressr::with_progress({
        # total tasks = sum of lengths of inner lists
        total <- sum(purrr::map_int(salsa_clean_new_list, length))
        p <- progressr::progressor(steps = total)
        
        purrr::map(salsa_clean_new_list, function(dflist) {
          furrr::future_map(
            dflist,
            ~ { p(); process_one(.x) },
            .options = furrr::furrr_options(seed = TRUE)
          )
        })
      })
    
  })
  
    
  air_toxicants_avg_ztrans_list <- air_toxicants_avg_list |>
    purrr::map(function(dflist){
      dflist |> 
        purrr::map(function(data){
          data |> 
            dplyr::mutate(across(starts_with("exp_"), 
                                 ~ scale(.x, center = TRUE) %>% as.vector()))
        })
    })

}

# make a heatmap for exposure correlations
pheatmap::pheatmap(
  cor(air_toxicants_avg_ztrans_list[[1]] %>%
        select(-c(rand_id, blood_date, exp_o3, exp_zinc)),
      use = "complete.obs"),
  display_numbers = TRUE
)

# Define covariates for adjustment -------------------------------------------

# Covariates to adjust for in the metabolomics data
# These should NOT include the exposure of interest (NOx, air toxicants)
# It would be good to separate the analysis by demcind status and also include
# demcind in the total population analysis

myvars_covar <- quote_all(age_at_blooddraw, gender, edu_year, mh62,
                          bmi_at_blooddraw, diab_at_blooddraw, demcind)


# Prepare covariate matrices -------------------------------------------------


list(
  list(salsa_clean_long_list_c18, salsa_clean_long_list_hilic),
  list(salsa_mapping_c18neg_surveylinked_03jun2025,
       salsa_mapping_hilicpos_surveylinked_03jun2025),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |> 
  purrr::pmap(function(datalist, data1, data2){
    list(salsa_clean_new_list, datalist) |> 
      purrr::pmap(function(dflist, salsa_data){
        dflist |> 
          purrr::map(function(data){
            data |> 
              dplyr::left_join(salsa_data |> 
                                 dplyr::select(rand_id, id, 
                                               blood_date, file.name_new), 
                        by = c("rand_id", "id", "blood_date")) |>
              dplyr::inner_join(
                data1 |>
                  dplyr::rename_all(str_to_lower),
                by = c("id", "file.name_new")) |>
                  dplyr::select(file.name_new, rand_id, blood_date,
                                all_of(myvars_covar)) |>
                  dplyr::filter(file.name_new %in% colnames(data2)) |>
                  # dplyr::select(-rand_id) |>
                  tibble::column_to_rownames(var = "file.name_new")
          })
      })
  }) %>% 
  set_names("covar_list_c18", "covar_list_hilic") |>
  list2env(.GlobalEnv)





# QGCOMP MIXTURE MODELS FOR AIR TOXICANTS ---------------------------------

source("scripts/gqcomp_modified.R")

## Identify air toxicant exposure variable names
exp_vars <- air_toxicants_avg_list[["total"]][["all"]] |>
  dplyr::select(starts_with("exp_")) |>
  names()

test_data <- covar_list_c18[["total"]][["all"]] |> 
  left_join(air_toxicants_avg_list[["total"]][["all"]], 
            by = c("rand_id", "blood_date"))

test_out <- run_qgcomp_boot_parallel(
  data = test_data,
  exposures_list = exp_vars,
  outcomes_list = "demcind",
  covariates_list = myvars_covar |> 
    discard(~str_detect(.x, "demcind")),
  q = 4, B = 200, seed = 42,
  workers = 6,                # set based on your machine
  progress_handler = "txtprogressbar"
)

test_out$results

skimr::skim(test_data)

test_model <- qgcomp::qgcomp.glm.boot(
  f      = as.formula(paste("demcind", "~", 
                            paste(c(exp_vars, myvars_covar |> 
                                      discard(~str_detect(.x, "demcind"))), 
                                  collapse = " + "))),
  data   = test_data,
  expnms = exp_vars,
  q      = 4,
  B      = 200,
  seed   = 42,   # IMPORTANT: distinct per outcome
  family = binomial(),
  rr = FALSE
)

summary(test_model)


test_out_weight <- run_qgcomp_noboot_parallel(
  data = test_data,
  exposures_list = exp_vars,
  outcomes_list = "demcind",
  covariates_list = myvars_covar |> 
    discard(~str_detect(.x, "demcind")),
  q = 4,
  seed = 42,
  workers = 6,                # set based on your machine
  id_cols = c("rand_id", "blood_date")
)

test_exp_df <- test_data |> 
  dplyr::left_join(test_out_weight$composites, by = c("rand_id", "blood_date"))


pheatmap::pheatmap(
  cor(test_exp_df %>%
        select(starts_with("exp_"), starts_with("comp_")),
      use = "complete.obs"),
  display_numbers = TRUE
)

### get the 1- 10-yr mean exposure prior to baseline for each air toxicant 
# {
#   lag_windows <- c(1, 3, 5, 10)
# 
#   air_toxicants_yearly <- lag_windows |> 
#     purrr::map(function(lag_yrs){
#       ls(envir = .GlobalEnv, all.names = TRUE) |> 
#     purrr::keep(~ str_detect(.x, regex(str_c(exp_data_names, collapse = "|"), 
#                                        ignore_case = TRUE))) |> 
#     purrr::discard(~ str_detect(.x, regex(str_c("no2|o3|pm2_5"),
#                                  ignore_case = TRUE))) |> 
#     purrr::keep(~ exists(.x, envir = .GlobalEnv, inherits = TRUE)) |> 
#     mget(envir = .GlobalEnv, inherits = TRUE) |> 
#     (\(x) x[order(names(x))])() |> 
#     purrr::map(function(data){
#       data |> 
#         dplyr::rename(value = all_of(setdiff(names(data), col_common)),
#                toxicant = .source_dir
#             ) |> 
#         dplyr::mutate(value = extreme_remove_percentile_win(value) # winsorize the extreme values
#       ) 
#     }) |> 
#     purrr::list_rbind() |> 
#     dplyr::mutate(year = lubridate::year(date)) |>
#     dplyr::select(-date) |> 
#     dplyr::relocate(year, .after = rand_id) |> 
#     dplyr::left_join(salsa_clean |> 
#                 dplyr::select(rand_id, bl_date, enrollment), by = "rand_id") |> 
#     dplyr::mutate(baseline_year = lubridate::year(bl_date),
#                   lag = lag_yrs) |>
#     dplyr::filter(
#       year <= baseline_year - lag
#       )    
#     }) |> 
#     purrr::list_rbind() |> 
#     dplyr::group_by(rand_id, toxicant, lag) |>
#     dplyr::summarise(
#       exp_mean = mean(value, na.rm = TRUE),
#       .groups = "drop"
#     ) |> 
#     # dplyr::mutate(var_name = stringr::str_c(toxicant, "_", lag, "y")) |>
#     # dplyr::select(rand_id, var_name, exp_mean) |>
#     tidyr::pivot_wider(
#       names_from  = c(toxicant, lag),
#       values_from = exp_mean
#     ) |> 
#     dplyr::rename_all(str_to_lower) |> 
#     dplyr::mutate(across(-rand_id,
#                          ~ scale(.x, center = TRUE) %>% as.vector())) 
# }

## get NOx IQR data

### first look at the cleaned caline nox data from Dr. Paul, i.e., salsa2_ap

# the nox_iqr variable is basically the rescaled version of nox (i.e., nox divided by IQR(nox)), with subtle difference
# I guess the we should use the IQR of nox in the total study population (n = 1789) to serve as the denominator for rescaling
# but what exactly is the nox variable here? 
# is it the average nox exposure during the enrollment year (https://pmc.ncbi.nlm.nih.gov/articles/PMC7591265/)? or some other period?
# we need to verify this with Dr. Paul

# test_nox <- salsa2_ap |> 
#   dplyr::select(rand_id, nox, nox_iqr) |> 
#   dplyr::mutate(
#     nox_iqr_check = nox / IQR(nox, na.rm = TRUE),
#     nox_iqr_check_origin = nox / 2.31 # pre-calculated IQR value from Dr. Paul's paper
#   )
# 
# test <- caline1789_nox_1998_2002 |> 
#   mutate(nox_avg = rowMeans(across(nox_1998:nox_2002), na.rm = TRUE))
# IQR(test$nox_avg, na.rm = TRUE)
# 
# ### check the 2002 caline nox data
# 
# caline_nox_2002 <- caline_2002 |> 
#   # merge with salsa_geocode_unique_ca to get rand_id
#   dplyr::left_join(salsa_geocode_unique_ca |> 
#     dplyr::select(rand_id, unique_id),
#     by = "unique_id"
#   ) |>
#   # filter(rand_id %in% salsa2_ap$rand_id) |> 
#   dplyr::group_by(rand_id) |>
#   dplyr::summarise(
#     nox_2002_avg = mean(nox, na.rm = TRUE),
#     .groups = "drop"
#   ) |> 
#   dplyr::mutate(
#     nox_iqr_check_2002 = nox_2002_avg / IQR(nox_2002_avg, na.rm = TRUE)
#   ) 
# 
# # obviously, the nox_iqr_check_2002 is different from the nox_iqr in salsa2_ap
# # and it doesn't make sense to use only 2002 data to calculate the IQR for the 
# # entire study period
# 
# ### calculate the yearly average nox exposure from caline data up to enrollment date
# caline_nox_long <- salsa2_ap |>
#   dplyr::select(rand_id, unique_id, starts_with("nox_19"),
#                 starts_with("nox_20")) |>
#   # make the dataset long format
#   tidyr::pivot_longer(
#     cols = starts_with("nox_"),
#     names_to = c("year", "month"),
#     names_pattern = "nox_(\\d{4})_(\\d{1,2})",
#     values_to = "value") |> 
#   dplyr::left_join(
#     salsa_data_04212016 |> 
#       dplyr::select(rand_id, enrollment),
#     by = "rand_id"
#   ) |> 
#   dplyr::mutate(
#     enroll_year = lubridate::year(enrollment),
#     enroll_month = lubridate::month(enrollment),
#     year = as.numeric(year),
#     month = as.numeric(month)
#   ) |> 
#   dplyr::group_by(rand_id, unique_id) |> 
#   dplyr::filter(
#     year == enroll_year,
#     month <= enroll_month
#   ) |>
#   dplyr::ungroup()
# 
# 
# 
# caline_nox_yearly <- caline_nox_long %>%
#   dplyr::group_by(rand_id, unique_id) %>%
#   dplyr::summarise(
#     yearly_avg = mean(value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# nox_iqr_value <- IQR(caline_nox_yearly$yearly_avg, na.rm = TRUE)
# nox_median_value <- median(caline_nox_yearly$yearly_avg, na.rm = TRUE)
# 
# caline_nox_iqr <- caline_nox_yearly %>%
#   dplyr::mutate(
#     monthly_std = yearly_avg/nox_iqr_value
#   ) |>
#   # group_by(rand_id, unique_id) %>%
#   # summarise(
#   #   monthly_avg_std = mean(monthly_std, na.rm = TRUE),
#   #   .groups = "drop"
#   # ) |>
#   dplyr::left_join(
#     salsa2_ap |>
#       dplyr::select(rand_id, unique_id, nox_iqr),
#         by = c("rand_id", "unique_id")
#   )



# =============================================================================
# SECTION: PREPARE FINAL NOx EXPOSURE DATA
# =============================================================================

## Create final NOx exposure dataset with IQR scaling


# nox_exposure_final <- salsa2_ap |>
#   dplyr::select(rand_id, nox, nox_iqr) |>
#   dplyr::rename(
#     exp_nox = nox,
#     exp_nox_iqr = nox_iqr
#   ) |>
#   dplyr::distinct(rand_id, .keep_all = TRUE)

# nox_exposure_final <- caline_nox_iqr |>
#   dplyr::select(rand_id, yearly_avg, monthly_std) |>
#   dplyr::rename(
#     exp_nox = yearly_avg,
#     exp_nox_iqr = monthly_std
#   ) |>
#   dplyr::distinct(rand_id, .keep_all = TRUE)


## Merge NOx with other air toxicants exposure data


# air_toxicants_all <- air_toxicants_avg_ztrans |>
#   dplyr::left_join(nox_exposure_final, by = "rand_id") |>
#   dplyr::mutate(
#     exp_nox_z = scale(exp_nox, center = TRUE)[,1]
#   )


# =============================================================================
# SECTION: CREATE FINAL ANALYSIS DATASET
# =============================================================================

## Merge SALSA covariates with all exposure data

# salsa_data_final <- salsa_clean |>
#   dplyr::mutate(rand_id = as.character(rand_id)) |>
#   dplyr::left_join(
#     air_toxicants_avg_ztrans |>
#       dplyr::mutate(rand_id = as.character(rand_id)),
#     by = "rand_id"
#   )
# 
# ## Summary of final dataset
# message("\nFinal analysis dataset summary:")
# message(paste0("  Total participants: ", nrow(salsa_data_final)))
# message(paste0("  Participants with air toxicants data: ",
#                sum(!is.na(salsa_data_final$exp_benzene))))
# message(paste0("  Participants with NOx data: ",
#                sum(!is.na(salsa_data_final$exp_nox))))




# =============================================================================
# SECTION: SAVE CLEANED DATA
# =============================================================================

## Create output directory
dir.create(here::here("data", "processed"), 
           showWarnings = FALSE, recursive = TRUE)

## Save cleaned datasets for downstream analysis

save(salsa_clean_total, salsa_clean_cox, salsa_clean_new_list,
     salsa_clean_long_c18, salsa_clean_long_hilic,
     file = here::here("data", "processed", "salsa_clean.RData"))

save(air_toxicants_avg_list, air_toxicants_avg_ztrans_list,
     file = here::here("data", "processed", "air_toxicants_exposure.RData"))

save(covar_list_c18, covar_list_hilic,
     file = here::here("data", "processed", "covar_matrices.RData"))

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
message("  - covar_matrices.RData")
message("  - metabolomics_met_link.RData")
# message("  - nox_exposure.RData")
# message("  - salsa_final_analysis.RData")

#--------------------------------End of the code--------------------------------