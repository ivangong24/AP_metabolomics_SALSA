## ---------------------------
##
## Script name: 3-clean_data.R
## Purpose of script: To clean and preprocess raw data for analysis
##
## Author: Yufan Gong
##
## Date Created: 2025-11-12
##
## Date Modified: 2026-04-02
##
## Copyright (c) Yufan Gong, 2025
## Email: ivangong@ucla.edu
##
## ---------------------------

# Load required packages -----------------------------------------------------

source(here::here("scripts", "1-functions.R"))

# Load raw data --------------------------------------------------------------

source(here::here("scripts", "2-load_data.R"))

# clean the original salsa data ----------------------------------------------
{
  salsa_clean_total <- salsa_data_04212016 |> 
    dplyr::left_join(pa_nses_08042023, by = "rand_id") |>
    dplyr::left_join(salsa_ruca_07102015, by = "rand_id") |>
    dplyr::select(rand_id, bl_date, enrollment, birth_date, 
                  gender, blage, ageatcind, ageatdem, ageatdc, 
                  ses3, dem, cind, demcind, 
                  mh62, mh65, mh66, mh67, finalapoe,
                  pa3_met_if_ca, quinyostct, county, ruca_metro) |> 
    dplyr::rename(
      edu_year = ses3
    ) |> 
    dplyr::mutate(
      bl_date = if_else(is.na(bl_date), enrollment, bl_date),
      bl_date = lubridate::ymd(bl_date),
      alcohol_drinking = case_when(
        mh65 == 1 | mh66 == 1 | mh67 == 1 ~ 1,
        mh65 == 0 & mh66 == 0 & mh67 == 0 ~ 0,
        TRUE ~ NA
      ),
      nses = case_when(
        quinyostct %in% c(1, "1") ~ 1,
        quinyostct %in% c(2, "2") ~ 2,
        quinyostct %in% c(3, "3") ~ 3,
        quinyostct %in% c(4, "4", 5, "5") ~ 4
      ),
      apoe = case_when(
        finalapoe %in% c(1, 3, "1", "3") ~ "Non Apoe4 carrier",
        finalapoe %in% c(2, 4, 5, "2", "4", "5") ~ "Apoe4 carrier",
        TRUE ~ NA_character_
      )) |> 
    dplyr::select(-c(mh65, mh66, mh67)) |> 
    # in 1789 participants, 3 missed blage, 10 missed edu_year,
    # 14 missed smoking status, 113 missed physical activity, 
    # 14 missed alcohol drinking
    # impute missing data using mice, method = predictive mean matching (pmm)
    mice::mice(m = 5, maxit = 50, method = "pmm", seed = 42) |>
    mice::complete(1) |>
    # join the linked id
    dplyr::left_join(salsa_id, by = "rand_id") |>
    # restrict to people who have metabolomics data
    dplyr::filter(
      id %in% salsa_mapping_c18neg_surveylinked_03jun2025$id |
        id %in% salsa_mapping_hilicpos_surveylinked_03jun2025$id) |> 
    # in the 952 participants with metabolomics, 2 missed blage,
    # 6 missed edu_year, 8 missed smoking status, 41 missed physical activity,
    # 8 missed alcohol drinking
    # join the follow up dates back
    dplyr::left_join(salsa_data_04212016 |>
                       dplyr::select(rand_id, dyear, dcyear, dcst, demst,
                                     contains("mse3"),
                                     ends_with("date"), ends_with("diab"), 
                                     ends_with("bmi"), ends_with("age"), 
                                     ends_with("weight"), ends_with("height")
                                     ) |> 
                       dplyr::select(-matches("sa|bl_date|birth_date|blage")),
                     by = "rand_id") |> 
    dplyr::mutate(mse3_ultimate = coalesce(fv6_mse3_new, fv5_mse3_new, 
                                           fv4_mse3_new, fv3_mse3_new, 
                                           av1_mse3_new, bl_mse3_new)
    ) |>
    # transform weight from lbs to kg, 
    # height from inches to cm, and calculate BMI for each visit
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
    set_variable_labels(
      mh62 = "Baseline Smoking status",
      bl_date = "Baseline visit date",
      county = "County of residence",
      alcohol_drinking = "Baseline Alcohol drinking status",
      nses = "Neighborhood socioeconomic status",
      apoe = "Apoe4 carrier status",
      pa3_met_if_ca = "Baseline Physical activity status",
      ruca_metro = "Urban Residence",
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
      alcohol_drinking = c("Never drinker" = 0, "Ever Drinker" = 1),
      pa3_met_if_ca = c("Low" = 0, "High" = 1),
      nses = c("Low" = 1, "Low-Middle" = 2, "Middle" = 3, "Middle-High" = 4),
      ruca_metro = c("Rural" = 0, "Urban" = 1),
      demcind = c("No Dementia/CIND" = 0, "Dementia/CIND" = 1),
      cind = c("No CIND" = 0, "CIND" = 1)
    ) |>
    modify_if(is.labelled, to_factor)
  
  
  salsa_clean_cox <- salsa_clean_total |> 
    # removed people with demcind at baseline, n = 42
    dplyr::filter(!(demcind == "Dementia/CIND" & dcyear == 0)) |>
    # removed people without follow-up, n = 2
    dplyr::filter(!(dplyr::if_all(av1_date:fv6_date, is.na)))
    # # removed people without 3MSE data during follow-ups, n = 5
    # dplyr::filter(!(dplyr::if_all(av1_mse3_new:fv6_mse3_new, is.na)))
  
}
# we only have 471 packyr info out of 952 with metabolomics data, so we will not use packyr as a covariate in the main analysis, but we can do sensitivity analysis with it
skimr::skim(salsa_clean_total)
skimr::skim(salsa_clean_cox)
look_for(salsa_data_04212016, "drink")
look_for(pa_nses_08042023 |> 
           filter(rand_id %in% salsa_clean_total$rand_id), "pa3_met") # 41 missing physical activity data
table(salsa_clean_total$nses)

table(salsa_clean_total$alcohol_drinking) 
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
   
   # merge the blood_date to salsa_clean and make it a long format dataset for later use
   
   list(salsa_blood_date_c18, salsa_blood_date_hilic) |>
     purrr::map(function(data){
       list(salsa_clean_total, salsa_clean_cox) |> 
         purrr::map(function(df){
           df |>
             dplyr::left_join(data |>
                                dplyr::select(rand_id, file.name_new,
                                              blood_date, batch, wave), 
                              by = "rand_id") |>
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
               mse3_at_blooddraw = dplyr::if_else(
                 is.na(mse3_at_blooddraw), 
                 mse3_ultimate, mse3_at_blooddraw),
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
   
   # check if the blood date for c18 and hilic are the same 
   # in the long format datasets after merging back to salsa_clean
   
   list(salsa_clean_long_list_c18, salsa_clean_long_list_hilic) |> 
     purrr::pmap(
       function(df_c18, df_hilic){
         all(df_c18$blood_date == df_hilic$blood_date, na.rm = TRUE)
       }
     )
   
   # create a new list with total, total for cox regression, 
   # and stratified datasets (no demcind and demcind), 
   # and check the data distribution after cleaning
   
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
  

  # ---- parallel plan + progress handlers for average exposure calculation ----
  future::plan(future::multisession, 
               workers = max(1, future::availableCores() - 1))
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
          salsa_data |> dplyr::select(rand_id, id, blood_date, wave),
          by = "rand_id"
        ) |>
        dplyr::mutate(
          blood_year = lubridate::year(blood_date),
          year = lubridate::year(date)
        ) |>
        dplyr::filter(year >= blood_year - 5, year < blood_year) |>
        dplyr::group_by(rand_id, year, blood_date, wave, toxicant) |>
        dplyr::summarize(value = mean(value, na.rm = TRUE), .groups = "drop") |>
        dplyr::group_by(rand_id, blood_date, wave, toxicant) |>
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
        dplyr::rename_all(stringr::str_to_lower) |> 
        dplyr::left_join(
          salsa_data |> dplyr::select(rand_id, blood_date, batch, wave),
          by = c("rand_id", "blood_date", "wave")
        ) |> 
        dplyr::relocate(batch, .after = wave)
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
  
    
  air_toxicants_avg_iqr_list <- air_toxicants_avg_list |>
    purrr::map(function(dflist){
      dflist |> 
        purrr::map(function(data){
          data |> 
            dplyr::mutate(
              across(starts_with("exp_"),
                     ~ as.numeric(.x) / IQR(.x, na.rm = TRUE),
                     .names = "{.col}_iqr")) 
        })
    })

}

# make a heatmap for exposure correlations
pheatmap::pheatmap(
  cor(air_toxicants_avg_iqr_list[["total"]][["all"]] %>%
        select(ends_with("_iqr"), -matches("o3|zinc")),
      use = "complete.obs"),
  display_numbers = TRUE
)

# Define covariates for adjustment -------------------------------------------

# Covariates to adjust for in the metabolomics data
# These should NOT include the exposure of interest (NOx, air toxicants)
# It would be good to separate the analysis by demcind status

covar_list <- list(
  covar = quote_all(age_at_blooddraw, gender, edu_year, mh62, 
                    ruca_metro, nses,
                    wave, batch, demcind),
  
  covar_sen = quote_all(age_at_blooddraw, gender, edu_year, mh62, 
                        ruca_metro, nses,
                        alcohol_drinking, pa3_met_if_ca,  
                        bmi_at_blooddraw, diab_at_blooddraw, 
                        wave, batch, demcind)
  
)


# Prepare covariate matrices -------------------------------------------------

# if using the cox dataset, 
# also include dcst as the survival time variable

list(
  list(salsa_clean_long_list_c18, salsa_clean_long_list_hilic),
  list(salsa_mapping_c18neg_surveylinked_03jun2025,
       salsa_mapping_hilicpos_surveylinked_03jun2025),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |> 
  purrr::pmap(function(datalist, data1, data2){
    list(salsa_clean_new_list, datalist, names(salsa_clean_new_list)) |> 
      purrr::pmap(function(dflist, salsa_data, list_name){
        dflist |> 
          purrr::map(function(data){
            covar_list |> 
              purrr::map(function(covars){
                if (list_name == "cox") {
                  covars <- c(covars, quote_all(dcst))
                }
                data |> 
                  dplyr::left_join(salsa_data |> 
                                     dplyr::select(rand_id, id, 
                                                   blood_date, file.name_new), 
                                   by = c("rand_id", "id", "blood_date")) |>
                  dplyr::inner_join(
                    data1 |>
                      dplyr::select(-c(wave, batch)) |> 
                      dplyr::rename_all(str_to_lower),
                    by = c("id", "file.name_new")) |>
                  dplyr::select(file.name_new, rand_id, blood_date,
                                all_of(covars)) |>
                  dplyr::filter(file.name_new %in% colnames(data2)) |>
                  # dplyr::select(-rand_id) |>
                  tibble::column_to_rownames(var = "file.name_new")
              })
          })
      })
  }) %>% 
  set_names("covar_list_c18", "covar_list_hilic") |>
  list2env(.GlobalEnv)

# WQS MIXTURE MODELS FOR AIR TOXICANTS ---------------------------------

library(gWQS)

source("scripts/wqs_modified.R")

## Identify air toxicant exposure variable names
exp_vars <- air_toxicants_avg_list[["total"]][["all"]] |>
  dplyr::select(starts_with("exp_")) |>
  dplyr::select(-matches("o3|zinc")) |> 
  names()

# For a sensitivity analysis, 
# we can also run the model with only traffic-related air toxicants/metals
exp_vars_traffic <- exp_vars |> 
  purrr::discard(~ str_detect(.x, regex("benzene|butadiene|lead", 
                                        ignore_case = TRUE)))

exp_vars_metal <- exp_vars |> 
  purrr::keep(~ str_detect(.x, regex("lead|chromium|nickel", 
                                        ignore_case = TRUE)))

# create combined dataset for analysis
combined_data_list <- list(covar_list_c18, air_toxicants_avg_list) |>
  purrr::pmap(function(covar_datalist, exp_datalist){
    list(covar_datalist, exp_datalist) |>
      purrr::pmap(function(covar_df_list, exp_df){
        covar_df_list |> 
           purrr::map(function(covar_df){
             covar_df |>
               dplyr::left_join(exp_df |>
                                  dplyr::select(rand_id, blood_date, 
                                                all_of(exp_vars)),
                                by = c("rand_id", "blood_date"))
           })
      })
  })

skim(combined_data_list[["total"]][["all"]][["covar"]])

system.time({
  list(
    exp_vars,
    exp_vars_traffic,
    exp_vars_metal
  ) |> 
    purrr::map(function(exp_list){
      covar_list_new <- covar_list |> 
        purrr::map(function(covars){
          covars |> 
            purrr::discard(~ stringr::str_detect(.x, "demcind|wave|batch"))
        })
      
      list(covar_list_new, combined_data_list[["total"]][["all"]]) |> 
        purrr::pmap(function(covars, data){
          run_wqs(
            data = data,
            outcome = "demcind",
            mix_name = exp_list,
            covariates = covars,
            id_cols = c("rand_id", "blood_date"),
            q = 4,
            validation = 0.6,
            b = 200,
            b1_pos = TRUE,
            b_constr = FALSE,
            rh = 5,
            family = "binomial",
            seed = 42
          )
        })
    }) |>
    set_names("wqs_model_weight_all", "wqs_model_weight_traffic", 
              "wqs_model_weight_metal") |>
    list2env(.GlobalEnv)
})

# extract the OR and 95% CI for the WQS index from each model
list(wqs_model_weight_all, wqs_model_weight_traffic, wqs_model_weight_metal) |> 
  purrr::map(function(modellist){
    modellist |> 
      purrr::map(function(model){
        exp(summary(model$model)$coefficients[,c("Estimate", "2.5 %", "97.5 %")])
      })
  }) |> 
  set_names("all", "traffic", "metal")

# extract the weights for each component from each model
list(wqs_model_weight_all, wqs_model_weight_traffic, wqs_model_weight_metal) |> 
  purrr::map(function(modellist){
    modellist |> 
      purrr::map(function(model){
        model$model$final_weights
      })
  }) |> 
  set_names("all", "traffic", "metal")

# prepare a datalist with the WQS index for each sample from each model 

list(wqs_model_weight_all, wqs_model_weight_traffic, wqs_model_weight_metal) |> 
  purrr::map(function(modellist){
    modellist |> 
      purrr::map(function(model){
        tibble::as_tibble(model$wqs_df) |> 
          rename(comp_wqs = wqs)
      })
  }) |> 
  purrr::pmap(function(df_all, df_traffic, df_metal){
    list(df_all, df_traffic, df_metal) |> 
      purrr::reduce(left_join, by = c("rand_id", "blood_date")) |> 
      dplyr::rename(comp_wqs_all = comp_wqs.x,
             comp_wqs_traffic = comp_wqs.y,
             comp_wqs_metal = comp_wqs)
  }) -> wqs_df_list



# QGCOMP MIXTURE MODELS FOR AIR TOXICANTS ---------------------------------

library(qgcomp)

source("scripts/qgcomp_modified.R")

# vignette("qgcomp-basic-vignette", package="qgcomp")
# vignette("qgcomp-advanced-vignette", package="qgcomp")

# run the qgcomp models with the same covariates as the WQS models, 
# to get the composite exposure index for each sample
list(
  exp_vars,
  exp_vars_traffic,
  exp_vars_metal
) |> 
  purrr::map(function(exp_list){
    covar_list_new <- covar_list |> 
      purrr::map(function(covars){
        covars |> 
          purrr::discard(~ stringr::str_detect(.x, "demcind|wave|batch"))
      })
    
    list(covar_list_new, combined_data_list[["total"]][["all"]]) |> 
      purrr::pmap(function(covars, data){
        run_qgcomp_noboot_parallel(
          data = data,
          exposures_list = exp_list,
          outcomes_list = "demcind",
          covariates_list = covars,
          q = 4,
          seed = 42,
          workers = 8,                # set based on your machine
          id_cols = c("rand_id", "blood_date")
        )
      })
  }) |> 
  set_names("qgcomp_model_weight_all", "qgcomp_model_weight_traffic",
            "qgcomp_model_weight_metal") |>
  list2env(.GlobalEnv)

# run the qgcomp cox models
list(
  exp_vars,
  exp_vars_traffic,
  exp_vars_metal
) |> 
  purrr::map(function(exp_list){
    covar_list_new <- covar_list |> 
      purrr::map(function(covars){
        covars |> 
          purrr::discard(~ stringr::str_detect(.x, "demcind|dcst|wave|batch"))
      })
    
    list(covar_list_new, combined_data_list[["cox"]][["all"]]) |> 
      purrr::pmap(function(covars, data){
        data_new <- data |> 
           dplyr::mutate(demcind = case_when(
             demcind == "Dementia/CIND" ~ 1,
             demcind == "No Dementia/CIND" ~ 0,
             TRUE ~ NA_real_
           ))
        
        run_qgcomp_cox_noboot_parallel(
          data = data_new,
          exposures_list = exp_list,
          time_var = "dcst",
          event_var = "demcind",
          covariates_list = covars,
          q = 4,
          seed = 42,
          workers = 8,                # set based on your machine
          id_var = "rand_id",
          id_cols = c("rand_id", "blood_date")
        )
      })
  }) |> 
  set_names("qgcomp_cox_model_weight_all", "qgcomp_cox_model_weight_traffic",
            "qgcomp_cox_model_weight_metal") |>
  list2env(.GlobalEnv)

# refit the models with the composite exposure index to get the OR and 95% CI 
# for the composite index, adjusting for the same covariates
list(qgcomp_model_weight_all, qgcomp_model_weight_traffic, 
     qgcomp_model_weight_metal) |> 
  purrr::map(function(modellist){
    covar_list_new <- covar_list |> 
      purrr::map(function(covars){
        covars |> 
          purrr::discard(~ stringr::str_detect(.x, "demcind|wave|batch"))
      })
    
    list(covar_list_new, combined_data_list[["total"]][["all"]], modellist) |> 
      purrr::pmap(function(covars, data, model){
        combined_data <- data |> 
          dplyr::left_join(model$composites, 
                           by = c("rand_id", "blood_date")) |> 
          rename(comp_qgcomp = comp_demcind) |> 
          dplyr::mutate(demcind = case_when(
            demcind == "Dementia/CIND" ~ 1,
            demcind == "No Dementia/CIND" ~ 0,
            TRUE ~ NA_real_
          ))
        
        test_model <- glm(as.formula(paste("demcind", "~ comp_qgcomp + ", 
                                           paste(covars, collapse = " + "))),
                          data = combined_data,
                          family = "binomial")
        
        exp(cbind(OR=coef(test_model),confint(test_model)))
      })
  }) |> 
  set_names("all", "traffic", "metal")

# do the same for the cox models to get the HR and 95% CI for the composite index
list(qgcomp_cox_model_weight_all, qgcomp_cox_model_weight_traffic,
     qgcomp_cox_model_weight_metal) |>
  purrr::map(function(modellist){
    covar_list_new <- covar_list |> 
      purrr::map(function(covars){
        covars |> 
          purrr::discard(~ stringr::str_detect(.x, "demcind|dcst|wave|batch"))
      })
    
    list(covar_list_new, combined_data_list[["cox"]][["all"]], modellist) |> 
      purrr::pmap(function(covars, data, model){
        combined_data <- data |> 
          dplyr::left_join(model$composites, 
                           by = c("rand_id", "blood_date")) |> 
          dplyr::rename(comp_qgcomp = comp_demcind) |>
          dplyr::mutate(demcind = case_when(
            demcind == "Dementia/CIND" ~ 1,
            demcind == "No Dementia/CIND" ~ 0,
            TRUE ~ NA_real_
          ))
        
        test_model <- coxph(as.formula(paste("Surv(dcst, demcind) ~ comp_qgcomp + ", 
                                           paste(covars, collapse = " + "))),
                            id = rand_id,
                            data = combined_data)
        
        exp(cbind(HR=coef(test_model),confint(test_model)))
        })
  }) |> 
  set_names("all", "traffic", "metal")
     


# prepare a datalist with the QGCOMP composite index for each sample from each model

list(qgcomp_model_weight_all, qgcomp_model_weight_traffic,
     qgcomp_model_weight_metal) |> 
  purrr::map(function(modellist){
    modellist |> 
      purrr::map(function(model){
        model$composites
      })
  }) |> 
  purrr::pmap(function(df_all, df_traffic, df_metal){
    list(df_all, df_traffic, df_metal) |> 
      purrr::reduce(left_join, by = c("rand_id", "blood_date")) |> 
      dplyr::rename(comp_qgcomp_all = comp_demcind.x,
             comp_qgcomp_traffic = comp_demcind.y,
             comp_qgcomp_metal = comp_demcind)
  }) -> qgcomp_df_list


list(qgcomp_cox_model_weight_all, qgcomp_cox_model_weight_traffic,
     qgcomp_cox_model_weight_metal) |> 
  purrr::map(function(modellist){
    modellist |> 
      purrr::map(function(model){
        model$composites
      })
  }) |> 
  purrr::pmap(function(df_all, df_traffic, df_metal){
    list(df_all, df_traffic, df_metal) |> 
      purrr::reduce(left_join, by = c("rand_id", "blood_date")) |> 
      dplyr::rename(comp_qgcomp_cox_all = comp_demcind.x,
                    comp_qgcomp_cox_traffic = comp_demcind.y,
                    comp_qgcomp_cox_metal = comp_demcind)
  }) -> qgcomp_cox_df_list




# PCA MIXTURE MODELS FOR AIR TOXICANTS ---------------------------------
#
# Note: PCA finds the direction of maximum *variance* in the exposure data,
# whereas WQS/QGcomp find the direction most predictive of *dementia/CIND*.
# These are fundamentally different objectives, so moderate-to-low correlation
# between PCA-PC1 and the outcome-driven composites is expected unless the
# highest-variance exposure pattern happens to coincide with the most
# disease-relevant one.
#
# To make the PCA composite as comparable as possible to WQS/QGcomp, we apply
# PCA to the same quartile-scored (0-3) exposure matrix that those methods use.

pca_composite_list <- list(
  all     = exp_vars,
  traffic = exp_vars_traffic,
  metal   = exp_vars_metal
) |>
  purrr::map(function(exp_list) {
    # use the same combined dataset as WQS/QGcomp
    dat <- combined_data_list[["total"]][["all"]][["covar"]] |>
      dplyr::select(rand_id, blood_date, all_of(exp_list)) |>
      tidyr::drop_na()

    # quartile-score each pollutant to match WQS/QGcomp input (0, 1, 2, 3)
    dat_q <- dat |>
      dplyr::mutate(
        across(all_of(exp_list),
               ~ as.integer(cut(.x,
                                breaks = quantile(.x, probs = 0:4 / 4, na.rm = TRUE),
                                include.lowest = TRUE, labels = FALSE)) - 1L)
      )

    # run PCA on the quartile-scored matrix (center + scale within prcomp)
    pca_fit <- prcomp(dat_q[exp_list], center = TRUE, scale. = TRUE)

    # proportion of variance explained by PC1
    pve <- summary(pca_fit)$importance[2, 1]
    message("PCA on ", length(exp_list), " pollutants: PC1 explains ",
            round(pve * 100, 1), "% of variance")

    # PC1 loadings — flip sign so that the majority of loadings are positive
    # (ensures higher score = higher overall exposure, matching WQS convention)
    loadings_pc1 <- pca_fit$rotation[, 1]
    if (sum(loadings_pc1 < 0) > sum(loadings_pc1 > 0)) {
      pca_fit$x[, 1] <- -pca_fit$x[, 1]
      loadings_pc1   <- -loadings_pc1
      message("  -> PC1 sign flipped so majority of loadings are positive")
    }
    message("  PC1 loadings: ",
            paste(names(loadings_pc1), round(loadings_pc1, 3),
                  sep = "=", collapse = ", "))

    # participant-level PC1 score
    comp_tbl <- dat |>
      dplyr::select(rand_id, blood_date) |>
      dplyr::mutate(comp_pca = as.numeric(pca_fit$x[, 1]))

    list(pca_fit = pca_fit, composites = comp_tbl,
         loadings = loadings_pc1, pve = pve)
  })

# Compare PCA composites with WQS and QGcomp composites
message("\n--- Correlation between PCA-PC1 and WQS/QGcomp composites ---")
purrr::iwalk(pca_composite_list, function(pca_obj, group_name) {
  comp <- pca_obj$composites |>
    dplyr::left_join(wqs_df_list[["covar"]], by = c("rand_id", "blood_date")) |>
    dplyr::left_join(qgcomp_df_list[["covar"]], by = c("rand_id", "blood_date"))

  wqs_col    <- paste0("comp_wqs_", group_name)
  qgcomp_col <- paste0("comp_qgcomp_", group_name)

  if (wqs_col %in% names(comp)) {
    r_wqs <- cor(comp$comp_pca, comp[[wqs_col]], use = "complete.obs")
    message(group_name, ": PCA vs WQS r = ", round(r_wqs, 3))
  }
  if (qgcomp_col %in% names(comp)) {
    r_qg <- cor(comp$comp_pca, comp[[qgcomp_col]], use = "complete.obs")
    message(group_name, ": PCA vs QGcomp r = ", round(r_qg, 3))
  }
})

# --- Scatter plot matrix: PCA vs WQS vs QGcomp for each mixture grouping ---

scatter_plot_list <- purrr::imap(pca_composite_list, function(pca_obj, group_name) {

  wqs_col    <- paste0("comp_wqs_", group_name)
  qgcomp_col <- paste0("comp_qgcomp_", group_name)

  comp <- pca_obj$composites |>
    dplyr::left_join(wqs_df_list[["covar"]], by = c("rand_id", "blood_date")) |>
    dplyr::left_join(qgcomp_df_list[["covar"]], by = c("rand_id", "blood_date")) |>
    dplyr::select(PCA = comp_pca,
                  WQS = dplyr::all_of(wqs_col),
                  QGcomp = dplyr::all_of(qgcomp_col)) |>
    tidyr::drop_na()

  # all pairwise combinations
  pairs <- list(
    c("WQS",    "QGcomp"),
    c("PCA",    "WQS"),
    c("PCA",    "QGcomp")
  )

  panels <- purrr::map(pairs, function(p) {
    r <- cor(comp[[p[1]]], comp[[p[2]]])
    ggplot(comp, aes(x = .data[[p[1]]], y = .data[[p[2]]])) +
      geom_point(alpha = 0.25, size = 0.8, color = "steelblue") +
      geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 0.7) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
               label = paste0("r = ", round(r, 3)),
               size = 3.5, fontface = "bold") +
      labs(x = p[1], y = p[2]) +
      theme_bw(base_size = 10)
  })

  patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(
      title = paste0("Composite index correlations — ",
                     tools::toTitleCase(group_name), " mixture"),
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )
})

# display all three groupings stacked
composite_corr_plot <- patchwork::wrap_plots(scatter_plot_list, ncol = 1)
print(composite_corr_plot)

ggsave(here::here("figures", "composite_index_correlations.png"),
       composite_corr_plot,
       width = 10, height = 9, dpi = 300, bg = "white")

# Prepare PCA composite df_list in the same structure as wqs_df_list /
# qgcomp_df_list (list keyed by covariate set). PCA is unsupervised so the
# same scores apply to every covariate set — replicate across sets.
pca_df <- pca_composite_list |>
  purrr::map(~ .x$composites |>
               dplyr::select(rand_id, blood_date, comp_pca)) |>
  purrr::reduce(dplyr::left_join, by = c("rand_id", "blood_date")) |>
  dplyr::rename(comp_pca_all     = comp_pca.x,
                comp_pca_traffic = comp_pca.y,
                comp_pca_metal   = comp_pca)

# replicate across covariate sets to match the structure of wqs_df_list
covar_set_names <- names(wqs_df_list)
pca_df_list <- purrr::set_names(
  purrr::map(covar_set_names, ~ pca_df),
  covar_set_names
)


# Merge the composite exposure back ---------------------------------------

# merge the WQS, QGCOMP, and PCA composite indices back to the combined
# dataset for later use in MWAS
combined_data_list |>
  purrr::map(function(dflist){
    dflist |>
      purrr::map(function(datalist){
        list(datalist, wqs_df_list, qgcomp_df_list,
             qgcomp_cox_df_list, pca_df_list) |>
          purrr::pmap(function(data, wqs_df, qgcomp_df,
                               qgcomp_cox_df, pca_df){
            data |>
              dplyr::mutate(
                across(all_of(exp_vars),
                       ~ as.numeric(.x) / IQR(.x, na.rm = TRUE),
                       .names = "{.col}_iqr")) |>
              dplyr::left_join(wqs_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(qgcomp_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(qgcomp_cox_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(pca_df,
                               by = c("rand_id", "blood_date"))
          })
      })
  }) -> combined_data_list_new

# merge the WQS, QGCOMP, and PCA composite indices back to the air toxicant
# exposure dataset for future loading
air_toxicants_avg_list |>
  map(function(dflist){
    dflist |>
      map(function(data){
         list(wqs_df_list, qgcomp_df_list,
              qgcomp_cox_df_list, pca_df_list) |>
          purrr::pmap(function(wqs_df, qgcomp_df,
                               qgcomp_cox_df, pca_df){
            data |>
              dplyr::mutate(
                across(all_of(exp_vars),
                       ~ as.numeric(.x) / IQR(.x, na.rm = TRUE),
                       .names = "{.col}_iqr")) |>
              dplyr::left_join(wqs_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(qgcomp_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(qgcomp_cox_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(pca_df,
                               by = c("rand_id", "blood_date"))
          })
      })
  }) -> air_toxicants_avg_list_new

# fit logistic regression models with the composite indices 
# to get the OR and 95% CI
test_combined_data <- combined_data_list_new[["total"]][["all"]][["covar"]]

test_combined_data |> 
  dplyr::select(starts_with("comp"), ends_with("_iqr")) |> 
  names() |>
  purrr::map(function(exp){
    model <- glm(as.formula(paste("demcind", "~", exp, "+", 
                       paste(covar_list[["covar"]] |> 
                               discard(~str_detect(.x, "demcind|wave|batch")), 
                             collapse = " + "))),
        data = test_combined_data,
        family = "binomial")
    
    exp(cbind(OR=coef(model),confint(model)))
  })

# check the correlation between the composite indices and visualize with a heatmap
tbl_composite_cor <- air_toxicants_avg_list_new[["total"]][["all"]][["covar"]] %>%
  select(starts_with("comp")) %>%
  rename(`Air toxicant composite (QGCOMP all)` = comp_qgcomp_all,
         `Air toxicant composite (QGCOMP traffic-related)` = comp_qgcomp_traffic,
         `Air toxicant composite (QGCOMP metal)` = comp_qgcomp_metal,
         `Air toxicant composite (QGCOMP cox all)` = comp_qgcomp_cox_all,
         `Air toxicant composite (QGCOMP cox traffic-related)` = comp_qgcomp_cox_traffic,
         `Air toxicant composite (QGCOMP cox metal)` = comp_qgcomp_cox_metal,
         `Air toxicant composite (WQS all)` = comp_wqs_all,
         `Air toxicant composite (WQS traffic-related)` = comp_wqs_traffic,
         `Air toxicant composite (WQS metal)` = comp_wqs_metal,
         `Air toxicant composite (PCA all)` = comp_pca_all,
         `Air toxicant composite (PCA traffic-related)` = comp_pca_traffic,
         `Air toxicant composite (PCA metal)` = comp_pca_metal) %>%
  cor(use = "pairwise.complete.obs") %>% 
  as_tibble(rownames = 'var_x') |> 
  pivot_longer(
    -var_x,
    names_to = "var_y", 
    values_to = "correlation"
  ) 

tbl_composite_cor |> 
  mutate(correlation = if_else(var_x > var_y, correlation, NA)) |> 
  ggplot(aes(var_x, var_y)) +
  geom_tile(fill = 'white', col = 'grey80') +
  geom_point(
    aes(fill = correlation, size = abs(correlation)), 
    color = 'black',
    shape = 21
  ) +
  geom_text(
    data = tbl_composite_cor |> 
      mutate(correlation = if_else(var_y > var_x, correlation, NA)),
    aes(label = round(correlation, 2)),
    color = 'black',
    size = 5,
    fontface = 'bold'
  ) +
  theme_minimal(
    base_size = 16
  ) +
  labs(
    x = element_blank(),
    y = element_blank(),
    fill = 'Correlation'
  ) +
  scale_fill_gradient2_tableau(
    trans = "reverse"
  ) + 
  # scale_fill_gradient2(
  #   high = 'firebrick2',
  #   mid = 'white',
  #   low = 'dodgerblue4',
  #   limits = c(0.8, 1),
  #   midpoint = 0.9
  # ) +
  scale_size_area(
    limits = c(0, 1),
    max_size = 18
  ) +
  coord_cartesian(expand = FALSE) +
  theme(legend.position = 'top',
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15, vjust = 0.5, angle = 30)) +
  guides(
    fill = guide_colorbar(
      barwidth = unit(10, 'cm')
    ),
    size = guide_none()
  )

pheatmap::pheatmap(
  cor(air_toxicants_avg_list_new[[1]][[1]][[1]] %>%
        dplyr::select(ends_with("_iqr"), starts_with("comp_")),
      use = "complete.obs"),
  display_numbers = TRUE
)



# =============================================================================
# SECTION: SAVE CLEANED DATA
# =============================================================================

## Create output directory
dir.create(here::here("data", "processed"), 
           showWarnings = FALSE, recursive = TRUE)

## Save cleaned datasets for downstream analysis

save(salsa_clean_total, salsa_clean_cox, salsa_clean_new_list,
     salsa_clean_long_list_c18, salsa_clean_long_list_hilic,
     salsa_blood_date_c18, salsa_blood_date_hilic,
     file = here::here("data", "processed", "salsa_clean.RData"))

save(air_toxicants_avg_list, air_toxicants_avg_iqr_list,
     air_toxicants_avg_list_new,
     file = here::here("data", "processed", "air_toxicants_exposure.RData"))

save(combined_data_list_new, file = here::here("data", "processed", 
                                               "combined_data_list_new.RData"))

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
message("  - combined_data_list_new.RData")
message("  - covar_matrices.RData")
message("  - metabolomics_met_link.RData")
# message("  - nox_exposure.RData")
# message("  - salsa_final_analysis.RData")

#--------------------------------End of the code--------------------------------