## ---------------------------
##
## Script name: 4-analyze_data.R
## Purpose of script: To analyze the cleaned data and generate results
##
## Author: Yufan Gong
##
## Date Created: 2025-11-17
##
## Date Modified: 2025-11-17
##
## Copyright (c) Yufan Gong, 2025
## Email: ivangong@ucla.edu
##
## ---------------------------

# logistic regression models (cross-sectional)
{
  ## Example: air toxicants and cind/demcind

  ### merge air toxicants with salsa data
  salsa_data_final <- salsa_clean |> 
    dplyr::left_join(air_toxicants_avg_ztrans, by = "rand_id") |> 
    dplyr::mutate(rand_id = as.character(rand_id))

  ### fit logistic regression models
  list(
    quote_all(cind, demcind),
    quote_all(ageatcind, ageatdc)
  ) |> 
    purrr::pmap(function(outcome, age_var){
      salsa_data_final |> 
        dplyr::select(starts_with("exp_")) |> 
        purrr::map(~{
        # build a formula like: cind ~ .x + ageatcind + gender
        f <- as.formula(str_c(outcome, "~ .x +", age_var, "+ gender + edu_year + smoking_status"))
        
        glm(
          formula = f,
          family  = binomial(link = "logit"),
          data    = salsa_data_final
        )
      })
    }) |> 
    purrr::set_names(c("logit_cind", "logit_demcind")) |>
    list2env(.GlobalEnv)

  list(logit_cind, logit_demcind) |> 
    purrr::map(function(models){
      models |> 
        purrr::map(tidy) |> 
        purrr::map(filter, term == ".x") |> 
        dplyr::bind_rows(.id = "term") |> 
        dplyr::rename(toxicant = term) |>
        dplyr::mutate(
          OR = exp(estimate),
          lowCL = exp(estimate - 1.96 * std.error),
          upCL = exp(estimate + 1.96 * std.error),
          dplyr::across(where(is.numeric) & !p.value, ~ round(.x, digits = 3))
        ) |>
        dplyr::arrange(p.value)
    }) |> 
    purrr::set_names(
      "logit_cind_results",
      "logit_demcind_results"
    )  

}

# cox-proportional regression models

## Fixed exposure Cox model
{
  # Example: air toxicants and cind/demcind
  salsa_data_final_cox <- salsa_clean |> 
    dplyr::left_join(air_toxicants_yearly, by = "rand_id") |> 
    dplyr::mutate(timediff_cind = ageatcind - blage,
                  timediff_demcind = ageatdem - blage)


  # fit Cox model for demcind
  # library(survival)

  list(
    quote_all(cind, demcind),
    quote_all(ageatcind, ageatdem),
    quote_all(timediff_cind, timediff_demcind)
  ) |> 
    purrr::pmap(function(outcome, age_var, time_diff_var){
      salsa_data_final_cox |> 
        dplyr::filter(!!sym(time_diff_var) > 0) |>
        dplyr::select(benzene_1:zinc_10) |>
        purrr::map(~{
          # build a formula like: Surv(blage, ageatcind, cind) ~ exp + gender + edu_year + smoking_status
          f <- as.formula(str_c("Surv(blage, ", age_var, ", ", outcome, 
          ") ~ .x + gender + edu_year + smoking_status"))
          survival::coxph(
            formula = f,
            data    = salsa_data_final_cox |> 
              dplyr::filter(!!sym(time_diff_var) > 0) # only include participants with valid time difference
            )
        })
      }) |> 
    purrr::set_names(c("cox_cind", "cox_demcind")) |>
    list2env(.GlobalEnv)

  options(tibble.print_max = 50) 

  list(cox_cind, cox_demcind) |> 
    purrr::map(function(models){
      models |> 
        purrr::map(tidy) |> 
        purrr::map(filter, term == ".x") |> 
        dplyr::bind_rows(.id = "term") |> 
        dplyr::rename(toxicant = term) |>
        dplyr::mutate(
          HR = exp(estimate),
          lowCL = exp(estimate - 1.96 * std.error),
          upCL = exp(estimate + 1.96 * std.error),
          dplyr::across(where(is.numeric) & !p.value, ~ round(.x, digits = 3))
        )
    }) |>
    purrr::set_names(
      "cox_cind_results",
      "cox_demcind_results"
    )

}

# replication of previous findings from Dr. Paul's paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC7591265/

{
  # create table 1 for baseline characteristics
  # exposure: nox_iqr, outcome: demcind
  # Covariates: baseline age(blage), sex(gender), years of education(ses3),
  # primary language spoken(language), occupation for most of the participant’s life(ses7_occgrp)
  # smoking status(mh62), baseline cognitive function(bl_mse3_new)
  salsa_clean_nox <- salsa2_ap |> 
    dplyr::select(rand_id, nox_iqr, nox_cal_q3, nox, demcind, dcst, blage, gender, ses3,
    language, ses7, mh62, bl_mse3_new, smoke_cigarettes_30day, age_start_smoke) |> 
    dplyr::mutate(
      smoke_status = case_when(
        smoke_cigarettes_30day == 1 ~ 3,
        smoke_cigarettes_30day == 0 & !is.na(age_start_smoke) | mh62 == 2 ~ 2,
        smoke_cigarettes_30day == 0 & is.na(age_start_smoke) | mh62 == 1  ~ 1,
        TRUE ~ 1
      )
    ) |> 
    select(-c(smoke_cigarettes_30day, age_start_smoke, mh62)) |> 
    rename(
      edu_year = ses3,
      occupation = ses7
    ) |> 
    # left_join(salsa_geocode_unique_ca |> 
    #   select(rand_id, county) |> 
    #   distinct(),
    #   by = "rand_id"
    # ) |>
    set_variable_labels(
      rand_id = "Participant ID",
      nox_iqr = "NOx (IQR-scaled)",
      nox = "NOx (ppb)",
      demcind = "Dementia or CIND",
      dcst = "Time to demcind diagnosis, last visit, or death",
      blage = "Baseline Age (years)",
      gender = "Gender",
      edu_year = "Years of Education",
      language = "Primary Language Spoken",
      occupation = "Occupation Group",
      smoke_status = "Smoking Status",
      bl_mse3_new = "Baseline Cognitive Function (MSE Score)"
    ) |> 
    set_value_labels(
      demcind = c("No Demintia/CIND" = 0, "Demintia/CIND" = 1),
      gender = c("Male" = 1, "Female" = 2),
      language = c("English" = 2, "Spanish" = 1),
      smoke_status = c("Never" = 1, "Former" = 2, "Current" = 3)
    ) |> 
    modify_if(is.labelled, to_factor) |> 
    mutate(language = fct_rev(language),
           gender = fct_rev(gender))
  
  vars_to_keep <- quote_all(demcind, blage, gender, edu_year, language, bl_mse3_new, smoke_status, nox)
  
  salsa_clean_nox |> 
    select(all_of(vars_to_keep)) |> 
    tbl_summary(
      missing = "no",
      by = demcind,
      statistic = list(all_continuous() ~ "{mean} ({sd})",
                       all_categorical() ~ "{n} ({p}%)"),
      digits = list(all_continuous() ~ 1,
                    all_categorical() ~ c(0, 1))
    ) |> 
    add_overall() |> 
    modify_header(label = "**Characteristics**") |> 
    modify_caption("Table 1. Demographic characteristics of SALSA participants (N = {N})") |> 
    bold_labels() |> 
    table1()


  # fit Cox model for demcind with NOx_iqr
  cox_nox_model <- survival::coxph(
    formula = Surv(dcst, demcind) ~ nox_iqr + blage + gender + edu_year + language + bl_mse3_new,
    id      = rand_id,
    data    = salsa_clean_nox
  )
  
  summary(cox_nox_model)  
}

caline_2002_avg <- caline_2002 |> 
  # merge with salsa_geocode_unique_ca to get rand_id
  dplyr::left_join(salsa_geocode_unique_ca |> 
    dplyr::select(rand_id, unique_id),
    by = "unique_id"
  ) |>
  dplyr::group_by(rand_id) |>
  dplyr::summarise(
    nox_caline_2002 = mean(nox, na.rm = TRUE),
    .groups = "drop"
  )

vars_to_drop <- quote_all(iqr, new, cat)

nox_caline_new <- complete2_together_lur2 |> 
  dplyr::select(rand_id, starts_with("nox_caline")) |> 
  dplyr::select(-matches(str_c(vars_to_drop, collapse = "|"))) |>
  dplyr::left_join(caline_2002_avg, by = "rand_id") |> 
  dplyr::relocate(nox_caline_2002, .after = nox_caline_2001) |> 
  # calculate the yearly average nox from baseline to 2002
  dplyr::mutate(
    nox_avg = rowMeans(across(starts_with("nox_caline_")), na.rm = TRUE)
  ) |> 
  # dplyr::mutate(nox_avg_iqr = nox_avg / IQR(nox_avg, na.rm = TRUE)) |> 
  dplyr::filter(rand_id %in% salsa2_ap$rand_id) |> 
  dplyr::mutate(nox_avg_iqr = nox_caline_bl / 2.31)


IQR(salsa2_ap$nox)
IQR(complete2_together_lur2$nox_caline_bl)

IQR(nox_caline_new$nox_avg)
