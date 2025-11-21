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
    quote_all(ageatcind, ageatdem)
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