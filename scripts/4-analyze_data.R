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
  salsa_data_final <- salsa_data_04212016 |> 
    dplyr::left_join(air_toxicants, by = "rand_id")

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
        f <- as.formula(str_c(outcome, "~ .x +", age_var, "+ gender"))
        
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
    map(function(models){
      models |> 
        map(tidy) |> 
        map(filter, term == ".x") |> 
        bind_rows(.id = "term") |> 
        rename(toxicant = term) |>
        mutate(
          OR = exp(estimate),
          lowCL = exp(estimate - 1.96 * std.error),
          upCL = exp(estimate + 1.96 * std.error),
          across(where(is.numeric) & !p.value, ~ round(.x, digits = 3))
        ) |>
        arrange(p.value)
    }) |> 
    purrr::set_names(
      "logit_cind_results",
      "logit_demcind_results"
    )  

}

# cox-proportional regression models
{
  # Example: air toxicants and cind/demcind
  
}