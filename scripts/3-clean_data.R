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

# clean air toxicants exposure data

## get the average exposure for each air toxicant
{
  col_common <- quote_all(rand_id, date, .source_dir)

  air_toxicants <- ls(envir = .GlobalEnv, all.names = TRUE) |> 
    purrr::keep(~ str_detect(.x, regex(str_c(exp_data_names, collapse = "|"), ignore_case = TRUE))) |> 
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
        dplyr::mutate(value = extreme_remove_percentile_win(value)) # winsorize the extreme values
    }) |> 
    purrr::list_rbind() |> 
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
    mutate(across(starts_with("exp_"), ~ scale(.x, center = TRUE)))
}

# merge air toxicants with salsa data
