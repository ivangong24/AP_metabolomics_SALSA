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

salsa_clean <- salsa_data_04212016 |> 
  dplyr::select(rand_id, bl_date, enrollment, blage, birth_date, gender, ageatcind, ageatdem,
  ses3, cind, demcind, contains("smoke")) |> 
  # dplyr::filter()
  # check all variables contains "smoke", if any of them is not NA, then classify as "ever smoker"
  dplyr::mutate(smoking_status = if_else(
    rowSums(across(contains("smoke"), ~ !is.na(.x))) > 0,
    "ever smoker",
    "never smoker"
  )) |> 
  dplyr::rename(
    edu_year = ses3
  ) |> 
  dplyr::select(-contains("smoke")) |> 
  dplyr::mutate(
    bl_date = if_else(is.na(bl_date), enrollment, bl_date),
    bl_date = lubridate::ymd(bl_date)
  ) |>
  # impoute missing data using mice, method = predictive mean matching (pmm)
  mice::mice(m = 5, maxit = 50, method = "pmm", seed = 42) |> 
  mice::complete(1) |> 
  dplyr::mutate(timediff_cind = ageatcind - blage,
    timediff_demcind = ageatdem - blage
    # index_cind = 
  )

## clean air toxicants exposure data

### get the average exposure for each air toxicant
{
  col_common <- quote_all(rand_id, date, .source_dir)

  air_toxicants_avg_ztrans <- ls(envir = .GlobalEnv, all.names = TRUE) |> 
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
    dplyr::mutate(across(starts_with("exp_"), ~ scale(.x, center = TRUE)))

}

### get the lagged mean exposure for each air toxicant
{
  lag_windows <- c(1, 3, 5, 10)

  air_toxicants_yearly <- lag_windows |> 
    purrr::map(function(lag_yrs){
      ls(envir = .GlobalEnv, all.names = TRUE) |> 
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
    dplyr::mutate(across(-rand_id, ~ scale(.x, center = TRUE))) |> 
    dplyr::rename_all(str_to_lower)

}

