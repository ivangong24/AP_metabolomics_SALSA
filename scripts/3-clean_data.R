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
    # dplyr::filter()
    # check all variables contains "smoke", if any of them is not NA, then classify as "ever smoker"
    dplyr::mutate(smoking_status = case_when(
      mh62 == 1 ~ "Never smoker",
      mh62 == 2 ~ "Former smoker",
      mh62 == 3 ~ "Current smoker",
      TRUE ~ NA_character_
    )) |> 
    dplyr::rename(
      edu_year = ses3
    ) |> 
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
}


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


### check the 2002 caline nox data

caline_nox_2002 <- caline_2002 |> 
  # merge with salsa_geocode_unique_ca to get rand_id
  dplyr::left_join(salsa_geocode_unique_ca |> 
    dplyr::select(rand_id, unique_id),
    by = "unique_id"
  ) |>
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


caline_nox_long <- salsa2_ap |> 
  dplyr::select(rand_id, unique_id, starts_with("nox_19"), starts_with("nox_20")) |> 
  # make the dataset long format
  tidyr::pivot_longer(
    cols = matches("^nox_\\d{4}_\\d{1,2}$"),
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