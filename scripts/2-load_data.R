## ---------------------------
##
## Script name: 2-load_data.R
## Purpose of script: To load data for analysis
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

# 1. Load the data --------------------------------------------------

{
  # set exposure data names
  exp_data_names <- list.dirs(here("data", "exposure", "AP_DATA"), recursive = TRUE) |> 
    purrr::discard(~ str_detect(basename(.x), regex("ap_data", ignore_case = TRUE))) |> 
    purrr::keep(~ length(list.files(.x,
      pattern = "\\.(csv|xlsx|txt|sas7bdat)$",
      ignore.case = TRUE)) > 0) |>
    basename() |>
    janitor::make_clean_names() |> 
    make.unique(sep = "_")
  
  # set up parallel processing
  future::plan(multisession)
  doFuture::registerDoFuture()

  # load exposure data in parallel
  system.time({
    list.dirs(here("data", "exposure", "AP_DATA"), recursive = TRUE) |> 
    purrr::discard(~ str_detect(basename(.x), regex("ap_data", ignore_case = TRUE))) |> 
    purrr::keep(~ length(list.files(.x,
      pattern = "\\.(csv|xlsx|txt|sas7bdat)$",
      ignore.case = TRUE)) > 0) |>
    furrr::future_map(sdir_merge, .progress = TRUE) |>
    purrr::set_names(exp_data_names) |>
    list2env(.GlobalEnv)
  })

  # set salsa data names
  salsa_data_names <- list.dirs(here("data", "salsa"), recursive = FALSE) |>
    list.files(pattern = "\\.(csv|xlsx|txt|sas7bdat)$", 
      full.names = TRUE, recursive = TRUE) |>
    (\(files) {
      files |>
        stringr::str_extract("[^/]+$") |>         # extract file name from full path
        stringr::str_remove("\\.[^.]+$") |>       # remove file extension
        stringr::str_to_lower()                   # convert to lower case
    })()

  # load salsa data
  list.dirs(here("data", "salsa"),recursive = FALSE) |> 
    list.files("\\.(csv|xlsx|txt|sas7bdat)$", full.names = TRUE, recursive = T) |> 
    purrr::map(read_file) |> 
    purrr::map(~rename_all(.x, str_to_lower)) |> 
    purrr::set_names(salsa_data_names) |> 
    list2env(.GlobalEnv)
}

#--------------------------------End of the code--------------------------------