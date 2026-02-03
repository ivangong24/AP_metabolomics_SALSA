## ---------------------------
##
## Script name: 2-load_data.R
## Purpose of script: To load data for analysis
##
## Author: Yufan Gong
##
## Date Created: 2025-11-12
##
## Date Modified: 2025-11-24
##
## Copyright (c) Yufan Gong, 2025
## Email: ivangong@ucla.edu
##
## ---------------------------

# 1. Load the data --------------------------------------------------

{
  # set AP exposure data names
  exp_data_names <- list.dirs(here("data", "exposure", "AP_DATA"), 
                              recursive = TRUE) |> 
    purrr::discard(~ str_detect(basename(.x), 
                                regex("ap_data", ignore_case = TRUE))) |> 
    purrr::keep(~ length(list.files(.x,
      pattern = "\\.(csv|xlsx|txt|sas7bdat)$",
      ignore.case = TRUE)) > 0) |>
    basename() |>
    janitor::make_clean_names() |> 
    make.unique(sep = "_")
  
  # set up parallel processing
  future::plan(multisession)
  doFuture::registerDoFuture()

  # load AP exposure data in parallel
  system.time({
    list.dirs(here("data", "exposure", "AP_DATA"), recursive = TRUE) |> 
      purrr::discard(~ str_detect(basename(.x), 
        regex("ap_data", ignore_case = TRUE))) |> 
      purrr::keep(~ length(list.files(.x,
        pattern = "\\.(csv|xlsx|txt|sas7bdat)$",
        ignore.case = TRUE)) > 0) |>
      furrr::future_map(sdir_merge, .progress = TRUE) |>
      purrr::set_names(exp_data_names) |>
      list2env(.GlobalEnv)
  })

  # set up caline data names
  caline_data_names <- list.dirs(here("data", "caline"), recursive = TRUE) |> 
    list.files(pattern = "\\.(csv|xlsx|sas7bdat)$", 
      full.names = TRUE, recursive = FALSE) |> 
    (\(files) {
      files |>
        stringr::str_extract("[^/]+$") |>         # extract file name from full path
        stringr::str_remove("\\.[^.]+$") |>       # remove file extension
        stringr::str_to_lower()                   # convert to lower case
    })()

  # load caline data
  list.dirs(here("data", "caline"), recursive = TRUE) |> 
    list.files("\\.(csv|xlsx|sas7bdat)$", full.names = TRUE, recursive = F) |> 
    purrr::map(read_file) |> 
    purrr::map(~rename_all(.x, str_to_lower)) |> 
    purrr::set_names(caline_data_names) |> 
    list2env(.GlobalEnv)

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
    list.files("\\.(csv|xlsx|txt|sas7bdat)$",
               full.names = TRUE, recursive = T) |> 
    purrr::map(read_file) |> 
    purrr::map(~rename_all(.x, str_to_lower)) |> 
    purrr::set_names(salsa_data_names) |> 
    list2env(.GlobalEnv)

  # set link datanames
  
  link_data_names <- list.dirs(here("data", "links"), recursive = FALSE) |>
    list.files(pattern = "\\.(csv|txt|Rdata)$", 
      full.names = TRUE, recursive = TRUE) |>
    (\(files) {
      files |>
        stringr::str_extract("[^/]+$") |>         # extract file name from full path
        stringr::str_remove("\\.[^.]+$") |>       # remove file extension
        stringr::str_to_lower()                   # convert to lower case
    })()
  
  # load link data
  list.dirs(here("data", "links"),recursive = FALSE) |> 
    list.files("\\.(csv|txt|Rdata)$", full.names = TRUE, recursive = T) |> 
    purrr::map(read_file) |> 
    purrr::map(~rename_all(.x, str_to_lower)) |> 
    purrr::set_names(link_data_names) |> 
    list2env(.GlobalEnv)


  # set metabolomics data names
  metabolomics_data_names <- list.dirs(here("data", "metabolomics"), 
                                       recursive = FALSE) |> 
    # purrr::discard(~ str_detect(.x, regex("processed", ignore_case = TRUE))) |> 
    purrr::map(function(dir){
      dir |> 
        list.files(pattern = "\\.(csv|txt|Rdata)$", 
          full.names = TRUE, recursive = TRUE) |>
        (\(files) {
          files |>
          stringr::str_extract("[^/]+$") |>         # extract file name from full path
          stringr::str_remove("\\.[^.]+$") |>       # remove file extension
          stringr::str_to_lower()                   # convert to lower case
        })()
    })


  # load raw metabolomics data
  list(
    list.dirs(here::here("data", "metabolomics"),recursive = FALSE),
    metabolomics_data_names,
    list.dirs(here::here("data", "metabolomics"),recursive = FALSE) |> 
      basename()
  ) |> 
    purrr::pmap(function(dir, dataname, dirname){
      dir |> 
        list.files("\\.(csv|txt|Rdata)$", full.names = TRUE, recursive = T) |> 
        purrr::map(read_file) |>
        # purrr::map(as.data.frame) |> 
        purrr::map(~rename_all(.x, str_to_lower)) |>
        purrr::set_names(
          str_c(dataname, dirname, sep = "_")
        ) |> 
        list2env(.GlobalEnv)
    })
  
}

#--------------------------------End of the code--------------------------------