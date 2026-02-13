## ---------------------------
##
## Script name: 6-annotation.R
## Purpose of script: To annotate metabolic features
##
## Author: Yufan Gong
##
## Date Created: 2026-01-29
##
## Date Modified: 2026-01-29
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: This script performs annotation using
##       1. Emory Inhouse library
##       2. HMDB
##
##       Dependencies: Run scripts 1-5 before this script.
## ---------------------------

# Load required packages -----------------------------------------------------

library(tidyverse)
library(writexl)



# Load required datasets --------------------------------------------------


# Load annotation files (m/z and retention time)
# These should contain: chemical_ID, mz, time (retention time)
annotation_c18 <- read_csv(
  here::here("data", "metabolomics", "c18neg", "xmsannotator_c18neg.csv"),
  show_col_types = FALSE
)

annotation_hilic <- read_csv(
  here::here("data", "metabolomics", "hilicpos", "xmsannotator_hilicpos.csv"),
  show_col_types = FALSE
)

inhouse_c18 <- read_delim(
  here::here("data", "metabolomics", "annotation", 
             "c18neg_inhouse_library.txt"),
  delim = "\t",
  show_col_types = FALSE) %>% 
  filter(!is.na(C18_name))

inhouse_hilic <- read_delim(
  here::here("data", "metabolomics", "annotation", 
             "hilicpos_inhouse_library.txt"),
  delim = "\t",
  show_col_types = FALSE) %>% 
  filter(!is.na(HILIC_name))

load(here::here("data", "metabolomics", "mz_links", 
                "c18_mz_rt_link.Rdata"))

load(here::here("data", "metabolomics", "mz_links", 
                "hil_mz_rt_link.Rdata"))

load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))

load(here::here("data", "processed",  "metabolomics_met_link.RData"))


# xMSannotator ------------------------------------------------------------
library(xMSannotator)

list(
  list(c18_mz_rt_link_mz_links, hil_mz_rt_link_mz_links),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |> 
  purrr::pmap(function(mz_rt_link, med_data){
    med_data |> 
      tibble::rownames_to_column("mz_rt") |> 
      dplyr::inner_join(mz_rt_link, by = "mz_rt") |>
      dplyr::relocate(mz, time, .after = mz_rt) |>
      dplyr::select(-mz_rt)
  }) |> 
  purrr::set_names("tbl_feature_c18", "tbl_feature_hilic") |>
  list2env(.GlobalEnv)

data(adduct_table)
data(adduct_weights)

# create directories for annotation from xMSannotator output
list("c18", "hilic") |> 
  purrr::map(function(mode){
    list("HMDB", "KEGG", "LipidMaps") |> 
      purrr::map(function(db){
        dir <- here::here("annotation", mode, db)
        if (!dir.exists(dir)) {
          dir.create(dir, recursive = TRUE)
          message(paste0("Directory created: ", dir))
        } else {
          message(paste0("Directory already exists: ", dir))
        }
      })
  }) |> 
  invisible()

# need to fix namespace for this package
# 1. replace "if (queryadductlist == "all" & mode == "pos") {"
# with "if (is.character(queryadductlist) && length(queryadductlist) == 1 &&
# queryadductlist == "all" && mode == "pos") {" also do it for "neg"
# 2. replace "get_peak_blocks_modulesvhclust" with "xMSannotator:::get_peak_blocks_modulesvhclust"

get("multilevelannotation", asNamespace("xMSannotator"))

num_nodes = 10
# takes ~ 4 hours to run (M3 pro chip)
list(
  list(tbl_feature_c18, tbl_feature_hilic),
  list("neg", "pos"),
  list("c18", "hilic"),
  list(
    c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H"),
    c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H",
      "M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O",
      "M+H-2H2O")
  ),
  list(c("M-H"), c("M+H"))
) |> 
  purrr::pmap(function(feature_tbl, mode, mode_name, adductlist, filter_adduct){
    list("HMDB", "KEGG", "LipidMaps") |> 
      purrr::map(function(db){
        xMSannotator::multilevelannotation(
          feature_tbl,
          max.mz.diff = 10, max.rt.diff = 37, 
          num_nodes = 10,
          # queryadductlist = c("M-H", "M-2H", "M-H2O-H"),
          queryadductlist = adductlist, 
          filter.by = filter_adduct,
          adduct_weights = adduct_weights,
          mode = mode,
          # DB to search
          db_name = db,
          # biofluid.location = "Blood",
          # other parameters
          num_sets = 300,
          # output directory
          outloc = here::here("annotation", mode_name, db)
        )
      })
  })


# Load stage 5 tables from xMSannotator output ----------------------------
annotation_names <- list.dirs(here::here("annotation"), recursive = FALSE) |>
  list.files(pattern = "\\.(csv)$", 
             full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Stage5")) |>
  (\(files) {
    files |>
      stringr::str_extract("[^/]+$") |>         # extract file name from full path
      stringr::str_remove("\\.[^.]+$") |>       # remove file extension
      stringr::str_to_lower()                   # convert to lower case
  })()


# check if there are overlap metabolite features within the inhouse library

list(
  list(c18_mz_rt_link, hil_mz_rt_link),
  list(inhouse_c18, inhouse_hilic)
) %>% 
  purrr::pmap(function(mz_rt_link, inhouse_df){
    list(mz_rt_link, inhouse_df) %>% 
      purrr::map(function(df){
        df %>% 
          dplyr::mutate(mz3 = round(mz, 3))
      }) %>% 
      purrr::reduce(inner_join, by = "mz3")
  }) %>% 
  set_names("inhouse_c18_overlap", "inhouse_hilic_overlap") %>%
  list2env(.GlobalEnv)

# create final annotation dataframes

list(
  list(annotation_c18, annotation_hilic),
  list(inhouse_c18, inhouse_hilic),
  list(met_c18, met_hilic)
) |>
  purrr::pmap(function(annot_df, inhouse_df, met_df){
    list(
      list(annot_df |> 
             dplyr::rename(confidence_level = Annotation.confidence.score) |>
             dplyr::filter(confidence_level >= 2), 
           inhouse_df),
      list("xmsannotator", "inhouse")
    ) |>
      purrr::pmap(function(df, ref){
        df |>
          dplyr::rename_all(str_to_lower) |>
          dplyr::rename(
            chemical_id = any_of(matches("hmdbid")),
            chemical_name = any_of(matches("name"))
          ) |>
          dplyr::mutate(met = str_c("mz_rt", round(mz, 4), 
                                    round(time, 4), sep = "_"),
                        reference = ref) |>
          dplyr::select(met, mz, time, chemical_id, chemical_name, reference)
      }) |>
      dplyr::bind_rows() |>
      dplyr::left_join(met_df, by = "met") |>
      dplyr::distinct()
  }) |>
  purrr::set_names("met_c18_all", "met_hilic_all") |>
  list2env(.GlobalEnv)

list(
  list(met_c18_all, met_hilic_all),
  list(inhouse_c18_overlap, inhouse_hilic_overlap)
) |> 
  purrr::pmap(function(data1, data2){
    data1 |> 
      dplyr::mutate(reference = if_else(met_new %in% data2$mz_rt, 
                                 "inhouse", reference)) |>
      dplyr::select(-met) |>
      dplyr::rename(met = met_new) |>
      dplyr::relocate(met) |>
      dplyr::filter(!is.na(met))
  }) %>%
  purrr::set_names("met_c18_final", "met_hilic_final") |>
  list2env(.GlobalEnv)

# Link to MWAS results

list(
  list(met_c18_final, met_hilic_final),
  list(significant_c18, significant_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(met_df, sig_mwas_df_list, mode){
    sig_mwas_df_list |>
      purrr::map(function(df){
        df |>
          dplyr::left_join(met_df |> 
                             dplyr::select(-c(mz, time)),
                           by = "met")
      })
  }) |>
  set_names("mwas_c18_annotated", "mwas_hilic_annotated") |>
  list2env(.GlobalEnv)


# Save annotated MWAS results -------------------------------------------------

list(
  list(mwas_c18_annotated, mwas_hilic_annotated),
  list("c18", "hilic")
) |>
  purrr::pmap(function(df_list, mode){
    df_list %>% 
      purrr::imap(function(df, exp_name){
        df %>%
          writexl::write_xlsx(
            here::here("tables", "mwas_results", 
                       glue::glue("mwas_{mode}_{exp_name}_sig_annotated.xlsx"))
          )
        message(paste0("MWAS results with annotation for ", 
                       exp_name, " ", mode, 
                       " saved successfully!"))
      })
  })


# Save R objects for downstream analysis -------------------------------------

save(mwas_c18_annotated, mwas_hilic_annotated, 
     file = here::here("data", "metabolomics", "processed",
                       "mwas_annotation.RData"))



#--------------------------------End of the code--------------------------------