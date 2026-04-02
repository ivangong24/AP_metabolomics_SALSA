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
  # extract the mode and database from the file path
  (\(files) {
    files |>
      #remove the directory and file extension
      str_extract("[^/]+/[^/]+/Stage5") |>
      #replace "/" to "_"
      str_replace_all("/", "_")
  })()
  

list.dirs(here::here("annotation"), recursive = FALSE) |>
  list.files(pattern = "\\.(csv)$", full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Stage5")) |>
  purrr::map(read_file) |> 
  purrr::map(~rename_all(.x, str_to_lower)) |> 
  purrr::set_names(annotation_names) |> 
  list2env(.GlobalEnv)


# Matched feature tables -------------------------------------------------------

# Function to calculate ppm mass error
ppm <- function(theo_mz, obs_mz) {
  ((obs_mz - theo_mz) / theo_mz) * 1e6
}

# Match features against inhouse library
matchInhouse <- function(id, mzr, rt, inhouse_tbl, name_col) {
  tbl <- inhouse_tbl |>
    dplyr::filter(abs(ppm(round(mz, 3), round(mzr, 3))) <= 10) |>
    dplyr::filter(abs(time - rt) <= 30)

  if (nrow(tbl) == 0) {
    tibble(id = id, mz = mzr, rt = rt,
           match_chemical = "", hmdbid = "", keggid = "")
  } else {
    tibble(id = id, mz = tbl$mz, rt = tbl$time,
           match_chemical = tbl[[name_col]],
           hmdbid = tbl$ID, keggid = tbl$KEGGID)
  }
}

# Match features against xMSannotator Stage 5 results
matchAnnotator <- function(id, mzr, rt, stage5_tbl) {
  tbl <- stage5_tbl |>
    dplyr::filter(abs(ppm(round(mz, 3), round(mzr, 3))) <= 10) |>
    dplyr::filter(abs(time - rt) <= 30) |>
    dplyr::rename(match_chemical = name, rt = time) |>
    dplyr::select(mz, rt, match_chemical, chemical_id, confidence, score)

  if (nrow(tbl) == 0) {
    tibble(id = id, mz = mzr, rt = rt,
           match_chemical = "", chemical_id = "",
           confidence = NA_real_, score = NA_real_)
  } else {
    tbl |> dplyr::mutate(id = id) |> dplyr::relocate(id)
  }
}

# Inhouse library matched tables
list(
  list(c18_mz_rt_link_mz_links, hil_mz_rt_link_mz_links),
  list(inhouse_c18, inhouse_hilic),
  list("C18_name", "HILIC_name")
) |>
  purrr::pmap(function(mz_link, inhouse_tbl, name_col) {
    lapply(seq_len(nrow(mz_link)), function(i) {
      matchInhouse(
        id = mz_link$mz_rt[i],
        mzr = mz_link$mz[i],
        rt = mz_link$time[i],
        inhouse_tbl = inhouse_tbl,
        name_col = name_col
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::filter(match_chemical != "")
  }) |>
  purrr::set_names("inhouse_c18_matched", "inhouse_hilic_matched") |>
  list2env(.GlobalEnv)

# xMSannotator Stage 5 matched tables
list(
  list(c18_mz_rt_link_mz_links, hil_mz_rt_link_mz_links),
  list("c18", "hilic")
) |>
  purrr::pmap(function(mz_link, mode) {
    list("HMDB", "KEGG", "LipidMaps") |>
      purrr::map(function(db) {
        stage5_tbl <- get(paste0(mode, "_", db, "_Stage5"))

        lapply(seq_len(nrow(mz_link)), function(i) {
          matchAnnotator(
            id = mz_link$mz_rt[i],
            mzr = mz_link$mz[i],
            rt = mz_link$time[i],
            stage5_tbl = stage5_tbl
          )
        }) |>
          dplyr::bind_rows() |>
          dplyr::filter(match_chemical != "")
      }) |>
      purrr::set_names(paste0(mode, "_", c("HMDB", "KEGG", "LipidMaps"), "_matched"))
  }) |>
  purrr::list_flatten() |>
  list2env(.GlobalEnv)

# Save full matched tables
save(
  inhouse_c18_matched, inhouse_hilic_matched,
  c18_HMDB_matched, c18_KEGG_matched, c18_LipidMaps_matched,
  hilic_HMDB_matched, hilic_KEGG_matched, hilic_LipidMaps_matched,
  file = here::here("annotation", "full_matched_tables.RData")
)


# Clean annotation -------------------------------------------------------------

load(here::here("data", "metabolomics", 
                "annotation", "full_matched_tables.RData"))

# Rank annotation: select best match per feature by score, then confidence,
# then database priority (HMDB > KEGG > LIPID MAPS)
rank_annotation <- function(tbl, feature_id) {
  tbl_subset <- tbl |> dplyr::filter(id == feature_id)
  max_score <- max(tbl_subset$score)
  if (sum(tbl_subset$score == max_score) == 1) {
    return(tbl_subset |> dplyr::filter(score == max_score))
  }

  max_confidence <- max(tbl_subset$confidence)
  tbl_subset <- tbl_subset |> dplyr::filter(confidence == max_confidence)

  if (any(tbl_subset$reference == "HMDB")) {
    return(tbl_subset |> dplyr::filter(reference == "HMDB"))
  }
  if (any(tbl_subset$reference == "KEGG")) {
    return(tbl_subset |> dplyr::filter(reference == "KEGG"))
  }
  return(tbl_subset)
}

# Clean and rank annotations for each mode (c18, hilic)
cleaned <- list(
  list("c18", "hilic"),
  list(inhouse_c18_matched, inhouse_hilic_matched),
  list(c18_HMDB_matched, hilic_HMDB_matched),
  list(c18_KEGG_matched, hilic_KEGG_matched),
  list(c18_LipidMaps_matched, hilic_LipidMaps_matched)
) |>
  purrr::pmap(function(mode, inhouse_matched, hmdb_matched,
                        kegg_matched, lipidmaps_matched) {
    # Filter xMSannotator matches by confidence >= 2
    hmdb_filtered <- hmdb_matched |> dplyr::filter(confidence >= 2)
    kegg_filtered <- kegg_matched |> dplyr::filter(confidence >= 2)
    lipidmaps_filtered <- lipidmaps_matched |> dplyr::filter(confidence >= 2)

    # Inhouse wide table: collapse matches per feature
    inhouse_wide <- inhouse_matched |>
      dplyr::select(id, match_chemical, hmdbid) |>
      dplyr::group_by(id) |>
      dplyr::reframe(
        compound = paste(match_chemical, collapse = "; "),
        chemical_id = paste(hmdbid, collapse = "; "),
        multiple_match = (dplyr::n() > 1)
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(reference = "In House Library",
                    confidence = NA)

    # Combine xMSannotator matches with database labels
    xms_long <- dplyr::bind_rows(
      hmdb_filtered |> dplyr::mutate(reference = "HMDB"),
      kegg_filtered |> dplyr::mutate(reference = "KEGG"),
      lipidmaps_filtered |> dplyr::mutate(reference = "LIPID MAPS")
    )

    # Rank and deduplicate annotations
    xms_long <- lapply(unique(xms_long$id), function(x) {
      rank_annotation(xms_long, x)
    }) |>
      dplyr::bind_rows() |>
      dplyr::distinct()

    # xMSannotator wide table: collapse per feature
    xms_wide <- xms_long |>
      dplyr::select(id, match_chemical, chemical_id, 
                    reference, confidence) |>
      dplyr::group_by(id) |>
      dplyr::reframe(
        compound = paste(match_chemical, collapse = "; "),
        chemical_id = paste(chemical_id, collapse = "; "),
        multiple_match = (dplyr::n() > 1),
        reference = unique(reference),
        confidence = unique(confidence)
      )

    # Combined wide: inhouse takes priority over xMSannotator
    wide <- dplyr::bind_rows(
      inhouse_wide,
      xms_wide |> dplyr::filter(!id %in% inhouse_wide$id)
    )

    list(xms_long = xms_long, wide = wide)
  }) |>
  purrr::set_names("c18", "hilic")

# Extract cleaned results
xms_c18_long <- cleaned$c18$xms_long
xms_hilic_long <- cleaned$hilic$xms_long
annotation_c18_wide <- cleaned$c18$wide
annotation_hilic_wide <- cleaned$hilic$wide

# Save cleaned annotation tables
save(
  inhouse_c18_matched, inhouse_hilic_matched,
  xms_c18_long, xms_hilic_long,
  file = here::here("data", "metabolomics", 
                    "annotation", "annotation_cleaned_long.RData")
)

save(
  annotation_c18_wide, annotation_hilic_wide,
  file = here::here("data", "metabolomics", 
                    "annotation", "annotation_cleaned_wide.RData")
)


# Link to MWAS results ---------------------------------------------------------

load(here::here("data", "metabolomics", 
                "annotation", "annotation_cleaned_long.RData"))

load(here::here("data", "metabolomics", 
                "annotation", "annotation_cleaned_wide.RData"))

list(
  list(annotation_c18_wide, annotation_hilic_wide),
  list(significant_results_list_c18, significant_results_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(annot_wide, sig_mwas_df_list, mode) {
    sig_mwas_df_list |>
      purrr::map(function(datalist) {
        datalist |> 
          purrr::map(function(dflist) {
            dflist |> 
              purrr::map(function(dfls){
                dfls |>
                  purrr::map(function(df){
                    df |> 
                      dplyr::left_join(annot_wide, by = c("met" = "id"))
                  })
              })
          })
      })
  }) |>
  purrr::set_names("mwas_annotated_list_c18", "mwas_annotated_list_hilic") |>
  list2env(.GlobalEnv)







# Save annotated MWAS results -------------------------------------------------

list(
  list(mwas_annotated_list_c18, mwas_annotated_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(datalist, mode) {
    datalist |> 
      purrr::imap(function(data_list, study){
        data_list |> 
          purrr::imap(function(dflist, population){
            dflist |>
              purrr::imap(function(dfls, covar_name){
                dfls |> 
                  purrr::imap(function(df, exp_name) {
                    df |>
                      writexl::write_xlsx(
                        here::here(
                          "tables", "mwas_results", study, 
                          population, covar_name,
                          glue::glue("mwas_{mode}_{exp_name}_{study}_{population}_{covar_name}_sig_annotated.xlsx"))
                      )
                    message(paste0("MWAS results with annotation for ",
                                   exp_name, " ", mode, " in ", 
                                   study, "_", population, 
                                   "with covariate ", covar_name,
                                   " saved successfully!"))
                  })
              })
          })
      })
  })


# Save R objects for downstream analysis -------------------------------------

save(mwas_annotated_list_c18, mwas_annotated_list_hilic,
     file = here::here("data", "metabolomics", "results",
                       "mwas_annotation.RData"))



#--------------------------------End of the code--------------------------------