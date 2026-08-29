## ---------------------------
##
## Script name: R5-annotation_revision.R
## Purpose of script: To attach compound annotations to the R1 revision MWAS
##                    results
##
## Author: Yufan Gong
##
## Date Created: 2026-08-27
##
## Date Modified: 2026-08-27
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: This is the revision analogue of scripts/5-annotation.R.
##
##        The annotation itself does NOT depend on the exposure. Only the
##        mixture index changed between the submitted analysis and the
##        revision, so the same feature tables, the same Emory in-house
##        library and the same xMSannotator Stage 5 output apply. Re-running
##        multilevelannotation() would take ~4 hours and reproduce identical
##        tables, so this script reuses the cleaned annotation from
##        data/metabolomics/annotation/annotation_cleaned_wide.RData and joins
##        it to the revision MWAS results.
##
##        Re-run scripts/5-annotation.R (Sections "xMSannotator" through
##        "Clean annotation") only if the feature tables themselves change.
##
##        Outputs -> revision_output/tables/mwas_results/... (_sig_annotated)
##                   revision_output/data/metabolomics/results/
##        Downstream: R7-visualization_revision.R
##
## ---------------------------

# Load required packages -----------------------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R5-annotation_revision.R")


# Load required datasets -----------------------------------------------------

load(here::here("data", "metabolomics",
                "annotation", "annotation_cleaned_long.RData"))

load(here::here("data", "metabolomics",
                "annotation", "annotation_cleaned_wide.RData"))

load(rev_here("data", "metabolomics", "results",
              "mwas_results_all_revision.RData"))

message("Annotated features available: c18 = ", nrow(annotation_c18_wide),
        ", hilic = ", nrow(annotation_hilic_wide))


# Link to MWAS results -------------------------------------------------------

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


# Annotate the full (unfiltered) result tables as well ------------------------

## The submitted pipeline only annotated the _sig tables. Reviewer 1 major
## comment 8 asks for identification confidence to be traceable, so the full
## tables carry the annotation columns here too.
list(
  list(annotation_c18_wide, annotation_hilic_wide),
  list(combined_results_list_c18, combined_results_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(annot_wide, mwas_df_list, mode) {
    mwas_df_list |>
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
  purrr::set_names("mwas_full_annotated_list_c18",
                   "mwas_full_annotated_list_hilic") |>
  list2env(.GlobalEnv)


# Save annotated MWAS results ------------------------------------------------

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
                rev_dir("tables", "mwas_results", study, population, covar_name)
                dfls |>
                  purrr::imap(function(df, exp_name) {
                    df |>
                      writexl::write_xlsx(
                        rev_here(
                          "tables", "mwas_results", study,
                          population, covar_name,
                          glue::glue("mwas_{mode}_{exp_name}_{study}_{population}_{covar_name}_sig_annotated.xlsx"))
                      )
                    message(paste0("MWAS results with annotation for ",
                                   exp_name, " ", mode, " in ",
                                   study, "_", population,
                                   " with covariate ", covar_name,
                                   " saved successfully!"))
                  })
              })
          })
      })
  })


# Summarize annotation coverage of the significant features -------------------

annotation_coverage <- list(
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
                    tibble::tibble(
                      platform    = mode,
                      study       = study,
                      population  = population,
                      covar_set   = covar_name,
                      exposure    = exp_name,
                      label       = unname(rev_label(exp_name)),
                      n_features  = nrow(df),
                      n_fdr05     = sum(df$adj.P.Val < 0.05, na.rm = TRUE),
                      n_vip_gt2   = sum(df$VIP_comp1 > 2, na.rm = TRUE),
                      n_annotated = sum(!is.na(df$compound) & df$compound != "",
                                        na.rm = TRUE),
                      n_inhouse   = sum(df$reference == "In House Library",
                                        na.rm = TRUE),
                      n_multiple_match = sum(df$multiple_match, na.rm = TRUE)
                    )
                  }) |>
                  purrr::list_rbind()
              }) |>
              purrr::list_rbind()
          }) |>
          purrr::list_rbind()
      }) |>
      purrr::list_rbind()
  }) |>
  purrr::list_rbind()

rev_save_table(annotation_coverage, "annotation_coverage", "annotation")
print(annotation_coverage |>
        dplyr::filter(population == "all", covar_set == "covar"), n = 40)


# Save R objects for downstream analysis -------------------------------------

save(mwas_annotated_list_c18, mwas_annotated_list_hilic,
     mwas_full_annotated_list_c18, mwas_full_annotated_list_hilic,
     annotation_coverage,
     file = rev_here("data", "metabolomics", "results",
                     "mwas_annotation_revision.RData"))

message("\nAnnotation completed!")
message("Results saved to:")
message("  - ", rev_here("tables", "mwas_results"), " (_sig_annotated.xlsx)")
message("  - ", rev_here("tables", "annotation"))

#--------------------------------End of the code--------------------------------
