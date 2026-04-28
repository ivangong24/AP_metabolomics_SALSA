## ---------------------------
##
## Script name: 6-pathway_analysis.R
## Purpose of script: To perform pathway analysis using Mummichog
##                    for significant metabolites from MWAS
##
## Author: Yufan Gong
##
## Date Created: 2026-01-29
##
## Date Modified: 2026-04-01
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: This script performs pathway enrichment analysis using:
##        1. Mummichog (via MetaboAnalystR)
##        2. Creates input files for external pathway tools
##
## ---------------------------


# Load MWAS results and annotation data --------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "2-load_data.R"))
load(here::here("data", "processed", "combined_data_list_new.RData"))
load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))

# Create output directories --------------------------------------------------

exp_toexclude_total <- quote_all(cox)
exp_toexclude_cox <- quote_all(exp, wqs, qgcomp_all, qgcomp_traffic, 
                               qgcomp_metal, comp_pca)

# clean up the combined datalist
combined_data_list_new <- combined_data_list_new |> 
  purrr::imap(function(datalist, study){
    datalist |> 
      purrr::imap(function(data_list, population){
        data_list |> 
          purrr::imap(function(data, covar_name){
            if (study == "total") {
              data |> 
                dplyr::select(
                  -matches(str_c(exp_toexclude_total, collapse = "|"))
                )
            } else {
              data |> 
                dplyr::select(
                  -matches(str_c(exp_toexclude_cox, collapse = "|"))
                )
            }
          })
      })
  })


covar_list <- list(
  covar = quote_all(age_at_blooddraw, gender, edu_year, mh62, 
                    ruca_metro, nses,
                    wave, batch, demcind),
  
  covar_sen = quote_all(age_at_blooddraw, gender, edu_year, mh62, 
                        ruca_metro, nses,
                        alcohol_drinking, pa3_met_if_ca,  
                        bmi_at_blooddraw, diab_at_blooddraw, 
                        wave, batch, demcind)
  
)

exposure_vars_list <- combined_data_list_new |> 
  purrr::map(function(datalist){
    datalist[["all"]][["covar"]] |>
      dplyr::select(starts_with("comp_"), ends_with("iqr")) |>
      colnames()
  })

list("Input", "Output") |> 
  purrr::map(function(dir){
    combined_data_list_new |> 
      purrr::imap(function(data, study){
        names(data) |> 
          purrr::walk(function(population){
            names(covar_list) |> 
              purrr::walk(function(covar_name){
                dir.create(here::here("metaboAnalyst", dir, study, 
                                      population, covar_name), 
                           showWarnings = FALSE, recursive = TRUE)
              })
          })
      })
  })

# =============================================================================
# SECTION 1: PREPARE MUMMICHOG INPUT FILES
# =============================================================================

# Function to create Mummichog input format ----------------------------------

create_mummichog_input <- function(mwas_result, mz_rt_link_df, mode) {
  # Mummichog input format:
  # m.z | rt | p.value | t.score | mode

  # Get chemical_ID column name (may vary)
  id_col <- intersect(
    c("met", "mz_rt"),
    tolower(colnames(mz_rt_link_df))
  )[1]

  # Get m/z and retention time columns
  mz_col <- intersect(c("mz", "m.z", "mass"), 
                      tolower(colnames(mz_rt_link_df)))[1]
  rt_col <- intersect(c("time", "rt", "retention_time"), 
                      tolower(colnames(mz_rt_link_df)))[1]

  mwas_result |>
    tibble::rownames_to_column("met") |>
    dplyr::left_join(
      mz_rt_link_df |>
        dplyr::rename(met = !!sym(id_col)),
      by = "met"
    ) |>
    dplyr::transmute(
      `m.z` = .data[[mz_col]],
      `rt` = .data[[rt_col]],
      `p.value` = P.Value,
      `t.score` = t,
      mode = mode
    ) |>
    dplyr::filter(!is.na(`m.z`)) |>
    dplyr::arrange(`p.value`)
}

# Create input files for each exposure ---------------------------------------

# exposure_vars <- names(combined_results_list_c18[["all"]][["covar"]])

list(
  list(mwas_results_list_c18, mwas_results_list_hilic),
  list(c18_mz_rt_link_mz_links, hil_mz_rt_link_mz_links),
  list("negative", "positive")
) |> 
  purrr::pmap(function(mwas_results_data_list, mz_rt_link_df, mode){
    list(mwas_results_data_list, exposure_vars_list, 
         names(exposure_vars_list)) |> 
      purrr::pmap(function(mwas_results_list, exposure_vars, study){
        mwas_results_list |> 
          purrr::imap(function(mwas_results_ls, population){
            mwas_results_ls |> 
              purrr::imap(function(mwas_results, covar_name){
                message(paste0("Creating Mummichog input for: ", 
                               study, "_", population, 
                               " - ", covar_name, " (", mode, ")"))
                exposure_vars |>
                  purrr::set_names() |>
                  purrr::map(function(exp) {
                    create_mummichog_input(
                      mwas_results[[exp]],
                      mz_rt_link_df,
                      mode = mode
                    )
                  })
              })
          })
      })
  }) |> 
  purrr::set_names("mummichog_input_list_c18", "mummichog_input_list_hilic") |>
  list2env(.GlobalEnv)


# Combine C18 and HILIC for each exposure ------------------------------------

list(mummichog_input_list_c18, mummichog_input_list_hilic, 
     exposure_vars_list) |> 
  purrr::pmap(function(c18_datalist, hilic_datalist, exposure_vars){
    list(c18_datalist, hilic_datalist) |> 
      purrr::pmap(function(c18_list_ls, hilic_list_ls){
        list(c18_list_ls, hilic_list_ls) |> 
          purrr::pmap(function(c18_list, hilic_list){
            exposure_vars |>
              purrr::set_names() |>
              purrr::map(function(exp) {
                dplyr::bind_rows(
                  c18_list[[exp]],
                  hilic_list[[exp]]
                ) |>
                  dplyr::arrange(`p.value`)
              })
          })
      }) 
  }) -> mummichog_input_list_combined


# mummichog_input_combined <- exposure_vars |>
#   purrr::set_names() |>
#   purrr::map(function(exp) {
#     dplyr::bind_rows(
#       mummichog_input_c18[[exp]],
#       mummichog_input_hilic[[exp]]
#     ) |>
#       dplyr::arrange(`p.value`)
#   })


# Write input files ----------------------------------------------------------

mummichog_input_list_combined |> 
  purrr::imap(function(datalist, study){
    datalist |> 
      purrr::imap(function(dflist, population) {
        dflist |> 
          purrr::imap(function(dfls, covar_name){
            dfls |> 
              purrr::imap(function(df, exp) {
                message(paste0("Creating Mummichog input for: ", 
                               study, "_", population, 
                               " - ", covar_name, " - ", exp))
                write.table(
                  df,
                  file = here::here("metaboAnalyst", "Input", study, population, 
                                    covar_name,
                                    paste0("mwas_", exp, "_", 
                                           study, "_", population, 
                                           "_", covar_name, ".txt")),
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE,
                  sep = "\t"
                )
              })
          })
      })
  })
  

message("Mummichog input files created in metaboAnalyst/Input/")


# =============================================================================
# SECTION 2: RUN MUMMICHOG PATHWAY ANALYSIS
# =============================================================================

source(here::here("scripts", "mummichog_pathway.R"))

# Create output directories for each exposure ----------------------------------

# exposure_vars <- names(combined_results_list_c18[["all"]][["covar"]])

list(combined_results_list_c18, exposure_vars_list, 
     names(exposure_vars_list)) |>
  purrr::pmap(function(data, exposure_vars, study){
    names(data) |> 
      purrr::walk(function(population) {
        names(covar_list) |>
          purrr::walk(function(covar_name) {
            exposure_vars |>
              purrr::walk(function(exp_name) {
                dir.create(here::here("metaboAnalyst", "Output", study,
                                      population, covar_name, exp_name),
                           showWarnings = FALSE, recursive = TRUE)
              })
          })
      })
  })


# Run Mummichog for each exposure and population --------------------------------

message("Running Mummichog pathway analysis...")

system.time({
  combined_results_list_c18 |> 
    purrr::imap(function(data, study){
      names(data) |>
        purrr::set_names() |>
        purrr::map(function(population) {
          names(covar_list) |>
            purrr::set_names() |>
            purrr::map(function(covar_name) {
              input_dir <- here::here("metaboAnalyst", "Input", study,
                                      population, covar_name)
              input_files <- list.files(input_dir, pattern = "\\.txt$",
                                        full.names = TRUE)
              input_files |>
                purrr::set_names(
                  basename(input_files) |>
                    stringr::str_remove("\\.txt$") |>
                    stringr::str_remove("^mwas_") |>
                    stringr::str_remove(paste0("_", study, "_", population,
                                               "_", covar_name, "$"))
                ) |>
                purrr::imap(function(input_file, exp_name) {
                  message(paste0("\n--- Running Mummichog for: ",
                                 exp_name, " (", study, "_", population, " - ",
                                 covar_name, ") ---"))
                  
                  output_dir <- here::here("metaboAnalyst", "Output",
                                            study, population, covar_name,
                                            exp_name)
                  result <- tryCatch(
                    run_mummichog(
                      input_file = input_file,
                      output_dir = output_dir,
                      p_cutoff = 0.1,
                      organism = "hsa_mfn",
                      instrument_ppm = 10.0,
                      ion_mode = "mixed",
                      adducts = c("M-H [1-]", "M-2H [2-]",
                                  "M-H2O-H [1-]", "M [1+]",
                                  "M+H [1+]", "M+Na [1+]"),
                      min_hits = 3,
                      num_permutations = 100
                    ),
                    error = function(e) {
                      warning(paste0("Mummichog failed for ", exp_name,
                                     " (", study, "_", population, " - ",
                                     covar_name, "): ", e$message))
                      return(NULL)
                    }
                  )
                  # Remove large mum.RData to free disk space
                  mum_rdata <- file.path(output_dir, "mum.RData")
                  if (file.exists(mum_rdata)) file.remove(mum_rdata)
                  result
                })
            })
        }) 
    }) -> mummichog_results_combined
})

# Save Mummichog R objects for downstream analysis
save(mummichog_results_combined,
     file = here::here("data", "metabolomics", "results",
                       "mummichog_results_all.RData"))

# Extract and save Mummichog result tables -------------------------------------

mummichog_results_combined |>
  purrr::imap(function(pop_list, study) {
    pop_list |>
      purrr::imap(function(covar_list_res, population) {
        covar_list_res |>
          purrr::imap(function(result_list, covar_name) {
            # Save result tables to Excel
            result_tables <- result_list |>
              purrr::compact() |>
              purrr::map(~ .x$result_table)

            if (length(result_tables) > 0) {
              dir.create(here::here("tables", "mummichog_results", study,
                                    population, covar_name),
                         showWarnings = FALSE, recursive = TRUE)

              result_tables |>
                purrr::imap(function(tbl, exp_name) {
                  writexl::write_xlsx(
                    tbl,
                    path = here::here(
                      "tables", "mummichog_results", study, population, covar_name,
                      glue::glue("mummichog_{exp_name}_{study}_{population}_{covar_name}.xlsx"))
                  )
                })
            }

            # Save plots
            dir.create(here::here("figures", "mummichog", study,
                                  population, covar_name),
                       showWarnings = FALSE, recursive = TRUE)

            result_list |>
              purrr::compact() |>
              purrr::imap(function(res, exp_name) {
                if (!is.null(res$plot)) {
                  ggsave(
                    filename = here::here(
                      "figures", "mummichog", study, population, covar_name,
                      glue::glue("mummichog_{exp_name}_{study}_{population}_{covar_name}.png")),
                    plot = res$plot +
                      ggtitle(paste0("Mummichog: ", exp_name)),
                    width = 10, height = 8, dpi = 300
                  )
                }
              })
          })
      })
  })

message("Mummichog pathway analysis completed!")
message("Results saved to:")
message("  - metaboAnalyst/Output/ (MetaboAnalystR output)")
message("  - tables/mummichog_results/ (result tables)")
message("  - figures/mummichog/ (bubble plots)")




# =============================================================================
# SECTION 3: RUN METAPONE PATHWAY ANALYSIS (ALTERNATIVE)
# =============================================================================

source(here::here("scripts", "metapone_pathway.R"))

# Write Metapone input files (reuse combined feature tables from Section 1) ----

combined_data_list_new |> 
  purrr::imap(function(data, study){
    names(data) |> 
      purrr::walk(function(population) {
        names(covar_list) |> 
          purrr::walk(function(covar_name) {
            dir.create(here::here("Metapone", "Input", study, 
                                  population, covar_name),
                       showWarnings = FALSE, recursive = TRUE)
          })
      })
  })


mummichog_input_list_combined |>
  purrr::imap(function(datalist, study) {
    datalist |> 
      purrr::imap(function(dflist, population) {
        dflist |>
          purrr::imap(function(dfls, covar_name){
            dfls |> 
              purrr::imap(function(df, exp) {
                message(paste0("Creating Metapone input for: ", 
                               study, "_", population, 
                               " - ", covar_name, " - ", exp))
                write.table(
                  df,
                  file = here::here("Metapone", "Input", study, 
                                    population, covar_name,
                                    paste0("mwas_", exp, "_", study, "_",
                                           population, "_", 
                                           covar_name, ".txt")),
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE,
                  sep = "\t"
                )
              })
          })
      })
  })


message("Metapone input files created in Metapone/Input/")

# Run metapone for each exposure and population --------------------------------

message("Running metapone pathway analysis...")

system.time({
  combined_results_list_c18 |> 
    purrr::imap(function(data, study){
      names(data) |>
        purrr::set_names() |>
        purrr::map(function(population) {
          names(covar_list) |>
            purrr::set_names() |>
            purrr::map(function(covar_name) {
              input_dir <- here::here("Metapone", "Input", 
                                      study ,population, covar_name)
              input_files <- list.files(input_dir, pattern = "\\.txt$",
                                        full.names = TRUE)
              input_files |>
                purrr::set_names(
                  basename(input_files) |>
                    stringr::str_remove("\\.txt$") |>
                    stringr::str_remove(paste0("^mwas_")) |>
                    stringr::str_remove(paste0("_", study, "_", population, 
                                               "_", covar_name, "$"))
                ) |>
                purrr::imap(function(input_file, exp_name) {
                  message(paste0("\n---Running metapone for: ", 
                                 exp_name, " (", study, "_", population, " - ", 
                                 covar_name, ") ---"))
                  
                  tryCatch(
                    run_metapone(
                      input_file = input_file,
                      p_cutoff = 0.05,
                      num_permutations = 200,
                      match_tol_ppm = 10,
                      pos.adductlist = c("M+H", "M+Na", "M+"),
                      neg.adductlist = c("M-H", "M-2H", "M-H2O-H")
                    ),
                    error = function(e) {
                      warning(paste0("metapone failed for ", exp_name,
                                     " (", study, "_", population, " - ", 
                                     covar_name, "): ", e$message))
                      return(NULL)
                    }
                  )
                })
            })
        })
    }) -> metapone_results_combined
})

# Save metapone R objects for downstream analysis
save(metapone_results_combined,
     file = here::here("data", "metabolomics", "results",
                       "metapone_results_all.RData"))

# Extract and save metapone result tables --------------------------------------

metapone_results_combined |>
  purrr::imap(function(pop_list, study) {
    pop_list |>
      purrr::imap(function(covar_list_res, population) {
        covar_list_res |>
          purrr::imap(function(result_list, covar_name) {
            # Save result tables to Excel
            result_tables <- result_list |>
              purrr::compact() |>
              purrr::map(~ .x$result_table)

            if (length(result_tables) > 0) {
              dir.create(here::here("tables", "metapone_results", study,
                                    population, covar_name),
                         showWarnings = FALSE, recursive = TRUE)

              result_tables |>
                purrr::imap(function(tbl, exp_name) {
                  writexl::write_xlsx(
                    tbl,
                    path = here::here(
                      "tables", "metapone_results", study, population, covar_name,
                      glue::glue("metapone_{exp_name}_{study}_{population}_{covar_name}.xlsx"))
                  )
                })
            }

            # Save plots
            dir.create(here::here("figures", "metapone", study,
                                  population, covar_name),
                       showWarnings = FALSE, recursive = TRUE)

            result_list |>
              purrr::compact() |>
              purrr::imap(function(res, exp_name) {
                if (!is.null(res$plot)) {
                  ggsave(
                    filename = here::here(
                      "figures", "metapone", study, population, covar_name,
                      glue::glue("metapone_{exp_name}_{study}_{population}_{covar_name}.png")),
                    plot = res$plot +
                      ggtitle(paste0("metapone: ", exp_name)),
                    width = 10, height = 8, dpi = 300
                  )
                }
              })
          })
      })
  })



message("metapone pathway analysis completed!")
message("Results saved to:")
message("  - Metapone/Input/ (input files)")
message("  - tables/metapone_results/ (result tables)")
message("  - figures/metapone/ (bubble plots)")


# =============================================================================
# SECTION 4: ALTERNATIVE - MANUAL METABOANALYST INPUT
# =============================================================================

# If MetaboAnalystR fails, create input for web interface --------------------

# The input files created in Section 1 can be uploaded to:
# https://www.metaboanalyst.ca/MetaboAnalyst/ModuleView.xhtml

message("\nAlternative: Manual pathway analysis")
message("If MetaboAnalystR fails, upload input files to MetaboAnalyst web interface:")
message("https://www.metaboanalyst.ca/")
message("Input files are located in: metaboAnalyst/Input/")


# Output message ----------------------------------------------------------

message("\nPathway analysis completed!")
message("Results saved to:")
message("  - metaboAnalyst/Input/ (input files)")
message("  - metaboAnalyst/Output/ (analysis results)")
message("  - tables/mwas_results/ (summary tables)")

#--------------------------------End of the code--------------------------------
