## ---------------------------
##
## Script name: R6-pathway_revision.R
## Purpose of script: To perform pathway analysis for the R1 revision MWAS
##                    (unsupervised PCA and cross-fitted mixture indices)
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
## Notes: This is the revision analogue of scripts/6-pathway_analysis.R.
##        Same two algorithms, same parameters, revision exposures:
##        1. Mummichog (via MetaboAnalystR)
##        2. metapone
##
##        Feature-level input is the FULL MWAS result table (all ~20,000
##        features, C18 and HILIC stacked), exactly as in the submitted
##        analysis - both algorithms take the whole p-value distribution, not a
##        pre-filtered feature list.
##
##        PERMUTATIONS. Both algorithms run at 1000 permutations, answering
##        Reviewer 1 major comment 5 (the submitted analysis used 100 for
##        Mummichog and 200 for metapone). Runtime scales roughly linearly with
##        the permutation count, so this is the dominant cost of the script.
##        Note that the revision pathway p-values are therefore estimated at a
##        finer resolution than the submitted ones: a permutation p reported as
##        0 now means < 1e-3 rather than < 1e-2, which matters when comparing
##        the two sets of results side by side.
##
##        Outputs -> revision_output/metaboAnalyst/{Input,Output}/...
##                   revision_output/Metapone/Input/...
##                   revision_output/{tables,figures}/{mummichog,metapone}_results/
##        Downstream: R7-visualization_revision.R
##
##        RUNTIME: several hours. Validate first with
##        Sys.setenv(SALSA_REVISION_QUICK = "true").
##
## ---------------------------

# Load MWAS results and annotation data --------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R6-pathway_revision.R")

MUMMICHOG_PERM <- rev_n(1000, 10)
METAPONE_PERM  <- rev_n(1000, 10)

load(rev_here("data", "processed", "combined_data_list_revision.RData"))
load(rev_here("data", "metabolomics", "results",
              "mwas_results_all_revision.RData"))

## m/z - retention time links. 2-load_data.R renames these objects on load;
## reading the two .Rdata files directly avoids re-running the whole data load
## just for two 10,000-row lookup tables.
load(here::here("data", "metabolomics", "mz_links", "c18_mz_rt_link.Rdata"))
load(here::here("data", "metabolomics", "mz_links", "hil_mz_rt_link.Rdata"))

c18_mz_rt_link_mz_links <- c18_mz_rt_link
hil_mz_rt_link_mz_links <- hil_mz_rt_link


## In QUICK mode cut the analysis grid down to one population, one covariate
## set and two exposures per study, so the whole script can be validated end
## to end in minutes. Nothing produced this way is reportable.
if (REV_QUICK) {
  ## Keep the post-diagnosis sensitivity population as well as `all`: it is
  ## the newest slice of the grid and the one a smoke test most needs to
  ## exercise. It exists in `total` only, hence the intersect.
  keep_population <- c("all", "all predx")
  keep_covar_set  <- "covar"

  trim_grid <- function(x) {
    x |> purrr::map(function(pop_list) {
      pop_list[intersect(keep_population, names(pop_list))] |>
        purrr::map(~ .x[keep_covar_set])
    })
  }

  ## Two composites plus one single pollutant, so the smoke test exercises the
  ## per-population exposure set (the pollutants are `all` only, so `all predx`
  ## gets the composites alone -- which is exactly the case that used to be
  ## driven by a per-study vector and now goes through exposures_for()).
  exposure_vars_list <- exposure_vars_list |>
    purrr::map(~ c(head(.x[!is_single_pollutant(.x)], 2),
                   head(.x[is_single_pollutant(.x)], 1)))
  combined_results_list_c18   <- trim_grid(combined_results_list_c18)
  combined_results_list_hilic <- trim_grid(combined_results_list_hilic)
  mwas_results_list_c18       <- trim_grid(mwas_results_list_c18)
  mwas_results_list_hilic     <- trim_grid(mwas_results_list_hilic)
  combined_data_list_revision <- trim_grid(combined_data_list_revision)
  covar_list                  <- covar_list[keep_covar_set]

  message("QUICK: grid trimmed to populations '",
          paste(keep_population, collapse = "', '"),
          "', covariate set '", keep_covar_set, "', exposures ",
          paste(unlist(exposure_vars_list), collapse = ", "))
}


# Create output directories --------------------------------------------------

## The directory grid comes from the MWAS results, not from R3's data list.
## R4 adds the `all predx` post-diagnosis population as a row filter on `all`,
## so R3's populations are a strict subset of the ones that actually have
## results. Building directories from the smaller list left the input writer
## with nowhere to write ("cannot open the connection").
combined_data_list_new <- mwas_results_list_c18

list("Input", "Output") |>
  purrr::map(function(dir){
    combined_data_list_new |>
      purrr::imap(function(data, study){
        names(data) |>
          purrr::walk(function(population){
            names(covar_list) |>
              purrr::walk(function(covar_name){
                rev_dir("metaboAnalyst", dir, study, population, covar_name)
              })
          })
      })
  })

message("Exposure variables for the revision pathway analysis:")
print(exposure_vars_list)


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

list(
  list(mwas_results_list_c18, mwas_results_list_hilic),
  list(c18_mz_rt_link_mz_links, hil_mz_rt_link_mz_links),
  list("negative", "positive")
) |>
  purrr::pmap(function(mwas_results_data_list, mz_rt_link_df, mode){
    list(mwas_results_data_list, names(mwas_results_data_list)) |>
      purrr::pmap(function(mwas_results_list, study){
        mwas_results_list |>
          purrr::imap(function(mwas_results_ls, population){
            mwas_results_ls |>
              purrr::imap(function(mwas_results, covar_name){
                message(paste0("Creating Mummichog input for: ",
                               study, "_", population,
                               " - ", covar_name, " (", mode, ")"))
                ## exposures_for(): the single pollutants are fitted in `all`
                ## only, so the exposure set varies by population.
                exposures_for(study, population) |>
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

list(mummichog_input_list_c18, mummichog_input_list_hilic) |>
  purrr::pmap(function(c18_datalist, hilic_datalist){
    list(c18_datalist, hilic_datalist) |>
      purrr::pmap(function(c18_list_ls, hilic_list_ls){
        list(c18_list_ls, hilic_list_ls) |>
          purrr::pmap(function(c18_list, hilic_list){
            names(c18_list) |>
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


# Write input files ----------------------------------------------------------

mummichog_input_list_combined |>
  purrr::imap(function(datalist, study){
    datalist |>
      purrr::imap(function(dflist, population) {
        dflist |>
          purrr::imap(function(dfls, covar_name){
            dfls |>
              purrr::imap(function(df, exp) {
                message(paste0("Writing Mummichog input for: ",
                               study, "_", population,
                               " - ", covar_name, " - ", exp))
                write.table(
                  df,
                  file = rev_here("metaboAnalyst", "Input", study, population,
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


message("Mummichog input files created in ",
        rev_here("metaboAnalyst", "Input"))


# =============================================================================
# SECTION 2: RUN MUMMICHOG PATHWAY ANALYSIS
# =============================================================================

source(here::here("scripts", "mummichog_pathway.R"))

# Create output directories for each exposure ---------------------------------

list(combined_results_list_c18, names(combined_results_list_c18)) |>
  purrr::pmap(function(data, study){
    names(data) |>
      purrr::walk(function(population) {
        names(covar_list) |>
          purrr::walk(function(covar_name) {
            exposures_for(study, population) |>
              purrr::walk(function(exp_name) {
                rev_dir("metaboAnalyst", "Output", study,
                        population, covar_name, exp_name)
              })
          })
      })
  })


# Run Mummichog for each exposure and population ------------------------------

message("Running Mummichog pathway analysis (",
        MUMMICHOG_PERM, " permutations)...")

system.time({
  combined_results_list_c18 |>
    purrr::imap(function(data, study){
      names(data) |>
        purrr::set_names() |>
        purrr::map(function(population) {
          names(covar_list) |>
            purrr::set_names() |>
            purrr::map(function(covar_name) {
              input_dir <- rev_here("metaboAnalyst", "Input", study,
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

                  output_dir <- rev_here("metaboAnalyst", "Output",
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
                      num_permutations = MUMMICHOG_PERM
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

## The tryCatch above turns a failure into a NULL slot and a warning, and
## warnings are easy to lose in a log this size. MetaboAnalystR pulls
## `mixed_adduct.qs` and `hsa_mfn.qs` off metaboanalyst.ca on every run, so a
## transient server error can silently cost one stratum -- and it did once, for
## cox/all/covar. Name any empty slots loudly before anything downstream
## quietly skips them.
mummichog_failures <- mummichog_results_combined |>
  purrr::imap(function(pop_list, study) {
    pop_list |>
      purrr::imap(function(covar_list_res, population) {
        covar_list_res |>
          purrr::imap(function(result_list, covar_name) {
            names(result_list)[purrr::map_lgl(result_list, is.null)] |>
              purrr::map_chr(~ paste(study, population, covar_name, .x,
                                     sep = " / "))
          })
      })
  }) |>
  unlist(use.names = FALSE)

if (length(mummichog_failures) > 0) {
  message("\n", strrep("!", 78))
  message("Mummichog produced no result for ", length(mummichog_failures),
          " stratum/strata:")
  purrr::walk(mummichog_failures, ~ message("  - ", .x))
  message("Re-run those strata before reporting. If the log shows ",
          "'Loading files from server unsuccessful', delete the stale ",
          "mixed_adduct.qs in the output folder first -- a served error page ",
          "is cached there and reused.")
  message(strrep("!", 78), "\n")
}

# Save Mummichog R objects for downstream analysis
rev_dir("data", "metabolomics", "results")

save(mummichog_results_combined,
     file = rev_here("data", "metabolomics", "results",
                     "mummichog_results_all_revision.RData"))

# Extract and save Mummichog result tables ------------------------------------

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
              rev_dir("tables", "mummichog_results", study,
                      population, covar_name)

              result_tables |>
                purrr::imap(function(tbl, exp_name) {
                  writexl::write_xlsx(
                    tbl,
                    path = rev_here(
                      "tables", "mummichog_results", study, population,
                      covar_name,
                      glue::glue("mummichog_{exp_name}_{study}_{population}_{covar_name}.xlsx"))
                  )
                })
            }

            # Save plots
            rev_dir("figures", "mummichog", study, population, covar_name)

            result_list |>
              purrr::compact() |>
              purrr::imap(function(res, exp_name) {
                if (!is.null(res$plot)) {
                  ggsave(
                    filename = rev_here(
                      "figures", "mummichog", study, population, covar_name,
                      glue::glue("mummichog_{exp_name}_{study}_{population}_{covar_name}.png")),
                    plot = res$plot +
                      ggtitle(paste0("Mummichog: ", rev_label(exp_name))),
                    width = 10, height = 8, dpi = 300
                  )
                }
              })
          })
      })
  })

message("Mummichog pathway analysis completed!")


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
            rev_dir("Metapone", "Input", study, population, covar_name)
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
                message(paste0("Writing Metapone input for: ",
                               study, "_", population,
                               " - ", covar_name, " - ", exp))
                write.table(
                  df,
                  file = rev_here("Metapone", "Input", study,
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


message("Metapone input files created in ", rev_here("Metapone", "Input"))

# Run metapone for each exposure and population -------------------------------

## PARALLEL ACROSS CELLS (changed 2026-09-03).
##
## metapone is the long pole of the whole revision: ~430 s per cell at 1000
## permutations, and the 3- and 10-year exposure arms roughly tripled the cell
## count. Run one cell per worker instead of one after another.
##
## Each call is independent -- run_metapone() reads one input file, loads the
## HMDB and pathway databases into its own environment and returns a result
## object, holding no state between calls -- so this changes the wall clock and
## nothing else. Workers are capped well below the core count because each one
## holds its own copy of those databases; this is a memory knob, not a speed
## one, and raising it past the point where the machine swaps makes the whole
## stage slower.
##
## Seeding goes through furrr's L'Ecuyer streams, so the permutation p-values
## are reproducible per cell and independent of how many workers ran and in
## what order. The sequential version set no seed at all.
METAPONE_WORKERS <- as.integer(
  Sys.getenv("SALSA_METAPONE_WORKERS",
             unset = max(1, min(6, future::availableCores() - 1))))

message("Running metapone pathway analysis (", METAPONE_PERM,
        " permutations) on ", METAPONE_WORKERS, " workers...")

## Flatten the grid to one job per cell so every worker gets a whole cell, then
## re-nest. The nesting is rebuilt from the job keys rather than carried
## through the parallel call, so a failed cell still lands in the right slot.
metapone_jobs <- combined_results_list_c18 |>
  purrr::imap(function(data, study){
    names(data) |>
      purrr::map(function(population) {
        names(covar_list) |>
          purrr::map(function(covar_name) {
            input_dir <- rev_here("Metapone", "Input",
                                  study, population, covar_name)
            input_files <- list.files(input_dir, pattern = "\\.txt$",
                                      full.names = TRUE)
            exp_names <- basename(input_files) |>
              stringr::str_remove("\\.txt$") |>
              stringr::str_remove("^mwas_") |>
              stringr::str_remove(paste0("_", study, "_", population,
                                         "_", covar_name, "$"))

            list(input_files, exp_names) |>
              purrr::pmap(function(input_file, exp_name){
                list(study = study, population = population,
                     covar_set = covar_name, exposure = exp_name,
                     input_file = input_file)
              })
          })
      })
  }) |>
  unlist(recursive = FALSE) |> unlist(recursive = FALSE) |>
  unlist(recursive = FALSE)

names(metapone_jobs) <- metapone_jobs |>
  purrr::map_chr(~ paste(.x$study, .x$population, .x$covar_set, .x$exposure,
                         sep = "|"))

message("  ", length(metapone_jobs), " metapone cells queued")

set.seed(42)
future::plan(future::multisession, workers = METAPONE_WORKERS)

system.time({
  metapone_flat <- metapone_jobs |>
    furrr::future_map(function(job) {
      tryCatch(
        run_metapone(
          input_file = job$input_file,
          p_cutoff = 0.05,
          num_permutations = METAPONE_PERM,
          match_tol_ppm = 10,
          pos.adductlist = c("M+H", "M+Na", "M+"),
          neg.adductlist = c("M-H", "M-2H", "M-H2O-H")
        ),
        error = function(e) {
          warning(paste0("metapone failed for ", job$exposure,
                         " (", job$study, "_", job$population, " - ",
                         job$covar_set, "): ", e$message), call. = FALSE)
          NULL
        }
      )
    },
    ## one cell per chunk: they differ in how many features pass p < 0.05, so
    ## static chunking leaves workers idle at the end.
    ##
    ## `packages` is not optional. run_metapone() namespaces its metapone and
    ## dplyr calls, but create_metapone_plot() -- which it calls to build the
    ## bubble plot -- uses bare ggplot(), aes() and geom_text_repel(). Those
    ## resolve in the parent because scripts/metapone_pathway.R attaches
    ## ggplot2 and ggrepel at the top level, which a fresh worker does not
    ## inherit, and the failure would arrive as "could not find function
    ## 'ggplot'" after the permutations had already been paid for.
    .options = furrr::furrr_options(
      seed = TRUE, chunk_size = 1,
      packages = c("metapone", "ggplot2", "ggrepel", "dplyr", "tibble")),
    .progress = TRUE)
})

future::plan(future::sequential)

n_failed <- sum(purrr::map_lgl(metapone_flat, is.null))
if (n_failed > 0) {
  message("  ", n_failed, " of ", length(metapone_flat),
          " metapone cells returned NULL (see warnings above)")
}

## Re-nest to [study][population][covar_set][exposure], the shape everything
## downstream reads.
metapone_results_combined <- combined_results_list_c18 |>
  purrr::imap(function(data, study){
    names(data) |>
      purrr::set_names() |>
      purrr::map(function(population) {
        names(covar_list) |>
          purrr::set_names() |>
          purrr::map(function(covar_name) {
            keys <- names(metapone_jobs) |>
              purrr::keep(~ startsWith(.x, paste(study, population, covar_name,
                                                 "", sep = "|")))
            metapone_flat[keys] |>
              purrr::set_names(purrr::map_chr(metapone_jobs[keys],
                                              ~ .x$exposure))
          })
      })
  })

rm(metapone_flat, metapone_jobs)
gc()

# Save metapone R objects for downstream analysis
save(metapone_results_combined,
     file = rev_here("data", "metabolomics", "results",
                     "metapone_results_all_revision.RData"))

# Extract and save metapone result tables -------------------------------------

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
              rev_dir("tables", "metapone_results", study,
                      population, covar_name)

              result_tables |>
                purrr::imap(function(tbl, exp_name) {
                  writexl::write_xlsx(
                    tbl,
                    path = rev_here(
                      "tables", "metapone_results", study, population,
                      covar_name,
                      glue::glue("metapone_{exp_name}_{study}_{population}_{covar_name}.xlsx"))
                  )
                })
            }

            # Save plots
            rev_dir("figures", "metapone", study, population, covar_name)

            result_list |>
              purrr::compact() |>
              purrr::imap(function(res, exp_name) {
                if (!is.null(res$plot)) {
                  ggsave(
                    filename = rev_here(
                      "figures", "metapone", study, population, covar_name,
                      glue::glue("metapone_{exp_name}_{study}_{population}_{covar_name}.png")),
                    plot = res$plot +
                      ggtitle(paste0("metapone: ", rev_label(exp_name))),
                    width = 10, height = 8, dpi = 300
                  )
                }
              })
          })
      })
  })


message("metapone pathway analysis completed!")


# =============================================================================
# SECTION 4: PATHWAY SUMMARY ACROSS EXPOSURES
# =============================================================================

## One long table of every pathway result, so the response letter can quote
## how the enrichment picture changes when the outcome-informed index is
## replaced by the unsupervised one.

extract_pathway_rows <- function(results, algorithm) {
  results |>
    purrr::imap(function(pop_list, study) {
      pop_list |>
        purrr::imap(function(covar_res, population) {
          covar_res |>
            purrr::imap(function(result_list, covar_name) {
              result_list |>
                purrr::compact() |>
                purrr::imap(function(res, exp_name) {
                  tbl <- res$result_table
                  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)

                  ## Take the RAW enrichment p-value: "P(Fisher)"/"FET" for
                  ## Mummichog, "p_value" for metapone. metapone's `adjust.p`
                  ## and `lfdr` are NA for every pathway without hits (359 of
                  ## 372 in a typical run), so ranking on them silently throws
                  ## away the real signal. FDR is recomputed below instead.
                  p_col <- intersect(
                    c("P(Fisher)", "FET", "p_value", "p-value",
                      "P(Gamma)", "Gamma"),
                    names(tbl))[1]
                  if (is.na(p_col)) {
                    p_col <- intersect(c("adjust.p", "lfdr"), names(tbl))[1]
                  }
                  if (is.na(p_col)) return(NULL)

                  hit_col <- intersect(
                    c("Hits.sig", "n_significant_metabolites"), names(tbl))[1]
                  tot_col <- intersect(
                    c("Hits.total", "Pathway.total", "Hits.all",
                      "n_mapped_metabolites"), names(tbl))[1]

                  tibble::tibble(
                    algorithm  = algorithm,
                    study      = study,
                    population = population,
                    covar_set  = covar_name,
                    exposure   = exp_name,
                    label      = unname(rev_label(exp_name)),
                    pathway    = tbl$pathway,
                    p_value    = as.numeric(tbl[[p_col]]),
                    p_column   = p_col,
                    n_sig      = if (is.na(hit_col)) NA_real_
                                 else as.numeric(tbl[[hit_col]]),
                    n_total    = if (is.na(tot_col)) NA_real_
                                 else as.numeric(tbl[[tot_col]])
                  )
                }) |>
                purrr::list_rbind()
            }) |>
            purrr::list_rbind()
        }) |>
        purrr::list_rbind()
    }) |>
    purrr::list_rbind()
}

pathway_results_long <- dplyr::bind_rows(
  extract_pathway_rows(mummichog_results_combined, "mummichog"),
  extract_pathway_rows(metapone_results_combined,  "metapone")
) |>
  dplyr::filter(!is.na(p_value)) |>
  dplyr::group_by(algorithm, study, population, covar_set, exposure) |>
  dplyr::mutate(
    p_fdr = stats::p.adjust(p_value, method = "BH"),
    ## A permutation p of exactly 0 means "below the resolution of the
    ## permutation count", not "impossible". Keep the raw value but carry a
    ## floored copy so -log10() stays finite in the figures.
    p_plot = pmax(p_value,
                  min(c(p_value[p_value > 0], 1), na.rm = TRUE) / 2)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(algorithm, study, population, covar_set, exposure, p_value)

rev_save_table(pathway_results_long, "pathway_results_long", "pathway")

pathway_summary <- pathway_results_long |>
  dplyr::group_by(algorithm, study, population, covar_set, exposure, label) |>
  dplyr::summarise(
    n_pathways_tested = dplyr::n(),
    p_column          = dplyr::first(p_column),
    n_p05             = sum(p_value < 0.05, na.rm = TRUE),
    n_fdr05           = sum(p_fdr < 0.05, na.rm = TRUE),
    top_pathway       = paste(utils::head(pathway[order(p_value)], 3),
                              collapse = "; "),
    top_p             = min(p_value, na.rm = TRUE),
    .groups = "drop"
  )

rev_save_table(pathway_summary, "pathway_summary", "pathway")
print(pathway_summary |>
        dplyr::filter(population == "all", covar_set == "covar"), n = 40)

message("\nPathway analysis completed!")
message("Results saved to:")
message("  - ", rev_here("metaboAnalyst", "Input"), " / Output")
message("  - ", rev_here("Metapone", "Input"))
message("  - ", rev_here("tables", "mummichog_results"))
message("  - ", rev_here("tables", "metapone_results"))
message("  - ", rev_here("tables", "pathway"))
message("  - ", rev_here("figures", "mummichog"), " / metapone")

#--------------------------------End of the code--------------------------------
