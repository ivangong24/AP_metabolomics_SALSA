## ---------------------------
##
## Script name: 7-pathway_analysis.R
## Purpose of script: To perform pathway analysis using Mummichog
##                    for significant metabolites from MWAS
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
## Notes: This script performs pathway enrichment analysis using:
##        1. Mummichog (via MetaboAnalystR)
##        2. Creates input files for external pathway tools
##
##        Dependencies: Run scripts 1-6 before this script.
## ---------------------------


# Load MWAS results and annotation data --------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "2-load_data.R"))
load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))

# Create output directories --------------------------------------------------

covar_list <- list(
  covar = quote_all(age_at_blooddraw, gender, edu_year, mh62, 
                    wave, batch, demcind),
  
  covar_sen = quote_all(age_at_blooddraw, gender, edu_year, mh62,
                        alcohol_drinking, pa3_met_if_ca, nses, 
                        bmi_at_blooddraw, diab_at_blooddraw, 
                        wave, batch, demcind)
  
)

list("Input", "Output") |> 
  purrr::map(function(dir){
    names(combined_results_list_c18) |> 
      purrr::walk(function(population){
        names(covar_list) |> 
          purrr::walk(function(covar_name){
            dir.create(here::here("metaboAnalyst", dir, population, covar_name), 
                       showWarnings = FALSE, recursive = TRUE)
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

exposure_vars <- names(combined_results_list_c18[["all"]][["covar"]])

list(
  list(mwas_results_list_c18, mwas_results_list_hilic),
  list(c18_mz_rt_link_mz_links, hil_mz_rt_link_mz_links),
  list("negative", "positive")
) |> 
  purrr::pmap(function(mwas_results_list, mz_rt_link_df, mode){
    mwas_results_list |> 
      purrr::imap(function(mwas_results_ls, population){
        mwas_results_ls |> 
          purrr::imap(function(mwas_results, covar_name){
            message(paste0("Creating Mummichog input for: ", 
                           population, " - ", covar_name, " (", mode, ")"))
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
  }) |> 
  purrr::set_names("mummichog_input_list_c18", "mummichog_input_list_hilic") |>
  list2env(.GlobalEnv)


# Combine C18 and HILIC for each exposure ------------------------------------

list(mummichog_input_list_c18, mummichog_input_list_hilic) |> 
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
  purrr::imap(function(dflist, population) {
    dflist |> 
      purrr::imap(function(dfls, covar_name){
        dfls |> 
          purrr::imap(function(df, exp) {
            message(paste0("Creating Mummichog input for: ", 
                           population, " - ", covar_name, " - ", exp))
            write.table(
              df,
              file = here::here("metaboAnalyst", "Input", population, 
                                covar_name,
                                paste0("mwas_", exp, "_", population, 
                                       "_", covar_name, ".txt")),
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE,
              sep = "\t"
            )
          })
      })
  })

message("Mummichog input files created in metaboAnalyst/Input/")


# =============================================================================
# SECTION 2: RUN MUMMICHOG PATHWAY ANALYSIS
# =============================================================================

#### Create directories in "metaboAnalyst/" for each exposure

exposure_vars <- names(combined_results_list_c18[["all"]][["covar"]])

names(combined_results_list_c18) |> 
  purrr::map(function(population){
    names(covar_list) |> 
      purrr::map(function(covar_name){
        exposure_vars %>% 
          map(function(exp_name){
            dir.create(here::here("metaboAnalyst", "Output", 
                                  population, covar_name, exp_name), 
                       showWarnings = FALSE, recursive = TRUE)
          })
      })
  })

#### Get list of directories and input files for Mummichog analysis
wd_num_list <- names(combined_results_list_c18) |> 
  purrr::map(function(population){
    names(covar_list) |> 
      purrr::map(function(covar_name){
        here("metaboAnalyst", "Output", population, covar_name) |> 
          list.dirs(recursive = FALSE)
      })
  })

sub_dir_list1 <- names(combined_results_list_c18)

sub_dir_list2 <- names(covar_list)


input_list <- names(combined_results_list_c18) |> 
  purrr::map(function(population){
     names(covar_list) |> 
      purrr::map(function(covar_name){
        here("metaboAnalyst", "Input", population, covar_name)|> 
          list.files(full.names = TRUE, recursive = FALSE) |> 
          basename()
      })
  })

  
pathway_name_list <- input_list |>  
  purrr::map(function(files_list){
    files_list |> 
      purrr::map(function(files){
        files |> 
          str_remove(".txt") |>  
          map(function(name){
            str_c("p_pathway","_", name)
          })
      })
  })


# Function to run Mummichog using MetaboAnalystR -----------------------------

setwd(here::here())
source("scripts/function_mummichog.r")

# Run Mummichog for each exposure --------------------------------------------

# Run analysis (this may take some time)
message("Running Mummichog pathway analysis...")

system.time({
  list(wd_num_list,
       sub_dir_list1, 
       input_list,
       pathway_name_list) |> 
    purrr::pmap(function(wd_mum_ls, sub_dir1, input_ls, pathway_name_ls){
      list(wd_mum_ls, sub_dir_list2, input_ls, pathway_name_ls) |>  
        purrr::pmap(function(wd_mumls, sub_dir2, inputls, pathway_namels){
          list(wd_mumls, inputls) |> 
            purrr::pmap(function(wd_mum, input){
              message(paste0("Running Mummichog pathway analysis in ", 
                             sub_dir1, " ", sub_dir2, " using ", input, " ---"))
              mummichog(wd_mum, input, sub_dir1, sub_dir2)
            }) %>% 
            purrr::set_names(pathway_namels) |> 
            list2env(.GlobalEnv)
        })
    })
})


# Function to extract pathway results from mSet ------------------------------

extract_pathway_results <- function(mSet) {
  if (is.null(mSet)) return(NULL)
  
  # Extract Mummichog results
  mum_results <- tryCatch({
    mSet$mummi.resmat |>
      as.data.frame() |>
      tibble::rownames_to_column("pathway") |>
      dplyr::arrange(P.Value)
  }, error = function(e) NULL)
  
  # Extract GSEA results
  gsea_results <- tryCatch({
    mSet$gsea.resmat |>
      as.data.frame() |>
      tibble::rownames_to_column("pathway") |>
      dplyr::arrange(P.Value)
  }, error = function(e) NULL)
  
  return(list(
    mummichog = mum_results,
    gsea = gsea_results
  ))
}


# Extract results for all exposures ------------------------------------------

if (!is.null(mummichog_results)) {
  pathway_results <- mummichog_results |>
    purrr::map(extract_pathway_results)
  
  # Save pathway results
  save(pathway_results,
       file = here::here("data", "metabolomics", "results",
                         "pathway_results_all.RData"))
}

# Function to create pathway summary table -----------------------------------

create_pathway_summary <- function(pathway_results, p_threshold = 0.05) {
  if (is.null(pathway_results)) return(NULL)
  
  # Summarize Mummichog results
  mum_summary <- pathway_results |>
    purrr::imap(function(res, exp) {
      if (is.null(res$mummichog)) return(NULL)
      res$mummichog |>
        dplyr::filter(P.Value < p_threshold) |>
        dplyr::mutate(exposure = exp) |>
        dplyr::select(exposure, pathway, P.Value, everything())
    }) |>
    purrr::compact() |>
    purrr::list_rbind()
  
  # Summarize GSEA results
  gsea_summary <- pathway_results |>
    purrr::imap(function(res, exp) {
      if (is.null(res$gsea)) return(NULL)
      res$gsea |>
        dplyr::filter(P.Value < p_threshold) |>
        dplyr::mutate(exposure = exp) |>
        dplyr::select(exposure, pathway, P.Value, everything())
    }) |>
    purrr::compact() |>
    purrr::list_rbind()
  
  return(list(
    mummichog = mum_summary,
    gsea = gsea_summary
  ))
}


# Create and save summary tables ---------------------------------------------

if (exists("pathway_results") && !is.null(pathway_results)) {
  pathway_summary <- create_pathway_summary(pathway_results)
  
  # Save to Excel
  if (!is.null(pathway_summary$mummichog) && nrow(pathway_summary$mummichog) > 0) {
    writexl::write_xlsx(
      pathway_summary$mummichog,
      here::here("tables", "mwas_results", "pathway_summary_mummichog.xlsx")
    )
  }
  
  if (!is.null(pathway_summary$gsea) && nrow(pathway_summary$gsea) > 0) {
    writexl::write_xlsx(
      pathway_summary$gsea,
      here::here("tables", "mwas_results", "pathway_summary_gsea.xlsx")
    )
  }
}



# mummichog_results <- tryCatch({
#   exposure_vars |>
#     purrr::set_names() |>
#     purrr::map(function(exp) {
#       message(paste0("Processing: ", exp))
#       run_mummichog(
#         input_file = here::here("metaboAnalyst", "Input",
#                                 paste0("mwas_", exp, ".txt")),
#         output_dir = here::here("metaboAnalyst", "Output", exp),
#         exposure_name = exp,
#         p_cutoff = 0.05,
#         organism = "hsa"
#       )
#     })
# }, error = function(e) {
#   message("Mummichog analysis failed. Error: ", e$message)
#   message("Please run pathway analysis manually using MetaboAnalyst web interface.")
#   return(NULL)
# })




# =============================================================================
# SECTION 3: RUN METAPONE PATHWAY ANALYSIS (ALTERNATIVE)
# =============================================================================

source(here::here("scripts", "metapone_pathway.R"))

# Write Metapone input files (reuse combined feature tables from Section 1) ----

names(combined_results_list_c18) |>
  purrr::walk(function(population) {
    names(covar_list) |> 
      purrr::walk(function(covar_name) {
        dir.create(here::here("Metapone", "Input", population, covar_name),
                   showWarnings = FALSE, recursive = TRUE)
      })
  })

mummichog_input_list_combined |>
  purrr::imap(function(dflist, population) {
    dflist |>
      purrr::imap(function(dfls, covar_name){
        dfls |> 
          purrr::imap(function(df, exp) {
            message(paste0("Creating Metapone input for: ", 
                           population, " - ", covar_name, " - ", exp))
            write.table(
              df,
              file = here::here("Metapone", "Input", population, covar_name,
                                paste0("mwas_", exp, "_", 
                                       population, "_", covar_name, ".txt")),
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE,
              sep = "\t"
            )
          })
      })
  })

message("Metapone input files created in Metapone/Input/")

# Run metapone for each exposure and population --------------------------------

message("Running metapone pathway analysis...")

system.time({
  names(combined_results_list_c18) |>
    purrr::set_names() |>
    purrr::map(function(population) {
      names(covar_list) |> 
        purrr::map(function(covar_name) {
          input_dir <- here::here("Metapone", "Input", population, covar_name)
          input_files <- list.files(input_dir, pattern = "\\.txt$",
                                    full.names = TRUE)
          input_files |>
            purrr::set_names(
              basename(input_files) |>
                stringr::str_remove("\\.txt$") |>
                stringr::str_remove(paste0("^mwas_")) |>
                stringr::str_remove(paste0("_", population, 
                                           "_", covar_name, "$"))
            ) |>
            purrr::imap(function(input_file, exp_name) {
              message(paste0("\n---Running metapone for: ", 
                             exp_name, " (", population, " - ", 
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
                                 " (", population, " - ", 
                                 covar_name, "): ", e$message))
                  return(NULL)
                }
              )
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
  purrr::imap(function(result_list, population) {
    # Save result tables to Excel
    result_tables <- result_list |>
      purrr::compact() |>
      purrr::map(~ .x$result_table)

    if (length(result_tables) > 0) {
      dir.create(here::here("tables", "metapone_results", population),
                 showWarnings = FALSE, recursive = TRUE)

      result_tables |>
        purrr::imap(function(tbl, exp_name) {
          writexl::write_xlsx(
            tbl,
            path = here::here(
              "tables", "metapone_results", population,
              glue::glue("metapone_{exp_name}_{population}.xlsx"))
          )
        })
    }

    # Save plots
    dir.create(here::here("figures", "metapone", population),
               showWarnings = FALSE, recursive = TRUE)

    result_list |>
      purrr::compact() |>
      purrr::imap(function(res, exp_name) {
        if (!is.null(res$plot)) {
          ggsave(
            filename = here::here(
              "figures", "metapone", population,
              glue::glue("metapone_{exp_name}_{population}.png")),
            plot = res$plot +
              ggtitle(paste0("metapone: ", exp_name)),
            width = 10, height = 8, dpi = 300
          )
        }
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


# =============================================================================
# SECTION 5: PATHWAY VISUALIZATION
# =============================================================================

# Function to create pathway enrichment plot ---------------------------------

create_pathway_barplot <- function(pathway_df, exposure_name, top_n = 20) {
  if (is.null(pathway_df) || nrow(pathway_df) == 0) return(NULL)

  plot_data <- pathway_df |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      pathway = forcats::fct_reorder(pathway, -P.Value),
      neg_log10_p = -log10(P.Value)
    )

  p <- ggplot(plot_data, aes(x = pathway, y = neg_log10_p)) +
    geom_bar(stat = "identity", fill = "#3B4CC0", alpha = 0.8) +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(
      title = paste0("Top Enriched Pathways - ", gsub("exp_", "", exposure_name)),
      x = "",
      y = expression(-log[10](P-value))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(face = "bold", size = 12)
    )

  return(p)
}


# Create pathway plots if results exist --------------------------------------

if (exists("pathway_results") && !is.null(pathway_results)) {
  dir.create(here::here("figures", "pathway"), showWarnings = FALSE, recursive = TRUE)

  purrr::iwalk(pathway_results, function(res, exp) {
    if (!is.null(res$mummichog) && nrow(res$mummichog) > 0) {
      p <- create_pathway_barplot(res$mummichog, exp)
      if (!is.null(p)) {
        ggsave(
          filename = here::here("figures", "pathway",
                                glue::glue("pathway_mummichog_{exp}.png")),
          plot = p,
          width = 10, height = 8, dpi = 300
        )
      }
    }
  })
}


message("\nPathway analysis completed!")
message("Results saved to:")
message("  - metaboAnalyst/Input/ (input files)")
message("  - metaboAnalyst/Output/ (analysis results)")
message("  - tables/mwas_results/ (summary tables)")
message("  - figures/pathway/ (visualizations)")

#--------------------------------End of the code--------------------------------
