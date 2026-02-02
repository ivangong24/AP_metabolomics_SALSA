## ---------------------------
##
## Script name: 8-pathway_analysis.R
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

# Load required packages -----------------------------------------------------

library(tidyverse)
library(MetaboAnalystR)
library(writexl)

# Load MWAS results and annotation data --------------------------------------

load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))

# Load annotation files (m/z and retention time)
# These should contain: chemical_ID, mz, time (retention time)
annotation_c18 <- read_csv(
  here::here("data", "metabolomics", "annotation", "xmsannotator_c18neg.csv"),
  show_col_types = FALSE
) |>
  rename_all(str_to_lower)

annotation_hilic <- read_csv(
  here::here("data", "metabolomics", "annotation", "xmsannotator_hilicpos.csv"),
  show_col_types = FALSE
) |>
  rename_all(str_to_lower)

# Create output directories --------------------------------------------------

dir.create(here::here("metaboAnalyst", "Input"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("metaboAnalyst", "Output"), showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# SECTION 1: PREPARE MUMMICHOG INPUT FILES
# =============================================================================

# Function to create Mummichog input format ----------------------------------

create_mummichog_input <- function(mwas_result, annotation_df, mode) {
  # Mummichog input format:
  # m.z | rt | p.value | t.score | mode

  # Get chemical_ID column name (may vary)
  id_col <- intersect(
    c("chemical_id", "id", "met", "feature"),
    tolower(colnames(annotation_df))
  )[1]

  # Get m/z and retention time columns
  mz_col <- intersect(c("mz", "m.z", "mass"), tolower(colnames(annotation_df)))[1]
  rt_col <- intersect(c("time", "rt", "retention_time"), tolower(colnames(annotation_df)))[1]

  mwas_result |>
    tibble::rownames_to_column("met") |>
    dplyr::left_join(
      annotation_df |>
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

exposure_vars <- names(mwas_results_c18)

## Create C18 input files (negative mode)
mummichog_input_c18 <- exposure_vars |>
  purrr::set_names() |>
  purrr::map(function(exp) {
    create_mummichog_input(
      mwas_results_c18[[exp]],
      annotation_c18,
      mode = "negative"
    )
  })

## Create HILIC input files (positive mode)
mummichog_input_hilic <- exposure_vars |>
  purrr::set_names() |>
  purrr::map(function(exp) {
    create_mummichog_input(
      mwas_results_hilic[[exp]],
      annotation_hilic,
      mode = "positive"
    )
  })


# Combine C18 and HILIC for each exposure ------------------------------------

mummichog_input_combined <- exposure_vars |>
  purrr::set_names() |>
  purrr::map(function(exp) {
    dplyr::bind_rows(
      mummichog_input_c18[[exp]],
      mummichog_input_hilic[[exp]]
    ) |>
      dplyr::arrange(`p.value`)
  })


# Write input files ----------------------------------------------------------

purrr::iwalk(mummichog_input_combined, function(df, exp) {
  write.table(
    df,
    file = here::here("metaboAnalyst", "Input",
                      paste0("mwas_", exp, "_combined.txt")),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )
})

message("Mummichog input files created in metaboAnalyst/Input/")


# =============================================================================
# SECTION 2: RUN MUMMICHOG PATHWAY ANALYSIS
# =============================================================================

# Function to run Mummichog using MetaboAnalystR -----------------------------

run_mummichog <- function(input_file, output_dir, exposure_name,
                           p_cutoff = 0.05,
                           organism = "hsa") {
  # Initialize MetaboAnalystR object
  mSet <- MetaboAnalystR::InitDataObjects("mass_table", "mummichog", FALSE)

  # Set peak format
  mSet <- MetaboAnalystR::SetPeakFormat(mSet, "mprt")

  # Read peak data
  mSet <- MetaboAnalystR::Read.PeakListData(mSet, input_file)

  # Set analysis parameters
  mSet <- MetaboAnalystR::UpdateMummichogParameters(
    mSet,
    "0.05",  # p-value cutoff
    "5"      # minimum number of hits
  )

  # Set organism
  mSet <- MetaboAnalystR::SetOrganism(mSet, organism)

  # Run Mummichog
  mSet <- MetaboAnalystR::PerformMummichog(mSet, "hsa_mfn")

  # Run GSEA
  mSet <- MetaboAnalystR::PerformGSEA(mSet, "hsa_mfn")

  # Save results
  save(mSet, file = file.path(output_dir, paste0("mummichog_", exposure_name, ".RData")))

  return(mSet)
}


# Run Mummichog for each exposure --------------------------------------------

# Create output directories for each exposure
purrr::walk(exposure_vars, function(exp) {
  dir.create(
    here::here("metaboAnalyst", "Output", exp),
    showWarnings = FALSE,
    recursive = TRUE
  )
})

# Run analysis (this may take some time)
message("Running Mummichog pathway analysis...")

mummichog_results <- tryCatch({
  exposure_vars |>
    purrr::set_names() |>
    purrr::map(function(exp) {
      message(paste0("Processing: ", exp))
      run_mummichog(
        input_file = here::here("metaboAnalyst", "Input",
                                paste0("mwas_", exp, "_combined.txt")),
        output_dir = here::here("metaboAnalyst", "Output", exp),
        exposure_name = exp,
        p_cutoff = 0.05,
        organism = "hsa"
      )
    })
}, error = function(e) {
  message("Mummichog analysis failed. Error: ", e$message)
  message("Please run pathway analysis manually using MetaboAnalyst web interface.")
  return(NULL)
})


# =============================================================================
# SECTION 3: EXTRACT AND SUMMARIZE PATHWAY RESULTS
# =============================================================================

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


# =============================================================================
# SECTION 4: CREATE SUMMARY TABLES
# =============================================================================

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


# =============================================================================
# SECTION 5: ALTERNATIVE - MANUAL METABOANALYST INPUT
# =============================================================================

# If MetaboAnalystR fails, create input for web interface --------------------

# The input files created in Section 1 can be uploaded to:
# https://www.metaboanalyst.ca/MetaboAnalyst/ModuleView.xhtml

message("\nAlternative: Manual pathway analysis")
message("If MetaboAnalystR fails, upload input files to MetaboAnalyst web interface:")
message("https://www.metaboanalyst.ca/")
message("Input files are located in: metaboAnalyst/Input/")


# =============================================================================
# SECTION 6: PATHWAY VISUALIZATION
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
