## ---------------------------
##
## Script name: mummichog_pathway.R
## Purpose of script: Function to run pathway analysis using Mummichog
##                    (via MetaboAnalystR)
##
## Author: Yufan Gong
##
## Date Created: 2026-01-29
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: Mummichog performs metabolomics pathway analysis using
##        a probabilistic approach on m/z features.
##        It does not require prior metabolite identification.
##
##        Input: Combined feature table with columns:
##               m.z, rt, p.value, t.score, mode
##        Output: mSet result object + bubble plot
##
##        Dependencies: MetaboAnalystR, ggrepel
##
##        2026-08-26: result extraction and the plot helper now resolve column
##        names by inspection instead of assuming `P.Value` / `Hits.all`.
##        MetaboAnalystR 4.3.0 returns FET / EASE / Gamma in `mummi.resmat`, so
##        the previous hard-coded `arrange(P.Value)` errored after the
##        permutations had already run. PerformPSEA also writes
##        `mummichog_pathway_enrichment_*.csv` with richer, better-named columns
##        (P(Fisher), P(EASE), P(Gamma), adjusted versions, Hits.sig,
##        Hits.total, cpd.hits), so that file is preferred when present.
##        Ported from the working version in ~/github/sleep_metabolomics_salsa.
## ---------------------------

library(MetaboAnalystR)
library(ggrepel)

#' Run Mummichog pathway analysis from a combined feature table file
#'
#' @param input_file Path to the tab-delimited input file. Must contain
#'   columns: m.z, rt, p.value, t.score, mode (where mode is "negative"
#'   or "positive"). This is the same combined C18 + HILIC feature table.
#' @param output_dir Path to the output directory for MetaboAnalystR results.
#'   MetaboAnalystR writes files to the working directory, so this function
#'   temporarily sets the working directory to output_dir.
#' @param p_cutoff Numeric, p-value cutoff for significant features
#'   (default 0.1).
#' @param organism Character, organism code for pathway database
#'   (default "hsa_mfn").
#' @param instrument_ppm Numeric, mass tolerance in ppm (default 10.0).
#' @param ion_mode Character, ionization mode: "mixed", "positive", or
#'   "negative" (default "mixed").
#' @param adducts Character vector of adduct types to consider.
#' @param min_hits Integer, minimum pathway size (default 3).
#' @param num_permutations Integer, number of permutations (default 100).
#' @param lib_cache Optional path to a directory holding the MetaboAnalystR
#'   reference libraries. MetaboAnalystR caches them in the WORKING directory,
#'   and this function gives every analysis its own working directory, so
#'   without a shared cache each run re-downloads them. See the note below.
#'
#' @return A list containing:
#'   - mSet: MetaboAnalystR mSet object
#'   - result_table: data frame of pathway results
#'   - plot: ggplot bubble plot of enriched pathways
run_mummichog <- function(input_file,
                          output_dir,
                          p_cutoff = 0.1,
                          organism = "hsa_mfn",
                          instrument_ppm = 10.0,
                          ion_mode = "mixed",
                          adducts = c("M-H [1-]", "M-2H [2-]",
                                      "M-H2O-H [1-]", "M [1+]",
                                      "M+H [1+]", "M+Na [1+]"),
                          min_hits = 3,
                          num_permutations = 100,
                          lib_cache = NULL) {

  # --- 1. Run MetaboAnalystR in the output directory ---
  wd_orig <- getwd()
  on.exit(setwd(wd_orig), add = TRUE)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ## REFERENCE LIBRARIES: seed the output directory from a shared cache.
  ##
  ## MetaboAnalystR's .get.my.lib() looks for the library file in the CURRENT
  ## WORKING DIRECTORY and downloads it from metaboanalyst.ca only when it is
  ## absent. Because every analysis is given its own working directory (it has
  ## to be -- PerformPSEA writes a dozen result files with fixed names), each
  ## one starts with an empty cache and re-downloads the same files. Across the
  ## 236-cell revision grid that was 237 copies of hsa_mfn.qs, one per
  ## directory, and it dominated the runtime of the whole Mummichog stage.
  ##
  ## Copying rather than symlinking is deliberate: MetaboAnalystR treats these
  ## as its own working files, and a symlink would let it write back through to
  ## the shared cache. They are under a megabyte, so a local copy costs
  ## milliseconds against a network round trip.
  lib_files <- c(paste0(organism, ".qs"),
                 paste0(ion_mode, "_adduct.qs"))

  if (!is.null(lib_cache)) {
    dir.create(lib_cache, showWarnings = FALSE, recursive = TRUE)
    for (f in lib_files) {
      from <- file.path(lib_cache, f)
      to   <- file.path(output_dir, f)
      if (file.exists(from) && !file.exists(to)) {
        file.copy(from, to, overwrite = FALSE)
      }
    }
  }

  setwd(output_dir)

  mSet <- InitDataObjects("mass_all", "mummichog", FALSE, default.dpi = 300)
  mSet <- SetPeakFormat(mSet, "rmp")
  mSet <- UpdateInstrumentParameters(mSet, instrument_ppm, ion_mode, "yes", 0.02)
  mSet <- Read.PeakListData(mSet, input_file)
  mSet <- SanityCheckMummichogData(mSet)
  mSet <- Setup.AdductData(mSet, adducts)
  mSet <- PerformAdductMapping(mSet, ion_mode)
  mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
  mSet <- SetMummichogPval(mSet, p_cutoff)
  mSet <- PerformPSEA(mSet, organism, "current", min_hits, num_permutations)
  mSet <- PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 300, width = 10)

  setwd(wd_orig)

  ## Populate the cache from whatever this run had to fetch, so the next cell
  ## finds it locally. The first analysis of a session pays the download once;
  ## every later one copies.
  if (!is.null(lib_cache)) {
    for (f in lib_files) {
      from <- file.path(output_dir, f)
      to   <- file.path(lib_cache, f)
      if (file.exists(from) && !file.exists(to)) {
        file.copy(from, to, overwrite = FALSE)
      }
    }
  }

  # --- 2. Read the enrichment result table. PerformPSEA writes
  #        `mummichog_pathway_enrichment_mummichog.csv` with the full, nicely
  #        named columns -- P(Fisher), P(EASE), P(Gamma) and their adjusted
  #        versions, Hits.sig / Hits.total, cpd.hits -- which is richer than the
  #        in-memory resmat (FET/EASE/Gamma). Rank by the Fisher p-value. ---
  csv_file <- list.files(
    output_dir, pattern = "^mummichog_pathway_enrichment.*\\.csv$",
    full.names = TRUE)
  if (length(csv_file) > 0) {
    result_table <- readr::read_csv(csv_file[1], show_col_types = FALSE) |>
      dplyr::rename(pathway = 1)
  } else {
    result_table <- mSet$mummi.resmat |>
      as.data.frame() |>
      tibble::rownames_to_column("pathway")
  }
  p_col <- intersect(c("P(Fisher)", "FET", "P(Gamma)", "Gamma"),
                     names(result_table))[1]
  if (!is.na(p_col)) {
    result_table <- dplyr::arrange(result_table, .data[[p_col]])
  }

  # --- 3. Create bubble plot ---
  plot <- create_mummichog_plot(result_table)

  return(list(
    mSet = mSet,
    result_table = result_table,
    plot = plot
  ))
}


#' Create bubble plot for Mummichog results
#'
#' @param result_table Data frame from Mummichog result extraction.
#'   Expected columns: pathway, Hits.sig, Hits.all, P.Value, and the
#'   enrichment column (column 8 in mummi.resmat).
#' @param p_threshold Numeric, p-value threshold for labeling (default 0.05).
#' @param min_hits Integer, minimum significant hits for labeling (default 3).
#' @param top_n Integer, maximum number of pathways to display (default 30).
#'
#' @return A ggplot object.
create_mummichog_plot <- function(result_table,
                                  p_threshold = 0.05,
                                  min_hits = 3,
                                  top_n = 30) {

  if (is.null(result_table) || nrow(result_table) == 0) return(NULL)

  # Column names differ between the PerformPSEA CSV ("P(Fisher)", "Hits.sig",
  # "Hits.total") and the in-memory resmat ("FET", "Hits.sig", "Hits.all"),
  # and between MetaboAnalystR versions. Resolve them by inspection.
  p_col   <- intersect(c("P(Fisher)", "FET", "P(Gamma)", "Gamma"),
                       names(result_table))[1]
  hit_col <- intersect(c("Hits.sig", "Hits_sig"), names(result_table))[1]
  tot_col <- intersect(c("Hits.total", "Pathway.total", "Hits.all"),
                       names(result_table))[1]
  if (is.na(p_col) || is.na(hit_col) || is.na(tot_col)) return(NULL)

  plot_data <- result_table |>
    dplyr::mutate(
      p_val = as.numeric(.data[[p_col]]),
      n_hit = as.numeric(.data[[hit_col]]),
      n_tot = as.numeric(.data[[tot_col]])
    ) |>
    dplyr::filter(p_val > 0, is.finite(p_val)) |>
    dplyr::arrange(p_val) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(neg_log10_p = -log10(p_val),
                  enrichment_factor = n_hit / n_tot) |>
    dplyr::filter(is.finite(enrichment_factor), is.finite(neg_log10_p))


  if (nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = enrichment_factor, y = neg_log10_p)) +
    geom_point(aes(size = sqrt(abs(enrichment_factor)),
                   color = neg_log10_p),
               stroke = 0.5) +
    scale_size_continuous(range = c(1, 6)) +
    scale_color_gradient(low = "yellow", high = "red",
                         name = expression(-log[10](p))) +
    xlab("Enrichment Factor") +
    ylab(expression(-log[10](p))) +
    theme_minimal() +
    theme(legend.position = "none")

  # Label significant pathways
  top_indices <- plot_data$p_val < p_threshold &
    plot_data$n_hit >= min_hits

  if (any(top_indices)) {
    p <- p +
      geom_text_repel(
        aes(label = pathway),
        data = plot_data[top_indices, ],
        size = 3, max.overlaps = 15
      )
  }

  return(p)
}
