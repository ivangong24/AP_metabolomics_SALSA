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
                          num_permutations = 100) {

  # --- 1. Run MetaboAnalystR in the output directory ---
  wd_orig <- getwd()
  on.exit(setwd(wd_orig), add = TRUE)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
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

  # --- 2. Extract result table ---
  result_table <- mSet$mummi.resmat |>
    as.data.frame() |>
    tibble::rownames_to_column("pathway") |>
    dplyr::arrange(P.Value)

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

  plot_data <- result_table |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      neg_log10_p = -log10(P.Value),
      enrichment_factor = Hits.sig / Hits.all
    ) |>
    dplyr::filter(is.finite(enrichment_factor), is.finite(neg_log10_p)) |>
    dplyr::arrange(dplyr::desc(neg_log10_p))

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
  top_indices <- plot_data$P.Value < p_threshold &
    plot_data$Hits.sig >= min_hits

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
