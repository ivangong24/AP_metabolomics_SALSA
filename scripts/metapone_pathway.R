## ---------------------------
##
## Script name: metapone_pathway.R
## Purpose of script: Function to run pathway analysis using metapone
##
## Author: Yufan Gong
##
## Date Created: 2026-03-10
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: metapone performs metabolomics pathway analysis using
##        a permutation-based approach on m/z features.
##        It does not require prior metabolite identification.
##
##        Input: Combined feature table (same format as Mummichog input)
##               with columns: m.z, rt, p.value, t.score, mode
##        Output: metapone result object + bubble plot
##
##        Dependencies: BiocManager::install("metapone")
## ---------------------------

library(metapone)
library(ggrepel)

#' Run metapone pathway analysis from a combined feature table file
#'
#' @param input_file Path to the tab-delimited input file. Must contain
#'   columns: m.z, rt, p.value, t.score, mode (where mode is "negative"
#'   or "positive"). This is the same combined C18 + HILIC feature table
#'   used for Mummichog.
#' @param p_cutoff Numeric, p-value cutoff for significant features
#'   (default 0.05).
#' @param num_permutations Integer, number of permutations (default 200).
#' @param match_tol_ppm Numeric, mass tolerance in ppm (default 10).
#'
#' @return A list containing:
#'   - result: metapone result object (metaponeResult S4 class)
#'   - result_table: data frame of pathway results
#'   - plot: ggplot bubble plot of enriched pathways
run_metapone <- function(input_file,
                         p_cutoff = 0.05,
                         num_permutations = 200,
                         match_tol_ppm = 10,
                         pos.adductlist = c("M+H", "M+Na", "M+"),
                         neg.adductlist = c("M-H", "M-2H", "M-H2O-H")
                         ) {

  # --- 1. Read and split by ionization mode ---
  feat_df <- read.delim(input_file, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)

  # metapone expects a list of matrices with columns: mz, time, p-value, t-stat
  # Split by mode and convert to matrices
  dat_list <- list()
  type_list <- list()

  mode_map <- c("negative" = "neg", "positive" = "pos")

  for (m in intersect(names(mode_map), unique(feat_df$mode))) {
    mat <- feat_df |>
      dplyr::filter(mode == m) |>
      dplyr::transmute(
        mz = m.z,
        time = rt,
        `p-value` = p.value,
        `t-statistics` = t.score
      ) |>
      dplyr::filter(!is.na(mz)) |>
      as.matrix()

    if (nrow(mat) > 0) {
      dat_list <- c(dat_list, list(mat))
      type_list <- c(type_list, list(mode_map[[m]]))
      message(paste0("  ", m, " features: ", nrow(mat)))
    }
  }

  if (length(dat_list) == 0) {
    warning("No valid features found in input file: ", input_file)
    return(NULL)
  }

  # --- 2. Load pathway and compound databases ---
  data(hmdbCompMZ, package = "metapone", envir = environment())
  data(pa, package = "metapone", envir = environment())

  # --- 3. Run metapone ---
  message(paste0("  Running metapone (", num_permutations, " permutations)..."))

  result <- metapone::metapone(
    dat = dat_list,
    type = type_list,
    pa = pa,
    hmdbCompMZ = hmdbCompMZ,
    pos.adductlist = pos.adductlist, 
    neg.adductlist = neg.adductlist,
    p.threshold = p_cutoff,
    n.permu = num_permutations,
    fractional.count.power = 0.5,
    max.match.count=10, 
    use.fgsea = FALSE,
    match.tol.ppm = match_tol_ppm
  )

  # --- 4. Extract result table ---
  result_table <- result@test.result |>
    as.data.frame() |>
    tibble::rownames_to_column("pathway")

  # Convert numeric columns
  num_cols <- c("n_significant_metabolites", "n_mapped_metabolites",
                "n_metabolites", "lfdr", "adjust.p")
  existing_num_cols <- intersect(num_cols, colnames(result_table))
  if (length(existing_num_cols) > 0) {
    result_table <- result_table |>
      dplyr::mutate(
        dplyr::across(dplyr::all_of(existing_num_cols), as.numeric)
      )
  }

  # Sort by adjusted p-value
  if ("adjust.p" %in% colnames(result_table)) {
    result_table <- result_table |> dplyr::arrange(adjust.p)
  } else if ("lfdr" %in% colnames(result_table)) {
    result_table <- result_table |> dplyr::arrange(lfdr)
  }

  # --- 5. Create bubble plot ---
  plot <- create_metapone_plot(result_table)

  return(list(
    result = result,
    result_table = result_table,
    plot = plot
  ))
}


#' Create bubble plot for metapone results
#'
#' @param result_table Data frame from metapone result extraction.
#' @param p_threshold Numeric, p-value threshold for labeling (default 0.05).
#' @param top_n Integer, maximum number of pathways to display (default 30).
#'
#' @return A ggplot object.
create_metapone_plot <- function(result_table,
                                 p_threshold = 0.05,
                                 top_n = 30) {

  if (is.null(result_table) || nrow(result_table) == 0) return(NULL)

  # Identify the p-value column
  p_col <- intersect(c("adjust.p", "lfdr"), colnames(result_table))[1]
  if (is.na(p_col)) return(NULL)

  # Identify count columns
  n_sig_col <- intersect(
    c("n_significant_metabolites", "n_significant"), colnames(result_table)
  )[1]
  n_map_col <- intersect(
    c("n_mapped_metabolites", "n_mapped"), colnames(result_table)
  )[1]

  if (is.na(n_sig_col) || is.na(n_map_col)) return(NULL)

  plot_data <- result_table |>
    dplyr::mutate(
      p_val = as.numeric(.data[[p_col]]),
      n_sig = as.numeric(.data[[n_sig_col]]),
      n_map = as.numeric(.data[[n_map_col]])
    ) |>
    dplyr::filter(!is.na(p_val), p_val > 0, is.finite(p_val)) |>
    dplyr::arrange(p_val) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      neg_log10_p = -log10(p_val),
      enrichment_factor = n_sig / n_map
    ) |>
    dplyr::filter(is.finite(enrichment_factor))

  if (nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = enrichment_factor, y = neg_log10_p)) +
    geom_point(aes(size = n_sig, color = neg_log10_p),
               stroke = 0.5) +
    scale_size_continuous(range = c(1, 6), name = "Significant hits") +
    scale_color_gradient(low = "yellow", high = "red",
                         name = expression(-log[10](p))) +
    xlab("Enrichment Factor") +
    ylab(expression(-log[10](p))) +
    theme_minimal() +
    theme(legend.position = "right")

  # Label significant pathways
  top_indices <- plot_data$p_val < p_threshold & plot_data$n_sig >= 2

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
