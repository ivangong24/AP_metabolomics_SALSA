## ---------------------------
##
## Script name: R7-visualization_revision.R
## Purpose of script: To create visualizations for the R1 revision MWAS and
##                    pathway results
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
## Notes: This is the revision analogue of scripts/7-visualization.R.
##
##        Sections:
##          1. Volcano plots (combined C18 + HILIC)
##          2. Manhattan and VIP Manhattan plots
##          3. Unsupervised versus cross-fitted effect-size scatter plots
##             (the figure that answers "do the findings hold?")
##          4. Revision versus submitted effect-size scatter plots
##          5. Heatmaps of the top features across exposures
##          6. Pathway enrichment bubble plots
##          7A. Biological panels for the combined figure (Reviewer 1
##              comment 14): single-pollutant effect comparison, effect
##              estimates with 95% CIs for the confirmed (Level 1)
##              identifications, and exposure-response plots for those
##              metabolites
##          7. Combined main-text panel and the SI coefficient scatters
##
##        Figure tree mirrors the primary pipeline:
##        revision_output/figures/mwas/{study}/{population}/{covar_set}/{exposure}/
##
## ---------------------------

# Load required packages -----------------------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

pacman::p_load("ggh4x", "ggtext")

rev_announce("R7-visualization_revision.R")

# Load results ---------------------------------------------------------------

load(rev_here("data", "metabolomics", "results",
              "mwas_results_all_revision.RData"))

load(rev_here("data", "metabolomics", "results",
              "mwas_annotation_revision.RData"))

load(rev_here("data", "processed", "combined_data_list_revision.RData"))

## In QUICK mode cut the figure grid down so the script can be validated end
## to end in minutes. Nothing produced this way is reportable.
if (REV_QUICK) {
  ## Two composites plus one single pollutant -- see the same trim in R6.
  exposure_vars_list <- exposure_vars_list |>
    purrr::map(~ c(head(.x[!is_single_pollutant(.x)], 2),
                   head(.x[is_single_pollutant(.x)], 1)))
}

study_names      <- names(mwas_results_list_c18)
## Populations are per-study, not global: the post-diagnosis sensitivity
## population `all predx` is built for `total` only (see R4), so indexing the
## cox arm with the total arm's population list would yield NULL results.
populations_for <- function(study) {
  pops <- names(mwas_results_list_c18[[study]])
  ## QUICK keeps `all predx` too -- it is the newest slice of the grid and the
  ## one a smoke test most needs to exercise. It exists in `total` only.
  if (REV_QUICK) intersect(c("all", "all predx"), pops) else pops
}
population_names <- populations_for("total")
covar_names      <- if (REV_QUICK) "covar" else names(covar_list)

message("Studies: ",     paste(study_names, collapse = ", "))
message("Populations: ", paste(population_names, collapse = ", "))
message("Covariate sets: ", paste(covar_names, collapse = ", "))


# Create output directories --------------------------------------------------

study_names |>
  purrr::walk(function(study){
    populations_for(study) |>
      purrr::walk(function(population){
        covar_names |>
          purrr::walk(function(covar_name){
            exposures_for(study, population) |>
              purrr::walk(function(exposure){
                rev_dir("figures", "mwas", study, population,
                        covar_name, exposure)
              })
          })
      })
  })


## ggsave wrapper that applies REV_TEXT. patchwork objects take a theme the
## same way a ggplot does, so one wrapper covers both.
rev_ggsave <- function(filename, plot, ...) {
  ggsave(filename = filename, plot = plot & REV_TEXT, ...)
}


## The shared significance palette, SIG_LEVELS and sig_class() now live in
## R1-revision_functions.R, so that R11's meta-analysis panel classifies and
## colours its pooled estimates by exactly the same rule as the MWAS figures.


## Identification confidence --------------------------------------------------
##
## annotation_confidence_level(), n_candidate_compounds() and
## candidate_display_name() live in R1-revision_functions.R, because R5 uses
## them to write the confidence_level and n_candidates columns into the
## annotated result tables and this figure has to agree with those tables.
##
## Only Level 1 identifications get a name printed. Reviewer 1 asks that
## implausible candidate names be removed; on these data Level 1 is the
## in-house library matches that resolve to a single compound, and nothing
## else reaches it -- Level 2 requires an MS/MS spectral match and this study
## acquired MS1 full scan only.
VOLCANO_LABEL_MAX_LEVEL <- 1L


# =============================================================================
# SECTION 1: VOLCANO PLOTS (Combined C18 + HILIC)
# =============================================================================

# Function to create combined volcano plot ------------------------------------

create_volcano_plot <- function(mwas_result_c18, mwas_result_hilic,
                                annotation_result_c18, annotation_result_hilic,
                                exposure_name,
                                fdr_threshold = 0.05,
                                n_labels = 10,
                                label_size = 4.5) {

  # Helper to prepare data for one column type
  prep_data <- function(mwas_result, annotation_result, column_type) {
    mwas_result |>
      tibble::rownames_to_column("met") |>
      dplyr::left_join(
        annotation_result |>
          dplyr::select(met, chemical_id, compound, multiple_match, reference,
                        dplyr::any_of("VIP_comp1")),
        by = "met"
      ) |>
      dplyr::group_by(met) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::mutate(column_type = column_type)
  }

  # Combine C18 and HILIC data
  plot_data <- dplyr::bind_rows(
    prep_data(mwas_result_c18, annotation_result_c18, "C18/neg-"),
    prep_data(mwas_result_hilic, annotation_result_hilic, "HILIC/pos+")
  ) |>
    dplyr::mutate(
      neg_log10_p = -log10(P.Value),
      ## Reviewer 1 minor comment 14 asks for legends distinguishing nominal,
      ## nominal, suggestive and FDR-significant features. All three come
      ## from sig_class() so the volcano, the Manhattan and the effect-estimate
      ## panel cannot disagree about what a colour means.
      significant = sig_class(P.Value, adj.P.Val)
    )

  # Get top metabolites for labeling (by p-value)
  ##
  ## Restricted to identifications at VOLCANO_LABEL_MAX_LEVEL or better. An
  ## accurate-mass database hit supplies a plausible-sounding name with no
  ## evidence behind it, and over half the annotated features here match more
  ## than one compound (multiple_match is TRUE for 1,494 of 2,814 C18 and
  ## 1,464 of 2,815 HILIC annotations), so labelling by p-value alone printed
  ## names the data cannot support. See annotation_confidence_level().
  plot_data <- plot_data |>
    dplyr::mutate(conf_level = annotation_confidence_level(reference, compound))

  labelable <- plot_data |>
    dplyr::filter(!is.na(compound), compound != "",
                  !is.na(conf_level), conf_level <= VOLCANO_LABEL_MAX_LEVEL)

  ## Two rules, unioned.
  ##
  ## The cap is a display limit, so on its own it can silently drop a genuine
  ## finding: across the grid it left exactly one identification at FDR < 0.10
  ## unnamed (Indole under comp_pca_all_w10, p = 0.003, q = 0.093 -- eleventh
  ## by raw p because ten other Level 1 features sat above it). An FDR-
  ## significant identification is the one thing this panel must never omit,
  ## so those are named unconditionally and the cap only governs the rest.
  ##
  ## The gate on the rest is P < 0.05, NOT the FDR tiers. Gating labels on FDR
  ## would leave 41 of the 60 cells with no named compound at all -- Level 1 is
  ## 119 of ~20,125 features and is not enriched at the top of the p-value
  ## distribution, so Level 1 AND FDR is a double filter that usually empties.
  ## Three of the four main-text cells would go completely unlabelled,
  ## including 1,3-butadiene, whose 34 FDR-significant features all happen to be
  ## Level 3. The label answers "which compound is this point"; the COLOUR
  ## answers how strong it is, and since 2026-09-02 those colours are two
  ## calibrated FDR tiers, so nothing is oversold by naming a nominal hit.
  always_label <- labelable |>
    dplyr::filter(adj.P.Val < 0.10) |>
    dplyr::pull(met)

  top_mets <- labelable |>
    dplyr::filter(P.Value < 0.05) |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = n_labels) |>
    dplyr::pull(met) |>
    union(x = always_label)

  plot_data <- plot_data |>
    dplyr::mutate(
      ## Not str_split_i(compound, ";", 1): taking the first of several
      ## candidates prints whichever name happens to sort first, which is how
      ## "Pectin" ended up labelling a pentose feature. Only features that
      ## collapse to a single compound are labelled at all (see
      ## annotation_confidence_level), and this prints that compound.
      compound_display = candidate_display_name(compound),
      label = ifelse(met %in% top_mets, compound_display, "")
    )

  # Create volcano plot
  ggplot(plot_data, aes(x = logFC, y = neg_log10_p)) +
    geom_point(
      data = . %>% filter(significant == "NS"),
      aes(shape = column_type),
      color = "grey70", alpha = 0.5, size = 1.5
    ) +
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant, shape = column_type),
      alpha = 0.7, size = 2.5
    ) +
    scale_color_manual(values = SIG_COLORS, name = "Significance",
                       drop = FALSE) +
    scale_shape_manual(
      values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
      name = "Column"
    ) +
    ## Only the labelled features are handed to ggrepel. Passing the full
    ## ~20,000-point frame with label = "" makes every point compete for
    ## space, and max.overlaps then silently drops labels in exactly the
    ## dense high-significance region the labels are meant to mark (3 of 10
    ## survived before this change). max.overlaps = Inf so a requested label
    ## is never dropped without us knowing.
    ## max.overlaps = Inf guarantees a requested label is never silently
    ## dropped, but it also means ggrepel will place labels on top of one
    ## another rather than give up, which is what caused the overlap. The fix
    ## is to give the solver room and time instead of letting it drop labels:
    ## a stronger repulsive force, many more iterations, and axis expansion
    ## below so there is white space to push into. seed makes the result
    ## reproducible across re-runs.
    geom_label_repel(
      data = ~ dplyr::filter(.x, label != ""),
      aes(label = label),
      size = label_size, max.overlaps = Inf,
      box.padding = 1.1, point.padding = 0.4,
      label.padding = 0.18, label.size = 0.25,
      force = 12, force_pull = 0.4,
      max.iter = 50000, max.time = 2,
      min.segment.length = 0, segment.color = "grey50",
      segment.size = 0.3, seed = 42,
      show.legend = FALSE
    ) +
    ## Head-room for the repelled labels, which otherwise pile into the top of
    ## the panel where the most significant features sit.
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    scale_x_continuous(expand = expansion(mult = 0.14)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40") +
    labs(
      title = rev_label(exposure_name),
      x = expression(bold(Delta*" "*log[2]*" abundance per 1-SD increase in index")),
      y = expression(bold(-log[10](P-value)))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 13),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12)
    )
}


# Create and save volcano plots ----------------------------------------------

study_names |>
  purrr::walk(function(study){
    populations_for(study) |>
      purrr::walk(function(pop) {
        covar_names |>
          purrr::walk(function(cov) {
            exposures_for(study, pop) |>
              purrr::walk(function(exp) {
                message("Volcano: ", study, "_", pop, " [", cov, "] - ", exp)
                p <- create_volcano_plot(
                  mwas_result_c18 =
                    mwas_results_list_c18[[study]][[pop]][[cov]][[exp]],
                  mwas_result_hilic =
                    mwas_results_list_hilic[[study]][[pop]][[cov]][[exp]],
                  annotation_result_c18 =
                    mwas_full_annotated_list_c18[[study]][[pop]][[cov]][[exp]],
                  annotation_result_hilic =
                    mwas_full_annotated_list_hilic[[study]][[pop]][[cov]][[exp]],
                  exposure_name = exp
                )
                rev_ggsave(
                  filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                      glue::glue("volcano_{exp}.png")),
                  plot = p, width = 10, height = 6, dpi = 300, bg = "white"
                )
              })
          })
      })
  })


# =============================================================================
# SECTION 2: MANHATTAN AND VIP MANHATTAN PLOTS
# =============================================================================

# Function to create Manhattan plot (faceted by platform) --------------------

create_manhattan <- function(mwas_result_c18, mwas_result_hilic,
                             exposure_name) {

  prep_platform <- function(mwas_result, platform_label) {
    if (!"met" %in% colnames(mwas_result)) {
      mwas_result <- mwas_result |> tibble::rownames_to_column("met")
    }
    if (!"VIP_comp1" %in% colnames(mwas_result)) {
      mwas_result$VIP_comp1 <- NA_real_
    }

    mwas_result |>
      dplyr::mutate(
        # Parse m/z from feature name: mz_rt_[mz]_[rt]
        mz = as.numeric(stringr::str_extract(met, "(?<=mz_rt_)[0-9.]+")),
        neg_log10_p = -log10(P.Value),
        significant = sig_class(P.Value, adj.P.Val),
        platform = platform_label
      )
  }

  plot_data <- dplyr::bind_rows(
    prep_platform(mwas_result_c18, "C18/neg−"),
    prep_platform(mwas_result_hilic, "HILIC/pos+")
  )

  count_labels <- plot_data |>
    dplyr::group_by(platform) |>
    dplyr::summarize(
      ## Cumulative, matching the nesting of the classes: a feature at
      ## FDR < 0.05 is also at FDR < 0.10 and also at P < 0.05, and a count
      ## that said otherwise would read as three disjoint findings.
      n_fdr = sum(adj.P.Val < 0.05),
      n_fdr10 = sum(adj.P.Val < 0.10),
      n_nom = sum(P.Value < 0.05),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = paste0("FDR < 0.05: ", n_fdr,
                     "\nFDR < 0.10: ", n_fdr10,
                     "\nP < 0.05: ", n_nom),
      mz = Inf,
      neg_log10_p = Inf
    )

  ggplot(plot_data, aes(x = mz, y = neg_log10_p)) +
    geom_point(
      data = \(d) dplyr::filter(d, significant == "NS"),
      color = "grey70", alpha = 0.4, size = 1.2
    ) +
    geom_point(
      data = \(d) dplyr::filter(d, significant != "NS"),
      aes(color = significant), alpha = 0.7, size = 1.8
    ) +
    scale_color_manual(values = SIG_COLORS, name = "Significance",
                       drop = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_text(
      data = count_labels,
      aes(x = mz, y = neg_log10_p, label = label),
      hjust = 1.05, vjust = 1.2, size = 3.5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ platform, scales = "free_x") +
    labs(
      title = rev_label(exposure_name),
      x = expression(bold("Mass-to-charge ratio (" * italic(m/z) * ")")),
      y = expression(bold(-log[10](P-value)))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.background = element_rect(colour = "grey80", fill = "white",
                                       linewidth = 0.5),
      legend.margin = margin(4, 6, 4, 6),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.key.size = unit(0.5, "cm")
    )
}


# Function to create VIP Manhattan plot --------------------------------------

create_vip_manhattan <- function(combined_result_c18, combined_result_hilic,
                                 exposure_name, vip_threshold = 2) {

  prep_platform <- function(df, platform_label) {
    df |>
      dplyr::mutate(
        mz = as.numeric(stringr::str_extract(met, "(?<=mz_rt_)[0-9.]+")),
        selected = case_when(
          VIP_comp1 > vip_threshold & adj.P.Val < 0.05 ~ "VIP > 2 & FDR < 0.05",
          VIP_comp1 > vip_threshold ~ "VIP > 2",
          TRUE ~ "NS"
        ),
        selected = factor(selected,
                          levels = c("VIP > 2 & FDR < 0.05", "VIP > 2", "NS")),
        platform = platform_label
      )
  }

  plot_data <- dplyr::bind_rows(
    prep_platform(combined_result_c18, "C18/neg−"),
    prep_platform(combined_result_hilic, "HILIC/pos+")
  )

  ggplot(plot_data, aes(x = mz, y = VIP_comp1)) +
    geom_point(
      data = \(d) dplyr::filter(d, selected == "NS"),
      color = "grey70", alpha = 0.4, size = 1.2
    ) +
    geom_point(
      data = \(d) dplyr::filter(d, selected != "NS"),
      aes(color = selected), alpha = 0.75, size = 1.8
    ) +
    scale_color_manual(
      values = c("VIP > 2 & FDR < 0.05" = "#B73F42", "VIP > 2" = "#436C85"),
      name = NULL, drop = FALSE
    ) +
    geom_hline(yintercept = vip_threshold, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    facet_wrap(~ platform, scales = "free_x") +
    labs(
      title = rev_label(exposure_name),
      x = expression(bold("Mass-to-charge ratio (" * italic(m/z) * ")")),
      y = "PLS VIP (component 1)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.position = "bottom"
    )
}


# Create and save Manhattan plots --------------------------------------------

study_names |>
  purrr::walk(function(study) {
    populations_for(study) |>
      purrr::walk(function(pop) {
        covar_names |>
          purrr::walk(function(cov) {
            exposures_for(study, pop) |>
              purrr::walk(function(exp) {
                message("Manhattan: ", study, "_", pop, " [", cov, "] - ", exp)

                rev_ggsave(
                  filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                      glue::glue("manhattan_{exp}.png")),
                  ## combined_results_list_*, not mwas_results_list_*: only
                  ## the former carries VIP_comp1. The Manhattan no longer
                  ## shades by VIP -- its middle band is FDR < 0.10 now, see
                  ## SIG_COLORS -- but the combined frames are still what the
                  ## rest of the grid is indexed on, so this stays as the
                  ## single input to both Manhattan variants.
                  plot = create_manhattan(
                    combined_results_list_c18[[study]][[pop]][[cov]][[exp]],
                    combined_results_list_hilic[[study]][[pop]][[cov]][[exp]],
                    exposure_name = exp),
                  width = 14, height = 6, dpi = 300, bg = "white"
                )

                ## comp_qgcomp_fw_all is a contrast and has no PLS model, so
                ## its VIP columns are all NA -- plotting them yields an axis
                ## of nothing. Skip rather than emit an empty figure.
                vip_c18 <- combined_results_list_c18[[study]][[pop]][[cov]][[exp]]
                if (!all(is.na(vip_c18$VIP_comp1))) {
                  rev_ggsave(
                    filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                        glue::glue("vip_manhattan_{exp}.png")),
                    plot = create_vip_manhattan(
                      vip_c18,
                      combined_results_list_hilic[[study]][[pop]][[cov]][[exp]],
                      exposure_name = exp),
                    width = 14, height = 6, dpi = 300, bg = "white"
                  )
                }
              })
          })
      })
  })


# =============================================================================
# SECTION 3: UNSUPERVISED VERSUS CROSS-FITTED EFFECT SIZES
# =============================================================================

## The figure that answers Reviewer 1's "show whether the findings hold". Each
## point is a feature; axes are the logFC from two exposure definitions in the
## same model, same sample. Points are restricted to features reaching P < 0.05
## under the reference (x-axis) exposure.

compare_logfc <- function(result_c18_x, result_hilic_x,
                          result_c18_y, result_hilic_y,
                          xlab, ylab, p_filter = 0.05) {

  prep <- function(res, platform) {
    res |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC, P.Value) |>
      dplyr::mutate(platform = platform)
  }

  dat_x <- dplyr::bind_rows(prep(result_c18_x, "C18/neg−"),
                            prep(result_hilic_x, "HILIC/pos+"))
  dat_y <- dplyr::bind_rows(prep(result_c18_y, "C18/neg−"),
                            prep(result_hilic_y, "HILIC/pos+"))

  plot_data <- dat_x |>
    dplyr::rename(logFC_x = logFC, p_x = P.Value) |>
    dplyr::inner_join(
      dat_y |> dplyr::rename(logFC_y = logFC, p_y = P.Value),
      by = c("met", "platform")
    ) |>
    dplyr::filter(p_x < p_filter)

  if (nrow(plot_data) < 3) return(NULL)

  stats_lab <- plot_data |>
    dplyr::group_by(platform) |>
    dplyr::summarise(
      r = stats::cor(logFC_x, logFC_y),
      n = dplyr::n(),
      concordant = mean(sign(logFC_x) == sign(logFC_y)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = paste0("r = ", round(r, 3),
                     "\nsame direction: ", round(100 * concordant), "%",
                     "\nn = ", n),
      logFC_x = -Inf, logFC_y = Inf
    )

  ggplot(plot_data, aes(x = logFC_x, y = logFC_y)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(alpha = 0.35, size = 1.3, colour = "#436C85") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#B73F42", linewidth = 0.7) +
    geom_text(data = stats_lab, aes(label = label),
              hjust = -0.1, vjust = 1.2, size = 3.4, fontface = "bold") +
    facet_wrap(~ platform, scales = "free") +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", colour = NA))
}


## Pairs to compare, by study
logfc_pairs <- list(
  total = list(
    c("comp_pca_all", "comp_qgcomp_fw_all"),
    c("comp_pca_all", "comp_wqs_cf_all"),
    c("comp_pca_all", "comp_qgcomp_cf_all"),
    ## Exposure-window sensitivity (Reviewer 1 comment 13). This is the
    ## "concordance of effect estimates between windows" the response letter
    ## promises: the same exposure, the same model, only the averaging window
    ## differs.
    c("comp_pca_all", "comp_pca_all_w1"),
    c("comp_pca_all", "comp_pca_all_w3"),
    c("comp_pca_all", "comp_pca_all_w10")
  ),
  ## The cox arm carries no unsupervised index (comp_pca_* is fitted in the
  ## total cohort only), so there is nothing to compare it against here. Its
  ## cross-fitted-versus-naive contrast comes from Section 4 instead.
  cox = list()
)

## Every other exposure against its own 5-year self, derived rather than listed
## so that adding a window in R1 reaches this figure. The pairing is the same
## contrast the three PCA pairs above draw -- one exposure, one model, the
## window as the only difference -- extended to the cross-fitted indices, the
## feature-wise contrast and the eight pollutants now that those carry their
## own 3- and 10-year arms.
##
## The pairs are filtered against exposures_for() at draw time, so a stratum
## that does not carry a window simply skips it.
logfc_pairs$total <- c(
  logfc_pairs$total,
  MWAS_WINDOWS |>
    purrr::map(function(w){
      windowed <- c(qgcomp_fw_exposure_for(w),
                    crossfit_exposure_for("comp_wqs_cf_", w),
                    crossfit_exposure_for("comp_qgcomp_cf_", w),
                    single_pollutant_exposures_for(w))
      purrr::map(windowed, ~ c(untagged_exposure(.x), .x))
    }) |>
    purrr::list_flatten()
)

message("logFC window-concordance pairs: ", length(logfc_pairs$total))

study_names |>
  purrr::walk(function(study){
    populations_for(study) |>
      purrr::walk(function(pop){
        covar_names |>
          purrr::walk(function(cov){
            logfc_pairs[[study]] |>
              purrr::walk(function(pair){
                x <- pair[1]; y <- pair[2]
                if (!all(c(x, y) %in% exposures_for(study, pop))) return(NULL)
                message("logFC comparison: ", study, "_", pop, " [", cov, "] ",
                        x, " vs ", y)

                p <- compare_logfc(
                  mwas_results_list_c18[[study]][[pop]][[cov]][[x]],
                  mwas_results_list_hilic[[study]][[pop]][[cov]][[x]],
                  mwas_results_list_c18[[study]][[pop]][[cov]][[y]],
                  mwas_results_list_hilic[[study]][[pop]][[cov]][[y]],
                  xlab = paste0("logFC - ", rev_label(x)),
                  ylab = paste0("logFC - ", rev_label(y))
                )
                if (is.null(p)) return(NULL)

                rev_dir("figures", "mwas", study, pop, cov, x)
                rev_ggsave(
                  filename = rev_here(
                    "figures", "mwas", study, pop, cov, x,
                    glue::glue("logfc_{x}_vs_{y}.png")),
                  plot = p + labs(
                    title = "Unsupervised versus cross-fitted exposure index",
                    subtitle = paste0(study, " / ", pop, " / ", cov,
                                      "; features with P < 0.05 under ",
                                      rev_label(x))),
                  width = 11, height = 5, dpi = 300, bg = "white"
                )
              })
          })
      })
  })


# =============================================================================
# SECTION 4: REVISION VERSUS SUBMITTED EFFECT SIZES
# =============================================================================

## How much does replacing the outcome-informed index change the MWAS? Each
## revision exposure is plotted against its naive counterpart from the
## submitted analysis (tables/mwas_results/... via data/metabolomics/results).

submitted_path <- here::here("data", "metabolomics", "results",
                             "mwas_results_all.RData")

if (file.exists(submitted_path)) {

  submitted <- new.env()
  load(submitted_path, envir = submitted)

  ## Revision exposure -> submitted counterpart
  ## The single pollutants are their own counterpart: the same IQR-scaled
  ## exposure was modelled in the submitted analysis, so this scatter isolates
  ## what the revision changed for them -- the covariate set (batch dropped)
  ## and the imputation -- with the exposure held fixed. For the composites the
  ## exposure changed too, so those panels confound the two.
  submitted_counterpart <- c(
    comp_pca_all           = "comp_wqs_all",
    comp_qgcomp_fw_all     = "comp_qgcomp_all",
    comp_wqs_cf_all        = "comp_wqs_all",
    comp_qgcomp_cf_all     = "comp_qgcomp_all",
    comp_qgcomp_cox_cf_all = "comp_qgcomp_cox_all",
    purrr::set_names(single_pollutant_exposures, single_pollutant_exposures)
  )

  study_names |>
    purrr::walk(function(study){
      covar_names |>
        purrr::walk(function(cov){
          exposures_for(study, "all") |>
            purrr::walk(function(exp){
              if (!exp %in% names(submitted_counterpart)) return(NULL)
              old <- submitted_counterpart[[exp]]
              old_c18 <-
                submitted$mwas_results_list_c18[[study]][["all"]][[cov]][[old]]
              old_hil <-
                submitted$mwas_results_list_hilic[[study]][["all"]][[cov]][[old]]
              if (is.null(old_c18) || is.null(old_hil)) return(NULL)

              message("Revision vs submitted: ", study, " [", cov, "] ",
                      exp, " vs ", old)

              p <- compare_logfc(
                old_c18, old_hil,
                mwas_results_list_c18[[study]][["all"]][[cov]][[exp]],
                mwas_results_list_hilic[[study]][["all"]][[cov]][[exp]],
                xlab = paste0("logFC - submitted ", old),
                ylab = paste0("logFC - revision ", exp)
              )
              if (is.null(p)) return(NULL)

              rev_ggsave(
                filename = rev_here(
                  "figures", "mwas", study, "all", cov, exp,
                  glue::glue("logfc_submitted_{old}_vs_revision_{exp}.png")),
                plot = p + labs(
                  title = "Submitted versus revision exposure definition",
                  subtitle = paste0(study, " / all / ", cov,
                                    "; features with P < 0.05 in the ",
                                    "submitted analysis")),
                width = 11, height = 5, dpi = 300, bg = "white"
              )
            })
        })
    })

} else {
  message("Submitted MWAS results not found at ", submitted_path,
          " - skipping Section 4.")
}


# =============================================================================
# SECTION 5: HEATMAP OF TOP METABOLITES
# =============================================================================

## NOTE ON SCALE: the columns of this heatmap are not on a common scale. The
## composites were divided by their SD in R4, so their logFC is per one-SD
## increase in the index; the single pollutants are IQR-scaled, so theirs is
## per interquartile range. SD/IQR runs from 0.73 (benzene) to 1.53
## (butadiene), so the shading compares SIGN and PATTERN across exposures, not
## magnitude. Any figure that needs comparable magnitudes has to rescale first.
##
## The columns also span three averaging windows now. Each is scaled within its
## own window -- the SD and the IQR are taken from that window's column in the
## study's `all` frame -- so a 3-year and a 10-year column are each "per one
## unit of themselves", not per a shared physical increment.
create_sig_heatmap <- function(mwas_results_exposure_list, column_type = "C18",
                               population = "all", covar_set = "covar",
                               study = "total", top_n = 50) {

  # Get top metabolites across all exposures
  all_sig <- mwas_results_exposure_list |>
    purrr::imap(function(df, exp) {
      df |>
        tibble::rownames_to_column("met") |>
        dplyr::filter(adj.P.Val < 0.1) |>
        dplyr::mutate(exposure = exp)
    }) |>
    purrr::list_rbind()

  if (nrow(all_sig) == 0) {
    message("No significant metabolites found for heatmap: ",
            study, "_", population, " [", covar_set, "] ", column_type)
    return(NULL)
  }

  top_mets <- all_sig |>
    dplyr::group_by(met) |>
    dplyr::summarize(min_p = min(adj.P.Val), .groups = "drop") |>
    dplyr::arrange(min_p) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(met)

  heatmap_data <- mwas_results_exposure_list |>
    purrr::imap(function(df, exp) {
      df |>
        tibble::rownames_to_column("met") |>
        dplyr::filter(met %in% top_mets) |>
        dplyr::select(met, logFC) |>
        dplyr::rename(!!exp := logFC)
    }) |>
    purrr::reduce(dplyr::full_join, by = "met") |>
    tibble::column_to_rownames("met") |>
    as.matrix()

  colnames(heatmap_data) <- gsub("comp_", "", colnames(heatmap_data))

  ## Drop features with a non-finite logFC anywhere, then bail out if nothing
  ## is left -- hclust() needs a complete numeric matrix.
  heatmap_data <- heatmap_data[
    apply(heatmap_data, 1, function(r) all(is.finite(r))), , drop = FALSE]

  if (nrow(heatmap_data) < 1) {
    message("No finite logFC values for heatmap: ",
            study, "_", population, " [", covar_set, "] ", column_type)
    return(NULL)
  }

  ## hclust() needs >= 2 objects on the axis it is clustering. A stratum can
  ## easily yield a single feature at FDR < 0.1 (the small demcind strata do),
  ## and a study arm can carry few exposures -- cluster only where it is
  ## possible rather than letting pheatmap abort the whole script.
  do_cluster_rows <- nrow(heatmap_data) >= 2
  do_cluster_cols <- ncol(heatmap_data) >= 2
  if (!do_cluster_rows || !do_cluster_cols) {
    message("Heatmap for ", study, "_", population, " [", covar_set, "] ",
            column_type, " is ", nrow(heatmap_data), " x ",
            ncol(heatmap_data), "; clustering disabled on the short axis.")
  }

  ## pheatmap builds its colour breaks with seq(min, max, length.out = 101) and
  ## then cut()s the values into them. If every cell holds the same value --
  ## which a 1 x 1 matrix always does, and a single-exposure arm can easily
  ## produce -- the breaks are not unique and cut() aborts the whole script.
  if (diff(range(heatmap_data)) == 0) {
    message("Heatmap for ", study, "_", population, " [", covar_set, "] ",
            column_type, " has a single distinct logFC value (",
            nrow(heatmap_data), " x ", ncol(heatmap_data),
            "); nothing to shade, skipping.")
    return(NULL)
  }

  pheatmap::pheatmap(
    heatmap_data,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    cluster_rows = do_cluster_rows,
    cluster_cols = do_cluster_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = paste0("Top ", top_n, " features - ", column_type,
                  " (", study, ", ", population, ", ", covar_set, ")"),
    filename = rev_here("figures", "mwas", study, population, covar_set,
                        glue::glue("heatmap_top{top_n}_{column_type}.png")),
    ## Width follows the exposure count. `total`/`all` carries 36 columns now
    ## that the 3- and 10-year arms are in, where it carried 15; at a fixed 10
    ## inches the column labels overlap into an unreadable band.
    width = max(10, 3.5 + 0.30 * ncol(heatmap_data)),
    height = 12
  )
}

study_names |>
  purrr::walk(function(study) {
    populations_for(study) |>
      purrr::walk(function(pop) {
        covar_names |>
          purrr::walk(function(cov) {
            create_sig_heatmap(mwas_results_list_c18[[study]][[pop]][[cov]],
                               "C18", population = pop, covar_set = cov,
                               study = study, top_n = 50)
            create_sig_heatmap(mwas_results_list_hilic[[study]][[pop]][[cov]],
                               "HILIC", population = pop, covar_set = cov,
                               study = study, top_n = 50)
          })
      })
  })


# =============================================================================
# SECTION 6: PATHWAY ENRICHMENT PLOTS
# =============================================================================

# Define pathway categories ---------------------------------------------------

amino_acid_metabolism <- c("Alanine", "Aspartate", "Asparagine",
                           "Arginine", "Histidine", "Lysine",
                           "Methionine", "Tryptophan", "Tyrosine",
                           "amino", "Dibasic", "Valine", "Glutamate",
                           "Glycine", "Serine", "Threonine", "Proline",
                           "Cysteine", "Phenylalanine", "Leucine",
                           "Glutathione")

carbohydrate_metabolism <- c("Fructose", "Galactose", "Starch", "Hexose",
                             "Blood", "Glycan", "Keratan", "Sialic",
                             "Hyaluronan", "Glucose", "Mannose", "Sucrose",
                             "Chondroitin", "Heparan")

lipid_metabolism <- c("Bile", "Fatty", "lipid", "Phytanic",
                      "Cholesterol", "neuroprostanes", "steroid",
                      "Sphingolipid", "Glycerophospholipid", "Phospholipid",
                      "Triglyceride", "Ceramide")

energy_metabolism <- c("Butanoate", "Carnitine", "Glycolysis",
                       "Pyruvate", "Pentose", "octadecatrienoate",
                       "TCA", "Citrate", "Oxidative", "Glyoxylate",
                       "Propanoate")

inflammation_metabolism <- c("Arachidonic", "Leukotriene", "Prostaglandin",
                             "linoleic", "Linoleate", "Eicosanoid")

vitamin_cofactor_metabolism <- c("Vitamin", "Biopterin", "Lipoate",
                                 "Porphyrin", "Catabolism", "Purine",
                                 "Caffeine", "Folate", "Riboflavin",
                                 "Thiamine", "Biotin", "Pantothenate")

signaling_metabolism <- c("Dynorphin", "Dopamine", "Serotonin",
                          "Catecholamine", "Neurotransmitter")

secondary_metabolite_metabolism <- c("Alkaloid")

xenobiotic_metabolism <- c("Xenobiotic", "Drug", "Benzoate")

## Ten possible classes, so a fixed palette rather than a Brewer scale
## (Set2 tops out at eight colours).
pathway_category_colours <- c(
  "amino acid"       = "#4C72B0",
  "carbohydrate"     = "#DD8452",
  "lipid"            = "#55A868",
  "energy"           = "#C44E52",
  "inflammation"     = "#8172B3",
  "vitamin/cofactor" = "#937860",
  "signaling"        = "#DA8BC3",
  "secondary"        = "#8C8C8C",
  "xenobiotic"       = "#CCB974",
  "other"            = "#BFBFBF"
)

categorize_pathway <- function(pathway_name) {
  case_when(
    str_detect(pathway_name, regex(str_c(amino_acid_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "amino acid",
    str_detect(pathway_name, regex(str_c(carbohydrate_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "carbohydrate",
    str_detect(pathway_name, regex(str_c(lipid_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "lipid",
    str_detect(pathway_name, regex(str_c(energy_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "energy",
    str_detect(pathway_name, regex(str_c(inflammation_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "inflammation",
    str_detect(pathway_name, regex(str_c(vitamin_cofactor_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "vitamin/cofactor",
    str_detect(pathway_name, regex(str_c(signaling_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "signaling",
    str_detect(pathway_name, regex(str_c(secondary_metabolite_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "secondary",
    str_detect(pathway_name, regex(str_c(xenobiotic_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "xenobiotic",
    TRUE ~ "other"
  )
}


# Read the pathway summary written by R6-pathway_revision.R -------------------

pathway_long_path <- rev_here("tables", "pathway", "pathway_results_long.xlsx")

if (file.exists(pathway_long_path)) {

  pathway_long <- readxl::read_xlsx(pathway_long_path) |>
    dplyr::mutate(category = categorize_pathway(pathway))

  ## Older R6 outputs predate the floored p column; derive it if absent.
  if (!"p_plot" %in% names(pathway_long)) {
    pathway_long <- pathway_long |>
      dplyr::group_by(algorithm, study, population, covar_set, exposure) |>
      dplyr::mutate(p_plot = pmax(p_value,
                                  min(c(p_value[p_value > 0], 1),
                                      na.rm = TRUE) / 2)) |>
      dplyr::ungroup()
  }

  ## Bubble plot: top pathways per exposure, primary covariate set,
  ## full cohort, both algorithms side by side.
  make_pathway_bubble <- function(data, study, algorithm, top_n = 15) {

    ## p_plot is p_value floored away from zero by R6 -- a permutation p of
    ## exactly 0 is the STRONGEST result, so filtering on p_value > 0 would
    ## drop precisely the pathways worth plotting.
    plot_data <- data |>
      dplyr::filter(study == !!study, algorithm == !!algorithm,
                    population == "all", covar_set == "covar",
                    !is.na(p_value), is.finite(p_plot), p_plot > 0) |>
      dplyr::group_by(exposure) |>
      dplyr::slice_min(p_value, n = top_n, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        neg_log10_p = -log10(p_plot),
        label = unname(rev_label(exposure)),
        pathway = stringr::str_trunc(pathway, 45)
      )

    if (nrow(plot_data) == 0) return(NULL)

    ## Size encodes the number of significant hits, not the p-value -- the
    ## x-axis already carries significance, and hit count is what separates a
    ## well-supported pathway from one resting on a single feature.
    if (all(is.na(plot_data$n_sig))) plot_data$n_sig <- 1

    ggplot(plot_data,
           aes(x = neg_log10_p, y = reorder(pathway, neg_log10_p))) +
      geom_point(aes(size = n_sig, colour = category), alpha = 0.85) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                 colour = "grey50") +
      scale_size_continuous(range = c(2, 7), name = "Significant hits") +
      scale_colour_manual(values = pathway_category_colours,
                          name = "Pathway class") +
      facet_wrap(~ label, scales = "free_y", ncol = 3) +
      labs(
        title = paste0(algorithm, " pathway enrichment - ", study,
                       " / all / covar"),
        subtitle = paste0("Top ", top_n,
                          " pathways per exposure. Dashed line: p = 0.05."),
        x = expression(-log[10](p)), y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 9),
            strip.background = element_rect(fill = "grey95", colour = NA),
            axis.text.y = element_text(size = 7))
  }

  tidyr::expand_grid(study = study_names,
                     algorithm = c("mummichog", "metapone")) |>
    purrr::pwalk(function(study, algorithm){
      p <- make_pathway_bubble(pathway_long, study, algorithm)
      if (is.null(p)) {
        message("No pathway rows for ", study, " / ", algorithm)
        return(NULL)
      }
      rev_save_plot(p, glue::glue("pathway_bubble_{algorithm}_{study}"),
                    "pathway", width = 14, height = 12)
    })

  # Pathway heatmap -----------------------------------------------------------
  ##
  ## Modelled on figures/pathway/pathway_heatmap_covar.png in the
  ## sleep_metabolomics_salsa project, redrawn on this project's palette.
  ##
  ## Two plots joined by patchwork: a left column that prints the pathway names
  ## on a class-coloured block, and a right heatmap whose squares are sized by
  ## the enrichment factor and filled by -log10 p, diverging at p = 0.05 on the
  ## same blue-white-red the MWAS figures use. Asterisks carry the SAME three
  ## classes as those figures (see sig_class): *** FDR < 0.05, ** FDR < 0.10,
  ## * p < 0.05.
  ##
  ## ONE FIGURE PER ALGORITHM. Mummichog and metapone test different pathway
  ## libraries -- 45 and 225 pathways respectively reach p < 0.05 in the full
  ## cohort -- so a shared row axis would be blank down one side and would
  ## imply a pathway had been tested and missed when it was never in the
  ## library at all.

  ## Rows: any pathway reaching p < 0.05 in ANY displayed cell. The cap is a
  ## safety limit on figure height, not a selection rule -- metapone tests 372
  ## pathways and 141 of them clear p < 0.05 somewhere in the window x
  ## population grid, which at one row per pathway is a poster rather than a
  ## figure. Where the cap binds the log says so, and the full table is
  ## pathway_results_long.xlsx.
  PATHWAY_MIN_COLS <- 1L
  PATHWAY_MAX_ROWS <- 45L

  ## Columns: the three exposure windows, and within each the cognitive strata.
  ##
  ## The 5-year index IS comp_pca_all -- the primary exposure -- so the middle
  ## column is the headline result and the outer two are its shorter and longer
  ## windows.
  ##
  ## All three windows carry all three strata (the window indices were widened
  ## to the cognitive strata on 2026-09-02, and the rest of the exposure set
  ## followed on 2026-09-03). `scales = "free_x"` and `space = "free_x"` are
  ## kept so that a cell missing for any other reason collapses rather than
  ## drawing an empty column.
  PATHWAY_WINDOW_COLUMNS <- tibble::tribble(
    ~exposure,          ~window,
    "comp_pca_all_w3",  "3-year",
    "comp_pca_all",     "5-year",
    "comp_pca_all_w10", "10-year")

  PATHWAY_POPULATION_LABELS <- c("all"        = "All",
                                 "no demcind" = "No dem/CIND",
                                 "demcind"    = "Dem/CIND")

  ## Eight hues chosen for separation from each other AND from the diverging
  ## fill, so the pollutant bars cannot be read as significance. Ordered so
  ## that neighbouring pollutants never take neighbouring hues.
  POLLUTANT_COLOURS <- c("#82B29B", "#DE476A", "#A8C3D1", "#A57E74",
                         "#7D518A", "#E29F34", "#3F3A39", "#E9B693")

  PATHWAY_CATEGORY_LABELS <- c(
    "amino acid" = "Amino acid", "carbohydrate" = "Carbohydrate",
    "lipid" = "Lipid", "energy" = "Energy", "inflammation" = "Inflammation",
    "vitamin/cofactor" = "Vitamin/cofactor", "signaling" = "Signalling",
    "secondary" = "Secondary", "xenobiotic" = "Xenobiotic", "other" = "Other")

  ## Class blocks are drawn from the SAME three colours the significance scale
  ## uses, at low alpha so they read as background rather than as data.
  rev_ramp <- function(n) {
    grDevices::colorRampPalette(
      c(SIG_COLORS[["FDR < 0.10"]], SIG_COLORS[["P < 0.05"]],
        SIG_COLORS[["FDR < 0.05"]]))(n)
  }

  rev_tints <- function(n, alpha = 0.22) scales::alpha(rev_ramp(n), alpha)

  ## Text drawn ON those blocks has to be darker than the ramp colour, not the
  ## ramp colour itself: at the tint's own hue the class name disappeared into
  ## its block.
  rev_tint_text <- function(n, mix = 0.5) {
    grDevices::rgb(t(grDevices::col2rgb(rev_ramp(n)) * mix),
                   maxColorValue = 255)
  }

  ## Enrichment factor: the pathway's own hit rate against the hit rate of the
  ## whole list that cell was tested on. n_sig / n_total alone is a proportion,
  ## not an enrichment, and would rank a large pathway with an average hit rate
  ## above a small one that is genuinely concentrated.
  add_enrichment <- function(df) {
    df |>
      dplyr::mutate(hit_rate = dplyr::if_else(n_total > 0, n_sig / n_total,
                                              NA_real_)) |>
      dplyr::group_by(exposure, population) |>
      dplyr::mutate(ef = hit_rate /
                      (sum(n_sig, na.rm = TRUE) / sum(n_total, na.rm = TRUE))) |>
      dplyr::ungroup()
  }

  create_pathway_heatmap <- function(data, algorithm_name,
                                     populations = names(PATHWAY_POPULATION_LABELS),
                                     min_cols = PATHWAY_MIN_COLS,
                                     max_rows = PATHWAY_MAX_ROWS) {

    base <- data |>
      dplyr::filter(algorithm == algorithm_name, study == "total",
                    covar_set == "covar",
                    !is.na(p_value), is.finite(p_plot), p_plot > 0)
    if (nrow(base) == 0) return(NULL)

    ## Middle panel: the three windows compared WITHIN each population.
    ## Populations are the facet columns and the windows the x axis, so a
    ## reader compares 3-, 5- and 10-year side by side in one stratum rather
    ## than comparing strata within one window.
    heat_src <- base |>
      dplyr::filter(exposure %in% PATHWAY_WINDOW_COLUMNS$exposure,
                    population %in% populations) |>
      add_enrichment() |>
      dplyr::left_join(PATHWAY_WINDOW_COLUMNS, by = "exposure") |>
      dplyr::mutate(
        neglogp = -log10(p_plot),
        stars = dplyr::case_when(!is.na(p_fdr) & p_fdr < 0.05 ~ "***",
                                 !is.na(p_fdr) & p_fdr < 0.10 ~ "**",
                                 p_value < 0.05               ~ "*",
                                 TRUE                         ~ ""),
        window = factor(window, levels = PATHWAY_WINDOW_COLUMNS$window),
        pop = factor(unname(PATHWAY_POPULATION_LABELS[population]),
                     levels = unname(PATHWAY_POPULATION_LABELS[populations])),
        category = unname(PATHWAY_CATEGORY_LABELS[categorize_pathway(pathway)]))
    if (nrow(heat_src) == 0) return(NULL)

    ## Rows: significant in at least min_cols of those columns, then the
    ## strongest max_rows by minimum p. Display limits only; the full table is
    ## pathway_results_long.xlsx.
    keep <- heat_src |>
      dplyr::filter(p_value < 0.05) |>
      dplyr::count(pathway, name = "n_cols") |>
      dplyr::filter(n_cols >= min_cols) |>
      dplyr::left_join(
        heat_src |>
          dplyr::group_by(pathway) |>
          dplyr::summarise(min_p = min(p_value, na.rm = TRUE), .groups = "drop"),
        by = "pathway") |>
      dplyr::arrange(min_p, dplyr::desc(n_cols))

    available <- nrow(keep)
    keep <- keep |> dplyr::slice_head(n = max_rows) |> dplyr::pull(pathway)

    if (length(keep) == 0) return(NULL)
    pd <- heat_src |> dplyr::filter(pathway %in% keep)

    ## Classes ordered by their best result, "Other" always last; pathways
    ## ordered within class so the strongest sits at the top of its block.
    cat_ord <- pd |>
      dplyr::group_by(category) |>
      dplyr::summarise(m = min(p_value, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(other = category == "Other") |>
      dplyr::arrange(other, m) |>
      dplyr::pull(category)

    path_ord <- pd |>
      dplyr::group_by(category, pathway) |>
      dplyr::summarise(m = min(p_value, na.rm = TRUE),
                       s = sum(neglogp, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(category = factor(category, levels = cat_ord)) |>
      dplyr::arrange(category, dplyr::desc(m), s) |>
      dplyr::pull(pathway)

    ## Truncate from the CENTRE, not the end. These are systematic names whose
    ## distinguishing part is the tail -- "Glycosphingolipid biosynthesis -
    ## ganglioseries" and "... - neolactoseries" share their first 34
    ## characters -- so trimming the end alone collapsed two pathways onto one
    ## label and duplicated a factor level. make.unique is a guard so a future
    ## collision degrades to a suffix rather than an error, and the mapping is
    ## built once so every panel labels a pathway the same way.
    ## End-truncation by default -- it reads better ("Aspartate and asparagine
    ## metab...") -- and centre-truncation ONLY for the names it would collide,
    ## which is where the tail is the distinguishing part. make.unique is the
    ## last guard so a collision degrades to a suffix rather than an error.
    lab_end  <- stringr::str_trunc(path_ord, 38)
    collides <- lab_end %in% lab_end[duplicated(lab_end)]
    path_lab <- purrr::set_names(
      make.unique(ifelse(collides,
                         stringr::str_trunc(path_ord, 38, side = "center"),
                         lab_end)),
      path_ord)
    short <- function(x) unname(path_lab[as.character(x)])

    pd <- pd |>
      dplyr::mutate(category = factor(category, levels = cat_ord),
                    pathway = factor(short(pathway), levels = unname(path_lab)))

    cat_fill <- purrr::set_names(rev_tints(length(cat_ord)), cat_ord)
    cat_text <- purrr::set_names(rev_tint_text(length(cat_ord)), cat_ord)
    pop_bg <- rev_tints(max(2, length(populations)), alpha = 0.30) |>
      utils::head(length(populations)) |>
      purrr::map(~ element_rect(fill = .x, colour = "grey80"))

    ## LEFT: class name and pathway names sharing one tinted block.
    ##
    ## The class sits left, the pathway names right. They collide if the block
    ## is too narrow for both, which is why the pathway names are truncated
    ## harder here than the underlying table needs -- the full names are in
    ## pathway_results_long.xlsx.
    row_df <- pd |> dplyr::distinct(pathway, category)
    catlab <- row_df |>
      dplyr::arrange(category, pathway) |>
      dplyr::group_by(category) |>
      dplyr::summarise(ymid = pathway[ceiling(dplyr::n() / 2)], .groups = "drop")

    ## ONE rect per class, not one per pathway. geom_rect draws a rect for
    ## every row it is given, and eight identical semi-transparent rects
    ## stacked in the amino-acid panel compounded 0.22 alpha into a solid
    ## block.
    block_df <- row_df |> dplyr::distinct(category)

    lab <- ggplot(row_df, aes(x = 1, y = pathway)) +
      geom_rect(data = block_df, inherit.aes = FALSE, aes(fill = category),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      geom_text(aes(label = pathway), x = 0.985, hjust = 1, size = 2.85,
                colour = "grey15") +
      geom_text(data = catlab,
                aes(y = ymid, label = category, colour = category),
                x = 0.02, hjust = 0, fontface = "bold", size = 2.9) +
      ggh4x::facet_grid2(category ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(values = cat_fill, guide = "none") +
      scale_colour_manual(values = cat_text, guide = "none") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(colour = "grey80", fill = NA),
            panel.spacing = unit(3, "pt"),
            axis.text = element_blank(), axis.ticks = element_blank(),
            plot.margin = margin(2, 0, 2, 2),
            strip.text = element_blank(), strip.background = element_blank())

    ## MIDDLE: the heatmap ---------------------------------------------------
    ##
    ## gradientn, not gradient2: gradient2 rescales symmetrically about the
    ## midpoint, so with p reaching 1e-5 on one side the whole non-significant
    ## half compressed into the palest tenth of the blue and never showed the
    ## colour at all. Anchoring the three stops at the data minimum, p = 0.05
    ## and the data maximum puts full blue back on the least significant cell.
    rng <- range(pd$neglogp, na.rm = TRUE)
    mid <- -log10(0.05)
    fill_stops <- sort(unique(c(min(rng[1], mid - 1e-6), mid, rng[2])))

    heat <- ggplot(pd, aes(x = window, y = pathway)) +
      geom_point(aes(size = ef, fill = neglogp), shape = 22,
                 colour = "grey75", stroke = 0.25) +
      geom_text(aes(label = stars), size = 3.2, vjust = 0.78,
                colour = "grey15", fontface = "bold") +
      ggh4x::facet_grid2(
        category ~ pop, scales = "free", space = "free",
        strip = ggh4x::strip_themed(
          background_x = pop_bg,
          background_y = rep(list(element_blank()), length(cat_ord)))) +
      scale_fill_gradientn(
        colours = c(SIG_COLORS[["FDR < 0.10"]], "white",
                    SIG_COLORS[["FDR < 0.05"]]),
        values = scales::rescale(fill_stops),
        limits = range(fill_stops),
        name = expression(bold(-log[10]~italic(p)))) +
      scale_size_area(max_size = 8, name = "Enrichment factor") +
      labs(x = NULL, y = NULL) +
      guides(
        fill = guide_colorbar(
          order = 1, title.position = "top", title.hjust = 0.5,
          theme = theme(legend.key.width = unit(3, "cm"),
                        legend.key.height = unit(0.5, "cm"))),
        size = guide_legend(order = 2, title.position = "top",
                            title.hjust = 0.5, nrow = 1,
                            override.aes = list(fill = "grey70",
                                                colour = "grey45"))) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            panel.spacing.x = unit(4, "pt"), panel.spacing.y = unit(3, "pt"),
            plot.margin = margin(2, 2, 2, 0),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            strip.text.y = element_blank(),
            ## A one-population figure names its stratum in the file name and
            ## the caption; repeating it in a strip over a single column is
            ## noise.
            strip.text.x = if (length(populations) > 1)
              element_text(face = "bold", size = 11) else element_blank(),
            strip.background.x = if (length(populations) > 1)
              element_rect() else element_blank(),
            legend.position = "bottom", legend.box = "horizontal",
            legend.title = element_text(face = "bold", size = 11),
            legend.background = element_blank(), legend.key = element_blank())

    ## RIGHT: how many single pollutants reach the same pathway --------------
    ##
    ## One segment per pollutant that reaches p < 0.05 for that pathway in the
    ## full cohort, so the bar length is the count and the colours say which
    ## ones. This is the attribution question beside the mixture result, the
    ## same role the tile half plays in the MWAS panel.
    ##
    ## The skeleton is completed to every pathway x pollutant with n = 0. A
    ## pathway that no pollutant reaches would otherwise be absent from this
    ## plot's data, and `space = "free_y"` sizes a facet from the rows its own
    ## data has -- the three panels would stop lining up.
    ## The 5-YEAR pollutants only. exposures_for() now returns the 3- and
    ## 10-year ones as well (23 in total), and this half of the figure is a
    ## count of how many pollutants reach a pathway at the primary window --
    ## mixing three windows into one bar would count the same pollutant up to
    ## three times, and POLLUTANT_COLOURS carries eight hues, one per species.
    pollutants <- intersect(exposures_for("total", "all"),
                            single_pollutant_exposures)

    bar_counts <- base |>
      dplyr::filter(population == "all", exposure %in% pollutants,
                    p_value < 0.05) |>
      dplyr::mutate(pathway = factor(short(pathway),
                                     levels = levels(pd$pathway)),
                    pollutant = rev_label_md(exposure)) |>
      dplyr::filter(!is.na(pathway)) |>
      dplyr::count(pathway, pollutant)

    bar_dat <- tidyr::crossing(
      row_df, pollutant = rev_label_md(pollutants)) |>
      dplyr::left_join(bar_counts, by = c("pathway", "pollutant")) |>
      dplyr::mutate(n = tidyr::replace_na(n, 0),
                    pollutant = factor(pollutant,
                                       levels = rev_label_md(pollutants)))

    bar_max <- max(1, max(
      bar_dat |> dplyr::group_by(pathway) |>
        dplyr::summarise(t = sum(n), .groups = "drop") |> dplyr::pull(t)))

    bar <- ggplot(bar_dat, aes(x = n, y = pathway, fill = pollutant)) +
      geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
      ggh4x::facet_grid2(category ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(values = POLLUTANT_COLOURS, name = "Single pollutant",
                        drop = FALSE) +
      scale_x_continuous(breaks = scales::breaks_width(1),
                         limits = c(0, bar_max), expand = expansion(mult = c(0, 0.04))) +
      labs(x = "Single pollutants at p < 0.05\n(full cohort)", y = NULL) +
      guides(fill = guide_legend(order = 3, title.position = "top",
                                 title.hjust = 0.5, nrow = 2)) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            panel.spacing = unit(3, "pt"),
            plot.margin = margin(2, 2, 2, 4),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(face = "bold", size = 10),
            strip.text = element_blank(), strip.background = element_blank(),
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 11),
            legend.text = ggtext::element_markdown(size = 10),
            legend.background = element_blank(), legend.key = element_blank())

    ## The label block was taking nearly a third of the figure. It only has to
    ## fit the pathway names, so it is narrowed and the heatmap given the room.
    layout <- lab + heat + bar +
      patchwork::plot_layout(widths = c(1.7, 2.1, 1.7), guides = "collect") &
      theme(legend.position = "bottom")

    list(plot = layout, n = dplyr::n_distinct(pd$pathway),
         n_available = available)
  }

  ## One figure over all three populations, plus one per population on the
  ## same layout. The per-population figures are not a different analysis --
  ## they are the same cells, re-drawn so a stratum can be read without the
  ## other two competing for the eye, and their row sets differ only because
  ## the "significant somewhere" rule is applied within what is shown.
  pathway_heatmap_targets <- c(
    list(list(pops = names(PATHWAY_POPULATION_LABELS), slug = "")),
    purrr::map(names(PATHWAY_POPULATION_LABELS), function(pop){
      list(pops = pop,
           slug = paste0("_", gsub("[^a-z0-9]+", "_", tolower(pop))))
    }))

  tidyr::expand_grid(alg = c("mummichog", "metapone"),
                     target = pathway_heatmap_targets) |>
    purrr::pwalk(function(alg, target){
      res <- create_pathway_heatmap(pathway_long, alg,
                                    populations = target$pops)
      if (is.null(res)) {
        message("No pathway rows for the ", alg, " heatmap",
                if (nzchar(target$slug)) paste0(" (", target$slug, ")") else "")
        return(invisible(NULL))
      }
      capped <- if (res$n < res$n_available)
        glue::glue(" (capped from {res$n_available})") else ""
      message("Pathway heatmap: ", alg,
              if (nzchar(target$slug)) target$slug else " [all populations]",
              " - ", res$n, " pathways", capped)

      ## Width scales with the number of population facets; the label and bar
      ## panels are a fixed cost, so a one-population figure must not be as
      ## wide as the three-population one or its squares stretch.
      width <- 7.5 + 2.6 * length(target$pops)
      ggsave(
        filename = file.path(rev_dir("figures", "pathway"),
                             glue::glue("pathway_heatmap_{alg}{target$slug}.png")),
        plot = res$plot, width = width, height = 0.30 * res$n + 3.2,
        dpi = 300, bg = "white", limitsize = FALSE)
    })

} else {
  message("Pathway summary not found at ", pathway_long_path,
          " - run R6-pathway_revision.R first. Skipping Section 6.")
}


# =============================================================================
# SECTION 7: COMBINED PANEL FIGURE
# =============================================================================

## The submitted figure (figures/mwas/.../combined_panel_*.png) was a single
## row of three panels -- a volcano, the subgroup-versus-full-cohort betas and
## the sensitivity-versus-primary betas. Reviewer 1 comment 14 calls the last
## two "largely descriptive displays of p-value distributions and coefficient
## agreement" and asks that they move to the Supporting Information and that
## the freed space carry biological content instead.
##
## The revision figure is therefore two rows:
##
##   A  volcano for this exposure                  \
##   B  Manhattan over m/z                          |  what the MWAS found
##   C  the same metabolites under each individual pollutant
##   D  effect estimates with 95% CIs for the confirmed identifications
##   E  exposure-response for the top Level 1 metabolites
##
## Panels C, D and E are built in Section 7A above. The two coefficient
## scatters keep their code here -- create_population_scatter() and
## create_covariate_scatter() are still taken verbatim from
## scripts/7-visualization.R, so the SI figure stays directly comparable with
## the submitted one.

create_population_scatter <- function(mwas_list_c18, mwas_list_hilic,
                                      covar_set, exposure_name) {

  # Helper to prep data for one column type
  prep_pop_data <- function(mwas_list, column_type) {
    df_all <- mwas_list[["all"]][[covar_set]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_all = logFC, P.Value_all = P.Value)

    df_demcind <- mwas_list[["demcind"]][[covar_set]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_demcind = logFC, P.Value_demcind = P.Value)

    df_no_demcind <- mwas_list[["no demcind"]][[covar_set]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_no_demcind = logFC, P.Value_no_demcind = P.Value)

    # Join all three populations
    df_all |>
      dplyr::inner_join(df_demcind, by = "met") |>
      dplyr::inner_join(df_no_demcind, by = "met") |>
      # Classify significance based on demcind and no demcind
      dplyr::mutate(
        significant = case_when(
          P.Value_demcind < 0.05 & P.Value_no_demcind < 0.05 ~ "Both P < 0.05",
          P.Value_demcind < 0.05 ~ "demcind P < 0.05",
          P.Value_no_demcind < 0.05 ~ "no demcind P < 0.05",
          TRUE ~ "NS"
        )
      ) |>
      # Pivot to long for faceting
      tidyr::pivot_longer(
        cols = c(logFC_demcind, logFC_no_demcind),
        names_to = "comparison",
        values_to = "logFC_sub",
        names_prefix = "logFC_"
      ) |>
      dplyr::mutate(
        comparison = dplyr::recode(comparison,
                                   "demcind" = "demcind vs all",
                                   "no_demcind" = "no demcind vs all"),
        column_type = column_type
      )
  }

  plot_data <- dplyr::bind_rows(
    prep_pop_data(mwas_list_c18, "C18/neg-"),
    prep_pop_data(mwas_list_hilic, "HILIC/pos+")
  ) |>
    # Only keep points where P < 0.05 in "all"
    dplyr::filter(P.Value_all < 0.05) |>
    dplyr::mutate(
      significant = factor(significant,
                           levels = c("Both P < 0.05",
                                      "demcind P < 0.05",
                                      "no demcind P < 0.05",
                                      "NS"))
    )

  # Per-facet correlation on P < 0.05 points (all points shown are P<0.05 in all)
  # cor_labels <- plot_data |>
  #   dplyr::group_by(comparison) |>
  #   dplyr::summarise(
  #     cor_test = list(cor.test(logFC_all, logFC_sub, use = "complete.obs")),
  #     .groups = "drop"
  #   ) |>
  #   dplyr::mutate(
  #     r = purrr::map_dbl(cor_test, ~ .x$estimate),
  #     p_val = purrr::map_dbl(cor_test, ~ .x$p.value),
  #     label = paste0("r = ", round(r, 3), "\np = ",
  #                    ifelse(p_val < 0.001,
  #                           formatC(p_val, format = "e", digits = 2),
  #                           round(p_val, 3)))
  #   )

  p <- ggplot(plot_data, aes(x = logFC_all, y = logFC_sub)) +
    # NS points
    geom_point(
      data = . %>% filter(significant == "NS"),
      aes(shape = column_type),
      color = "grey70", alpha = 0.5, size = 1.5
    ) +
    # Significant points
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant, shape = column_type),
      alpha = 0.7, size = 2
    ) +
    # Fit line on all displayed points (all P < 0.05 in "all")
    geom_smooth(
      method = "lm", se = TRUE,
      color = "black", linewidth = 0.8, linetype = "solid",
      alpha = 0.2
    ) +
    scale_color_manual(
      values = c(
        "Both P < 0.05" = "#B73F42",
        "demcind P < 0.05" = "#DE9960",
        "no demcind P < 0.05" = "#436C85"
      ),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
      name = "Column"
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    # Per-facet correlation labels
    ggpubr::stat_cor(
      method = "pearson",
      label.x.npc = "left", label.y.npc = "top",
      size = 5.5, fontface = "italic"
    ) +
    # geom_text(
    #   data = cor_labels,
    #   aes(x = Inf, y = -Inf, label = label),
    #   hjust = 1.1, vjust = -0.3, size = 4, fontface = "italic",
    #   inherit.aes = FALSE
    # ) +
    facet_wrap(~ comparison) +
    labs(
      title = gsub("exp_", "", exposure_name),
      x = "MWAS beta coefficients among all participants",
      y = "MWAS beta coefficients among the subgroup population"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 13),
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12)
    )

  return(p)
}

create_covariate_scatter <- function(mwas_list_c18, mwas_list_hilic,
                                     population, cov1, cov2,
                                     exposure_name) {

  # Helper to prep and join data for one column type
  prep_cov_data <- function(mwas_list, column_type) {
    df1 <- mwas_list[[population]][[cov1]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_cov1 = logFC, P.Value_cov1 = P.Value)

    df2 <- mwas_list[[population]][[cov2]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_cov2 = logFC, P.Value_cov2 = P.Value)

    dplyr::inner_join(df1, df2, by = "met") |>
      dplyr::mutate(column_type = column_type)
  }

  plot_data <- dplyr::bind_rows(
    prep_cov_data(mwas_list_c18, "C18/neg-"),
    prep_cov_data(mwas_list_hilic, "HILIC/pos+")
  ) |>
    dplyr::mutate(
      significant = case_when(
        P.Value_cov1 < 0.05 & P.Value_cov2 < 0.05 ~ "Both P < 0.05",
        P.Value_cov1 < 0.05 ~ "Primary \u2013 P < 0.05",
        P.Value_cov2 < 0.05 ~ "Sensitivity \u2013 P < 0.05",
        TRUE ~ "NS"
      ),
      significant = factor(significant,
                           levels = c("Both P < 0.05",
                                      "Primary \u2013 P < 0.05",
                                      "Sensitivity \u2013 P < 0.05",
                                      "NS"))
    )

  # Correlation on P < 0.05 points only
  sig_data <- plot_data |>
    dplyr::filter(significant != "NS")

  # cor_test <- cor.test(sig_data$logFC_cov1, sig_data$logFC_cov2,
  #                      use = "complete.obs")
  # r <- cor_test$estimate
  # p_val <- cor_test$p.value
  # cor_label <- paste0("r = ", round(r, 3), "\np = ",
  #                     ifelse(p_val < 0.001, formatC(p_val, format = "e", digits = 2),
  #                            round(p_val, 3)))

  p <- ggplot(plot_data, aes(x = logFC_cov1, y = logFC_cov2)) +
    # NS points
    geom_point(
      data = . %>% filter(significant == "NS"),
      aes(shape = column_type),
      color = "grey70", alpha = 0.4, size = 1.5
    ) +
    # Significant points
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant, shape = column_type),
      alpha = 0.7, size = 2
    ) +
    # Fit line on P < 0.05 points only
    geom_smooth(
      data = sig_data,
      aes(x = logFC_cov1, y = logFC_cov2),
      method = "lm", se = TRUE,
      color = "black", linewidth = 0.8, linetype = "solid",
      alpha = 0.2, inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Both P < 0.05" = "#B73F42",
        setNames("#DE9960", "Primary \u2013 P < 0.05"),
        setNames("#436C85", "Sensitivity \u2013 P < 0.05")
      ),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
      name = "Column"
    ) +
    # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    ggpubr::stat_cor(
      method = "pearson",
      label.x.npc = "left", label.y.npc = "top",
      size = 5.5, fontface = "italic"
    ) +
    # annotate("text", x = Inf, y = -Inf,
    #          label = cor_label,
    #          hjust = 1.1, vjust = -0.3, size = 4, fontface = "italic") +
    labs(
      title = paste0(gsub("exp_", "", exposure_name), " (", population, ")"),
      x = "MWAS beta coefficients in the primary analysis",
      y = "MWAS beta coefficients in the sensitivity analysis (additional covariates)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 13),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12)
    )

  return(p)
}


# =============================================================================
# SECTION 7A: BIOLOGICAL PANELS FOR THE COMBINED FIGURE
# =============================================================================

## Reviewer 1 comment 14: "Consider consolidating them, moving the
## coefficient-comparison panels to the Supporting Information, and replacing
## the freed space with content that carries biological information:
## single-pollutant effect comparison, effect estimates for confirmed
## identifications with confidence intervals, or exposure-response plots for
## the Level 1 metabolites."
##
## All three are built here and all three go into the freed space:
##
##   C  ONE panel with two halves sharing the metabolite axis --
##        left   effect estimates with 95% confidence intervals for the Level 1
##               (confirmed) identifications under the exposure this figure is
##               about
##        right  the same metabolites, same order, under each individual
##               pollutant the index is built from, so the mixture signal can
##               be read against its components
##   D  exposure-response plots for the Level 1 metabolites -- covariate-
##      adjusted partial residuals against the exposure, with the fitted linear
##      slope and a loess curve so a departure from linearity is visible
##
## The metabolites are the PANEL_N_FEATURES lowest-P Level 1 identifications
## in that cell; panel D takes the PANEL_N_RESPONSE lowest-P of those. Both
## subtitles say so, because neither is guessable from the axes.
##
## Effect sizes are reported as percent differences throughout, which is what
## add_percent_difference() writes into the annotated tables (R1 minor comment
## 6) and what makes these panels comparable with Hu et al. and Qi et al.

## Level 1 only. VOLCANO_LABEL_MAX_LEVEL is the same constant the volcano
## labels with, so a compound named in panel A is a compound plotted in
## panels C and D.
PANEL_MAX_LEVEL     <- VOLCANO_LABEL_MAX_LEVEL
## How many identified metabolites the two effect-size panels carry. Twelve
## fills the panel height at REV_TEXT_PANEL sizes without crowding.
PANEL_N_FEATURES    <- 12
## Exposure-response is a scatter per metabolite, so it takes the top few only.
PANEL_N_RESPONSE    <- 4
## Below this many complete observations a per-feature fit is not worth
## drawing; the small cognitive strata can fall under it.
PANEL_MIN_N         <- 50


# Exposure-response extract ---------------------------------------------------

## Written by R4 (see the "Exposure-response extract" section there): the
## in-house-library feature abundances plus the standardized analysis frames.
## Panel D needs the abundances themselves, which no other object in this
## script carries. If the extract has not been rebuilt since that section was
## added, panel D is skipped and the figure falls back to three panels rather
## than failing.
er_extract_path <- rev_here("data", "processed",
                            "exposure_response_extract_revision.RData")

ER_AVAILABLE <- file.exists(er_extract_path)

if (ER_AVAILABLE) {
  load(er_extract_path)
  message("Exposure-response extract loaded: ",
          ncol(er_abundance_c18) - 3, " C18 and ",
          ncol(er_abundance_hilic) - 3, " HILIC features")
} else {
  message("No exposure-response extract at ", er_extract_path,
          " - re-run R4-mwas_revision.R to build it. ",
          "Panel D will be omitted.")
}


# Which pollutants a composite is actually built from -------------------------

## comp_pca_all_w1 and comp_pca_all_w3 are SEVEN-pollutant indices. CALINE4
## NOx exists for model years 1998-2002 only, so R2 drops it from the short
## windows on a coverage threshold (see window_pollutant_coverage there);
## comp_pca_all_w10 carries all eight. Showing NOx beside a seven-pollutant
## index would invite the reader to attribute that index to a component that
## is not in it, so the single-pollutant half of the effect-estimate panel is
## restricted to the components the index was actually built from.
##
## Read from R2's own loadings table rather than hardcoded here: R2 decides
## the composition from a coverage threshold, and a change there has to reach
## this figure.
window_extract_path <- rev_here("data", "processed",
                                "exposure_windows_revision.RData")

## WHICH SINGLE-POLLUTANT MODELS BELONG BESIDE WHICH INDEX
##
## Two things have to match, and until 2026-09-03 neither did for the window
## arms:
##
##   COMPOSITION. comp_pca_all_w1 and comp_pca_all_w3 are seven-pollutant
##     indices (CALINE4 NOx exists for 1998-2002 only), so showing eight
##     columns beside them implies NOx contributed to an index it is not in.
##
##   WINDOW. A ten-year index sat beside five-year pollutant models. Nothing in
##     the figure said so, and the whole point of the half is to ask whether
##     the index's signal tracks one of its components -- which is unanswerable
##     if the components were averaged over a different period.
##
## Both are now derived from the exposure name. An exposure at window w takes
## that window's single pollutants when they were fitted (w3 and w10 carry
## their own; see MWAS_WINDOWS in R1); anything else takes the 5-year set,
## intersected with the index's composition when R2's loadings table gives one.
## comp_pca_all_w1 is the case that still falls back: there is no 1-year
## single-pollutant arm, so it shows the 5-year models for its seven
## components and create_single_pollutant_panel() says so on the panel.
WINDOW_INDEX_COMPOSITION <- if (file.exists(window_extract_path)) {
  env <- new.env()
  load(window_extract_path, envir = env)
  env$window_pca_loadings |>
    dplyr::distinct(window, pollutant) |>
    dplyr::mutate(exposure = paste0("comp_pca_all_", window),
                  var      = paste0("exp_", pollutant, "_iqr")) |>
    dplyr::filter(exposure %in% window_pca_exposures) |>
    (\(d) split(d$var, d$exposure))()
} else {
  message("No exposure-window extract at ", window_extract_path,
          " - the single-pollutant half will show all eight pollutants ",
          "even for the seven-pollutant window indices.")
  list()
}

if (length(WINDOW_INDEX_COMPOSITION) > 0) {
  purrr::iwalk(WINDOW_INDEX_COMPOSITION, function(vars, exp){
    message("  ", exp, " is built from ", length(vars), " pollutants")
  })
}

## The pollutant exposures to show beside `exposure_name`, and the window they
## were averaged over. Returns a list so the caller can say which it got.
components_for_exposure <- function(exposure_name, available) {
  w <- exposure_window(exposure_name)

  if (!identical(w, "")) {
    matched <- intersect(available, single_pollutant_exposures_for(w))
    if (length(matched) > 0) {
      return(list(pollutants = matched, window = w, matched = TRUE))
    }
  }

  ## Fall back to the 5-year models, narrowed to the index's own components
  ## where those are known.
  base <- intersect(available, single_pollutant_exposures)
  comp <- WINDOW_INDEX_COMPOSITION[[exposure_name]]
  if (!is.null(comp)) base <- intersect(base, comp)

  ## `matched` is about the index's own averaging window, not about whether a
  ## tag was parsed: comp_pca_all_w1 lands here with w == "" but is a one-year
  ## index, and calling that a match would be exactly the silent mismatch this
  ## function exists to remove.
  list(pollutants = base,
       window     = "",
       matched    = isTRUE(exposure_window_years(exposure_name) == 5))
}


# Identified-metabolite effect estimates --------------------------------------

PLATFORM_TAGS <- c(c18 = "C18", hilic = "HILIC")

## Distinct axis labels. Two features can carry the same compound name -- the
## same metabolite seen on both columns, or twice on one column at different
## retention times -- and a duplicated label silently collapses two rows of a
## forest plot into one. The platform tag separates the first case and the
## observed m/z the second.
unique_feature_labels <- function(df) {
  base <- paste0(df$compound_display, " (", df$platform_tag, ")")
  dup  <- base %in% base[duplicated(base)]
  mz   <- stringr::str_extract(df$met, "(?<=mz_rt_)[0-9.]+")
  ifelse(dup,
         paste0(df$compound_display, " (", df$platform_tag, " ",
                round(as.numeric(mz), 2), ")"),
         base)
}

## Effect estimates for the identifications in one cell of the grid.
##
## Reads the annotated tables rather than the raw topTables, because those are
## what carry compound, confidence_level and the percent-difference interval.
## pct_diff is recomputed defensively for any cell written before
## add_percent_difference() was introduced.
identified_effects <- function(study, population, covar_set, exposure_name,
                               max_level = PANEL_MAX_LEVEL) {

  grab <- function(annot_list, platform) {
    tag <- unname(PLATFORM_TAGS[[platform]])
    df  <- annot_list[[study]][[population]][[covar_set]][[exposure_name]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!"pct_diff" %in% names(df)) df <- add_percent_difference(df)
    df |>
      dplyr::filter(!is.na(confidence_level), confidence_level <= max_level,
                    !is.na(compound), compound != "") |>
      ## One row per feature: the annotated tables are long over candidate
      ## compounds, and a Level 1 feature resolves to one compound anyway.
      dplyr::group_by(met) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::mutate(platform = platform, platform_tag = tag)
  }

  out <- dplyr::bind_rows(grab(mwas_full_annotated_list_c18,   "c18"),
                          grab(mwas_full_annotated_list_hilic, "hilic"))

  if (nrow(out) == 0) return(out)

  out <- out |>
    dplyr::mutate(
      compound_display = candidate_display_name(compound),
      ## Same three display classes the volcano and the Manhattan use. VIP is
      ## a DISPLAY class here as it is there -- see the note on
      ## filter_significant() in create_volcano_plot above.
      significant = sig_class(P.Value, adj.P.Val)
    )

  out$label <- unique_feature_labels(out)
  out
}

## The metabolites the whole bottom row is about: the strongest identifications
## under THIS exposure, ranked by p-value. Both halves of panel C and panel D
## draw from this one set, so the row reads as views of the same metabolites
## rather than unrelated displays.
##
## No p-value gate here, deliberately. These are effect ESTIMATES, and the
## reviewer asked for the estimates for the confirmed identifications, not for
## the ones that passed a threshold; filtering the forest on its own p-value
## would be a selection on the outcome. The rank-and-cap is a display limit --
## how many rows fit -- and the subtitle says so.
##
## The volcano labels (create_volcano_plot) apply the SAME ranking to the SAME
## Level 1 set, gated at P < 0.05 and capped lower, plus every identification
## at FDR < 0.10 regardless of rank. A label there marks a point a reader is
## being pointed at rather than reporting an estimate, which is why it carries
## a gate at all.
##
## The labelled compounds in panel A are a subset of the rows in panel C: the
## gate can only remove higher-p features, and no FDR < 0.10 identification
## anywhere in the grid falls outside C's twelve rows (checked). If
## PANEL_N_FEATURES is ever lowered, re-check that.
panel_feature_set <- function(study, population, covar_set, exposure_name,
                              n = PANEL_N_FEATURES) {
  eff <- identified_effects(study, population, covar_set, exposure_name)
  if (nrow(eff) == 0) return(eff)
  eff |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = n) |>
    ## Ordered by effect size, not by p-value: panels C and D share this order,
    ## which is what lets a reader carry a row across from one to the other.
    dplyr::mutate(label = forcats::fct_reorder(label, pct_diff))
}


# Panel C: single-pollutant effect comparison ---------------------------------

## Percent difference per IQR for each individual pollutant, over the
## metabolites panel D shows. This is the attribution question the reviewer
## raises -- whether the mixture signal tracks one component -- put next to the
## mixture result rather than in a separate table.
##
## MATCHED TO THE INDEX IT SITS BESIDE (changed 2026-09-03). The panel used to
## be hardwired to population `all` and to the 5-year pollutants, whatever the
## index in the rest of the figure was. Both mismatches were invisible in the
## figure and both undercut what the half is for:
##
##   POPULATION. A dementia/CIND index sat beside full-cohort pollutant
##     estimates, so a difference between the two halves could be the stratum
##     rather than the pollutant. The single pollutants are now fitted in all
##     three cognitive strata (STRATA_POPULATIONS in R1), so the panel reads
##     the SAME population as the rest of the figure whenever those models
##     exist. `all predx` does not carry them and falls back to `all`.
##
##   WINDOW. A 10-year index sat beside 5-year pollutant estimates. The 3- and
##     10-year windows now carry their own single-pollutant models, so a
##     windowed index gets windowed pollutants.
##
## Whatever it ends up reading is stated in the subtitle, so a fallback is
## legible rather than silent. For the cox arm no single-pollutant models exist
## at all and the function returns NULL.
create_single_pollutant_panel <- function(study, population, covar_set,
                                          feature_set, exposure_name) {

  if (nrow(feature_set) == 0) return(NULL)

  ## The population to read from: this figure's own, when it carries the
  ## single-pollutant models, and the full cohort otherwise.
  pop_use <- population
  available <- exposures_for(study, pop_use) |> purrr::keep(is_single_pollutant)
  if (length(available) == 0) {
    pop_use   <- "all"
    available <- exposures_for(study, pop_use) |>
      purrr::keep(is_single_pollutant)
  }
  if (length(available) == 0) return(NULL)

  comp       <- components_for_exposure(exposure_name, available)
  pollutants <- comp$pollutants
  if (length(pollutants) == 0) return(NULL)

  plot_data <- pollutants |>
    purrr::set_names() |>
    purrr::map(function(pollutant){
      identified_effects(study, pop_use, covar_set, pollutant,
                         max_level = PANEL_MAX_LEVEL) |>
        dplyr::select(dplyr::any_of(c("met", "pct_diff", "P.Value",
                                      "adj.P.Val")))
    }) |>
    purrr::list_rbind(names_to = "pollutant") |>
    dplyr::inner_join(feature_set |> dplyr::select(met, label), by = "met")

  if (nrow(plot_data) == 0) return(NULL)

  ## What the panel actually shows, said plainly. Two clauses, each present
  ## only when it differs from the rest of the figure -- an unqualified
  ## subtitle means the pollutants are matched on both counts.
  pop_note <- if (identical(pop_use, population)) {
    unname(rev_population_label(pop_use))
  } else {
    paste0("full cohort (", unname(rev_population_label(population)),
           " not fitted)")
  }

  window_note <- if (identical(comp$window, "")) {
    if (comp$matched) {
      "5-year window"
    } else {
      paste0("5-year window (index is ",
             exposure_window_years(exposure_name), "-year)")
    }
  } else {
    paste0(WINDOW_YEARS[[comp$window]], "-year window")
  }

  plot_data <- plot_data |>
    dplyr::mutate(
      ## Markdown, rendered by element_markdown() below: NO2 and PM2.5 carry
      ## subscripts in their names and are wrong written flat.
      pollutant = factor(rev_label_md(pollutant),
                         levels = rev_label_md(pollutants)),
      ## Symmetric fill limits so that zero is white and a positive and a
      ## negative difference of the same size read as equally strong.
      ## The same three classes the colours use elsewhere in the figure,
      ## encoded as a count of asterisks -- the fill is already spent on the
      ## effect size here, so significance has to ride on the glyph.
      stars = dplyr::case_when(adj.P.Val < 0.05 ~ "***",
                               adj.P.Val < 0.10 ~ "**",
                               P.Value   < 0.05 ~ "*",
                               TRUE             ~ "")
    )

  fill_limit <- max(abs(plot_data$pct_diff), na.rm = TRUE)

  ggplot(plot_data, aes(x = pollutant, y = label, fill = pct_diff)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = stars), size = 5, vjust = 0.78, colour = "grey15") +
    scale_fill_gradient2(
      low = SIG_COLORS[["FDR < 0.10"]], mid = "white",
      high = SIG_COLORS[["FDR < 0.05"]], midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      ## Four breaks, not the default five: at the collected legend's text size
      ## five labels ran into each other on the bar ("-10-5 0 5 10").
      n.breaks = 4,
      name = "% difference per IQR"
    ) +
    ## The title fits on one line beside the bar in the collected bottom
    ## strip, and the bar is given room to actually read as a gradient -- the
    ## default key is the size of a discrete legend swatch.
    guides(fill = guide_colourbar(
      title.position = "left", title.vjust = 0.85,
      theme = theme(legend.key.width  = unit(5, "cm"),
                    legend.key.height = unit(0.55, "cm")))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = "Single-pollutant models",
      ## The unit is named here, not only on the colour bar: the bar lives in
      ## the collected legend strip at the foot of the figure, far from these
      ## tiles, while the forest declares "per 1 SD" directly beneath itself.
      ## The two halves use different increments -- an IQR is the same 25th-to
      ## -75th percentile shift for every pollutant whatever its skew, which is
      ## what makes the eight columns comparable with each other, whereas the
      ## index has no physical scale and only an SD to be expressed in -- so
      ## each half has to state its own.
      ## Two lines. The asterisk key made this longer than the tiles are wide,
      ## and it ran into the forest's subtitle to its left -- the two halves
      ## are one panel, so there is no margin between them to absorb it.
      subtitle = paste0(pop_note, ", ", window_note, ", per IQR\n",
                        "* P < 0.05, ** FDR < 0.10, *** FDR < 0.05"),
      x = NULL, y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.x   = element_text(angle = 45, hjust = 1),
      axis.line     = element_blank(),
      axis.ticks    = element_blank(),
      legend.position = "right"
    )
}


# Panel D: effect estimates for the confirmed identifications -----------------

## Forest plot of the percent difference and its 95% interval for the Level 1
## identifications. The interval comes from the coefficient's own standard
## error, recovered as beta / t and back-transformed on the log2 scale by
## add_percent_difference() in R1 -- so it is the moderated limma standard
## error, not one refitted here.
create_identified_forest <- function(feature_set, exposure_name) {

  if (nrow(feature_set) == 0) return(NULL)
  if (!all(c("pct_diff_lo", "pct_diff_hi") %in% names(feature_set))) return(NULL)

  unit <- if (is_single_pollutant(exposure_name)) "per IQR" else "per 1 SD"

  ggplot(feature_set, aes(x = pct_diff, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = pct_diff_lo, xmax = pct_diff_hi,
                      colour = significant),
                  orientation = "y", width = 0, linewidth = 0.9,
                  na.rm = TRUE) +
    geom_point(aes(colour = significant), size = 3.2, na.rm = TRUE) +
    ## Same palette object the volcano uses, so the Significance legend the
    ## combined panel collects covers this panel too and only one is drawn.
    ## NS is not in SIG_COLORS and falls through to na.value, exactly as it
    ## does on the volcano.
    scale_colour_manual(values = SIG_COLORS, name = "Significance",
                        drop = FALSE) +
    labs(
      title = "Confirmed (Level 1) identifications",
      ## The merged panel puts this beside the single-pollutant tiles, so the
      ## left half has to say which exposure it is the estimate for -- and how
      ## these twelve were picked out of the Level 1 set.
      ## Markdown (<br>, not \n) -- element_markdown() below renders the
      ## subscripts when this figure is about NO2 or PM2.5.
      subtitle = glue::glue(
        "{rev_label_md(exposure_name)}<br>",
        "{nrow(feature_set)} lowest-P Level 1 features"),
      ## Short enough to stay on one line at REV_TEXT_PANEL sizes:
      ## rev_wrap_axis_titles() breaks at 38 characters, and a two-line title
      ## here runs into the exposure-response panel's own axis title.
      x = glue::glue("% difference {unit} (95% CI)"),
      y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14),
      axis.text  = element_text(size = 13)
    )
}


# Panel D: exposure-response for the Level 1 metabolites ----------------------

## Covariate-adjusted partial residuals against the exposure, one facet per
## metabolite.
##
## The partial residual is  e + beta * x  from an ordinary least squares fit of
## the feature on the exposure and the SAME covariate set the MWAS used, so the
## red line's slope IS that model's exposure coefficient and the x axis is the
## exposure on the scale the coefficients are reported in. The loess curve is
## the point of the panel: it shows whether the linear term the MWAS fits is a
## fair summary of the relationship, which no volcano or Manhattan can.
##
## Two deliberate simplifications, both display-only. The fit is OLS rather
## than limma's duplicateCorrelation model, so participants with two draws
## contribute twice; and the moderated variance is not used. Neither changes
## the drawn slope materially -- the exposure coefficient is unweighted in both
## -- but the panel is a description of the relationship, not a second estimate
## of it, and the caption should say so.
partial_residual_data <- function(frame, abundance, met, exposure_name,
                                  covars) {

  if (!met %in% names(abundance)) return(NULL)
  if (!exposure_name %in% names(frame)) return(NULL)

  covars <- intersect(covars, names(frame))

  d <- frame |>
    dplyr::select(dplyr::all_of(c("file.name_new", exposure_name, covars))) |>
    dplyr::inner_join(
      abundance |> dplyr::select(dplyr::all_of(c("file.name_new", met))),
      by = "file.name_new")

  ## Copied into fixed names rather than renamed: a feature id is
  ## "mz_rt_373.8841_87.2473" and an exposure can be "exp_pm2.5_iqr", neither
  ## of which is a syntactic name, and both go into a formula below.
  d$.x <- d[[exposure_name]]
  d$.y <- d[[met]]

  d <- d |>
    dplyr::select(dplyr::all_of(c(".x", ".y", covars))) |>
    tidyr::drop_na()

  if (nrow(d) < PANEL_MIN_N) return(NULL)

  ## A covariate that is constant within a stratum (diab_at_blooddraw in a
  ## small subgroup, say) makes model.matrix rank-deficient; drop those rather
  ## than let lm() return an NA coefficient.
  keep <- covars |>
    purrr::keep(function(v) dplyr::n_distinct(d[[v]]) > 1)

  fit <- stats::lm(
    stats::as.formula(paste(".y ~ .x", paste(c("", keep), collapse = " + "))),
    data = d)

  beta <- stats::coef(fit)[[".x"]]
  if (!is.finite(beta)) return(NULL)

  tibble::tibble(x = d$.x, partial = stats::residuals(fit) + beta * d$.x)
}

create_exposure_response_panel <- function(study, population, covar_set,
                                           exposure_name, feature_set,
                                           n_features = PANEL_N_RESPONSE) {

  if (!ER_AVAILABLE || nrow(feature_set) == 0) return(NULL)
  ## comp_qgcomp_fw_all is a contrast over the eight quantized pollutants, not
  ## a scored column, so there is no exposure axis to plot it against.
  if (is_qgcomp_fw(exposure_name)) return(NULL)

  frames <- list(c18   = er_data_list_c18[[study]][[population]][[covar_set]],
                 hilic = er_data_list_hilic[[study]][[population]][[covar_set]])
  abundances <- list(c18 = er_abundance_c18, hilic = er_abundance_hilic)
  covars <- er_covars_list[[covar_set]]

  top <- feature_set |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = n_features) |>
    dplyr::mutate(facet_label = ifelse(duplicated(compound_display) |
                                         duplicated(compound_display,
                                                    fromLast = TRUE),
                                       as.character(label),
                                       compound_display))

  plot_data <- seq_len(nrow(top)) |>
    purrr::map(function(i){
      row <- top[i, ]
      pd <- partial_residual_data(frames[[row$platform]],
                                  abundances[[row$platform]],
                                  row$met, exposure_name, covars)
      if (is.null(pd)) return(NULL)
      ## Short facet labels. The full "Taurine (HILIC)" form panels C and D
      ## use is wider than a facet strip in this slot and gets clipped; the
      ## platform is already on the same metabolite one panel to the left.
      pd |> dplyr::mutate(label = as.character(row$facet_label))
    }) |>
    purrr::compact() |>
    purrr::list_rbind()

  if (nrow(plot_data) == 0) return(NULL)

  ## Facets in the order the metabolites were ranked, not alphabetically.
  plot_data <- plot_data |>
    dplyr::mutate(label = factor(label, levels = intersect(
      as.character(top$facet_label), unique(label))))

  unit <- if (is_single_pollutant(exposure_name)) "IQR units" else "SD units"

  ggplot(plot_data, aes(x = x, y = partial)) +
    geom_point(colour = "grey55", alpha = 0.25, size = 1.1) +
    geom_smooth(method = "loess", formula = y ~ x, span = 1,
                se = TRUE, colour = SIG_COLORS[["FDR < 0.10"]],
                fill = SIG_COLORS[["FDR < 0.10"]],
                alpha = 0.18, linewidth = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = SIG_COLORS[["FDR < 0.05"]], linewidth = 0.9) +
    ## Compound names run long ("5-Hydroxy-L-tryptophan") and a facet strip
    ## here is a quarter of a third of a row wide, so the label is wrapped
    ## rather than clipped.
    facet_wrap(~ label, scales = "free_y", nrow = 2,
               labeller = ggplot2::label_wrap_gen(16)) +
    labs(
      title = "Exposure-response",
      ## Says which four of the twelve in panel C these are, and on what
      ## quantity -- neither is guessable from the axes.
      subtitle = glue::glue(
        "{n_features} lowest-P Level 1 features; partial residuals\n",
        "Red: linear fit, blue: loess"),
      ## Comma, not a second parenthesis: several labels already carry one
      ## ("PCA, 1-year window (7 pollutants)") and stacking a second reads
      ## badly.
      x = glue::glue("{rev_label_md(exposure_name)}, {unit}"),
      ## NOT the same quantity as the volcano's x axis. That one is a
      ## coefficient -- the change in log2 abundance per 1 SD of the index.
      ## This is a level: the feature's own log2 abundance with the covariates'
      ## contribution removed, so it is centred near zero and the SLOPE across
      ## the panel is what corresponds to the volcano's x value. "Partial
      ## residual" is named in the subtitle rather than crowded into the axis.
      y = expression(bold("Adjusted "*log[2]*" abundance"))
    ) +
    theme_classic() +
    theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title    = element_text(face = "bold", size = 14),
      strip.text    = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}


# Level 1 effect-estimate table -----------------------------------------------

## The figure panels show the top identifications per cell; this writes every
## Level 1 estimate in every cell, which is what the response letter needs to
## cite and what a reader who wants a number rather than a picture will look
## for.
message("\nAssembling the Level 1 effect-estimate table ...")

level1_effect_table <- study_names |>
  purrr::set_names() |>
  purrr::map(function(study){
    ## set_names() at every level: list_rbind(names_to = ) writes the LIST
    ## names into the column, and mapping over a bare character vector
    ## produces an unnamed list, which would put 1, 2, 3 there instead of
    ## the population and covariate-set names.
    populations_for(study) |>
      purrr::set_names() |>
      purrr::map(function(pop){
        covar_names |>
          purrr::set_names() |>
          purrr::map(function(cov){
            exposures_for(study, pop) |>
              purrr::set_names() |>
              purrr::map(function(exp){
                identified_effects(study, pop, cov, exp) |>
                  dplyr::select(dplyr::any_of(c(
                    "met", "compound", "confidence_level", "n_candidates",
                    "adduct", "delta_ppm", "platform", "logFC",
                    "pct_diff", "pct_diff_lo", "pct_diff_hi",
                    "P.Value", "adj.P.Val", "VIP_comp1")))
              }) |>
              purrr::list_rbind(names_to = "exposure")
          }) |>
          purrr::list_rbind(names_to = "covar_set")
      }) |>
      purrr::list_rbind(names_to = "population")
  }) |>
  purrr::list_rbind(names_to = "study") |>
  dplyr::mutate(exposure_label = rev_label(exposure), .after = "exposure") |>
  dplyr::arrange(study, population, covar_set, exposure, P.Value)

rev_save_table(level1_effect_table, "level1_effect_estimates", "mwas")

message("  ", nrow(level1_effect_table), " Level 1 estimates across ",
        dplyr::n_distinct(level1_effect_table$exposure), " exposures")


# Build the panels ------------------------------------------------------------

## Boxed legends under each panel, as in the submitted figure.
panel_legend_theme <- theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.background = element_rect(colour = "grey80", fill = "white",
                                   linewidth = 0.5),
  legend.margin = margin(4, 6, 4, 6),
  legend.title = element_text(face = "bold", size = 10),
  legend.text = element_text(size = 9)
)

## Reviewer 1 comment 14: the volcano and the Manhattan are consolidated into
## the main-text panel, and the two coefficient-agreement scatters move to the
## Supporting Information. build_panel_parts() therefore builds the SI figure
## now; build_main_panel_parts() builds the main-text one.

build_main_panel_parts <- function(exposure_name, population, covar_set,
                                   study) {

  volcano <- create_volcano_plot(
    mwas_result_c18 = mwas_results_list_c18[[study]][[population]][[covar_set]][[exposure_name]],
    mwas_result_hilic = mwas_results_list_hilic[[study]][[population]][[covar_set]][[exposure_name]],
    annotation_result_c18 = mwas_full_annotated_list_c18[[study]][[population]][[covar_set]][[exposure_name]],
    annotation_result_hilic = mwas_full_annotated_list_hilic[[study]][[population]][[covar_set]][[exposure_name]],
    exposure_name = exposure_name,
    label_size = 6
  ) +
    labs(title = NULL) +
    guides(shape = guide_legend(order = 1), color = guide_legend(order = 2)) +
    panel_legend_theme

  ## ONE significance legend in the combined panel.
  ##
  ## guides = "collect" merges guides only when they compare EQUAL, and these
  ## two do not: the volcano maps `shape` as well as `colour`, so its colour
  ## guide draws its keys with a shape aesthetic attached and the Manhattan's
  ## does not. patchwork therefore kept both, and the second was clipped at the
  ## right edge of the figure. Suppressing the Manhattan's colour guide is
  ## deterministic where relying on guide equality is not; the volcano keeps
  ## the Significance legend (and the Column legend, which the Manhattan has no
  ## counterpart for, since it facets by platform instead).
  man <- create_manhattan(
    combined_results_list_c18[[study]][[population]][[covar_set]][[exposure_name]],
    combined_results_list_hilic[[study]][[population]][[covar_set]][[exposure_name]],
    exposure_name = exposure_name
  ) +
    labs(title = NULL) + guides(colour = "none") + panel_legend_theme

  ## Bottom row -- the biological content Reviewer 1 comment 14 asks for. All
  ## three read from ONE ranked set of identifications so the row is three
  ## views of the same metabolites (see panel_feature_set in Section 7A).
  ##
  ## Each can be absent for a good reason and the layout adapts rather than
  ## failing: the cox arm fits no single-pollutant models (panel C), a cell
  ## may carry no Level 1 identification at all (C), and panel D needs
  ## both the R4 exposure-response extract and a scored exposure column, which
  ## the comp_qgcomp_fw_all contrast does not have.
  feature_set <- panel_feature_set(study, population, covar_set, exposure_name)

  forest <- create_identified_forest(feature_set, exposure_name)
  ## The tiles exist only as the right half of the merged effect-estimate
  ## panel -- they have no y axis of their own -- so they are not built when
  ## there is no forest to attach them to.
  single <- if (is.null(forest)) NULL else
    create_single_pollutant_panel(study, population, covar_set, feature_set,
                                  exposure_name)
  er     <- create_exposure_response_panel(study, population, covar_set,
                                           exposure_name, feature_set)

  ## The volcano carries the ONE Significance legend for the whole figure, so
  ## the forest's identical colour guide is suppressed here for the same
  ## reason the Manhattan's is above -- guide equality is not something to
  ## rely on when the two scales differ in their other aesthetics.
  if (!is.null(forest)) forest <- forest + guides(colour = "none")

  ## `identified` before `single_pollutant`: the two are merged into one panel
  ## below, with the forest on the left carrying the metabolite names and the
  ## pollutant tiles on the right sharing that axis, and this is the order they
  ## are merged in.
  parts <- purrr::compact(list(volcano           = volcano,
                               manhattan         = man,
                               identified        = forest,
                               single_pollutant  = single,
                               exposure_response = er)) |>
    purrr::map(rev_wrap_axis_titles) |>
    purrr::map(~ .x + REV_TEXT + REV_TEXT_PANEL)

  ## Panel-title and strip sizes, applied AFTER the shared themes so they are
  ## not overwritten by them. REV_TEXT sets an 18 pt plot title and
  ## REV_TEXT_PANEL a 19 pt strip, both sized for a panel that spans half the
  ## figure. Panels C, D and E are a third of a row wide, and E puts four
  ## facets inside that, so at the shared sizes their titles and strip labels
  ## are wider than the space they have and ggplot clips them.
  bio <- intersect(c("single_pollutant", "identified", "exposure_response"),
                   names(parts))
  parts[bio] <- parts[bio] |>
    purrr::map(~ .x + theme(plot.title = element_text(face = "bold",
                                                      size = 15, hjust = 0.5),
                            plot.subtitle = element_text(size = 12,
                                                         hjust = 0.5)))

  ## ggtext, applied element by element rather than globally ----------------
  ##
  ## Only the three elements that actually carry a species name are markdown.
  ## The tiles' SUBTITLE deliberately is not: it is the asterisk key
  ## ("* P < 0.05, ** FDR < 0.10, ..."), and commonmark reads a leading
  ## asterisk as emphasis, so rendering it as markdown would eat the very
  ## glyphs it exists to explain.
  if (!is.null(parts$single_pollutant)) {
    parts$single_pollutant <- parts$single_pollutant +
      theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1,
                                                   size = 18))
  }
  if (!is.null(parts$identified)) {
    parts$identified <- parts$identified +
      theme(plot.subtitle = ggtext::element_markdown(size = 12, hjust = 0.5))
  }
  if (!is.null(parts$exposure_response)) {
    parts$exposure_response <- parts$exposure_response +
      theme(axis.title.x = ggtext::element_markdown(face = "bold", size = 19))
  }

  if (!is.null(parts$exposure_response)) {
    ## The x axis names the exposure, and the composite names are long
    ## ("PCA all toxicants (unsupervised)"), so rev_wrap_axis_titles() breaks
    ## it at 38 characters. Since the bottom row went from three panels to two
    ## this panel is wide enough to carry it on one line, which is what keeps
    ## its axis title the same weight as every other panel's -- so the wrap is
    ## undone here. It has to happen AFTER that call, not before, or the call
    ## simply re-wraps it.
    ##
    ## The strip text is the one size that does have to come down: four facets
    ## inside one panel slot leaves each strip a quarter of it.
    xlab <- parts$exposure_response$labels$x
    parts$exposure_response <- parts$exposure_response +
      theme(strip.text = element_text(face = "bold", size = 13))
    if (is.character(xlab) && length(xlab) == 1) {
      parts$exposure_response <- parts$exposure_response +
        labs(x = stringr::str_replace_all(xlab, "\n", " "))
    }
  }

  ## Tags -------------------------------------------------------------------
  ##
  ## Set by hand, and set HERE rather than in create_combined_panel(), because
  ## the forest and the tiles become one nested patchwork immediately below
  ## and a tag added to a patchwork has no single plot to attach to. The
  ## letters are assigned over the FINAL element list, so the merged panel
  ## takes one letter, not two, and the sequence has no gap in it.
  final_names <- setdiff(names(parts), "single_pollutant")
  tags <- purrr::set_names(LETTERS[seq_along(final_names)], final_names)

  parts[final_names] <- purrr::map2(
    parts[final_names], tags[final_names],
    ~ .x + labs(tag = .y) +
      theme(plot.tag = element_text(face = "bold", size = 30)))

  ## Merge the effect estimates and the single-pollutant tiles ---------------
  ##
  ## They carry the same metabolites in the same order, so the names belong on
  ## the figure once. The forest keeps the y axis; the tiles drop theirs and
  ## sit flush against it, which is also what lets the reader carry a row
  ## straight across from the mixture estimate to its eight components.
  if (!is.null(parts$identified) && !is.null(parts$single_pollutant)) {
    parts$single_pollutant <- parts$single_pollutant +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(5.5, 5.5, 5.5, 0))

    parts$identified <- patchwork::wrap_plots(
      list(parts$identified, parts$single_pollutant), nrow = 1,
      widths = IDENTIFIED_SUBPANEL_WIDTHS)

    parts$single_pollutant <- NULL
  }

  list(plots   = parts,
       weights = panel_weights(names(parts),
                               merged = inherits(parts$identified, "patchwork")))
}

build_panel_parts <- function(exposure_name, population, covar_set, study) {

  ## The volcano is not built here any more -- it is panel A of the main-text
  ## figure (build_main_panel_parts), and this function now returns only the
  ## two coefficient-agreement scatters that comment 14 moves to the SI.

  ## Both scatter panels compare one exposure across cells of the grid, so
  ## they can only be drawn where that exposure was actually fitted. The
  ## single pollutants are fitted in `all` only
  ## (SINGLE_POLLUTANT_POPULATIONS in R1), so they have no
  ## subgroup-versus-full-cohort contrast to draw -- panel B is dropped for
  ## them rather than erroring on a NULL topTable.
  has_exposure <- function(pop, cov) {
    !is.null(mwas_results_list_c18[[study]][[pop]][[cov]][[exposure_name]]) &&
      !is.null(mwas_results_list_hilic[[study]][[pop]][[cov]][[exposure_name]])
  }

  ## Panel B is the subgroup-versus-full-cohort contrast. As in
  ## scripts/7-visualization.R it is the same comparison regardless of which
  ## population the rest of the panel shows, so it is included throughout and
  ## dropped only if a study lacks the cognitive strata.
  p2 <- NULL
  if (all(c("demcind", "no demcind") %in% names(mwas_results_list_c18[[study]])) &&
      all(purrr::map_lgl(c("all", "demcind", "no demcind"),
                         has_exposure, cov = covar_set))) {
    p2 <- create_population_scatter(
      mwas_results_list_c18[[study]], mwas_results_list_hilic[[study]],
      covar_set = covar_set, exposure_name = exposure_name
    ) +
      labs(title = NULL) + guides(shape = "none") + panel_legend_theme
  }

  ## Panel C is the primary-versus-sensitivity covariate contrast, within the
  ## population being shown. Every exposure carries both covariate sets, so
  ## this normally draws; the guard is for a grid that has been trimmed.
  p3 <- NULL
  if (all(purrr::map_lgl(c("covar", "covar_sen"),
                         function(cs) has_exposure(population, cs)))) {
    p3 <- create_covariate_scatter(
      mwas_results_list_c18[[study]], mwas_results_list_hilic[[study]],
      population = population, cov1 = "covar", cov2 = "covar_sen",
      exposure_name = exposure_name
    ) +
      labs(title = NULL) + guides(shape = "none") + panel_legend_theme
  }

  purrr::compact(list(p2, p3)) |>
    purrr::map(rev_wrap_axis_titles) |>
    purrr::map(~ .x + REV_TEXT + REV_TEXT_PANEL)
}

## Annotation is applied once, to the finished layout. Tagging a patchwork and
## then nesting it inside another restarts the sequence, so the Manhattan
## variant is assembled from the raw panel list rather than from the annotated
## three-panel object -- that is what lets its Manhattan carry tag D.
annotate_panel <- function(layout) {
  layout +
    ## No figure title: the exposure, population and covariate set belong in
    ## the manuscript caption, and the title crowded the top of the panel.
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(face = "bold", size = 30))
    )
}

## Main-text panel geometry.
##
## The bottom row does not always carry the same panels (see
## build_main_panel_parts), so the figure is assembled row by row from
## explicit width weights rather than from one flat panel list. Each row is a
## patchwork; the rows are then stacked.
##
## The Manhattan takes the wider share of the top row because the m/z axis is
## what carries the retention-time / void-region argument, and a squeezed one
## cannot show it. On the bottom row the effect-estimate panel is widest
## because it is two sub-panels sharing one metabolite axis.
MAIN_PANEL_WEIGHTS <- c(volcano           = 5,
                        manhattan         = 8,
                        identified        = 5,
                        single_pollutant  = 4,
                        exposure_response = 5)

MAIN_PANEL_ROW <- c(volcano           = 1L,
                    manhattan         = 1L,
                    identified        = 2L,
                    exposure_response = 2L)

## Forest against tiles inside the merged panel. The forest takes slightly
## more because it also carries the metabolite names.
IDENTIFIED_SUBPANEL_WIDTHS <- c(5, 4)

## Row weights for the elements that survived. When the forest and the tiles
## are merged the merged panel claims both their weights, so the bottom row
## keeps the proportions it had when they were two panels.
panel_weights <- function(nms, merged) {
  w <- MAIN_PANEL_WEIGHTS[nms]
  if (merged && "identified" %in% nms) {
    w[["identified"]] <- w[["identified"]] +
      MAIN_PANEL_WEIGHTS[["single_pollutant"]]
  }
  w
}

create_combined_panel <- function(exposure_name, population, covar_set, study) {

  built   <- build_main_panel_parts(exposure_name, population, covar_set, study)
  panels  <- built$plots
  weights <- built$weights

  rows <- MAIN_PANEL_ROW[names(panels)]

  row_plots <- sort(unique(rows)) |>
    purrr::map(function(r){
      keep <- rows == r
      patchwork::wrap_plots(unname(panels[keep]), nrow = 1,
                            widths = unname(weights[names(panels)[keep]]))
    })

  ## patchwork aligns panel edges across EVERYTHING in one assembly, so the
  ## bottom row's metabolite names -- "N6,N6,N6-Trimethyl-L-lysine (HILIC)" --
  ## widened the shared left column and pushed the volcano's y-axis title an
  ## inch away from its own axis. free() releases the lower rows from that
  ## alignment: each row now sizes its own margins, and nothing is clipped
  ## (which `type = "space"` would do).
  ##
  ## BOTH SIDES, not just "l". Freeing only the left still left the two rows
  ## aligned on the right, which is a constraint on the same layout from the
  ## other end: the exposure-response panel had to end where the Manhattan
  ## ends, and since its width share is fixed the slack came out as an inch
  ## and a half of white space between its tag and its own y-axis title.
  ## Rendering the bottom row on its own showed no such gap, which is what
  ## identified the cross-row alignment as the cause rather than anything in
  ## the panel itself.
  if (length(row_plots) > 1) {
    row_plots[-1] <- row_plots[-1] |>
      purrr::map(~ patchwork::free(.x, type = "panel", side = "lr"))
  }

  ## guides = "collect" merges identical legends into one strip beneath the
  ## whole figure instead of drawing a boxed legend under each panel, which is
  ## what was colliding. It only merges guides that compare equal, hence the
  ## shared SIG_COLORS above. The volcano's Column (shape) guide has no
  ## counterpart in the Manhattan, which facets by platform instead, so it is
  ## carried alongside rather than merged; the single-pollutant tiles' fill
  ## gradient is its own legend and joins the same strip.
  layout <- purrr::reduce(row_plots, `/`) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.justification = "center")

  list(panel = layout,
       n     = length(panels),
       nrow  = dplyr::n_distinct(rows))
}

## Supporting Information: the two coefficient-agreement scatters, moved out of
## the main text per Reviewer 1 comment 14. Returns NULL when neither can be
## drawn -- the single-pollutant exposures are fitted in `all` only, so they
## have no subgroup contrast (see build_panel_parts).
create_si_scatter_panel <- function(exposure_name, population, covar_set,
                                    study) {
  panels <- build_panel_parts(exposure_name, population, covar_set, study)
  if (length(panels) == 0) return(NULL)
  layout <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.justification = "center")
  list(panel = annotate_panel(layout), n = length(panels))
}

study_names |>
  purrr::walk(function(study){
    populations_for(study) |>
      purrr::walk(function(pop){
        covar_names |>
          purrr::walk(function(cov){
            exposures_for(study, pop) |>
              purrr::walk(function(exp){
                message("Combined panel: ", study, "_", pop, " [", cov, "] - ", exp)
                rev_dir("figures", "mwas", study, pop, cov, exp)

                main <- create_combined_panel(exp, pop, cov, study)
                ## plain ggsave: the parts already carry REV_TEXT plus the
                ## panel-specific sizes, and rev_ggsave's `& REV_TEXT` would
                ## overwrite the latter.
                ##
                ## Height follows the number of rows the layout actually
                ## built. A cell with no Level 1 identification and no
                ## single-pollutant models falls back to the one-row
                ## volcano + Manhattan figure and keeps its old size, so
                ## nothing is stretched to fill space it does not use.
                ggsave(
                  filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                      glue::glue("combined_panel_{exp}.png")),
                  plot = main$panel,
                  width  = if (main$nrow > 1) 21 else 14,
                  height = if (main$nrow > 1) 14 else 6,
                  dpi = 300, bg = "white"
                )

                si <- create_si_scatter_panel(exp, pop, cov, study)
                if (!is.null(si)) {
                  ggsave(
                    filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                        glue::glue("si_scatters_{exp}.png")),
                    plot = si$panel, width = 10 * si$n, height = 9,
                    dpi = 300, bg = "white"
                  )
                }
              })
          })
      })
  })


message("\nVisualization completed!")
message("Figures saved to:")
message("  - ", rev_here("figures", "mwas"))
message("  - ", rev_here("figures", "pathway"))
message("  - ", rev_here("figures", "composites"), " (from R3)")

#--------------------------------End of the code--------------------------------
