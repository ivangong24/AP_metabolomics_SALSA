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
##          8. Mantel-style composite-toxicant network with the composite-
##             composite correlation inset, paired with the fold-to-fold
##             spread of the cross-fitted WQS / QGcomp weights (Reviewer 1
##             comment 1: correlation among the revision exposures, and
##             weight uncertainty rather than point estimates)
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
                                 "Porphyrin", "Catabolism",
                                 "Caffeine", "Folate", "Riboflavin",
                                 "Thiamine", "Biotin", "Pantothenate")

## Purine and pyrimidine metabolism are nucleotide metabolism. Purine used to
## sit under vitamin/cofactor, presumably because its degradation products run
## through the same cofactor chemistry, and pyrimidine matched nothing at all
## and fell through to "other" -- so the two halves of one class were drawn in
## two different blocks.
nucleotide_metabolism <- c("Purine", "Pyrimidine", "Nucleotide")

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
  "nucleotide"       = "#64B5CD",
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
    str_detect(pathway_name, regex(str_c(nucleotide_metabolism, collapse = "|"),
                                   ignore_case = TRUE)) ~ "nucleotide",
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

  ## The pathway figures are drawn for the CROSS-SECTIONAL ("total") study
  ## only. The time-to-event composite family is withdrawn under Reviewer 1
  ## comment 1, so a cox pathway panel would enrich an exposure the paper no
  ## longer reports.
  tidyr::expand_grid(study = intersect(study_names, "total"),
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

  ## The same rule the mummichog call applies (min_hits = 3 in
  ## R6-pathway_revision.R): a pathway needs at least three matched
  ## significant metabolites to be reported. Mummichog enforces it internally,
  ## so for its figures this is a no-op; metapone does not, and without it a
  ## single-metabolite pathway can carry the strongest enrichment factor in
  ## the figure. Applying it to both is what makes the two comparable.
  PATHWAY_MIN_HITS <- 3

  ## Columns: the cognitive strata, at the PRIMARY exposure only.
  ##
  ## This used to carry the 3-, 5- and 10-year windows side by side. It no
  ## longer does. The window comparison is a sensitivity analysis reported
  ## under Reviewer 1 comment 13, and putting three windows in the main pathway
  ## figure invited exactly the reading that comment 1 commits us to avoiding
  ## -- comparing exposures by how many pathways they light up, across indices
  ## correlated at r = 0.97. The figure shows comp_pca_all, the 5-year primary
  ## index, and the other windows stay in the Supporting Information table.
  ##
  ## `scales = "free_x"` and `space = "free_x"` are kept downstream so that a
  ## cell missing for any other reason collapses rather than drawing an empty
  ## column.
  PATHWAY_WINDOW_COLUMNS <- tibble::tribble(
    ~exposure,          ~window,
    "comp_pca_all",     "5-year")

  ## The three averaging windows, for the window-comparison figure only. The
  ## caution above still applies -- these indices correlate at r = 0.97, so
  ## the figure is a stability check, not three findings -- and the caption
  ## has to say so.
  PATHWAY_WINDOW_COLUMNS_3 <- tibble::tribble(
    ~exposure,             ~window,
    "comp_pca_all_w3",     "3-year",
    "comp_pca_all",        "5-year",
    "comp_pca_all_w10",    "10-year")

  ## Eight pollutants, so a fixed qualitative palette. Subscripts are drawn
  ## with ggtext in the legend, which is why the labels carry markup.
  PATHWAY_POLLUTANT_COLOURS <- c(
    "Benzene"       = "#55A868",
    "1,3-Butadiene" = "#C44E52",
    "Chromium"      = "#9DBBCB",
    "Lead"          = "#8C5A2B",
    "NO<sub>2</sub>"  = "#8172B3",
    "Nickel"        = "#DFA53A",
    "PM<sub>2.5</sub>" = "#3B3B3B",
    "NO<sub>x</sub>"  = "#DE9960")

  PATHWAY_POLLUTANT_LABELS <- c(
    "exp_benzene_iqr"   = "Benzene",
    "exp_butadiene_iqr" = "1,3-Butadiene",
    "exp_chromium_iqr"  = "Chromium",
    "exp_lead_iqr"      = "Lead",
    "exp_no2_iqr"       = "NO<sub>2</sub>",
    "exp_nickel_iqr"    = "Nickel",
    "exp_pm2.5_iqr"     = "PM<sub>2.5</sub>",
    "exp_nox_iqr"       = "NO<sub>x</sub>")

  PATHWAY_POPULATION_LABELS <- c("all"        = "All",
                                 "no demcind" = "No dem/CIND",
                                 "demcind"    = "Dem/CIND")

  PATHWAY_CATEGORY_LABELS <- c(
    "amino acid" = "Amino acid", "carbohydrate" = "Carbohydrate",
    "lipid" = "Lipid", "energy" = "Energy", "inflammation" = "Inflammation",
    "vitamin/cofactor" = "Vitamin/cofactor", "nucleotide" = "Nucleotide",
    "signaling" = "Signalling",
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
                                     max_rows = PATHWAY_MAX_ROWS,
                                     window_cols = PATHWAY_WINDOW_COLUMNS,
                                     single_pollutant_panel = FALSE,
                                     ref_window = NULL,
                                     pathway_whitelist = NULL,
                                     min_hits = PATHWAY_MIN_HITS) {

    base <- data |>
      dplyr::filter(algorithm == algorithm_name, study == "total",
                    covar_set == "covar",
                    !is.na(p_value), is.finite(p_plot), p_plot > 0)

    ## Optionally restrict to a common library. Metapone tests 372 pathways
    ## against mummichog's 94, and most of the extra ones are KEGG disease and
    ## signalling entries -- glucagon signalling, amoebiasis, central carbon
    ## metabolism in cancer -- which have no counterpart in the other tool and
    ## land in "Other". Matching on the normalised name is how the two
    ## libraries are compared, since metapone lower-cases its pathway names.
    if (!is.na(min_hits) && "n_sig" %in% names(base)) {
      dropped <- dplyr::n_distinct(base$pathway)
      base <- dplyr::filter(base, !is.na(.data$n_sig), .data$n_sig >= min_hits)
      dropped <- dropped - dplyr::n_distinct(base$pathway)
      if (dropped > 0)
        message("    min_hits = ", min_hits, ": dropped ", dropped,
                " pathway(s) with fewer matched metabolites")
    }

    if (!is.null(pathway_whitelist)) {
      norm_name <- function(x) tolower(trimws(gsub("\\s+", " ", x)))
      before <- dplyr::n_distinct(base$pathway)
      base <- dplyr::filter(base,
                            norm_name(.data$pathway) %in%
                              norm_name(pathway_whitelist))
      message("    library restriction: ", dplyr::n_distinct(base$pathway),
              " of ", before, " pathways retained")
    }
    if (nrow(base) == 0) return(NULL)

    ## Middle panel: the three windows compared WITHIN each population.
    ## Populations are the facet columns and the windows the x axis, so a
    ## reader compares 3-, 5- and 10-year side by side in one stratum rather
    ## than comparing strata within one window.
    heat_src <- base |>
      dplyr::filter(exposure %in% window_cols$exposure,
                    population %in% populations) |>
      add_enrichment() |>
      dplyr::left_join(window_cols, by = "exposure") |>
      dplyr::mutate(
        neglogp = -log10(p_plot),
        stars = dplyr::case_when(!is.na(p_fdr) & p_fdr < 0.05 ~ "***",
                                 !is.na(p_fdr) & p_fdr < 0.10 ~ "**",
                                 p_value < 0.05               ~ "*",
                                 TRUE                         ~ ""),
        window = factor(window, levels = window_cols$window),
        pop = factor(unname(PATHWAY_POPULATION_LABELS[population]),
                     levels = unname(PATHWAY_POPULATION_LABELS[populations])),
        category = unname(PATHWAY_CATEGORY_LABELS[categorize_pathway(pathway)]))
    if (nrow(heat_src) == 0) return(NULL)

    ## Rows: either significant in a NAMED reference column, or significant
    ## in at least min_cols of them.
    ##
    ## The reference rule is what the window figures want. A pathway earns its
    ## row by the primary analysis, and the other windows are then shown for
    ## it whether or not they reach p < 0.05 -- which is the point, since a
    ## window that misses is as informative as one that hits. Selecting on
    ## "significant somewhere" instead would let a pathway in on the strength
    ## of a sensitivity window alone.
    ##
    ## Either way max_rows caps the figure height. Display limits only; the
    ## full table is pathway_results_long.xlsx.
    keep <- heat_src |>
      dplyr::filter(p_value < 0.05,
                    if (is.null(ref_window)) TRUE
                    else as.character(window) == ref_window) |>
      dplyr::count(pathway, name = "n_cols") |>
      dplyr::filter(n_cols >= if (is.null(ref_window)) min_cols else 1L) |>
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
    cat_ord <- (if (is.null(ref_window)) pd else
                  dplyr::filter(pd, as.character(window) == ref_window)) |>
      dplyr::group_by(category) |>
      dplyr::summarise(m = min(p_value, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(other = category == "Other") |>
      dplyr::arrange(other, m) |>
      dplyr::pull(category)

    ## Rank by the REFERENCE window when there is one, not by the best p
    ## across the windows shown. With three windows the cross-window minimum
    ## reorders the rows against the primary analysis the text reports --
    ## glutathione metabolism, whose strongest result is at ten years, jumped
    ## above glutamate metabolism, which is stronger at five. The figure has
    ## to agree with the analysis it is a figure of.
    rank_src <- if (is.null(ref_window)) pd else
      dplyr::filter(pd, as.character(window) == ref_window)
    if (nrow(rank_src) == 0) rank_src <- pd

    ## Ties need a real tiebreaker. When ranking is on one window, `s` is that
    ## window's own -log10 p and carries no information beyond `m`, so two
    ## pathways with the same p fell back to whatever order the table happened
    ## to hold them in. Metapone makes this common rather than rare: its
    ## permutation p-values are granular at 1/1000, so arginine and proline
    ## and alanine and aspartate metabolism both read p = 0.001 and their
    ## order was arbitrary. The enrichment factor breaks the tie on the
    ## quantity the squares are already sized by, largest nearest the top.
    path_ord <- rank_src |>
      dplyr::group_by(category, pathway) |>
      dplyr::summarise(m = min(p_value, na.rm = TRUE),
                       s = sum(neglogp, na.rm = TRUE),
                       e = max(ef, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(category = factor(category, levels = cat_ord),
                    e = dplyr::if_else(is.finite(.data$e), .data$e, 0)) |>
      dplyr::arrange(category, dplyr::desc(m), s, e) |>
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
      geom_text(aes(label = pathway), x = 0.985, hjust = 1, size = 3.3,
                colour = "grey15") +
      geom_text(data = catlab,
                aes(y = ymid, label = category, colour = category),
                x = 0.02, hjust = 0, fontface = "bold", size = 3.4) +
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

    ## Light grey banding behind alternate rows, as in the submitted pathway
    ## figure. With a dozen rows and up to three columns the eye loses the row
    ## on the way across, and a band is cheaper than gridlines, which would
    ## compete with the squares.
    ##
    ## Drawn as one tile per CELL, not as a rect spanning the panel. A rect
    ## needs numeric ymin/ymax, which switches the y scale to continuous, and
    ## that breaks `space = "free"`: every facet then gets the height of the
    ## whole row set and the rows stop lining up with the label block.
    ##
    ## EVERY row is in the layer, banded ones opaque and the rest transparent.
    ## Passing only the banded rows looks right in a single panel and is wrong
    ## in a facetted one: with free scales the panel's positions come from the
    ## levels each layer actually carries, so four banded rows out of eight
    ## were re-indexed onto positions one to four and the band collapsed into
    ## a block at the bottom of the facet.
    ## Parity runs down the figure AS DRAWN, not along the factor. Categories
    ## stack top to bottom, but within a facet a discrete axis puts level one
    ## at the BOTTOM, so numbering the levels in order alternates correctly
    ## inside a block and can collide across one: a class with an odd row
    ## count leaves its last band touching the first band of the class below,
    ## which is what put the lipid and energy blocks against each other.
    ## Ordering by category and then by DESCENDING level walks the rows in the
    ## order a reader sees them, so no two adjacent rows share a band.
    stripe <- pd |>
      dplyr::distinct(pathway, category) |>
      dplyr::arrange(category, dplyr::desc(pathway)) |>
      dplyr::mutate(band = dplyr::row_number() %% 2 == 0) |>
      tidyr::expand_grid(dplyr::distinct(pd, window, pop))

    heat <- ggplot(pd, aes(x = window, y = pathway)) +
      geom_tile(data = stripe, inherit.aes = FALSE,
                aes(x = window, y = pathway),
                width = 1, height = 1,
                fill = ifelse(stripe$band, "grey95", NA)) +
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
            ## With a single exposure column the tick label repeats the
            ## caption in every facet; it earns its place only when there is
            ## more than one column to tell apart. Horizontal: the columns are
            ## short labels and there is room, so the 45-degree tilt bought
            ## nothing but a harder read.
            axis.text.x = if (nrow(window_cols) > 1)
              element_text(size = 11)
            else element_blank(),
            axis.ticks.x = if (nrow(window_cols) > 1)
              element_line() else element_blank(),
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

    ## RIGHT (optional): which single pollutants reach the same pathway.
    ##
    ## One stacked bar per pathway, one segment per pollutant at p < 0.05 in
    ## the full cohort. This is an attribution aid, not a ranking: the bar
    ## counts nominally significant pollutants and nothing more, and comment 1
    ## commits us to not comparing exposures by counts. The caption has to say
    ## so, and the composite columns to its left carry the actual result.
    sp <- NULL
    if (single_pollutant_panel) {
      sp_src <- base |>
        dplyr::filter(population == "all",
                      exposure %in% names(PATHWAY_POLLUTANT_LABELS),
                      p_value < 0.05, pathway %in% keep) |>
        dplyr::mutate(
          pollutant = factor(unname(PATHWAY_POLLUTANT_LABELS[exposure]),
                             levels = names(PATHWAY_POLLUTANT_COLOURS)),
          category = unname(PATHWAY_CATEGORY_LABELS[categorize_pathway(pathway)]))

      ## Every retained pathway needs a row even when no pollutant reaches it,
      ## or facet_grid2 drops it and the three panels stop lining up.
      sp_rows <- pd |> dplyr::distinct(pathway, category)
      sp_dat <- sp_src |>
        dplyr::mutate(category = factor(category, levels = cat_ord),
                      pathway = factor(short(pathway),
                                       levels = unname(path_lab))) |>
        dplyr::count(pathway, category, pollutant, name = "n")

      ## A zero-height bar for every pollutant, so that a pollutant reaching
      ## none of these pathways still gets a key. drop = FALSE alone keeps the
      ## level but leaves its swatch empty, because the glyph is drawn from
      ## the layer's data and there is none.
      sp_keys <- tidyr::expand_grid(
        pollutant = factor(names(PATHWAY_POLLUTANT_COLOURS),
                           levels = names(PATHWAY_POLLUTANT_COLOURS)),
        sp_rows |> dplyr::slice(1)) |>
        dplyr::mutate(n = 0)

      sp <- ggplot(sp_rows, aes(y = pathway)) +
        geom_blank() +
        geom_col(data = sp_keys, aes(x = n, fill = pollutant), width = 0.7) +
        geom_col(data = sp_dat, aes(x = n, fill = pollutant),
                 width = 0.7, colour = "white", linewidth = 0.2) +
        ggh4x::facet_grid2(category ~ ., scales = "free_y", space = "free_y") +
        scale_fill_manual(values = PATHWAY_POLLUTANT_COLOURS,
                          name = "Single pollutant", drop = FALSE) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.04)),
                           breaks = scales::breaks_width(1)) +
        guides(fill = guide_legend(order = 3, title.position = "top",
                                   title.hjust = 0.5, nrow = 2,
                                   byrow = TRUE)) +
        ## As the axis title, not a caption: patchwork collects the legends
        ## into a strip under the whole figure, and a caption sits below that,
        ## a long way from the panel it names.
        labs(x = "Single pollutants at p < 0.05\n(full cohort)", y = NULL) +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank(),
              plot.margin = margin(2, 2, 2, 2),
              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 10),
              axis.title.x = element_text(face = "bold", size = 10,
                                          lineheight = 1.1),
              strip.text.y = element_blank(), strip.background.y = element_blank(),
              legend.position = "bottom",
              legend.title = element_text(face = "bold", size = 11),
              legend.text = ggtext::element_markdown(size = 10),
              legend.background = element_blank(), legend.key = element_blank())
    }

    ## The label block was taking nearly a third of the figure. It only has to
    ## fit the pathway names, so it is narrowed and the heatmap given the room.
    ## The label block carries a class name on the left and a pathway name on
    ## the right of the same strip, so it needs width for both: at the larger
    ## type "Vitamin/cofactor" ran into "Vitamin B1 (thiamin) metabolism".
    layout <- if (is.null(sp)) {
      lab + heat +
        patchwork::plot_layout(widths = c(2.05, 2.1), guides = "collect")
    } else {
      lab + heat + sp +
        patchwork::plot_layout(widths = c(2.05, 2.1, 1.2), guides = "collect")
    }
    layout <- layout & theme(legend.position = "bottom")

    list(plot = layout, n = dplyr::n_distinct(pd$pathway),
         n_available = available, has_sp = !is.null(sp))
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
      ## The single-pollutant bars go on the full-cohort figure only. In the
      ## stratum figures they would be the same eight full-cohort bars beside
      ## a different stratum's heatmap, which invites reading them as that
      ## stratum's result.
      ## The full-cohort figure carries the three averaging windows and the
      ## single-pollutant bars; the rows are the pathways the 5-year primary
      ## analysis calls significant, so the other two windows are read as a
      ## stability check on that result rather than as results of their own.
      ## The stratum figures keep the single 5-year column: their question is
      ## the stratum, and three windows per stratum is a different figure.
      ## The single-pollutant bars are built (single_pollutant_panel) but not
      ## used: a count of nominally significant pollutants per pathway is not
      ## interpretable as attribution, and comment 1 commits us to not
      ## comparing exposures by counts. Kept behind the flag rather than
      ## deleted, since the panel itself took some getting right.
      ## Metapone is restricted to the pathways mummichog also tests, so the
      ## two figures compare like with like. Without it the metapone figure is
      ## dominated by library entries the primary analysis never saw.
      full_cohort <- identical(target$pops, "all")
      wl <- if (alg == "metapone")
        unique(pathway_long$pathway[pathway_long$algorithm == "mummichog"])
        else NULL
      res <- create_pathway_heatmap(
        pathway_long, alg, populations = target$pops,
        single_pollutant_panel = FALSE,
        pathway_whitelist = wl,
        window_cols = if (full_cohort) PATHWAY_WINDOW_COLUMNS_3
                      else PATHWAY_WINDOW_COLUMNS,
        ref_window  = if (full_cohort) "5-year" else NULL)
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

      ## Width scales with the number of population facets; the label panel is
      ## a fixed cost, so a one-population figure must not be as wide as the
      ## three-population one or its squares stretch. Narrower since the
      ## single-pollutant panel was dropped and each population now
      ## contributes one column rather than three.
      width <- 6.0 + 1.6 * length(target$pops) +
        if (isTRUE(res$has_sp)) 2.4 else 0
      ggsave(
        filename = file.path(rev_dir("figures", "pathway"),
                             glue::glue("pathway_heatmap_{alg}{target$slug}.png")),
        plot = res$plot, width = width, height = 0.30 * res$n + 3.2,
        dpi = 300, bg = "white", limitsize = FALSE)
    })

  ## Window-comparison figure: the full cohort only, the three averaging
  ## windows side by side, rows restricted to pathways reaching p < 0.05 under
  ## at least one of them. No single-pollutant panel -- the eight pollutants
  ## are a separate question from the stability of the composite, and carrying
  ## them here was what made the earlier version of this figure read as a
  ## league table of exposures.
  ## min_cols = every window, not just one. A row here is a pathway that
  ## clears p < 0.05 under the 3-, 5- AND 10-year index, which is the point of
  ## the figure: what survives the choice of averaging window rather than what
  ## any one window happens to return.
  purrr::walk(c("mummichog", "metapone"), function(alg){
    res <- create_pathway_heatmap(pathway_long, alg, populations = "all",
                                  window_cols = PATHWAY_WINDOW_COLUMNS_3,
                                  min_cols = nrow(PATHWAY_WINDOW_COLUMNS_3))
    if (is.null(res)) {
      message("No pathway rows for the ", alg, " window heatmap")
      return(invisible(NULL))
    }
    message("Pathway window heatmap: ", alg, " - ", res$n, " pathways",
            if (res$n < res$n_available)
              glue::glue(" (capped from {res$n_available})") else "")
    ggsave(
      filename = file.path(rev_dir("figures", "pathway"),
                           glue::glue("pathway_heatmap_{alg}_windows.png")),
      ## Requiring every window leaves few rows, so the figure is narrowed
      ## and its rows given more height: at the wider size the squares were
      ## lost in white space.
      plot = res$plot, width = if (res$n <= 6) 6.4 else 8.2,
      height = (if (res$n <= 6) 0.55 else 0.30) * res$n + 3.0,
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
##      adjusted partial residuals against the exposure, with the fitted
##      linear slope
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
##     exist. `all predx` carries the 5-year pollutants too as of 2026-09-09
##     (PREDX_EXTRA_EXPOSURES in R1), so it now reads its own models rather
##     than falling back to `all`; the fallback remains for any population
##     that carries no single-pollutant models.
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
      ##
      ## The WINDOW IS STRIPPED from the tick labels. Every column of this
      ## panel is at the same window and the subtitle already names it, so
      ## carrying it on all eight labels turned them into "Benzene, 10-year
      ## window" repeated eight times at a 45-degree rotation -- the exact
      ## crowding Reviewer 1 minor comment 14 is about. The species name is
      ## what distinguishes the columns; the window is panel-level and belongs
      ## in the subtitle, once.
      pollutant = factor(rev_label_md(untagged_exposure(pollutant)),
                         levels = rev_label_md(untagged_exposure(pollutants))),
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
      ## Four breaks, not the default five: five labels ran into each other on
      ## the bar ("-10-5 0 5 10").
      n.breaks = 4,
      name = "% difference per IQR"
    ) +
    ## The bar joins the collected legend strip at the foot of the figure,
    ## beside Column and Significance, and is boxed like them (see
    ## `panel_legend_theme` below) so the three read as three separate keys
    ## rather than one run-on row. What tells a reader the tiles are in IQR
    ## units is the panel's own x-axis title, not the position of the bar.
    ##
    ## Title BESIDE the bar, not above it: the strip is one row and a stacked
    ## title makes this legend taller than the two it sits next to.
    guides(fill = guide_colourbar(
      title.position = "left", title.vjust = 0.85,
      theme = theme(legend.key.width  = unit(4, "cm"),
                    legend.key.height = unit(0.5, "cm")))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = "Single-pollutant models",
      ## Two lines. The asterisk key made this longer than the tiles are wide,
      ## and it ran into the forest's subtitle to its left -- the two halves
      ## are one panel, so there is no margin between them to absorb it.
      subtitle = paste0(pop_note, ", ", window_note, "\n",
                        "* P < 0.05, ** FDR < 0.10, *** FDR < 0.05"),
      ## This half of panel C carries its OWN x-axis title. The forest beside
      ## it declares "% difference per 1 SD (95% CI)" and the two halves use
      ## different increments -- an IQR is the same 25th-to-75th percentile
      ## shift for every pollutant whatever its skew, which is what makes the
      ## eight columns comparable with each other, whereas the index has no
      ## physical scale and only an SD to be expressed in. Naming the unit
      ## under the tiles it applies to is what keeps the two readable side by
      ## side; REV_TEXT_PANEL gives both the same 19 pt bold as panel D's.
      x = "% difference per IQR", y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.x   = element_text(angle = 45, hjust = 1),
      axis.line     = element_blank(),
      axis.ticks    = element_blank()
    ) +
    ## panel_legend_theme for the position and spacing the other collected
    ## legends use, but WITHOUT its border: a box drawn round a continuous
    ## colour bar crowds the gradient in a way it does not crowd the discrete
    ## keys beside it.
    panel_legend_theme +
    theme(legend.background = element_blank())
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
## exposure on the scale the coefficients are reported in.
##
## LINEAR FIT ONLY. The panel used to carry a loess curve beside the line, to
## show whether the linear term the MWAS fits is a fair summary of the
## relationship. At a quarter of a third of a row per facet the two curves
## read as one thick band rather than as a comparison, so the loess was
## dropped; the linear fit is the quantity the manuscript actually reports.
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
        "Red: linear fit"),
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
      labs(title = NULL) + panel_legend_theme
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
      labs(title = NULL) +
      ## The Column key is carried by the subgroup panel above. Both panels
      ## draw the same two shapes, and patchwork collects the guides into one
      ## row, so a second identical key here only makes that row too wide for
      ## the figure -- which is what clipped both ends of it.
      guides(shape = "none") + panel_legend_theme
  }

  purrr::compact(list(p2, p3)) |>
    purrr::map(rev_wrap_axis_titles) |>
    purrr::map(~ .x + REV_TEXT + REV_TEXT_PANEL +
                 ## REV_TEXT sets legend text at 13/14, which is right for a
                 ## single key and too wide for three collected into one row
                 ## under the figure: the last one lost its final label off
                 ## the right edge. Applied after REV_TEXT so it wins.
                 theme(legend.title = element_text(face = "bold", size = 11),
                       legend.text  = element_text(size = 10)))
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
  ## gradient is its own legend and joins the same strip, boxed to match.
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


# =============================================================================
# SECTION 8: MANTEL-STYLE COMPOSITE-TOXICANT NETWORK AND WEIGHT STABILITY
# =============================================================================

## The revision analogue of the Mantel-style network in
## scripts/7-visualization.R, which anchored the SUBMITTED (naive, outcome-
## informed) composites.
##
## Reviewer 1 comment 1 asks for two things that belong in one figure: the
## correlation among the revision exposures (unsupervised PC1, cross-fitted
## WQS, cross-fitted QGcomp, cross-fitted Cox QGcomp), and the UNCERTAINTY of
## the WQS / QGcomp weights rather than point estimates alone.
##
## Panel A anchors the network on the revision composites instead:
## toxicant-toxicant Pearson
## correlations in the upper triangle, curves from each composite to each
## toxicant (width = |r|, colour = p-value tier), and the composite-composite
## correlation matrix inset at lower left.
##
## Panel B is the fold-to-fold weight spread, which is the stability evidence:
## every cross-fitting fold contributes one weight per pollutant, so the SD
## across folds and -- more tellingly -- whether a pollutant CHANGES SIGN
## between folds says how identified the weighting is. WQS is sign-constrained
## and cannot flip; QGcomp is not, and does.

message("\n=== Section 8: composite-toxicant network ===")

if (!requireNamespace("linkET", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("Hy4m/linkET", upgrade = "never")
}

## PLAIN TEXT, unlike every other figure in the revision, which carries
## NO<sub>2</sub> / NO<sub>x</sub> / PM<sub>2.5</sub> through ggtext. linkET
## draws the qcorrplot's species labels itself rather than as ggplot axis
## text, so an element_markdown() on this plot is ignored and the raw
## "<sub>" tags print verbatim. Flat names are the lesser evil here; the
## subscripts survive everywhere the labels are ordinary axis text.
tox_label_lookup <- c(exp_benzene_iqr   = "Benzene",
                      exp_butadiene_iqr = "1,3-Butadiene",
                      exp_chromium_iqr  = "Chromium",
                      exp_nickel_iqr    = "Nickel",
                      exp_lead_iqr      = "Lead",
                      exp_nox_iqr       = "NOx",
                      `exp_pm2.5_iqr`   = "PM2.5",
                      exp_no2_iqr       = "NO2")

## Short labels throughout -- the linkET anchor labels sit outside the panel
## and long ones are clipped, and the inset needs to match them anyway.
mix_label_lookup <- c(comp_pca_all           = "PC1",
                      comp_wqs_cf_all        = "WQS-CF",
                      comp_qgcomp_cf_all     = "QGcomp-CF",
                      comp_qgcomp_cox_cf_all = "QGcomp Cox-CF")

## The MAIN network anchors PC1 at all three averaging windows instead of the
## four index types. The toxicant block behind the curves stays the 5-year
## matrix throughout -- the primary exposure window, and the only one all
## three anchors can be compared against on common nodes.
##
## The 3-year index is built from SEVEN pollutants, not eight: NOx is
## unavailable inside a 3-year window for most specimens (see R1). Its curve
## to NOx is therefore a correlation with a pollutant that did not enter it,
## which is worth knowing when that curve reads thinner than its neighbours.
## Labels kept SHORT. linkET clips the topmost anchor's label at the panel
## edge no matter how large the left plot margin is -- the margin moves the
## panel, not the label's offset from its node -- so "PC1 (3-yr)" lost its
## leading character while the two below it rendered in full.
## The window alone, without the "PC1" prefix: all three anchors ARE PC1, and
## panels A and B name it either side of this one, so the prefix bought
## nothing and cost the topmost label its first character.
## The SUPPLEMENT network drops PC1: it is the main figure's anchor, and the
## question this one answers is how the OUTCOME-INFORMED indices relate to
## the pollutants once cross-fitted.
index_label_lookup <- c(comp_wqs_cf_all        = "WQS-CF",
                        comp_qgcomp_cf_all     = "QGcomp-CF",
                        comp_qgcomp_cox_cf_all = "QGcomp Cox-CF")

window_label_lookup <- c(comp_pca_all_w3  = "3-year",
                         comp_pca_all     = "5-year",
                         comp_pca_all_w10 = "10-year")

poll_label_lookup <- c(benzene = "Benzene", butadiene = "1,3-Butadiene",
                       chromium = "Chromium", nickel = "Nickel",
                       lead = "Lead", nox = "NO<sub>x</sub>",
                       `pm2.5` = "PM<sub>2.5</sub>", no2 = "NO<sub>2</sub>")

## Pearson r between a composite and one pollutant, with the p-value of the
## test of H0: rho = 0.
##
## TWO CAVEATS, both of which the figure caption should carry.
##
## 1. The test treats the 1,546 specimens as independent, but they come from
##    952 participants and exposure is assigned at the address, so repeated
##    draws from one person carry almost the same exposure value. The
##    effective n is nearer 952 than 1,546 and these p-values are therefore
##    anticonservative. Re-running one row per participant leaves every r
##    within 0.03 and every p still below 1e-12, so no curve changes tier --
##    but the p is smaller than the data support.
##
## 2. More fundamentally, a composite is a linear combination of the same
##    pollutants it is being correlated with, so rho = 0 is not a hypothesis
##    anyone holds. The p-value tier is close to tautological here; it is the
##    MAGNITUDE of r, encoded as curve width, that carries the information.
##
## `use = "pairwise.complete.obs"` used to be passed here. cor.test() has no
## such argument -- it went into `...` and was silently discarded -- and it
## was never needed, since cor.test() drops incomplete pairs itself.
cor_one <- function(x, y){
  tst <- suppressWarnings(stats::cor.test(x, y, method = "pearson"))
  tibble::tibble(r = unname(tst$estimate), p = tst$p.value)
}

## Composite-composite bubble matrix. `compact = TRUE` strips it down so the
## same plot can be dropped into the network figure as an inset.
build_mix_bubble <- function(data, compact = FALSE){
  cor_levels <- unname(mix_label_lookup)
  cor_long <- data |>
    dplyr::select(dplyr::all_of(names(mix_label_lookup))) |>
    dplyr::rename_with(~ unname(mix_label_lookup[.x])) |>
    cor(use = "pairwise.complete.obs") |>
    tibble::as_tibble(rownames = "var_x") |>
    tidyr::pivot_longer(-var_x, names_to = "var_y",
                        values_to = "correlation") |>
    dplyr::mutate(var_x = factor(var_x, levels = cor_levels),
                  var_y = factor(var_y, levels = cor_levels))

  bubble_size <- if (compact) 6.5 else 18
  text_size   <- if (compact) 2.4 else 5.5
  base_size   <- if (compact) 8   else 15

  ggplot(cor_long, aes(var_x, var_y)) +
    geom_tile(fill = "white", col = "grey85") +
    geom_point(data = dplyr::filter(cor_long,
                                    as.integer(var_x) > as.integer(var_y)),
               aes(fill = correlation, size = abs(correlation)),
               colour = "black", shape = 21, stroke = 0.4) +
    geom_text(data = dplyr::filter(cor_long,
                                   as.integer(var_y) > as.integer(var_x)),
              aes(label = sprintf("%.2f", correlation)),
              colour = "grey15", size = text_size, fontface = "bold") +
    ## `compact` is the inset form -- no colour bar, no title, tiny text.
    ## The standalone figure needs the bar, or its bubbles decode to nothing.
    scale_fill_gradient2(low = "#436C85", mid = "white", high = "#B73F42",
                         midpoint = 0, limits = c(-1, 1),
                         breaks = seq(-1, 1, 0.5), name = "Pearson r",
                         guide = if (compact) "none" else "colourbar") +
    scale_size_area(limits = c(0, 1), max_size = bubble_size) +
    coord_cartesian(expand = FALSE) +
    labs(x = NULL, y = NULL,
         title = if (compact) "Composite-composite r" else NULL) +
    guides(size = "none",
           fill = if (compact) "none" else
             guide_colourbar(title.position = "top", title.hjust = 0.5,
                             barwidth = unit(7, "cm"),
                             barheight = unit(0.5, "cm"),
                             frame.colour = "grey40",
                             ticks.colour = "grey40")) +
    theme_minimal(base_size = base_size) +
    theme(plot.title      = element_text(face = "bold",
                                         size = if (compact) 8.5 else 14,
                                         hjust = 0.5, margin = margin(b = 2)),
          axis.text.y     = element_text(size = if (compact) 6 else 12,
                                         colour = "grey25"),
          axis.text.x     = element_text(size = if (compact) 6 else 12,
                                         colour = "grey25",
                                         angle = 30, vjust = 1, hjust = 1),
          legend.position = if (compact) "none" else "bottom",
          panel.grid      = element_blank(),
          plot.background = element_rect(fill = "white", colour = "grey60",
                                         linewidth = 0.5),
          plot.margin     = margin(6, 8, 6, 8))
}

# PCA structure: loadings and the clustering they imply -----------------------

## Built ONCE, outside the covariate loop: the unsupervised index is a
## decomposition of the pollutant matrix alone and does not know about
## covariates.
##
## This figure is what the response letter's claim rests on -- that PC1
## "loads positively on all eight components ... so it is interpretable as a
## single axis of overall exposure burden". A sentence asserting that is not
## evidence; the loadings and the geometry are.
##
##   A  PC1 loadings, every bar on the same side of zero, with the variance
##      PC1 explains named on the axis.
##   B  the PC1-PC2 plane: specimen scores behind, one arrow per pollutant.
##      Pollutants pointing the same way co-vary, which is what makes a
##      single axis a fair summary; the spread along PC2 is what that axis
##      discards.
##
## The three colours are hierarchical clustering (Ward) of the pollutants in
## the PC1-PC2 loading plane, cut at k = 3. It is a READING AID for the
## biplot -- it says which arrows group together -- and not a claim about
## emission sources. Anything source-attributing would need a receptor model
## and is not what this figure is for.

message("\n=== Section 8b: PCA loadings and structure ===")

pca_fit_all <- pca_index_list[["all"]]$pca_fit
pca_pve_all <- pca_index_list[["all"]]$pve

pve_lab <- function(i) paste0("PC", i, " (", round(100 * pca_pve_all[i], 1), "%)")

## make_pca_index() flips PC1 so that most loadings come out positive, but it
## flips only the vectors it returns -- the stored prcomp object keeps the
## sign prcomp chose. Reapply the same rule, or this figure disagrees with
## pca_loadings and with the scored index every other figure uses.
pc1_sign <- if (sum(pca_fit_all$rotation[, 1] < 0) >
                sum(pca_fit_all$rotation[, 1] > 0)) -1 else 1

poll_label_md <- c(benzene   = "Benzene",
                   butadiene = "1,3-Butadiene",
                   chromium  = "Chromium",
                   nickel    = "Nickel",
                   lead      = "Lead",
                   nox       = "NO<sub>x</sub>",
                   `pm2.5`   = "PM<sub>2.5</sub>",
                   no2       = "NO<sub>2</sub>")

## The biplot labels are PLOTMATH, not markdown. Three of the arrows are
## nearly parallel (nickel, PM2.5 and benzene all point right and slightly
## down), so their labels collide and need ggrepel to separate them --
## and ggrepel cannot render ggtext markdown. plotmath gives the same
## subscripts through `parse = TRUE`; names are quoted so the comma and
## hyphen in "1,3-Butadiene" survive parsing.
poll_label_pm <- c(benzene   = "\"Benzene\"",
                   butadiene = "\"1,3-Butadiene\"",
                   chromium  = "\"Chromium\"",
                   nickel    = "\"Nickel\"",
                   lead      = "\"Lead\"",
                   nox       = "NO[x]",
                   `pm2.5`   = "PM[2.5]",
                   no2       = "NO[2]")

pca_load_df <- tibble::tibble(
    pollutant = rownames(pca_fit_all$rotation) |> stringr::str_remove("^exp_"),
    PC1       = unname(pca_fit_all$rotation[, 1]) * pc1_sign,
    PC2       = unname(pca_fit_all$rotation[, 2])) |>
  dplyr::mutate(label    = unname(poll_label_md[pollutant]),
                label_pm = unname(poll_label_pm[pollutant]))

pca_hc <- stats::hclust(
  stats::dist(as.matrix(pca_load_df[, c("PC1", "PC2")])), method = "ward.D2")
pca_load_df$cluster <- factor(stats::cutree(pca_hc, k = 3))

CLUSTER_COLS <- c(`1` = "#B73F42", `2` = "#436C85", `3` = "#DE9960")

## Scores carry their ids, so the quartile-scored exposures can be joined on
## rather than bound by position. make_pca_index() builds `composites` from
## the same drop_na()'d frame it fits prcomp to, in the same row order, which
## is what makes the mutate below safe -- and the check that follows is what
## proves it rather than assuming it.
pca_scores_df <- pca_index_list[["all"]]$composites |>
  dplyr::mutate(PC1 = unname(pca_fit_all$x[, 1]) * pc1_sign,
                PC2 = unname(pca_fit_all$x[, 2])) |>
  dplyr::select(rand_id, blood_date, PC1, PC2)

stopifnot(nrow(pca_scores_df) == nrow(pca_fit_all$x))

## Panel A of the network figure: the SAME PC1-PC2 plane, once per pollutant,
## with the points coloured by that pollutant's quartile score. The quartile
## score is the right colour quantity because it is the matrix PCA actually
## decomposed, and because it puts all eight pollutants on one 0-3 scale, so a
## single legend serves every facet -- raw concentrations would need eight.
##
## One facet per pollutant rather than one plane coloured by whichever
## pollutant dominates: the four strongest contributors (PM2.5, benzene,
## 1,3-butadiene, nickel) are so highly intercorrelated that a
## dominant-pollutant colouring overplots them into an undifferentiated
## middle, and the separation a reader sees is then mostly lead and NOx. The
## facets show each pollutant's gradient across the plane on its own.
pca_vars <- pca_index_list[["all"]]$vars

pca_quartiles <- combined_data_list_revision[["total"]][["all"]][["covar"]] |>
  dplyr::select(rand_id, blood_date, dplyr::all_of(pca_vars)) |>
  tidyr::drop_na() |>
  quantile_score_matrix(pca_vars)

## Strip labels wrap, panel B's axis labels do not. "1,3-Butadiene" is the
## longest of the eight and overruns a facet strip at this panel width, but
## it fits on panel B's y axis, so only the strip variant is broken. The break
## is "<br>", not "\n": ggtext renders these through gridtext, which reads
## markdown, and a bare newline is not a line break in markdown.
poll_label_facet <- poll_label_md
poll_label_facet[["butadiene"]] <- "1,3-<br>Butadiene"

## Facets ordered by PC1 loading, so the panel reads in the same order as the
## loadings plot beside it and the strongest contributors come first.
pca_facet_order <- pca_load_df |>
  dplyr::arrange(dplyr::desc(PC1)) |>
  dplyr::pull(pollutant)

pca_facet_df <- pca_scores_df |>
  dplyr::inner_join(pca_quartiles, by = c("rand_id", "blood_date")) |>
  tidyr::pivot_longer(dplyr::all_of(pca_vars),
                      names_to = "pollutant", values_to = "quartile") |>
  dplyr::mutate(
    pollutant = stringr::str_remove(pollutant, "^exp_"),
    label     = factor(unname(poll_label_facet[pollutant]),
                       levels = unname(poll_label_facet[pca_facet_order])))

p_pca_facets <- ggplot(pca_facet_df, aes(PC1, PC2, colour = quartile)) +
  geom_point(size = 0.75, alpha = 0.75) +
  facet_wrap(~ label, nrow = 2) +
  scale_colour_gradientn(
    colours = c("#436C85", "#9DBBCB", "grey92", "#E39B7B", "#B73F42"),
    limits = c(0, 3), breaks = 0:3,
    name = "Quartile of exposure") +
  guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                  barwidth = unit(6, "cm"),
                                  barheight = unit(0.45, "cm"),
                                  frame.colour = "grey40",
                                  ticks.colour = "grey40")) +
  labs(x = pve_lab(1), y = pve_lab(2)) +
  ## strip.text is left as PLAIN text here and re-declared as markdown in
  ## rev_save_plot()'s `post`. REV_TEXT is applied first and sets a plain
  ## element_text; ggplot2 will not merge that over an element_markdown, so
  ## setting markdown at this point aborts the save instead of styling it.
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92", colour = "grey60"),
        strip.text       = element_text(face = "bold", size = 11),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())


p_pca_load <- pca_load_df |>
  dplyr::mutate(label = forcats::fct_reorder(label, PC1)) |>
  ggplot(aes(x = label, y = PC1, fill = cluster)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", PC1)),
            hjust = -0.25, size = 4, colour = "grey20") +
  scale_fill_manual(values = CLUSTER_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  coord_flip() +
  labs(x = NULL,
       y = paste0("Loading on ", pve_lab(1)),
       title = "PC1 loadings",
       subtitle = paste0("All eight pollutants load positively, so PC1 is an ",
                         "overall-burden contrast")) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, colour = "grey25"),
        axis.text.y   = ggtext::element_markdown(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

## Arrows are drawn in score units so both can share one pair of axes; the
## scale factor is cosmetic and the arrow LENGTHS are therefore comparable
## with each other but not with the point cloud.
arrow_scale <- 0.85 * max(abs(c(pca_scores_df$PC1, pca_scores_df$PC2))) /
  max(abs(c(pca_load_df$PC1, pca_load_df$PC2)))

p_pca_biplot <- ggplot(pca_scores_df, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(colour = "grey65", alpha = 0.28, size = 1.1) +
  geom_segment(data = pca_load_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale, yend = PC2 * arrow_scale,
                   colour = cluster),
               arrow = grid::arrow(length = unit(0.22, "cm"), type = "closed"),
               linewidth = 1.1, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = pca_load_df,
    aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale,
        label = label_pm, colour = cluster),
    parse = TRUE, size = 4.6, fontface = "bold", show.legend = FALSE,
    seed = 42, box.padding = 0.55, point.padding = 0.3,
    min.segment.length = Inf, max.overlaps = Inf,
    ## Push labels away from the origin so a label never lands on the arrow
    ## it belongs to, or on the point cloud in the middle.
    nudge_x = pca_load_df$PC1 * arrow_scale * 0.18,
    nudge_y = pca_load_df$PC2 * arrow_scale * 0.18) +
  scale_colour_manual(values = CLUSTER_COLS) +
  labs(x = pve_lab(1), y = pve_lab(2),
       title = "Specimens and pollutants in the PC1-PC2 plane",
       subtitle = paste0("Grey points: ", nrow(pca_scores_df),
                         " specimens.\nArrows: pollutant loadings, coloured ",
                         "by clustering in this plane (k = 3, a reading aid).")) +
  ## coord_equal keeps the angles between arrows honest, which is the whole
  ## point of a biplot; the expansion is what stops a repelled label at the
  ## edge from being clipped.
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.10)) +
  scale_y_continuous(expand = expansion(mult = 0.10)) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, colour = "grey25"),
        plot.margin      = margin(6, 14, 6, 6),
        panel.grid.minor = element_blank())

p_pca_structure <- (p_pca_load | p_pca_biplot) +
  patchwork::plot_layout(widths = c(1, 1.25)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

rev_save_plot(p_pca_structure, "pca_loadings_structure", "composites",
              width = 16, height = 7.5,
              post = theme(axis.text.y = ggtext::element_markdown(size = 13)))

rev_save_table(
  pca_load_df |>
    dplyr::transmute(pollutant, label_md = label,
                     pc1_loading = round(PC1, 4),
                     pc2_loading = round(PC2, 4),
                     cluster = as.integer(cluster)),
  "pca_loadings_pc1_pc2", "composites")

message("  PC1 explains ", round(100 * pca_pve_all[1], 1),
        "% and PC2 ", round(100 * pca_pve_all[2], 1), "% of the variance")


network_cor_tables <- covar_names |>
  purrr::set_names() |>
  purrr::map(function(covar_set){

    ## `comp_qgcomp_cox_cf_all` lives in the `cox` arm, which drops the 42
    ## prevalent dementia/CIND cases and the 2 without follow-up, so it joins
    ## in with NA on the remaining specimens and every correlation involving
    ## it is taken on pairwise-complete observations.
    df <- dplyr::full_join(
      combined_data_list_revision[["total"]][["all"]][[covar_set]] |>
        dplyr::select(rand_id, blood_date,
                      dplyr::any_of(names(tox_label_lookup)),
                      ## Raw windowed pollutants, so each windowed PC1 can be
                      ## correlated against the pollutants of ITS OWN window.
                      dplyr::matches("^exp_.*_(w3|w10)$"),
                      dplyr::any_of(names(mix_label_lookup)),
                      dplyr::any_of(names(window_label_lookup))),
      combined_data_list_revision[["cox"]][["all"]][[covar_set]] |>
        dplyr::select(rand_id, blood_date,
                      dplyr::any_of("comp_qgcomp_cox_cf_all")),
      by = c("rand_id", "blood_date"))

    tox_cols    <- df |>
      dplyr::select(dplyr::any_of(names(tox_label_lookup))) |> colnames()
    tox_renamed <- df |>
      dplyr::select(dplyr::all_of(tox_cols)) |>
      dplyr::rename_with(~ unname(tox_label_lookup[.x]))

    ## Couple data for an arbitrary set of anchors, so the index network and
    ## the window network are built by the same code and cannot drift.
    ##
    ## `col_for` resolves which pollutant column an anchor should be
    ## correlated against. The nodes are SPECIES, but a windowed index has to
    ## be compared with that species at its own averaging window -- a 3-year
    ## index against 5-year benzene would confound the window with the
    ## weighting, which is the comparison the panel exists to make.
    couple_data_for <- function(lookup, col_for = function(s, e) e){
      specs <- names(lookup) |> purrr::keep(~ .x %in% names(df))
      tidyr::expand_grid(spec_orig = specs, env_orig = tox_cols) |>
        dplyr::mutate(poll_col = purrr::map2_chr(spec_orig, env_orig, col_for)) |>
        ## Drops the 3-year NOx curve: NOx has no 3-year average for most
        ## specimens, so the column does not exist and no correlation is
        ## defined. Dropping the edge is honest; drawing it against 5-year
        ## NOx would not be.
        dplyr::filter(poll_col %in% names(df)) |>
        dplyr::mutate(stats = purrr::map2(spec_orig, poll_col,
                                          ~ cor_one(df[[.x]], df[[.y]]))) |>
        tidyr::unnest(stats) |>
        dplyr::transmute(
          spec = factor(unname(lookup[spec_orig]),
                        levels = unname(lookup[specs])),
          env  = factor(unname(tox_label_lookup[env_orig]),
                        levels = unname(tox_label_lookup[tox_cols])),
          r, p,
          sign = factor(ifelse(r < 0, "negative", "positive"),
                        levels = c("positive", "negative")),
          rd = cut(abs(r), breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                   labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", "> 0.6")),
          pd = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05",
                              "> 0.05")))
    }

    ## Which pollutant column a windowed PC1 anchor is compared against.
    window_pollutant_col <- function(spec_orig, env_orig){
      species <- env_orig |> stringr::str_remove("^exp_") |>
        stringr::str_remove("_iqr$")
      suffix <- switch(spec_orig,
                       comp_pca_all_w3  = "_w3",
                       comp_pca_all_w10 = "_w10",
                       "_iqr")
      paste0("exp_", species, suffix)
    }

    mix_tox_cor <- couple_data_for(mix_label_lookup)
    win_tox_cor <- couple_data_for(window_label_lookup, window_pollutant_col)
    idx_tox_cor <- couple_data_for(index_label_lookup)

    ## The network is built from whichever composites are passed as anchors,
    ## because it is drawn twice.
    ##
    ## MAIN FIGURE: PC1 at the three averaging windows. PC1 is the
    ## manuscript's primary exposure, so the comparison a reader needs is
    ## across ITS windows, not across index types. Three anchors keep 24
    ## curves, which stays legible where the four-index version's 32
    ## overlapped into an unreadable band.
    ##
    ## SUPPLEMENT: the four index types, for the reviewer's question about how
    ## the cross-fitted indices relate to the individual pollutants. The
    ## agreement AMONG index types is quantified more precisely still in
    ## composite_composite_correlation_* and pc1_vs_crossfit_composites_*.
    ## Curves are coloured BY ANCHOR, not by p-value.
    ##
    ## The p-value is gone because it was never informative here: a composite
    ## is a linear combination of the pollutants it is correlated with, so
    ## rho = 0 is not a hypothesis anyone holds, and every curve in the main
    ## figure fell in the same "< 0.001" tier. Colouring by anchor puts the
    ## channel to work instead -- it is what lets a reader follow one window's
    ## eight curves through the bundle, which is the overlap problem the
    ## single-colour version had.
    ##
    ## Curvature is assigned `by = "from"`, so each anchor's curves bow by a
    ## different amount and the three families separate rather than tracing
    ## the same arcs.
    ## mark_size is the r printed in each cell, label_size the anchor node
    ## labels ("3-year", "PC1", ...). Both default to what the two supplement
    ## networks use; the main-text figure passes larger values, since it is
    ## printed at 21 inches wide and its type was reading small against the
    ## rest of the figure.
    make_network <- function(couple_data, anchor_cols, legend_name,
                             mark_size = 2.5, label_size = 3.88) {
      linkET::qcorrplot(linkET::correlate(tox_renamed),
                        type = "upper", diag = FALSE) +
        linkET::geom_square(colour = "grey85", size = 0.3) +
        ## Print r in each cell as well as encoding it by square size and
        ## fill. mark = "" suppresses the significance stars: with n = 1,546
        ## every pairwise correlation is significant, so the stars would mark
        ## nothing and only crowd the cell.
        linkET::geom_mark(size = mark_size, colour = "grey15",
                          digits = 2, nsmall = 2,
                          mark = c("", "", "")) +
        linkET::geom_couple(aes(colour = spec, size = rd, linetype = sign),
                            data = couple_data,
                            label.size = label_size,
                            curvature = linkET::nice_curvature(0.16,
                                                               by = "from"),
                            alpha = 0.8) +
        ## Sign has to be shown, not left to the width. Window-matching turns
        ## the 10-year PC1's correlation with 10-year 1,3-butadiene NEGATIVE
        ## (r = -0.25) while the 5-year pair is +0.64 -- the butadiene
        ## surfaces at those two windows are close to uncorrelated with each
        ## other. A width-only encoding would draw that curve identically to
        ## a positive one of the same magnitude.
        scale_linetype_manual(values = c(positive = "solid",
                                         negative = "22"),
                              name = "Sign of r") +
        scale_fill_gradient2(low = "#436C85", mid = "white", high = "#B73F42",
                             midpoint = 0, limits = c(-1, 1),
                             breaks = seq(-1, 1, 0.5),
                             name = "Pairwise r\n(toxicants)") +
        scale_size_manual(values = c("< 0.2" = 0.4, "0.2 - 0.4" = 1.2,
                                     "0.4 - 0.6" = 2.2, "> 0.6" = 3.5),
                          ## drop = TRUE (the default) so each figure's
                          ## legend lists only the tiers it actually draws.
                          ## The window network has no |r| < 0.2 edge -- its
                          ## weakest is 10-year lead at 0.21 -- and a legend
                          ## key for a width that appears nowhere on the page
                          ## is just something for a reader to hunt for. The
                          ## index network does have them, and keeps the key.
                          name = "|r| with pollutant\n(curve width)") +
        scale_colour_manual(values = anchor_cols, name = legend_name,
                            drop = FALSE) +
        ## All four guides in one stack on the right.
        guides(colour   = guide_legend(order = 1, ncol = 1,
                                       override.aes = list(linewidth = 2.5)),
               size     = guide_legend(order = 2, ncol = 1,
                                       override.aes = list(colour = "grey35")),
               linetype = guide_legend(order = 3, ncol = 1,
                                       override.aes = list(colour = "grey35",
                                                           linewidth = 1.2)),
               fill     = guide_colorbar(order = 4,
                                         barwidth  = unit(0.6, "cm"),
                                         barheight = unit(3.8, "cm"),
                                         frame.colour = "grey40",
                                         ticks.colour = "grey40")) +
        ## The pollutant names STAY ON THE RIGHT.
        ##
        ## Moving them with scale_y_discrete(position = "left") replaces the
        ## scale linkET built, and that scale carries the REVERSED row order
        ## the upper triangle depends on -- the tiles came back scrambled
        ## relative to the diagonal and the curve endpoints no longer landed
        ## on their own nodes.
        ##
        ## The reference figure gets left-hand names a different way: they are
        ## not axis text at all there, they are node labels drawn along the
        ## diagonal. That would need geom_couple to label the env side, not an
        ## axis move. Worth doing if the placement matters, but it is a
        ## different mechanism, not a position argument.
        ##
        ## There is also a reason the right edge is the natural side for an
        ## UPPER triangle: every row ends at the right margin, so a
        ## right-hand label sits next to its own tiles, whereas a left-hand
        ## label for NO2 -- whose only tile is the last column -- would sit an
        ## entire matrix away from the data it names.
        ## Some ggplot2 versions carry linkET's raw x/y aesthetic names
        ## through as axis titles; the axes here are pollutant names.
        labs(x = NULL, y = NULL) +
        theme(axis.title      = element_blank(),
              legend.title    = element_text(face = "bold", size = 11),
              legend.text     = element_text(size = 10),
              legend.key.size = unit(0.55, "cm"),

              ## Left margin buys room for the anchor labels, which linkET
              ## draws outside the panel.
              plot.margin     = margin(10, 10, 10, 62))
    }

    WINDOW_COLS <- stats::setNames(c("#DE9960", "#436C85", "#7C9A4E"),
                                   unname(window_label_lookup))
    INDEX_COLS  <- stats::setNames(c("#B73F42", "#436C85", "#DE9960"),
                                   unname(index_label_lookup))
    ## PC1 plus the three cross-fitted indices. The three keep exactly the
    ## colours INDEX_COLS gives them, so an index is the same colour in both
    ## supplement figures; only PC1's is new.
    MIX_COLS    <- c(stats::setNames("#5B4B8A",
                                     unname(mix_label_lookup[["comp_pca_all"]])),
                     INDEX_COLS)

    p_network <- make_network(win_tox_cor, WINDOW_COLS,
                              "Averaging window\n(curve colour)",
                              mark_size = 3.6, label_size = 5.2)

    p_network_all <- make_network(idx_tox_cor, INDEX_COLS,
                                  "Composite index\n(curve colour)")

    ## Same network with PC1 added, for the combined supplement figure. The
    ## standalone `..._allindices_` figure stays PC1-free: there the question
    ## is how the outcome-informed indices behave, while here PC1 is the
    ## unsupervised reference the other three are being read against.
    p_network_mix <- make_network(mix_tox_cor, MIX_COLS,
                                  "Composite index\n(curve colour)")

    ## --- Panel B: PC1 against each cross-fitted composite -------------------
    ##
    ## The scatter row is the quantitative claim the network only gestures at:
    ## how closely each outcome-informed index, once cross-fitted, tracks the
    ## unsupervised one. R-squared is annotated per facet because a reader
    ## comparing panels wants the number, not an eyeballed slope.
    cf_scatter_levels <- c("comp_wqs_cf_all", "comp_qgcomp_cf_all",
                           "comp_qgcomp_cox_cf_all")

    scatter_df <- cf_scatter_levels |>
      purrr::keep(~ .x %in% names(df)) |>
      purrr::map(function(cf){
        tibble::tibble(composite = unname(mix_label_lookup[cf]),
                       pc1 = df[["comp_pca_all"]], y = df[[cf]])
      }) |>
      purrr::list_rbind() |>
      tidyr::drop_na() |>
      dplyr::mutate(composite = factor(composite,
                                       levels = unname(mix_label_lookup[cf_scatter_levels])))

    scatter_r2 <- scatter_df |>
      dplyr::group_by(composite) |>
      dplyr::summarise(r = stats::cor(pc1, y), .groups = "drop") |>
      dplyr::mutate(label = paste0("italic(r) == ", sprintf("%.2f", r)))

    p_pc_scatter <- ggplot(scatter_df, aes(pc1, y)) +
      geom_point(colour = "grey55", alpha = 0.25, size = 0.9) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = "#B73F42", linewidth = 0.9) +
      geom_text(data = scatter_r2, parse = TRUE,
                aes(x = -Inf, y = Inf, label = label),
                hjust = -0.25, vjust = 1.5, size = 4.4, colour = "grey15") +
      ## One facet per COLUMN: each composite is on its own scale, and
      ## stacking them puts the three y axes one above another where a reader
      ## compares them, rather than side by side where the x axis repeats.
      facet_wrap(~ composite, ncol = 1, scales = "free_y") +
      labs(x = "PC1 (unsupervised index)",
           y = "Cross-fitted composite") +
      theme_bw(base_size = 12) +
      theme(strip.background = element_rect(fill = "grey92", colour = "grey60"),
            strip.text       = element_text(face = "bold", size = 11),
            panel.grid.minor = element_blank())

    ## --- Panel B: fold-to-fold weight spread --------------------------------
    cf_levels <- c("comp_wqs_cf_all", "comp_qgcomp_cf_all",
                   "comp_qgcomp_cox_cf_all")

    weights_5yr <- crossfit_weights |>
      dplyr::filter(covar_set == !!covar_set, window == "w5",
                    composite %in% cf_levels) |>
      dplyr::mutate(composite = factor(unname(mix_label_lookup[composite]),
                                       levels = unname(mix_label_lookup[cf_levels])),
                    pollutant = unname(poll_label_lookup[pollutant]))

    ## `flips` is the headline: TRUE when a pollutant's weight is positive in
    ## one fold and negative in another, i.e. the composite is not identified
    ## strongly enough to fix even the DIRECTION of that pollutant.
    weight_spread <- weights_5yr |>
      dplyr::group_by(composite, pollutant) |>
      dplyr::summarise(mean_w = mean(weight), sd_w = stats::sd(weight),
                       min_w = min(weight), max_w = max(weight),
                       flips = length(unique(sign(weight))) > 1,
                       .groups = "drop")

    poll_order <- weight_spread |>
      dplyr::group_by(pollutant) |>
      dplyr::summarise(m = mean(mean_w), .groups = "drop") |>
      dplyr::arrange(m)

    weights_5yr   <- dplyr::mutate(weights_5yr,
                                   pollutant = factor(pollutant,
                                                      levels = poll_order$pollutant))
    weight_spread <- dplyr::mutate(weight_spread,
                                   pollutant = factor(pollutant,
                                                      levels = poll_order$pollutant))

    flip_cols <- c(`FALSE` = "#436C85", `TRUE` = "#B73F42")
    flip_labs <- c("stable sign", "sign flips across folds")

    p_weights <- ggplot(weight_spread, aes(pollutant, mean_w)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
      geom_linerange(aes(ymin = mean_w - sd_w, ymax = mean_w + sd_w,
                         colour = flips), linewidth = 0.9) +
      geom_point(data = weights_5yr, aes(y = weight), colour = "grey35",
                 size = 1.1, alpha = 0.75,
                 position = position_nudge(x = 0.22)) +
      geom_point(aes(fill = flips), shape = 21, size = 3, stroke = 0.5,
                 colour = "black") +
      facet_wrap(~ composite, ncol = 1, scales = "free_x") +
      coord_flip() +
      scale_colour_manual(values = flip_cols, labels = flip_labs, name = NULL) +
      scale_fill_manual(values = flip_cols, labels = flip_labs, name = NULL) +
      labs(x = NULL,
           y = "Mixture weight (mean ± SD across folds)",
           title = "Fold-to-fold stability of the cross-fitted mixture weights",
           subtitle = paste0("All-toxicant groupings, ", covar_set,
                             ". Grey points: one per cross-fitting fold.\n",
                             "Red: the weight changes sign between folds, so ",
                             "the composite does not fix\neven the direction ",
                             "in which that pollutant enters it.")) +
      theme_bw(base_size = 12) +
      theme(plot.title       = element_text(face = "bold", size = 14),
            plot.subtitle    = element_text(size = 11, colour = "grey25"),
            strip.background = element_rect(fill = "grey92", colour = "grey60"),
            strip.text       = element_text(face = "bold", size = 11),
            legend.position  = "bottom",
            panel.grid.minor = element_blank())

    ## THE MAIN FIGURE: PCA coordinates over PC1-versus-composite scatters on
    ## the left, the toxicant network on the right. The left column reads top
    ## to bottom as one argument -- here is the unsupervised exposure space,
    ## and here is how closely each cross-fitted index reproduces it -- and
    ## the network then places both against the individual pollutants.
    ## The left column reads as one argument about the unsupervised index:
    ## the exposure plane, then the weight each pollutant carries in the axis
    ## that plane is built on. Two panels, not three -- the network is the
    ## figure's main display and takes close to two thirds of the width.
    ## The PC1-loadings bar chart that used to sit under the exposure plane is
    ## gone from this figure: the exposure-performance panel below plots the
    ## same loadings against the adjusted R2 of each surface, so the bar chart
    ## added nothing. It keeps its own standalone figure (pca_loadings.png).

    ## The distribution row of panel A and panel B answer the second half of
    ## Reviewer 1 comment 13, which
    ## asks for the pollutant distributions and the performance of each
    ## exposure surface. They flank the network because they describe the
    ## same eight components the network relates to each other: what each
    ## component looks like, and how well each is modelled.

    ## Distributions: the eight components span six orders of magnitude, so every facet
    ## takes its own axis. The panel is for the SHAPE and the relative spread
    ## -- which is where the narrow within-Sacramento contrast shows up -- not
    ## for comparing absolute levels between pollutants.
    ## Butadiene is stored in ppt; displayed in ppb so the axis reads 5-7
    ## rather than 5000-7000, and to match the unit the companion paper uses.
    dist_rescale <- c(butadiene = 1e-3)          # ppt -> ppb
    dist_units <- c(benzene = "ppb", butadiene = "ppb",
                    chromium = "&mu;g/m<sup>3</sup>",
                    nickel   = "&mu;g/m<sup>3</sup>",
                    lead     = "&mu;g/m<sup>3</sup>",
                    nox      = "ppb",
                    `pm2.5`  = "&mu;g/m<sup>3</sup>", no2 = "ppb")
    ## Name and unit together: this is now the only panel identifying each
    ## component, so the strip has to carry both.
    dist_labs <- stats::setNames(
      paste0(poll_label_md, "<br>(", dist_units[names(poll_label_md)], ")"),
      names(poll_label_md))

    dist_df <- combined_data_list_revision[["total"]][["all"]][[covar_set]] |>
      dplyr::select(dplyr::any_of(paste0("exp_", names(poll_label_md)))) |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "pollutant",
                          values_to = "value") |>
      dplyr::mutate(
        pollutant = stringr::str_remove(pollutant, "^exp_"),
        value     = value * dplyr::coalesce(unname(dist_rescale[pollutant]), 1),
        pollutant = factor(pollutant, levels = pca_facet_order))

    p_dist <- ggplot(dist_df, aes(x = value, y = "")) +
      geom_violin(fill = "#436C85", colour = NA, alpha = 0.45, width = 0.95) +
      geom_boxplot(width = 0.17, outlier.size = 0.25, outlier.alpha = 0.30,
                   linewidth = 0.32, fill = "white", colour = "grey25") +
      facet_wrap(~ pollutant, scales = "free_x", nrow = 2,
                 labeller = ggplot2::labeller(pollutant = dist_labs)) +
      scale_x_continuous(
        breaks = function(r) {
          n <- if (max(abs(r)) < 0.01) 2 else 3
          signif(seq(r[1] + 0.18 * diff(r), r[2] - 0.18 * diff(r),
                     length.out = n), 2)
        },
        labels = scales::label_number(drop0trailing = TRUE)) +
      ## C and D are deliberately small; keep their type proportionate

      ## Short enough for one line at this width. The window is per-draw and
      ## evaluated at whichever address was current then; that detail lives in
      ## the Methods and the figure caption rather than on the axis.
      labs(x = "Mean of the five calendar years preceding each blood draw",
           y = NULL) +
      theme_bw(base_size = 11) +
      theme(axis.text.x   = element_text(size = 7.5),
            strip.text    = element_text(size = 8.5, lineheight = 1.05),
            axis.title.x  = element_text(size = 12, face = "bold"),
            panel.spacing.x = unit(0.9, "lines"),
            axis.text.y   = element_blank(),
            axis.ticks.y  = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor   = element_blank())

    ## Model performance: adjusted R2 as published for these same surfaces -- Yu et al.,
    ## Environ Res 2025 (doi:10.1016/j.envres.2025.123105), Table S3 -- against
    ## the weight each pollutant carries in PC1. Plotting them together is the
    ## point: the Reviewer's closing argument is that measurement error differs
    ## across components and propagates into the index, and the two axes are
    ## exactly those two quantities. CALINE4 NOx is a Gaussian dispersion
    ## model, not a fitted regression, so it has no R2 and is annotated rather
    ## than plotted.
    lur_r2 <- tibble::tribble(
      ~pollutant,   ~adj_r2,
      "benzene",    0.81,
      "butadiene",  0.62,
      "chromium",   0.76,
      "lead",       0.59,
      "nickel",     0.70,
      "no2",        0.84,
      "pm2.5",      0.65)

    r2_df <- dplyr::inner_join(pca_load_df, lur_r2, by = "pollutant")

    p_r2 <- ggplot(r2_df, aes(adj_r2, PC1)) +
      geom_point(shape = 21, size = 4, stroke = 0.5,
                 fill = "#436C85", colour = "black") +
      ggrepel::geom_text_repel(aes(label = label_pm), parse = TRUE,
                               size = 5.4, colour = "grey20",
                               box.padding = 0.45, min.segment.length = 0,
                               segment.colour = "grey70", seed = 42) +
      annotate("text", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.45,
               size = 4.8, colour = "grey45", lineheight = 1.05,
               label = paste0("NOx (CALINE4) is a dispersion model, not a\n",
                              "fitted surface, so it has no R\u00b2")) +
      scale_x_continuous(limits = c(0.5, 0.95), breaks = seq(0.5, 0.9, 0.1)) +
      expand_limits(y = 0) +
      ## bold() has to be inside the expression: theme(face = "bold") styles a
      ## plain string but does not reach into plotmath, so this title rendered
      ## lighter than the y-axis title beside it.
      labs(x = expression(bold("Adjusted"~R^2~"of the exposure surface (published)")),
           y = "Loading on PC1") +
      theme_bw(base_size = 11) +
      theme(axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.text.x  = element_text(size = 9),
            axis.text.y  = element_text(size = 9),
            panel.grid.minor = element_blank())

    ## The PC1/PC2 exposure-plane grid is no longer part of this figure -- it
    ## restated what the network and the loading axis already show. It is still
    ## built above and still has its own figure (pca_loadings_structure.png).
    ##
    ## One trap worth recording: do NOT wrap the left column in
    ## wrap_elements(). rev_save_plot applies its themes with `&`, which does
    ## not descend into a wrapped element, and the ggtext markdown strips
    ## silently revert to raw "<sub>" text. Also parenthesise the composition:
    ## `a | b + plot_layout(...)` binds the layout to b alone, because `+` has
    ## higher precedence than `|` in R.
    ## TYPE SIZES ARE SET PER PANEL HERE, NOT GLOBALLY.
    ##
    ## rev_save_plot applies REV_TEXT (and any `post`) with `&`, which reaches
    ## every panel equally, and this figure needs them to differ: panels A and
    ## B were reading small at 21 inches, while panel C's pollutant names and
    ## its four-part legend are already at the right size and grow into the
    ## curves if enlarged. So each panel carries REV_TEXT plus its own sizes,
    ## and the composition adds only the tag, which nothing else sets.
    ##
    ## Panel A's facet strips are species names and need markdown. Adding
    ## element_markdown AFTER REV_TEXT is what makes it work: REV_TEXT sets
    ## strip.text as plain text and ggplot2 merges markdown over plain text
    ## but not the reverse.
    ## Panel A's type is bounded by its column width, not by taste: at the
    ## old 1 : 2.2 split its eight facets were about 1.6 inches apiece, where
    ## 18 pt tick labels run together, a 16 pt strip clips "1,3-Butadiene",
    ## and an 18 pt axis title is wider than the column and loses its first
    ## word off the left edge. The left column is therefore widened below,
    ## which panel C can afford — its curves sit in a lot of white space.
    pA <- p_dist + REV_TEXT +
      theme(## 14, not 15: panel B's larger y-axis numbers and title widen its
            ## label gutter, and patchwork aligns the two panels in this
            ## column, so panel A's facets lost the width that let a 15 pt
            ## strip hold "1,3-Butadiene" without clipping.
            strip.text      = ggtext::element_markdown(face = "bold",
                                                       size = 14),
            axis.text.x     = element_text(size = 15),
            ## Same size as panel B's axis titles, set just below.
            axis.title.x    = element_text(size = 18, face = "bold"),
            panel.spacing.x = unit(1.1, "lines"),
            ## The axis title is nearly as wide as the column; a little left
            ## margin keeps it off the edge of the page.
            plot.margin     = margin(5, 5, 5, 14))

    ## The .x / .y elements have to be named individually. p_r2 sets
    ## axis.text.x, axis.text.y, axis.title.x and axis.title.y explicitly, and
    ## a specific element beats a generic `axis.text` / `axis.title` added
    ## afterwards -- setting the generic ones here left the panel at its
    ## original 9 pt numbers and 12 pt titles.
    pB <- p_r2 + REV_TEXT +
      theme(axis.text.x  = element_text(size = 16),
            axis.text.y  = element_text(size = 16),
            axis.title.x = element_text(size = 18, face = "bold"),
            axis.title.y = element_text(size = 18, face = "bold"))

    ## Panel C keeps REV_TEXT's axis.text (the pollutant names) and its
    ## legend sizes unchanged; only the in-cell r values and the node labels
    ## grew, and those are geom sizes set on make_network above.
    pC <- p_network + REV_TEXT

    p_left <- (pA / pB) + patchwork::plot_layout(heights = c(1, 1))

    p_main <- (p_left | pC) +
      patchwork::plot_layout(widths = c(1.4, 2.2)) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 32))

    network_png <- file.path(rev_dir("figures", "composites"),
                             paste0("composite_toxicant_network_",
                                    covar_set, ".png"))
    ggplot2::ggsave(network_png, p_main, width = 21, height = 10.5,
                    dpi = 300, bg = "white")
    message("  figure -> ", network_png)

    rev_save_plot(p_network_all,
                  paste0("composite_toxicant_network_allindices_", covar_set),
                  "composites", width = 13, height = 9.5)

    ## The PC1-versus-composite scatters are no longer a panel of the main
    ## figure, but they are the only place the agreement is quantified, so
    ## they keep a figure of their own.
    rev_save_plot(p_pc_scatter,
                  paste0("pc1_vs_crossfit_composites_", covar_set),
                  "composites", width = 7, height = 10,
                  post = theme(strip.text = element_text(face = "bold",
                                                         size = 12)))

    ## The composite-composite matrix is its own figure now. It was an inset
    ## in the corner of the network, where it was too small to read and too
    ## easily taken for part of the network's own colour scale.
    ## No panel title. It named the covariate set with the internal label
    ## ("covar" / "covar_sen"), which means nothing to a reader; the caption
    ## says which set the figure shows.
    rev_save_plot(build_mix_bubble(df, compact = FALSE),
                  paste0("composite_composite_correlation_", covar_set),
                  "composites", width = 8, height = 8)

    rev_save_plot(p_weights,
                  paste0("crossfit_weight_uncertainty_", covar_set),
                  "composites", width = 9, height = 11,
                  post = theme(axis.text.y = ggtext::element_markdown(size = 14)))

    ## The two supplement displays as ONE figure. They answer the same
    ## question from two directions -- what the cross-fitted indices are
    ## correlated with, and how reproducible the weights behind them are --
    ## so a reader who has one wants the other on the same page.
    ##
    ## The weight panel's title and subtitle are dropped HERE ONLY. The
    ## subtitle explains the grey points and the red markers, which a
    ## standalone figure needs and a captioned two-panel figure does not; the
    ## title would be the only one in a two-panel layout, and at REV_TEXT's
    ## 18 pt it is wider than the panel and clips. The standalone version
    ## keeps both.
    p_supp <- (p_network_mix |
                 (p_weights + labs(title = NULL, subtitle = NULL))) +
      patchwork::plot_layout(widths = c(1.55, 1)) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 20))

    rev_save_plot(p_supp,
                  paste0("crossfit_network_and_weights_", covar_set),
                  "composites", width = 20, height = 10.5,
                  post = theme(axis.text.y = ggtext::element_markdown(size = 14)))

    ## The numbers behind the figure, so the response letter can quote them.
    mix_mix_cor <- df |>
      dplyr::select(dplyr::all_of(names(mix_label_lookup))) |>
      dplyr::rename_with(~ unname(mix_label_lookup[.x])) |>
      cor(use = "pairwise.complete.obs") |>
      tibble::as_tibble(rownames = "composite_a") |>
      tidyr::pivot_longer(-composite_a, names_to = "composite_b",
                          values_to = "r") |>
      dplyr::filter(composite_a != composite_b) |>
      dplyr::mutate(covar_set = covar_set, r = round(r, 3), .before = 1)

    list(mix_mix       = mix_mix_cor,
         mix_tox       = mix_tox_cor |>
           dplyr::transmute(covar_set = covar_set, composite = spec,
                            pollutant = env, r = round(r, 3), p = signif(p, 3)),
         weight_spread = weight_spread |>
           dplyr::mutate(covar_set = covar_set, .before = 1) |>
           dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                       ~ round(.x, 4))))
  })

composite_network_correlations <- network_cor_tables |>
  purrr::map("mix_mix") |> purrr::list_rbind()
composite_toxicant_correlations <- network_cor_tables |>
  purrr::map("mix_tox") |> purrr::list_rbind()
crossfit_weight_spread <- network_cor_tables |>
  purrr::map("weight_spread") |> purrr::list_rbind()

rev_save_table(composite_network_correlations,
               "composite_network_correlations", "composites")
rev_save_table(composite_toxicant_correlations,
               "composite_toxicant_correlations", "composites")
rev_save_table(crossfit_weight_spread, "crossfit_weight_spread", "composites")

message("Pollutants whose cross-fitted weight changes sign between folds:")
crossfit_weight_spread |>
  dplyr::filter(covar_set == "covar", flips) |>
  dplyr::select(composite, pollutant, mean_w, sd_w, min_w, max_w) |>
  print(n = 20)


# 9. Participant and specimen flow diagram (R1 comment 2) --------------------- ## FLOWSTART
##
## Reviewer 1 comment 2 asks for a participant and sample flow diagram, and for
## the specimen structure behind the pooled repeated-measures MWAS to be made
## explicit: how many specimens each participant contributed and from which
## visit wave, and how the post-diagnosis specimens are handled.
##
## DESIGN. A CONSORT/STROBE participant-flow layout: one vertical spine of
## cohort boxes, with exclusions hanging off to the right, so the reading path
## is a single top-to-bottom line and every side box is visibly a subtraction.
## Portrait canvas, narrow spine, restrained palette -- a neutral ground, one
## structural blue for the cohort boxes, warm grey for exclusions, and a single
## green accent reserved for the two analysis populations.
##
## This is a PARTICIPANT-INCLUSION chart, not an analysis chart: it says who
## contributed which specimens, not how the models were fitted. Model details
## belong in the Methods.
##
## The incident (Cox-eligible) cohort is deliberately absent. The Cox-weighted
## composite is a methodological diagnostic rather than a reported analysis
## (see the reply to comment 1), so giving it a branch here would imply a role
## in the findings it no longer has.
##
## Every count is recomputed from the analysis frames rather than transcribed,
## so the figure cannot drift away from the pipeline.

load(here::here("data", "processed", "salsa_clean.RData"))
load(rev_here("data", "processed", "combined_data_list_revision.RData"))

flow_total <- combined_data_list_revision[["total"]][["all"]][["covar"]]
flow_cox   <- combined_data_list_revision[["cox"]][["all"]][["covar"]]

## post-diagnosis specimens, on the same rule R4 uses for the `all predx`
## population: prevalent cases contribute no pre-diagnosis specimen at all,
## incident cases contribute the draws taken after baseline + dcst years.
flow_dx <- salsa_clean_cox |>
  dplyr::distinct(rand_id, .keep_all = TRUE) |>
  dplyr::select(rand_id, bl_date, dcst, demcind_incident = demcind) |>
  dplyr::mutate(dx_date = dplyr::if_else(
    demcind_incident == "Dementia/CIND",
    as.Date(bl_date) + as.numeric(dcst) * 365.25, as.Date(NA)))

flow_flagged <- flow_total |>
  dplyr::left_join(dplyr::select(flow_dx, rand_id, dx_date, demcind_incident),
                   by = "rand_id") |>
  dplyr::mutate(
    prevalent = is.na(demcind_incident) & demcind == "Dementia/CIND",
    postdx    = prevalent |
      (!is.na(dx_date) & as.Date(blood_date) > dx_date))

flow_n <- list(
  enrolled      = 1789L,
  n_part_total  = dplyr::n_distinct(flow_total$rand_id),
  n_spec_total  = nrow(flow_total),
  n_part_cox    = dplyr::n_distinct(flow_cox$rand_id),
  n_spec_cox    = nrow(flow_cox),
  n_events      = sum(salsa_clean_cox$demcind[!duplicated(salsa_clean_cox$rand_id)] ==
                        "Dementia/CIND"),
  n_spec_postdx = sum(flow_flagged$postdx),
  n_part_postdx = dplyr::n_distinct(flow_flagged$rand_id[flow_flagged$postdx]),
  n_spec_prev   = sum(flow_flagged$prevalent),
  n_part_prev   = dplyr::n_distinct(flow_flagged$rand_id[flow_flagged$prevalent])
)
flow_n$n_excluded    <- flow_n$enrolled - flow_n$n_part_total
flow_n$n_spec_predx  <- flow_n$n_spec_total - flow_n$n_spec_postdx
flow_n$n_part_predx  <- dplyr::n_distinct(flow_flagged$rand_id[!flow_flagged$postdx])
flow_n$n_spec_incpdx <- flow_n$n_spec_postdx - flow_n$n_spec_prev
flow_n$n_part_incpdx <- flow_n$n_part_postdx - flow_n$n_part_prev

## Cognitive status of the 952 analytic participants, for the flow figure:
## prevalent cases (dementia/CIND at the baseline visit), incident cases
## (diagnosed during follow-up) and those intact throughout. The three sum to
## the analytic sample, which is what reconciles Table 1's n = 141 with the
## 908-participant incident cohort (Reviewer 1 comment 3).
flow_n$n_part_prevalent <- flow_n$n_part_prev
flow_n$n_part_incident  <- flow_n$n_events
flow_n$n_part_intact    <- flow_n$n_part_total - flow_n$n_part_prevalent -
  flow_n$n_part_incident
stopifnot(flow_n$n_part_prevalent + flow_n$n_part_incident +
            flow_n$n_part_intact == flow_n$n_part_total)

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

draws_per_participant <- flow_total |>
  dplyr::count(rand_id, name = "n_draws") |>
  dplyr::count(n_draws, name = "n_participants")

spec_per_wave <- flow_total |>
  dplyr::mutate(year = lubridate::year(blood_date)) |>
  dplyr::group_by(wave) |>
  dplyr::summarise(n_specimens = dplyr::n(),
                   yr_min = min(year), yr_max = max(year), .groups = "drop")


## ---- palette -------------------------------------------------------------
## Neutral ground, one structural blue, warm grey for subtractions, a single
## green reserved for the analysis populations. Four hues, no more.
FLOW_PAL <- list(
  cohort_fill = "#E9F0F7", cohort_line = "#40658B",
  drop_fill   = "#F4F1ED", drop_line   = "#9C8F82",
  note_fill   = "#F8F8F6", note_line   = "#C6C6C0",
  main_fill   = "#DFEBE2", main_line   = "#3C6B4D",
  ink         = "#1F2933", ink_soft    = "#4A5560",
  rule        = "#7A8794"
)

## ---- geometry ------------------------------------------------------------
## Spine centred at x = 34 (width 50); exclusions at x = 79 (width 36).
SPINE_X <- 34; SPINE_W <- 50
SIDE_X  <- 79; SIDE_W  <- 36

flow_boxes <- tibble::tribble(
  ~id,       ~x,      ~y,   ~w,      ~h, ~kind, ~tag, ~head, ~body,

  "enrol",   SPINE_X, 95,   SPINE_W,  8, "cohort", NA,
  "SALSA cohort at enrolment",
  paste0(fmt_n(flow_n$enrolled), " participants\n",
         "Sacramento Valley, 1998-1999"),

  "excl1",   SIDE_X,  90.0, SIDE_W,   9, "drop", NA,
  "Excluded",
  paste0("No plasma specimen assayed\nby LC-HRMS\n",
         fmt_n(flow_n$n_excluded), " participants"),

  "metab",   SPINE_X, 80.5, SPINE_W, 12.6, "cohort", NA,
  "Plasma metabolomics analytic sample",
  paste0(fmt_n(flow_n$n_spec_total), " specimens from ",
         fmt_n(flow_n$n_part_total), " participants\n(dementia/CIND: ",
         fmt_n(flow_n$n_part_prevalent), " prevalent at baseline, ",
         fmt_n(flow_n$n_part_incident), " incident in follow-up;\n",
         fmt_n(flow_n$n_part_intact), " cognitively intact throughout)\n",
         "C18-negative and HILIC-positive run on every specimen"),

  "struct",  SPINE_X, 57,   SPINE_W, 28, "note", NA,
  NA, NA,

  "primary", SPINE_X, 34,   SPINE_W, 11, "main", "PRIMARY",
  "Pooled repeated-measures analysis",
  paste0(fmt_n(flow_n$n_spec_total), " specimens, ",
         fmt_n(flow_n$n_part_total), " participants\n",
         "Unsupervised PC1 air-toxicant index"),

  "excl2",   SIDE_X,  21,   SIDE_W,   9, "drop", NA,
  "Excluded",
  paste0(fmt_n(flow_n$n_spec_postdx), " post-diagnosis specimens\nfrom ",
         fmt_n(flow_n$n_part_postdx), " participants\n(",
         fmt_n(flow_n$n_spec_prev), " prevalent, ",
         fmt_n(flow_n$n_spec_incpdx), " incident)"),

  "sens",    SPINE_X, 9,    SPINE_W, 11, "main", "SENSITIVITY",
  "Pre-diagnosis specimens only",
  paste0(fmt_n(flow_n$n_spec_predx), " specimens, ",
         fmt_n(flow_n$n_part_predx), " participants\n",
         "Same exposure and model as the primary")
)

## ---- specimen-structure panel contents -----------------------------------
## Two labelled blocks of right-aligned labels and bold left-aligned counts,
## stacked rather than side by side so the spine stays narrow. This is the box
## that has to carry the most numbers, so it gets real alignment instead of a
## wrapped sentence.
STRUCT_TOP <- 57 + 28 / 2           # box top edge
SUB_X <- SPINE_X - SPINE_W / 2 + 4  # subheads, outdented
LAB_X <- SPINE_X - SPINE_W / 2 + 8  # row labels, left aligned
VAL_X <- SPINE_X + SPINE_W / 2 - 9  # counts, right aligned
STRUCT_TITLE_Y <- STRUCT_TOP - 2.6

## One row per printed line, with the vertical step that precedes it, so the
## panel is laid out by content rather than by hand-tuned coordinates. The
## wave subhead says the waves overlap in calendar time, because the timeline
## figure counts the same specimens by year of draw and the two groupings
## therefore do not match row for row.
struct_rows <- dplyr::bind_rows(
  tibble::tibble(kind = "subhead", label = "Draws per participant",
                 value = NA, gap = 3.1),
  draws_per_participant |>
    dplyr::transmute(kind = "row",
                     label = paste0(n_draws, dplyr::if_else(n_draws > 1,
                                                            " draws", " draw")),
                     value = fmt_n(n_participants), gap = 1.9),
  tibble::tibble(kind = "subhead",
                 label = "Specimens by visit wave (waves overlap in time)",
                 value = NA, gap = 3.1),
  spec_per_wave |>
    dplyr::transmute(kind = "row",
                     label = paste0("Wave ", wave, " (", yr_min,
                                    dplyr::if_else(yr_max > yr_min,
                                                   paste0("-", yr_max), ""), ")"),
                     value = fmt_n(n_specimens), gap = 1.9)
) |>
  dplyr::mutate(y = (STRUCT_TITLE_Y - 1.2) - cumsum(.data$gap),
                lab_x = LAB_X, val_x = VAL_X, sub_x = SUB_X)

## ---- drawing -------------------------------------------------------------
## Rounded boxes via grid, one grob per box: ggplot2 has no rounded rect and
## sharp corners read as heavier than this palette wants.
box_grob <- function(fill, colour) {
  grid::roundrectGrob(
    r  = grid::unit(2.4, "pt"),
    gp = grid::gpar(fill = fill, col = colour, lwd = 1.1))
}

flow_box_layers <- flow_boxes |>
  purrr::pmap(function(x, y, w, h, kind, ...) {
    ggplot2::annotation_custom(
      box_grob(FLOW_PAL[[paste0(kind, "_fill")]],
               FLOW_PAL[[paste0(kind, "_line")]]),
      xmin = x - w / 2, xmax = x + w / 2,
      ymin = y - h / 2, ymax = y + h / 2)
  })

## Connectors: spine segments end just above the next box, side connectors run
## from the spine out to the exclusion box.
flow_arrows <- tibble::tribble(
  ~x,       ~y,   ~xend,             ~yend,
  SPINE_X,  91.0, SPINE_X,           86.9,   # enrol  -> metab
  SPINE_X,  74.1, SPINE_X,           71.6,   # metab  -> structure
  SPINE_X,  43.0, SPINE_X,           40.1,   # struct -> primary
  SPINE_X,  28.5, SPINE_X,           15.1    # primary-> sensitivity
)

flow_side <- tibble::tribble(
  ~x,      ~y,   ~xend,               ~yend,
  SPINE_X, 90.0, SIDE_X - SIDE_W / 2, 90.0,
  SPINE_X, 21.0, SIDE_X - SIDE_W / 2, 21.0
)

has_body <- flow_boxes |> dplyr::filter(!is.na(body))
has_tag  <- flow_boxes |> dplyr::filter(!is.na(tag))

flow_plot <- ggplot2::ggplot() +
  flow_box_layers +
  ggplot2::geom_segment(
    data = flow_arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45, colour = FLOW_PAL$rule,
    arrow = grid::arrow(length = grid::unit(0.20, "cm"), type = "closed")) +
  ggplot2::geom_segment(
    data = flow_side,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45, colour = FLOW_PAL$rule,
    arrow = grid::arrow(length = grid::unit(0.20, "cm"), type = "closed")) +
  ## accent tag (PRIMARY / SENSITIVITY)
  ggplot2::geom_text(
    data = has_tag,
    ggplot2::aes(x = x, y = y + h / 2 - 2.2, label = tag),
    fontface = "bold", size = 2.75, colour = FLOW_PAL$main_line) +
  ## box heading (the structure panel draws its own title, so it has none here)
  ggplot2::geom_text(
    data = dplyr::filter(flow_boxes, !is.na(head)),
    ggplot2::aes(x = x,
                 y = dplyr::if_else(is.na(tag), y + h / 2 - 2.3,
                                    y + h / 2 - 5.0),
                 label = head),
    fontface = "bold", size = 3.5, colour = FLOW_PAL$ink) +
  ## box body
  ggplot2::geom_text(
    data = has_body,
    ggplot2::aes(x = x,
                 y = dplyr::if_else(is.na(tag), y + h / 2 - 4.1,
                                    y + h / 2 - 6.8),
                 label = body),
    vjust = 1, size = 3.0, colour = FLOW_PAL$ink_soft, lineheight = 1.25) +
  ## specimen-structure panel: centred title, outdented subheads, then a
  ## two-column table with the counts right aligned on a common edge
  ggplot2::annotate("text", x = SPINE_X, y = STRUCT_TITLE_Y,
                    label = "Specimen structure", fontface = "bold",
                    size = 3.5, colour = FLOW_PAL$ink) +
  ggplot2::geom_text(
    data = dplyr::filter(struct_rows, kind == "subhead"),
    ggplot2::aes(x = sub_x, y = y, label = label),
    hjust = 0, fontface = "bold", size = 2.95, colour = FLOW_PAL$ink) +
  ggplot2::geom_text(
    data = dplyr::filter(struct_rows, kind == "row"),
    ggplot2::aes(x = lab_x, y = y, label = label),
    hjust = 0, size = 2.9, colour = FLOW_PAL$ink_soft) +
  ggplot2::geom_text(
    data = dplyr::filter(struct_rows, kind == "row"),
    ggplot2::aes(x = val_x, y = y, label = value),
    hjust = 1, fontface = "bold", size = 2.9, colour = FLOW_PAL$ink) +
  ggplot2::coord_cartesian(xlim = c(1, 99), ylim = c(2.2, 100),
                           expand = FALSE, clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.margin     = ggplot2::margin(8, 8, 8, 8),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA))

rev_dir("figures", "cohort")
ggplot2::ggsave(rev_here("figures", "cohort", "participant_specimen_flow.png"),
                flow_plot, width = 7.5, height = 9.8, dpi = 400, bg = "white")
ggplot2::ggsave(rev_here("figures", "cohort", "participant_specimen_flow.pdf"),
                flow_plot, width = 7.5, height = 9.8, bg = "white")

## The counts behind the diagram, so the letter and the figure share a source.
## The incident-cohort rows are kept even though the figure no longer draws
## that arm: the response letter still reports them when answering the
## Reviewer's questions about the Cox weight derivation.
flow_counts_table <- tibble::tibble(
  quantity = c("Enrolled in SALSA",
               "Excluded: no LC-HRMS plasma specimen",
               "Participants with metabolomics",
               "Specimens with metabolomics",
               "Specimens: post-diagnosis (any)",
               "Participants contributing a post-diagnosis specimen",
               "Specimens: post-diagnosis, prevalent case",
               "Specimens: post-diagnosis, incident case",
               "Specimens: pre-diagnosis population",
               "Participants: pre-diagnosis population",
               "Specimens: incident (Cox-eligible) cohort",
               "Participants: incident (Cox-eligible) cohort",
               "Incident dementia/CIND events",
               "Rows entering the Cox weight model",
               "Participants: prevalent case at baseline",
               "Participants: incident case during follow-up",
               "Participants: cognitively intact throughout"),
  n = c(flow_n$enrolled, flow_n$n_excluded, flow_n$n_part_total,
        flow_n$n_spec_total, flow_n$n_spec_postdx, flow_n$n_part_postdx,
        flow_n$n_spec_prev, flow_n$n_spec_incpdx, flow_n$n_spec_predx,
        flow_n$n_part_predx, flow_n$n_spec_cox, flow_n$n_part_cox,
        flow_n$n_events, flow_n$n_part_cox,
        flow_n$n_part_prevalent, flow_n$n_part_incident,
        flow_n$n_part_intact))

rev_save_table(flow_counts_table, "participant_specimen_flow", "cohort")
rev_save_table(draws_per_participant, "draws_per_participant", "cohort")
rev_save_table(spec_per_wave, "specimens_per_wave", "cohort")

message("Participant and specimen flow diagram written to ",
        rev_here("figures", "cohort"))
## FLOWEND


# 10. Pre-diagnosis versus full-sample agreement (R1 comment 2) --------------- ## PREDXSTART
##
## Reviewer 1 comment 2 asks whether specimens drawn after a dementia/CIND
## diagnosis were retained, and asks for them to be excluded or separately
## evaluated. The `all predx` population is that separate evaluation, and this
## section quantifies what it shows.
##
## The count of features passing FDR is NOT the statistic to compare on -- the
## reply to comment 1 shows that a few per cent of attenuation moves the count
## by an order of magnitude on this data set's very dense Benjamini-Hochberg
## boundary. What is reported instead, per exposure and per platform, is the
## correlation of the feature-level coefficients, the sign concordance among
## the features the full sample calls significant, where those features land in
## the pre-diagnosis FDR ranking, and the attenuation of the t statistics.

predx_agreement_for <- function(results_list, platform, covar_set = "covar") {
  full  <- results_list[["total"]][["all"]][[covar_set]]
  predx <- results_list[["total"]][["all predx"]][[covar_set]]

  intersect(names(full), names(predx)) |>
    purrr::map(function(exp_name) {
      a <- full[[exp_name]]
      b <- predx[[exp_name]]
      if (is.null(a) || is.null(b)) return(NULL)
      a$feature <- rownames(a)
      b$feature <- rownames(b)
      m <- dplyr::inner_join(a, b, by = "feature",
                             suffix = c("_all", "_pre"))
      sig <- m |> dplyr::filter(adj.P.Val_all < 0.05)
      has_sig <- nrow(sig) > 0

      tibble::tibble(
        platform     = platform,
        covar_set    = covar_set,
        exposure     = exp_name,
        n_fdr05_all  = nrow(sig),
        n_fdr05_pre  = sum(m$adj.P.Val_pre < 0.05, na.rm = TRUE),
        coef_r       = round(stats::cor(m$logFC_all, m$logFC_pre), 3),
        t_r          = round(stats::cor(m$t_all, m$t_pre), 3),
        ## through-origin slope: the attenuation of the effect sizes
        slope        = round(unname(stats::coef(
                         stats::lm(logFC_pre ~ 0 + logFC_all, data = m))[1]), 3),
        pct_same_sign = if (has_sig) round(100 * mean(
                          sign(sig$logFC_all) == sign(sig$logFC_pre)), 1)
                        else NA_real_,
        n_still_fdr10 = if (has_sig) sum(sig$adj.P.Val_pre < 0.10) else NA_integer_,
        n_still_fdr25 = if (has_sig) sum(sig$adj.P.Val_pre < 0.25) else NA_integer_,
        median_t_ratio = if (has_sig) round(
                           stats::median(abs(sig$t_pre)) /
                           stats::median(abs(sig$t_all)), 3) else NA_real_
      )
    }) |>
    purrr::compact() |>
    purrr::list_rbind()
}

predx_agreement <- tidyr::expand_grid(
    platform  = c("C18", "HILIC"),
    covar_set = names(covar_list)
  ) |>
  purrr::pmap(function(platform, covar_set) {
    L <- if (identical(platform, "C18")) mwas_results_list_c18
         else mwas_results_list_hilic
    predx_agreement_for(L, platform, covar_set)
  }) |>
  purrr::list_rbind() |>
  dplyr::mutate(label = rev_label(exposure), .after = exposure)

rev_save_table(predx_agreement, "predx_vs_all_agreement", "mwas")

## ---- coefficient scatter: full sample versus pre-diagnosis -----------------
##
## The agreement table above is the evidence; this is the picture of it. Each
## feature is a point: its coefficient in the full sample against its
## coefficient in the pre-diagnosis population, with the two platforms shown
## side by side so a reader compares them without flipping between figures.
##
## FORM. A coefficient-versus-coefficient scatter is the right form for
## "do two estimates of the same quantity agree" -- the question is about the
## joint distribution of two continuous estimates, and the reference is a line
## rather than a level. The only reference drawn is the dashed 1:1 identity:
## the question the figure answers is whether the two estimates agree, and a
## fitted slope invites reading a small departure from 1 as a finding when on
## C18 it is fitted almost entirely on null coefficients. The slope is still
## computed and kept in predx_vs_all_agreement.xlsx for anyone who wants it.
##
## COLOUR. Two marks doing different jobs, not two peer categories. Features the
## full sample calls FDR < 0.05 are the subject and take the accent; every other
## feature is context and takes a recessive grey. The palette validator flags
## that grey on lightness and chroma, which is correct for a categorical palette
## and wrong for this one -- the grey is deliberately recessive. What matters is
## separability, and the pair clears it with room (normal-vision dE 30.3, worst
## CVD dE 27.8). The contrast warning is relieved by the printed statistics in
## each panel and by predx_vs_all_agreement.xlsx.
##
## Scales are free per panel: a C18 coefficient and a HILIC coefficient are on
## different scales, and forcing a shared axis would compress one platform into
## the middle of the other's range for no gain -- the claim here is about the
## slope within a panel, not about magnitudes across platforms.
PREDX_ACCENT <- "#2F6B8F"
PREDX_MUTED  <- "#B8BEC4"

predx_scatter_for <- function(results_by_platform, covar_set = "covar",
                              exposures = NULL, ncol = 2,
                              subtitle_extra = NULL) {

  dat <- results_by_platform |>
    purrr::imap(function(results_list, platform) {
      full  <- results_list[["total"]][["all"]][[covar_set]]
      predx <- results_list[["total"]][["all predx"]][[covar_set]]
      keep  <- intersect(names(full), names(predx))
      if (!is.null(exposures)) keep <- intersect(keep, exposures)
      if (length(keep) == 0) return(NULL)

      keep |>
        purrr::map(function(exp_name) {
          a <- full[[exp_name]]; b <- predx[[exp_name]]
          a$feature <- rownames(a); b$feature <- rownames(b)
          dplyr::inner_join(a, b, by = "feature",
                            suffix = c("_all", "_pre")) |>
            dplyr::transmute(exposure = exp_name, platform = platform,
                             feature, x = logFC_all, y = logFC_pre,
                             sig = adj.P.Val_all < 0.05)
        }) |>
        purrr::list_rbind()
    }) |>
    purrr::compact() |>
    purrr::list_rbind()

  if (nrow(dat) == 0) return(NULL)

  ## Panel order is exposure-major, platform-minor, so the two platforms of one
  ## exposure always sit next to each other.
  exp_levels <- if (is.null(exposures)) unique(dat$exposure) else
    intersect(exposures, unique(dat$exposure))
  plat_levels <- names(results_by_platform)

  ## The platform is named in the strip only when both are shown; with one
  ## platform the title already says which, and the suffix is just noise.
  panel_name <- function(exposure, platform) {
    lab <- unname(rev_label(as.character(exposure)))
    if (length(plat_levels) > 1) paste0(lab, "  |  ", platform) else lab
  }

  panel_levels <- tidyr::expand_grid(
      exposure = factor(exp_levels, levels = exp_levels),
      platform = factor(plat_levels, levels = plat_levels)) |>
    dplyr::mutate(panel = panel_name(exposure, platform)) |>
    dplyr::pull(panel)

  dat <- dat |>
    dplyr::mutate(panel = factor(panel_name(exposure, platform),
                                 levels = panel_levels))

  ## per-panel statistics, printed in the panel rather than in a legend
  stats <- dat |>
    dplyr::group_by(panel) |>
    dplyr::summarise(
      r     = stats::cor(x, y),
      n_sig = sum(sig),
      same  = if (sum(sig) > 0) mean(sign(x[sig]) == sign(y[sig])) else NA_real_,
      .groups = "drop") |>
    dplyr::mutate(label = paste0(
      "r = ", sprintf("%.3f", r),
      dplyr::if_else(n_sig > 0,
                     paste0("\n", n_sig, " at FDR < 0.05, ",
                            sprintf("%.0f%%", 100 * same), " same sign"),
                     "\nno feature at FDR < 0.05")))

  ## Label anchored to the panel's own corner, with a little padding, so it
  ## clears the point cloud in panels whose spread fills the frame.
  lab_pos <- dat |>
    dplyr::group_by(panel) |>
    dplyr::summarise(x = min(x) - 0.04 * diff(range(x)),
                     y = max(y) + 0.10 * diff(range(y)), .groups = "drop") |>
    dplyr::left_join(dplyr::select(stats, panel, label), by = "panel")

  ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey88") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey88") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22",
                         linewidth = 0.5, colour = "grey45") +
    ## context first, subject on top, so the accent is never buried
    ggplot2::geom_point(data = dplyr::filter(dat, !sig),
                        colour = PREDX_MUTED, size = 0.75, alpha = 0.45) +
    ggplot2::geom_point(data = dplyr::filter(dat, sig),
                        colour = PREDX_ACCENT, size = 1.5, alpha = 0.9) +
    ggplot2::geom_text(
      data = lab_pos,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = 3.2, lineheight = 1.2, colour = "grey25") +
    ggplot2::facet_wrap(~ panel, scales = "free", ncol = ncol) +
    ggplot2::labs(
      x = "Coefficient, full sample (1,546 specimens)",
      y = "Coefficient, pre-diagnosis only (1,408 specimens)",
      title = paste0("Full-sample versus pre-diagnosis coefficients",
                     if (length(results_by_platform) == 1)
                       paste0(": ", names(results_by_platform)) else ""),
      subtitle = paste0(
        "Blue = FDR < 0.05 in the full sample. Dashed line is 1:1.\n",
        "Covariate set: ", covar_set,
        if (!is.null(subtitle_extra)) paste0("\n", subtitle_extra) else "")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25,
                                               colour = "grey92"),
      strip.background = ggplot2::element_rect(fill = "grey96",
                                               colour = "grey85"),
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(colour = "grey35", size = 9,
                                               lineheight = 1.2))
}

rev_dir("figures", "mwas", "predx")

predx_platforms <- list(C18   = mwas_results_list_c18,
                        HILIC = mwas_results_list_hilic)

## THE FULL-SAMPLE ARM IS THE POOLED PRIMARY ANALYSIS.
##
## The manuscript's primary analysis is Rubin-pooled over m = 10 (R10), and at
## m = 10 it calls 106 HILIC features at FDR < 0.05. The single-imputation fit
## calls 107, so without this swap the panel annotation contradicts the text.
##
## Only the full-sample arm is swapped. R10 varies the exposure and the
## covariate set, never the population, so the pre-diagnosis subset has no
## pooled counterpart and stays single-imputation; the subtitle says so rather
## than leaving a reader to assume both sides were fitted the same way. The
## claim the panel supports -- that the features significant in the primary
## analysis keep their sign when post-diagnosis specimens are dropped -- is
## unaffected, since it is the primary analysis that defines the feature set.
swap_pooled_primary <- function(results_list, platform_key,
                                covar_set = "covar",
                                exposure = "comp_pca_all") {
  pooled_path <- rev_here("tables", "imputation", "mwas_pooled.xlsx")
  if (!file.exists(pooled_path)) {
    warning("No pooled estimates at ", pooled_path,
            " -- predx figure falls back to single imputation.")
    return(results_list)
  }
  pooled <- readxl::read_xlsx(pooled_path) |>
    dplyr::filter(.data$platform == platform_key) |>
    dplyr::transmute(feature = .data$met, p_logFC = .data$qbar,
                     p_t = .data$t_pooled, p_P = .data$p_pooled,
                     p_fdr = .data$fdr_pooled)

  cell <- results_list[["total"]][["all"]][[covar_set]][[exposure]]
  if (is.null(cell) || nrow(cell) == 0) return(results_list)

  out <- tibble::rownames_to_column(cell, "feature") |>
    dplyr::left_join(pooled, by = "feature")
  matched <- sum(!is.na(out$p_P))
  if (matched == 0) stop("predx: no pooled estimates matched for ", platform_key)
  out <- out |>
    dplyr::mutate(logFC     = dplyr::coalesce(.data$p_logFC, .data$logFC),
                  t         = dplyr::coalesce(.data$p_t, .data$t),
                  P.Value   = dplyr::coalesce(.data$p_P, .data$P.Value),
                  adj.P.Val = dplyr::coalesce(.data$p_fdr, .data$adj.P.Val)) |>
    dplyr::select(-p_logFC, -p_t, -p_P, -p_fdr) |>
    as.data.frame()
  rownames(out) <- out$feature
  out$feature <- NULL
  message("  predx: ", platform_key, " full-sample arm pooled (",
          matched, " features, ", sum(out$adj.P.Val < 0.05),
          " at FDR < 0.05)")
  results_list[["total"]][["all"]][[covar_set]][[exposure]] <- out
  results_list
}

predx_platforms_pooled <- list(
  C18   = swap_pooled_primary(mwas_results_list_c18,   "c18"),
  HILIC = swap_pooled_primary(mwas_results_list_hilic, "hilic"))

## The primary exposure, both platforms side by side, for the letter.
## Drawn through si_agreement_scatter so it matches the other SI coefficient
## scatters (R1-revision_functions.R); the older bespoke version is kept above
## for the multi-exposure grids, where a four-level legend in twenty panels
## would cost more than it explains.
predx_si_data <- function(results_by_platform, pooled_tbl,
                          covar_set = "covar", exposure = "comp_pca_all") {
  plat_key <- c("C18" = "c18", "HILIC" = "hilic")
  plat_lab <- c("C18" = "C18/neg-", "HILIC" = "HILIC/pos+")
  results_by_platform |>
    purrr::imap(function(results_list, platform) {
      a <- results_list[["total"]][["all"]][[covar_set]][[exposure]]
      b <- results_list[["total"]][["all predx"]][[covar_set]][[exposure]]
      if (is.null(a) || is.null(b)) return(NULL)
      a$feature <- rownames(a); b$feature <- rownames(b)
      n_sig <- sum(a$adj.P.Val < 0.05)
      same  <- a$feature[a$adj.P.Val < 0.05]
      dplyr::inner_join(a, b, by = "feature", suffix = c("_all", "_pre")) |>
        dplyr::transmute(
          feature,
          x = .data$logFC_all, y = .data$logFC_pre,
          column_type = unname(plat_lab[platform]),
          panel = unname(plat_lab[platform]),
          class = si_scatter_class(.data$P.Value_all, .data$P.Value_pre,
                                   "Full sample P < 0.05",
                                   "Pre-diagnosis P < 0.05"),
          fdr_all = .data$adj.P.Val_all,
          sign_agree = sign(.data$logFC_all) == sign(.data$logFC_pre))
    }) |>
    purrr::compact() |>
    purrr::list_rbind()
}

predx_dat <- predx_si_data(predx_platforms_pooled)

## Sign concordance among the features the primary analysis calls significant,
## which is the claim the Results paragraph makes; printed in the subtitle
## because stat_cor already occupies the corner of each facet.
predx_note <- predx_dat |>
  dplyr::filter(.data$fdr_all < 0.05) |>
  dplyr::group_by(.data$panel) |>
  dplyr::summarise(n = dplyr::n(), same = sum(.data$sign_agree),
                   .groups = "drop") |>
  dplyr::mutate(txt = paste0(.data$panel, ": ", .data$same, " of ", .data$n,
                             " FDR-significant features keep their sign")) |>
  dplyr::pull(.data$txt) |>
  paste(collapse = "  \u00b7  ")

si_agreement_scatter(
  predx_dat,
  xlab = "Coefficient, full sample (1,546 specimens)",
  ylab = "Coefficient, pre-diagnosis only (1,408 specimens)",
  title = "Full-sample versus pre-diagnosis coefficients",
  subtitle = paste0(
    predx_note, "\n",
    "Full sample: Rubin-pooled over m = 10 imputations.",
    " Pre-diagnosis: single completed dataset.")) |>
  ggplot2::ggsave(filename = rev_here("figures", "mwas", "predx",
                                      "predx_vs_all_pc1.png"),
                  width = 13, height = 7, dpi = 400, bg = "white")


## Cross-fitted WQS against the unsupervised PC1 index, ALL features.
##
## The generic logfc_* comparison figure restricts to features with P < 0.05
## under the reference exposure, which reports r = 0.987 on HILIC. The
## manuscript quotes r = 0.945, which is the correlation across the whole
## feature set, so this figure draws all of them and the two now agree.
wqs_si_data <- function(results_by_platform, covar_set = "covar") {
  plat_lab <- c("C18" = "C18/neg-", "HILIC" = "HILIC/pos+")
  results_by_platform |>
    purrr::imap(function(results_list, platform) {
      cell <- results_list[["total"]][["all"]][[covar_set]]
      a <- cell[["comp_pca_all"]]; b <- cell[["comp_wqs_cf_all"]]
      if (is.null(a) || is.null(b)) return(NULL)
      a$feature <- rownames(a); b$feature <- rownames(b)
      dplyr::inner_join(a, b, by = "feature", suffix = c("_pca", "_wqs")) |>
        dplyr::transmute(
          feature,
          x = .data$logFC_pca, y = .data$logFC_wqs,
          column_type = unname(plat_lab[platform]),
          panel = unname(plat_lab[platform]),
          class = si_scatter_class(.data$P.Value_pca, .data$P.Value_wqs,
                                   "PC1 P < 0.05", "WQS-CF P < 0.05"),
          fdr_pca = .data$adj.P.Val_pca,
          sign_agree = sign(.data$logFC_pca) == sign(.data$logFC_wqs))
    }) |>
    purrr::compact() |>
    purrr::list_rbind()
}

wqs_dat <- wqs_si_data(predx_platforms_pooled)

wqs_note <- wqs_dat |>
  dplyr::filter(.data$fdr_pca < 0.05) |>
  dplyr::group_by(.data$panel) |>
  dplyr::summarise(n = dplyr::n(), same = sum(.data$sign_agree),
                   .groups = "drop") |>
  dplyr::mutate(txt = paste0(.data$panel, ": ", .data$same, " of ", .data$n,
                             " FDR-significant features keep their sign")) |>
  dplyr::pull(.data$txt) |>
  paste(collapse = "  \u00b7  ")

si_agreement_scatter(
  wqs_dat,
  xlab = "Coefficient, unsupervised PC1 index",
  ylab = "Coefficient, cross-fitted WQS index",
  title = "Unsupervised versus cross-fitted outcome-informed index",
  subtitle = paste0(
    wqs_note, "\n",
    "All features. PC1: Rubin-pooled over m = 10 imputations.",
    " WQS-CF: single completed dataset.")) |>
  ggplot2::ggsave(filename = rev_here("figures", "mwas", "predx",
                                      "wqs_cf_vs_pc1_all_features.png"),
                  width = 13, height = 7, dpi = 400, bg = "white")

## The full exposure set: composites and the eight single pollutants. One
## figure per platform -- twenty panels on a single sheet would shrink each
## below the size at which a slope is readable.
predx_grid_exposures <- intersect(
  c("comp_pca_all", "comp_wqs_cf_all", single_pollutant_exposures),
  names(mwas_results_list_hilic[["total"]][["all predx"]][["covar"]]))

predx_platforms |>
  purrr::iwalk(function(L, platform) {
    p <- predx_scatter_for(stats::setNames(list(L), platform), "covar",
                           exposures = predx_grid_exposures, ncol = 4)
    if (is.null(p)) return(invisible(NULL))
    ggplot2::ggsave(
      rev_here("figures", "mwas", "predx",
               paste0("predx_vs_all_grid_", tolower(platform), ".png")),
      p, width = 12, height = 9.5, dpi = 400, bg = "white")
  })

message("Pre-diagnosis coefficient scatters written to ",
        rev_here("figures", "mwas", "predx"))

message("Pre-diagnosis versus full-sample agreement, HILIC / primary covariates:")
predx_agreement |>
  dplyr::filter(platform == "HILIC", covar_set == "covar") |>
  dplyr::select(exposure, n_fdr05_all, n_fdr05_pre, coef_r, pct_same_sign,
                n_still_fdr25, median_t_ratio) |>
  print(n = 30)
## PREDXEND


message("\nVisualization completed!")
message("Figures saved to:")
message("  - ", rev_here("figures", "mwas"))
message("  - ", rev_here("figures", "pathway"))
message("  - ", rev_here("figures", "composites"),
        " (composite scatters from R3; network + weight stability from Section 8)")

#--------------------------------End of the code--------------------------------
