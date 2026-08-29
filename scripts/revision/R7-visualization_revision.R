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
##
##        Figure tree mirrors the primary pipeline:
##        revision_output/figures/mwas/{study}/{population}/{covar_set}/{exposure}/
##
## ---------------------------

# Load required packages -----------------------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

pacman::p_load("ggh4x")

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
  exposure_vars_list <- exposure_vars_list |> purrr::map(~ head(.x, 2))
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
            exposure_vars_list[[study]] |>
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


# =============================================================================
# SECTION 1: VOLCANO PLOTS (Combined C18 + HILIC)
# =============================================================================

# Function to create combined volcano plot ------------------------------------

create_volcano_plot <- function(mwas_result_c18, mwas_result_hilic,
                                annotation_result_c18, annotation_result_hilic,
                                exposure_name,
                                fdr_threshold = 0.05,
                                n_labels = 10) {

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
      ## VIP-prioritised and FDR-significant features. The annotated frames
      ## carry VIP_comp1, so the volcano uses the same three classes as the
      ## Manhattan instead of collapsing the middle one.
      significant = case_when(
        adj.P.Val < fdr_threshold ~ "FDR < 0.05",
        P.Value < 0.05 & !is.na(VIP_comp1) & VIP_comp1 > 2 ~ "P < 0.05 & VIP > 2",
        P.Value < 0.05 ~ "P < 0.05",
        TRUE ~ "NS"
      ),
      significant = factor(
        significant,
        levels = c("FDR < 0.05", "P < 0.05 & VIP > 2", "P < 0.05", "NS"))
    )

  # Get top metabolites for labeling (by p-value)
  top_mets <- plot_data |>
    dplyr::filter(!is.na(compound), compound != "", P.Value < 0.05) |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = n_labels) |>
    dplyr::pull(met)

  plot_data <- plot_data |>
    dplyr::mutate(
      # Keep only the first compound when multiple matches are separated by ";"
      compound_first = stringr::str_squish(
        stringr::str_split_i(compound, ";", 1)
      ),
      label = ifelse(met %in% top_mets, compound_first, "")
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
    scale_color_manual(
      values = c("FDR < 0.05"         = "#BE3F42",
                 "P < 0.05 & VIP > 2" = "#436C85",
                 "P < 0.05"           = "#DE9960"),
      name = "Significance", drop = FALSE
    ) +
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
    geom_label_repel(
      data = ~ dplyr::filter(.x, label != ""),
      aes(label = label),
      size = 4.5, max.overlaps = Inf,
      box.padding = 0.6, point.padding = 0.3,
      min.segment.length = 0, segment.color = "grey50",
      show.legend = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40") +
    labs(
      title = rev_label(exposure_name),
      x = "Log Fold Change",
      y = expression(-log[10](P-value))
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
            exposure_vars_list[[study]] |>
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
        significant = case_when(
          adj.P.Val < 0.05 ~ "FDR < 0.05",
          P.Value < 0.05 & !is.na(VIP_comp1) & VIP_comp1 > 2 ~ "P < 0.05 & VIP > 2",
          P.Value < 0.05 ~ "P < 0.05",
          TRUE ~ "NS"
        ),
        significant = factor(significant,
                             levels = c("FDR < 0.05", "P < 0.05 & VIP > 2",
                                        "P < 0.05", "NS")),
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
      n_fdr = sum(adj.P.Val < 0.05),
      n_p_vip = sum(P.Value < 0.05 & !is.na(VIP_comp1) & VIP_comp1 > 2),
      n_nom = sum(P.Value < 0.05),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = paste0("FDR < 0.05: ", n_fdr,
                     "\nP < 0.05 & VIP > 2: ", n_p_vip,
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
    scale_color_manual(
      values = c("FDR < 0.05" = "#B73F42",
                 "P < 0.05 & VIP > 2" = "#436C85",
                 "P < 0.05" = "#DE9960"),
      name = "Significance", drop = FALSE
    ) +
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
      y = expression(-log[10](P-value))
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
            exposure_vars_list[[study]] |>
              purrr::walk(function(exp) {
                message("Manhattan: ", study, "_", pop, " [", cov, "] - ", exp)

                rev_ggsave(
                  filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                      glue::glue("manhattan_{exp}.png")),
                  ## combined_results_list_*, not mwas_results_list_*: only
                  ## the former carries VIP_comp1. Fed the topTable output the
                  ## Manhattan filled VIP with NA and could never draw the
                  ## "P < 0.05 & VIP > 2" class, which is 456 (C18) and 425
                  ## (HILIC) features under the primary exposure.
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
    c("comp_pca_all", "comp_qgcomp_cf_all")
  ),
  ## The cox arm carries no unsupervised index (comp_pca_* is fitted in the
  ## total cohort only), so there is nothing to compare it against here. Its
  ## cross-fitted-versus-naive contrast comes from Section 4 instead.
  cox = list()
)

study_names |>
  purrr::walk(function(study){
    populations_for(study) |>
      purrr::walk(function(pop){
        covar_names |>
          purrr::walk(function(cov){
            logfc_pairs[[study]] |>
              purrr::walk(function(pair){
                x <- pair[1]; y <- pair[2]
                if (!all(c(x, y) %in% exposure_vars_list[[study]])) return(NULL)
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
  submitted_counterpart <- c(
    comp_pca_all           = "comp_wqs_all",
    comp_qgcomp_fw_all     = "comp_qgcomp_all",
    comp_wqs_cf_all        = "comp_wqs_all",
    comp_qgcomp_cf_all     = "comp_qgcomp_all",
    comp_qgcomp_cox_cf_all = "comp_qgcomp_cox_all"
  )

  study_names |>
    purrr::walk(function(study){
      covar_names |>
        purrr::walk(function(cov){
          exposure_vars_list[[study]] |>
            purrr::walk(function(exp){
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
    width = 10,
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

  ## Heatmap-style overview: -log10(p) for the pathways that reach p < 0.05
  ## under at least one revision exposure.
  pathway_overview <- pathway_long |>
    dplyr::filter(population == "all", covar_set == "covar",
                  !is.na(p_value), is.finite(p_plot), p_plot > 0) |>
    dplyr::group_by(pathway) |>
    dplyr::filter(any(p_value < 0.05)) |>
    dplyr::ungroup() |>
    dplyr::mutate(neg_log10_p = -log10(p_plot),
                  label = unname(rev_label(exposure)),
                  category = categorize_pathway(pathway))

  if (nrow(pathway_overview) > 0) {
    overview_plot <- ggplot(
      pathway_overview,
      aes(x = label, y = reorder(stringr::str_trunc(pathway, 45),
                                 neg_log10_p))) +
      geom_tile(aes(fill = neg_log10_p), colour = "white", linewidth = 0.3) +
      scale_fill_gradient(low = "#F7F7F7", high = "#B73F42",
                          name = expression(-log[10](p))) +
      ggh4x::facet_nested(category ~ algorithm + study,
                          scales = "free", space = "free_y") +
      labs(title = "Revision pathway enrichment across exposure definitions",
           subtitle = "Full cohort, primary covariate set",
           x = NULL, y = NULL) +
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y = element_text(angle = 0, face = "bold", size = 7),
            strip.background = element_rect(fill = "grey95", colour = NA))

    rev_save_plot(overview_plot, "pathway_overview_heatmap", "pathway",
                  width = 15, height = 14)
  }

} else {
  message("Pathway summary not found at ", pathway_long_path,
          " - run R6-pathway_revision.R first. Skipping Section 6.")
}


# =============================================================================
# SECTION 7: COMBINED PANEL FIGURE
# =============================================================================

## Layout follows the submitted analysis (figures/mwas/.../combined_panel_*.png):
## a single row of three panels --
##   A  volcano for this exposure
##   B  stratified subgroup betas against the full-cohort betas
##   C  sensitivity-covariate betas against primary-covariate betas
##
## create_population_scatter() and create_covariate_scatter() are taken
## verbatim from scripts/7-visualization.R so the revision panels are directly
## comparable with the submitted ones.

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

build_panel_parts <- function(exposure_name, population, covar_set, study) {

  p1 <- create_volcano_plot(
    mwas_result_c18 = mwas_results_list_c18[[study]][[population]][[covar_set]][[exposure_name]],
    mwas_result_hilic = mwas_results_list_hilic[[study]][[population]][[covar_set]][[exposure_name]],
    annotation_result_c18 = mwas_full_annotated_list_c18[[study]][[population]][[covar_set]][[exposure_name]],
    annotation_result_hilic = mwas_full_annotated_list_hilic[[study]][[population]][[covar_set]][[exposure_name]],
    exposure_name = exposure_name
  ) +
    labs(title = NULL) +
    guides(shape = guide_legend(order = 1), color = guide_legend(order = 2)) +
    panel_legend_theme

  ## Panel B is the subgroup-versus-full-cohort contrast. As in
  ## scripts/7-visualization.R it is the same comparison regardless of which
  ## population the rest of the panel shows, so it is included throughout and
  ## dropped only if a study lacks the cognitive strata.
  p2 <- NULL
  if (all(c("demcind", "no demcind") %in% names(mwas_results_list_c18[[study]]))) {
    p2 <- create_population_scatter(
      mwas_results_list_c18[[study]], mwas_results_list_hilic[[study]],
      covar_set = covar_set, exposure_name = exposure_name
    ) +
      labs(title = NULL) + guides(shape = "none") + panel_legend_theme
  }

  p3 <- create_covariate_scatter(
    mwas_results_list_c18[[study]], mwas_results_list_hilic[[study]],
    population = population, cov1 = "covar", cov2 = "covar_sen",
    exposure_name = exposure_name
  ) +
    labs(title = NULL) + guides(shape = "none") + panel_legend_theme

  purrr::compact(list(p1, p2, p3)) |> purrr::map(rev_wrap_axis_titles)
}

## Annotation is applied once, to the finished layout. Tagging a patchwork and
## then nesting it inside another restarts the sequence, so the Manhattan
## variant is assembled from the raw panel list rather than from the annotated
## three-panel object -- that is what lets its Manhattan carry tag D.
annotate_panel <- function(layout, exposure_name, population, covar_set) {
  layout +
    patchwork::plot_annotation(
      title = paste0("MWAS Results: ", rev_label(exposure_name),
                     " (", population, ", ", covar_set, ")"),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
        plot.tag   = element_text(face = "bold", size = 30))
    )
}

create_combined_panel <- function(exposure_name, population, covar_set, study) {
  panels <- build_panel_parts(exposure_name, population, covar_set, study)
  list(panel = annotate_panel(patchwork::wrap_plots(panels, nrow = 1),
                              exposure_name, population, covar_set),
       n = length(panels))
}

## Second variant: the same row with the Manhattan stacked underneath, for
## when the m/z distribution of the hits matters as much as their volcano
## position. Saved alongside, not instead of, the three-panel version.
create_combined_panel_manhattan <- function(exposure_name, population,
                                            covar_set, study) {
  panels <- build_panel_parts(exposure_name, population, covar_set, study)

  man <- create_manhattan(
    combined_results_list_c18[[study]][[population]][[covar_set]][[exposure_name]],
    combined_results_list_hilic[[study]][[population]][[covar_set]][[exposure_name]],
    exposure_name = exposure_name
  ) +
    labs(title = NULL) + panel_legend_theme
  man <- rev_wrap_axis_titles(man)

  layout <- patchwork::wrap_plots(panels, nrow = 1) / man +
    patchwork::plot_layout(heights = c(1.6, 1))

  list(panel = annotate_panel(layout, exposure_name, population, covar_set),
       n = length(panels))
}

study_names |>
  purrr::walk(function(study){
    populations_for(study) |>
      purrr::walk(function(pop){
        covar_names |>
          purrr::walk(function(cov){
            exposure_vars_list[[study]] |>
              purrr::walk(function(exp){
                message("Combined panel: ", study, "_", pop, " [", cov, "] - ", exp)
                rev_dir("figures", "mwas", study, pop, cov, exp)

                three <- create_combined_panel(exp, pop, cov, study)
                rev_ggsave(
                  filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                      glue::glue("combined_panel_{exp}.png")),
                  plot = three$panel, width = 10 * three$n, height = 9,
                  dpi = 300, bg = "white"
                )

                withman <- create_combined_panel_manhattan(exp, pop, cov, study)
                rev_ggsave(
                  filename = rev_here("figures", "mwas", study, pop, cov, exp,
                                      glue::glue("combined_panel_manhattan_{exp}.png")),
                  plot = withman$panel, width = 10 * withman$n, height = 16,
                  dpi = 300, bg = "white"
                )
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
