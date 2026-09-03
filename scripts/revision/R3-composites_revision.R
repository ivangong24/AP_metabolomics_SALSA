## ---------------------------
##
## Script name: R3-composites_revision.R
## Purpose of script: To build the unsupervised PCA exposure index and the
##                    cross-fitted WQS / QGcomp indices used by the R1
##                    revision MWAS
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
## Notes: Addresses Reviewer 1 major comment 1 and Reviewer 2 comment 1 - the
##        objection that the WQS and QGcomp indices were estimated by
##        regressing dementia/CIND on the pollutant mixture and then reused as
##        exposures in the MWAS.
##
##        Two answers are built here:
##
##        A. UNSUPERVISED PCA (the new primary exposure). PC1 of the same
##           quartile-scored (0-3) pollutant matrix that WQS and QGcomp use.
##           Never sees the outcome, so there is no data reuse at all.
##
##        B. CROSS-FITTED WQS / QGCOMP. Weights are estimated in K-1 folds and
##           applied to the held-out fold, using the TRAINING fold's quartile
##           breakpoints. Folds are assigned to participants so repeated blood
##           draws stay together and no participant's exposure score depends on
##           their own outcome.
##
##        Two implementation choices worth stating in the response letter:
##          - Weight derivation uses ONE ROW PER PARTICIPANT (first blood
##            draw). The submitted analysis fit the weight models at the
##            specimen level, which double-counts participants with repeated
##            draws.
##          - The cross-fitted index keeps the weight convention of the naive
##            index: WQS uses final_weights$Estimate, QGcomp uses
##            coefficient / sum(|coefficient|), matching get_weights() in
##            scripts/qgcomp_modified.R.
##
##        Outputs -> revision_output/{data,tables,figures}/composites/
##        Downstream: R4-mwas_revision.R consumes
##        revision_output/data/processed/combined_data_list_revision.RData
##
##        RUNTIME: a few hours, dominated by 30 gWQS fits
##        (3 groupings x 2 covariate sets x 5 folds, each b = 200, rh = 5).
##        Validate first with Sys.setenv(SALSA_REVISION_QUICK = "true").
##
## ---------------------------

# Load required packages -----------------------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

pacman::p_load("gWQS", "qgcomp", "survival")

source(here::here("scripts", "wqs_modified.R"))

rev_announce("R3-composites_revision.R")

set.seed(42)

K_FOLDS <- rev_n(5,   3)
B_WQS   <- rev_n(200, 5)
RH_WQS  <- rev_n(5,   2)


# Load exposure data ---------------------------------------------------------

load(here::here("data", "processed", "combined_data_list_new.RData"))

dat_total <- combined_data_list_new[["total"]][["all"]][["covar"]]
dat_cox   <- combined_data_list_new[["cox"]][["all"]][["covar"]]

## covar_sen adds four covariates to the weight models, so the weight-derivation
## frame is built separately for each covariate set
dat_total_sen <- combined_data_list_new[["total"]][["all"]][["covar_sen"]]
dat_cox_sen   <- combined_data_list_new[["cox"]][["all"]][["covar_sen"]]

dat_total_list <- list(covar = dat_total, covar_sen = dat_total_sen)
dat_cox_list   <- list(covar = dat_cox,   covar_sen = dat_cox_sen)

exp_groups <- make_exp_groups(dat_total)


# Windowed pollutant columns -------------------------------------------------

## R2-exposure_windows_revision.R rebuilds each pollutant's average over a 3-
## and a 10-year window. Those columns are attached to the weight-derivation
## frames here so the cross-fitted WQS and QGcomp weights can be RE-DERIVED per
## window rather than transported from the five-year fit -- a weight is a
## property of the pollutant matrix it was estimated on, and carrying the
## five-year weights onto a ten-year matrix would confound the window with the
## weighting.
##
## Optional: R3 still runs without R2, and the window arm simply does not
## appear downstream.
window_path <- rev_here("data", "processed", "exposure_windows_revision.RData")

if (file.exists(window_path)) {
  load(window_path)
  if (!exists("window_exposure_df")) {
    stop("exposure_windows_revision.RData predates the windowed-pollutant ",
         "export. Re-run R2-exposure_windows_revision.R.")
  }
  message("Exposure-window indices loaded: ",
          paste(setdiff(names(window_pca_df), c("rand_id", "blood_date")),
                collapse = ", "))
  message("Windowed pollutant columns loaded: ",
          paste(setdiff(names(window_exposure_df), c("rand_id", "blood_date")),
                collapse = ", "))
} else {
  window_pca_df      <- tibble::tibble(rand_id = character(),
                                       blood_date = as.Date(character()))
  window_exposure_df <- window_pca_df
  warning("No exposure-window extract at ", window_path,
          " -- run R2-exposure_windows_revision.R first if the window ",
          "sensitivity analysis is wanted.", call. = FALSE)
}

## Which windows this run can actually build. A window whose pollutant columns
## did not arrive is dropped here rather than failing inside gWQS.
window_cols <- setdiff(names(window_exposure_df), c("rand_id", "blood_date"))

crossfit_windows <- MWAS_WINDOWS |>
  purrr::keep(function(w) all(pollutants_for_window(w) %in% window_cols))

if (!identical(sort(crossfit_windows), sort(MWAS_WINDOWS))) {
  warning("Cross-fitted indices will not be built for window(s) ",
          paste(setdiff(MWAS_WINDOWS, crossfit_windows), collapse = ", "),
          ": their pollutant columns are missing from R2's extract.",
          call. = FALSE)
}

attach_windows <- function(data) {
  data |> dplyr::left_join(window_exposure_df,
                           by = c("rand_id", "blood_date"))
}

dat_total_list <- dat_total_list |> purrr::map(attach_windows)
dat_cox_list   <- dat_cox_list   |> purrr::map(attach_windows)
dat_total      <- dat_total_list[["covar"]]

message("Mixture groupings:")
purrr::iwalk(exp_groups, ~ message("  ", .y, " (", length(.x), "): ",
                                   paste(.x, collapse = ", ")))

message("\nSpecimens: total = ", nrow(dat_total), " from ",
        dplyr::n_distinct(dat_total$rand_id), " participants")
message("Specimens: cox   = ", nrow(dat_cox), " from ",
        dplyr::n_distinct(dat_cox$rand_id), " participants")


# =============================================================================
# SECTION 1: UNSUPERVISED PCA EXPOSURE INDEX
# =============================================================================

message("\n=== Section 1: unsupervised PCA indices ===")

pca_index_list <- exp_groups |>
  purrr::imap(function(vars, grouping){
    message("  grouping: ", grouping)
    make_pca_index(dat_total, vars, name = paste0("comp_pca_", grouping))
  })

pca_df <- pca_index_list |>
  purrr::map("composites") |>
  purrr::reduce(dplyr::full_join, by = c("rand_id", "blood_date"))


# PC1 loadings and variance explained ----------------------------------------

pca_loadings <- pca_index_list |>
  purrr::imap(function(pca_obj, grouping){
    tibble::tibble(
      grouping    = grouping,
      pollutant   = names(pca_obj$loadings) |> str_remove("^exp_"),
      pc1_loading = round(unname(pca_obj$loadings), 4),
      pc1_pve     = round(unname(pca_obj$pve[1]), 4)
    )
  }) |>
  purrr::list_rbind()

pca_variance <- pca_index_list |>
  purrr::imap(function(pca_obj, grouping){
    tibble::tibble(
      grouping   = grouping,
      component  = names(pca_obj$pve),
      prop_var   = round(unname(pca_obj$pve), 4),
      cum_var    = round(cumsum(unname(pca_obj$pve)), 4)
    )
  }) |>
  purrr::list_rbind()

rev_save_table(pca_loadings, "pca_loadings", "composites")
rev_save_table(pca_variance, "pca_variance_explained", "composites")

print(pca_loadings, n = 30)


# Sanity check against the stored comp_pca_* columns -------------------------

## 3-clean_data.R already computes the same PC1. Recomputing here keeps this
## script self-contained; the check confirms nothing drifted.
pca_check <- exp_groups |>
  names() |>
  purrr::map(function(grouping){
    col <- paste0("comp_pca_", grouping)
    joined <- dat_total |>
      dplyr::select(rand_id, blood_date, stored = dplyr::all_of(col)) |>
      dplyr::inner_join(pca_df |>
                          dplyr::select(rand_id, blood_date,
                                        recomputed = dplyr::all_of(col)),
                        by = c("rand_id", "blood_date"))
    tibble::tibble(
      grouping = grouping,
      n        = nrow(joined),
      r        = round(stats::cor(joined$stored, joined$recomputed,
                                  use = "complete.obs"), 6)
    )
  }) |>
  purrr::list_rbind()

message("\nRecomputed vs stored comp_pca_* (expect r = 1):")
print(pca_check)


# Scree and loading plots ----------------------------------------------------

scree_plot <- pca_variance |>
  ggplot(aes(x = component, y = prop_var, group = grouping)) +
  geom_col(fill = "#4C72B0", alpha = 0.85) +
  geom_line(aes(y = cum_var), colour = "#BE3F42", linewidth = 0.7) +
  geom_point(aes(y = cum_var), colour = "#BE3F42", size = 1.8) +
  facet_wrap(~ grouping, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Principal component",
       y = "Proportion of variance",
       title = "PCA of the quartile-scored air toxicant matrix",
       subtitle = "Bars: variance per component. Line: cumulative variance.") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

rev_save_plot(scree_plot, "pca_scree", "composites", width = 10, height = 4)

loading_plot <- pca_loadings |>
  ggplot(aes(x = reorder(pollutant, pc1_loading), y = pc1_loading)) +
  geom_col(aes(fill = pc1_loading > 0), show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#4C72B0", `FALSE` = "#DE9960")) +
  facet_wrap(~ grouping, scales = "free_y") +
  coord_flip() +
  labs(x = NULL, y = "PC1 loading",
       title = "PC1 loadings of the unsupervised air toxicant index") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

rev_save_plot(loading_plot, "pca_loadings", "composites",
              width = 10, height = 4)


# =============================================================================
# SECTION 2: CROSS-FITTED WQS AND QGCOMP INDICES
# =============================================================================

message("\n=== Section 2: cross-fitted indices (K = ", K_FOLDS, ") ===")

## Every combination of covariate set x mixture grouping x method. WQS and
## QGcomp logistic are derived in the pooled cohort (study = total); QGcomp Cox
## is derived in the incident cohort (study = cox).
## The window dimension. "" is the five-year primary; w3 and w10 re-derive the
## same weights on their own pollutant matrix.
##
## The Cox arm stays at the five-year window: it is the incident cohort's
## single exposure, the window question is asked of the pooled cross-sectional
## analysis, and each extra Cox index would add three more populations to the
## pathway grid for an arm that is already the smallest.
crossfit_specs <- tidyr::expand_grid(
  covar_set = names(covars_weight_list),
  grouping  = names(exp_groups),
  method    = c("wqs", "qgcomp_glm", "qgcomp_cox"),
  window    = c("", crossfit_windows)
) |>
  dplyr::filter(!(method == "qgcomp_cox" & window != "")) |>
  dplyr::mutate(
    study    = dplyr::if_else(method == "qgcomp_cox", "cox", "total"),
    prefix   = dplyr::recode(method,
                             wqs        = "comp_wqs_cf_",
                             qgcomp_glm = "comp_qgcomp_cf_",
                             qgcomp_cox = "comp_qgcomp_cox_cf_"),
    comp_name = purrr::map2_chr(prefix, window,
                                ~ crossfit_exposure_for(.x, .y))
  )

## `grouping` is `all` throughout, so the mixture members are that window's
## pollutant columns.
crossfit_vars <- crossfit_specs$window |> purrr::map(pollutants_for_window)

print(crossfit_specs, n = 40)

system.time({
  crossfit_fits <- seq_len(nrow(crossfit_specs)) |>
    purrr::map(function(i){
      spec <- crossfit_specs[i, ]
      vars <- crossfit_vars[[i]]

      message("\n  ", spec$comp_name, " [", spec$covar_set, "] over ",
              length(vars), " pollutants")

      data <- if (spec$study == "cox") {
        dat_cox_list[[spec$covar_set]]
      } else {
        dat_total_list[[spec$covar_set]]
      }

      crossfit_composite(
        data          = data,
        vars          = vars,
        covars_weight = covars_weight_list[[spec$covar_set]],
        method        = spec$method,
        k             = K_FOLDS,
        b_wqs         = B_WQS,
        rh_wqs        = RH_WQS,
        seed          = 42
      )
    })
})

names(crossfit_fits) <- paste(crossfit_specs$covar_set,
                              crossfit_specs$comp_name, sep = "|")


# Assemble one composite table per covariate set -----------------------------

crossfit_df_list <- names(covars_weight_list) |>
  purrr::set_names() |>
  purrr::map(function(covar_set){
    idx <- which(crossfit_specs$covar_set == covar_set)

    idx |>
      purrr::map(function(i){
        crossfit_fits[[i]]$scores |>
          dplyr::select(rand_id, blood_date, score) |>
          dplyr::rename(!!crossfit_specs$comp_name[i] := score)
      }) |>
      purrr::reduce(dplyr::full_join, by = c("rand_id", "blood_date"))
  })


# Fold-level weights ---------------------------------------------------------

## `pollutant` is stripped back to the bare name -- exp_benzene_w3 -> benzene --
## so the same pollutant lines up across windows in the summary table and the
## stability plot. Which window a row belongs to is carried by `window`, not by
## the pollutant name.
crossfit_weights <- seq_len(nrow(crossfit_specs)) |>
  purrr::map(function(i){
    spec <- crossfit_specs[i, ]
    crossfit_fits[[i]]$weights |>
      dplyr::mutate(covar_set = spec$covar_set,
                    composite = spec$comp_name,
                    window    = dplyr::if_else(spec$window == "", "w5",
                                               spec$window),
                    pollutant = str_remove(pollutant, "^exp_") |>
                      str_remove(paste0("_", spec$window, "$")),
                    .before = 1)
  }) |>
  purrr::list_rbind()

crossfit_weight_summary <- crossfit_weights |>
  dplyr::group_by(covar_set, composite, window, pollutant) |>
  dplyr::summarise(
    mean_weight = round(mean(weight), 4),
    sd_weight   = round(stats::sd(weight), 4),
    min_weight  = round(min(weight), 4),
    max_weight  = round(max(weight), 4),
    .groups = "drop"
  ) |>
  dplyr::arrange(covar_set, composite, dplyr::desc(mean_weight))

rev_save_table(crossfit_weights, "crossfit_weights_by_fold", "composites")
rev_save_table(crossfit_weight_summary, "crossfit_weights_summary",
               "composites")

## Fold-to-fold spread in the weights is itself the stability evidence
## Reviewer 1 asked for at the end of major comment 1.
## Every window, so the figure answers "are the weights stable" and "do they
## move with the averaging window" in one place. Faceted by composite, which
## now carries the window in its name.
weight_stability_plot <- crossfit_weights |>
  dplyr::filter(covar_set == "covar") |>
  ggplot(aes(x = reorder(pollutant, weight), y = weight)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92", width = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 1.6, colour = "#4C72B0") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ composite, scales = "free_x") +
  coord_flip() +
  labs(x = NULL, y = "Mixture weight",
       title = "Fold-to-fold stability of the cross-fitted mixture weights",
       subtitle = paste0("All-toxicant groupings, primary covariate set, K = ",
                         K_FOLDS, " folds")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

rev_save_plot(weight_stability_plot, "crossfit_weight_stability", "composites",
              width = 14, height = 9)


# =============================================================================
# SECTION 3: AGREEMENT BETWEEN NAIVE, CROSS-FITTED AND UNSUPERVISED INDICES
# =============================================================================

message("\n=== Section 3: composite agreement ===")

comparison_df_list <- names(covars_weight_list) |>
  purrr::set_names() |>
  purrr::map(function(covar_set){
    dat_total |>
      dplyr::select(rand_id, blood_date,
                    dplyr::starts_with("comp_wqs_"),
                    dplyr::starts_with("comp_qgcomp_")) |>
      dplyr::left_join(pca_df, by = c("rand_id", "blood_date")) |>
      dplyr::left_join(crossfit_df_list[[covar_set]],
                       by = c("rand_id", "blood_date")) |>
      ## The PCA window indices, so the window-versus-window pairs below can be
      ## read off the same frame as everything else.
      dplyr::left_join(window_pca_df, by = c("rand_id", "blood_date"))
  })

comparison_pairs <- tibble::tribble(
  ~comparison,                                      ~a,                      ~b,
  "WQS all: naive vs cross-fitted",                 "comp_wqs_all",          "comp_wqs_cf_all",
  "QGcomp all: naive vs cross-fitted",              "comp_qgcomp_all",       "comp_qgcomp_cf_all",
  "QGcomp Cox all: naive vs cross-fitted",          "comp_qgcomp_cox_all",   "comp_qgcomp_cox_cf_all",
  "PCA all vs naive WQS all",                       "comp_pca_all",          "comp_wqs_all",
  "PCA all vs naive QGcomp all",                    "comp_pca_all",          "comp_qgcomp_all",
  "PCA all vs naive QGcomp Cox all",                "comp_pca_all",          "comp_qgcomp_cox_all",
  "PCA all vs cross-fitted WQS all",                "comp_pca_all",          "comp_wqs_cf_all",
  "PCA all vs cross-fitted QGcomp all",             "comp_pca_all",          "comp_qgcomp_cf_all",
  "PCA all vs cross-fitted QGcomp Cox all",         "comp_pca_all",          "comp_qgcomp_cox_cf_all"
)

## The window contrasts (R1 comment 13). Each index against its own 5-year
## self, so the correlation isolates the averaging window: same method, same
## weighting rule, same cohort. Reported alongside the PC1 window correlations
## R2 already writes, which cover the unsupervised index only.
window_comparison_pairs <- crossfit_windows |>
  purrr::map(function(w){
    yrs <- WINDOW_YEARS[[w]]
    tibble::tibble(
      comparison = c(
        paste0("PCA all: 5-year vs ", yrs, "-year"),
        paste0("Cross-fitted WQS all: 5-year vs ", yrs, "-year"),
        paste0("Cross-fitted QGcomp all: 5-year vs ", yrs, "-year")),
      a = c("comp_pca_all", "comp_wqs_cf_all", "comp_qgcomp_cf_all"),
      b = c(paste0("comp_pca_all_", w),
            crossfit_exposure_for("comp_wqs_cf_", w),
            crossfit_exposure_for("comp_qgcomp_cf_", w))
    )
  }) |>
  purrr::list_rbind()

comparison_pairs <- dplyr::bind_rows(comparison_pairs, window_comparison_pairs)

cor_pair <- function(data, a, b) {
  if (!all(c(a, b) %in% names(data))) return(NA_real_)
  ok <- stats::complete.cases(data[[a]], data[[b]])
  if (sum(ok) < 3) return(NA_real_)
  round(stats::cor(data[[a]][ok], data[[b]][ok]), 3)
}

composite_correlations <- comparison_df_list |>
  purrr::imap(function(data, covar_set){
    comparison_pairs |>
      dplyr::mutate(
        covar_set = covar_set,
        n = purrr::map2_int(a, b, function(x, y){
          if (!all(c(x, y) %in% names(data))) return(NA_integer_)
          sum(stats::complete.cases(data[[x]], data[[y]]))
        }),
        r = purrr::map2_dbl(a, b, ~ cor_pair(data, .x, .y)),
        .before = 1
      )
  }) |>
  purrr::list_rbind()

rev_save_table(composite_correlations, "composite_correlations", "composites")
print(composite_correlations |> dplyr::filter(covar_set == "covar"), n = 20)


# Descriptive statistics of every composite ----------------------------------

composite_descriptives <- comparison_df_list |>
  purrr::imap(function(data, covar_set){
    data |>
      dplyr::select(dplyr::starts_with("comp_")) |>
      tidyr::pivot_longer(dplyr::everything(),
                          names_to = "composite", values_to = "value") |>
      dplyr::group_by(composite) |>
      dplyr::summarise(
        covar_set = covar_set,
        n_nonmiss = sum(!is.na(value)),
        mean = round(mean(value, na.rm = TRUE), 4),
        sd   = round(stats::sd(value, na.rm = TRUE), 4),
        iqr  = round(stats::IQR(value, na.rm = TRUE), 4),
        min  = round(min(value, na.rm = TRUE), 4),
        max  = round(max(value, na.rm = TRUE), 4),
        .groups = "drop"
      ) |>
      dplyr::relocate(covar_set)
  }) |>
  purrr::list_rbind()

rev_save_table(composite_descriptives, "composite_descriptives", "composites")


# Scatter plots: naive vs cross-fitted vs PCA --------------------------------

make_scatter <- function(data, x, y, xlab, ylab) {
  d <- data |>
    dplyr::select(x = dplyr::all_of(x), y = dplyr::all_of(y)) |>
    tidyr::drop_na()
  r <- round(stats::cor(d$x, d$y), 3)

  ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.25, size = 0.8, colour = "steelblue") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "firebrick", linewidth = 0.7) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = paste0("r = ", r), size = 3.5, fontface = "bold") +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 10)
}

composite_scatter <- list(
  c("comp_wqs_all",        "comp_wqs_cf_all",        "WQS (naive)",        "WQS (cross-fitted)"),
  c("comp_qgcomp_all",     "comp_qgcomp_cf_all",     "QGcomp (naive)",     "QGcomp (cross-fitted)"),
  c("comp_qgcomp_cox_all", "comp_qgcomp_cox_cf_all", "QGcomp Cox (naive)", "QGcomp Cox (cross-fitted)"),
  c("comp_pca_all",        "comp_wqs_cf_all",        "PCA (unsupervised)", "WQS (cross-fitted)"),
  c("comp_pca_all",        "comp_qgcomp_cf_all",     "PCA (unsupervised)", "QGcomp (cross-fitted)"),
  c("comp_pca_all",        "comp_qgcomp_cox_cf_all", "PCA (unsupervised)", "QGcomp Cox (cross-fitted)")
) |>
  purrr::map(function(spec){
    make_scatter(comparison_df_list[["covar"]], spec[1], spec[2],
                 spec[3], spec[4])
  }) |>
  patchwork::wrap_plots(nrow = 2) +
  patchwork::plot_annotation(
    title = "Air toxicant composite indices, all-toxicant grouping",
    subtitle = paste0("Primary covariate set. Top row: outcome-informed ",
                      "indices before and after cross-fitting. ",
                      "Bottom row: unsupervised PC1 against each."),
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

rev_save_plot(composite_scatter, "composite_correlations", "composites",
              width = 12, height = 8)


# =============================================================================
# SECTION 4: MERGE COMPOSITES BACK INTO THE ANALYSIS DATA
# =============================================================================

message("\n=== Section 4: assembling combined_data_list_revision ===")

## Drop every composite from the submitted analysis and attach the revision
## exposures, so the downstream scripts cannot accidentally pick up an
## outcome-informed naive index. The window objects were loaded at the top of
## the script, where the cross-fitted window indices needed them.

## The RAW windowed pollutant columns travel too, not just the indices built
## from them. R4 derives the windowed single-pollutant exposures
## (exp_*_w3_iqr) and the windowed quartile scores (exp_*_w3_q) for the
## feature-wise QGcomp contrast from these, exactly as it derives the 5-year
## ones from exp_benzene and friends -- so the scaling convention is the same
## at every window and is set in one place.
window_raw_cols <- crossfit_windows |>
  purrr::map(pollutants_for_window) |>
  unlist(use.names = FALSE)

combined_data_list_revision <- combined_data_list_new |>
  purrr::imap(function(datalist, study){
    datalist |>
      purrr::imap(function(data_list, population){
        data_list |>
          purrr::imap(function(data, covar_name){
            message("  ", study, "_", population, " [", covar_name, "]")

            data |>
              dplyr::select(-dplyr::starts_with("comp_")) |>
              dplyr::left_join(pca_df, by = c("rand_id", "blood_date")) |>
              dplyr::left_join(crossfit_df_list[[covar_name]],
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(window_pca_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::left_join(window_exposure_df,
                               by = c("rand_id", "blood_date")) |>
              dplyr::select(dplyr::any_of(c(
                names(data)[!str_detect(names(data), "^comp_")],
                window_raw_cols,
                rev_exposure_vars_list[[study]]
              )))
          })
      })
  })

## Confirm every study carries exactly the exposures it should.
##
## Three kinds of exposure, checked three ways:
##
##   scored composites  built here, so checked against what was just attached
##   feature-wise QGcomp a CONTRAST, never a column -- checked by the presence
##                       of the raw pollutant columns its quartile scores are
##                       derived from in R4
##   single pollutants   not built here either. The 5-year _iqr columns come
##                       through untouched from 3-clean_data.R; the windowed
##                       ones are derived in R4 from the raw windowed columns
##                       attached above, so it is the RAW column that has to be
##                       present at this point.
message("\nExposures available by study:")
combined_data_list_revision |>
  purrr::iwalk(function(datalist, study){
    frame   <- datalist[["all"]][["covar"]]
    present <- frame |> dplyr::select(dplyr::starts_with("comp_")) |> names()
    message("  ", study, ": ", paste(present, collapse = ", "))

    wanted <- rev_exposure_vars_list[[study]]

    composites <- wanted |>
      purrr::discard(is_single_pollutant) |>
      purrr::discard(is_qgcomp_fw)
    missing <- setdiff(composites, present)
    if (length(missing) > 0) {
      warning("Missing revision composites for ", study, ": ",
              paste(missing, collapse = ", "))
    }

    ## What each exposure needs to be derivable in R4: an _iqr column for the
    ## 5-year pollutants, a raw column for everything windowed and for every
    ## feature-wise contrast.
    needed <- c(
      intersect(wanted, single_pollutant_exposures),
      wanted |> purrr::keep(is_single_pollutant) |>
        purrr::keep(~ exposure_window(.x) != "") |>
        (\(x) stringr::str_remove(x, "_iqr$"))(),
      wanted |> purrr::keep(is_qgcomp_fw) |>
        purrr::map(~ pollutants_for_window(exposure_window(.x))) |>
        unlist(use.names = FALSE)
    ) |> unique()

    missing_p <- setdiff(needed, names(frame))
    if (length(missing_p) > 0) {
      warning("Missing pollutant columns for ", study, ": ",
              paste(missing_p, collapse = ", "))
    } else if (length(needed) > 0) {
      message("  ", study, " pollutant columns present: ", length(needed))
    }
  })

## Completeness of each revision exposure
composite_missingness <- combined_data_list_revision |>
  purrr::imap(function(datalist, study){
    datalist |>
      purrr::imap(function(data_list, population){
        data_list |>
          purrr::imap(function(data, covar_name){
            data |>
              dplyr::select(dplyr::starts_with("comp_")) |>
              tidyr::pivot_longer(dplyr::everything(),
                                  names_to = "composite",
                                  values_to = "value") |>
              dplyr::group_by(composite) |>
              dplyr::summarise(n_total = dplyr::n(),
                               n_nonmiss = sum(!is.na(value)),
                               .groups = "drop") |>
              dplyr::mutate(study = study, population = population,
                            covar_set = covar_name, .before = 1)
          }) |>
          purrr::list_rbind()
      }) |>
      purrr::list_rbind()
  }) |>
  purrr::list_rbind()

rev_save_table(composite_missingness, "composite_missingness", "composites")
print(composite_missingness |>
        dplyr::filter(population == "all", covar_set == "covar"), n = 20)


# Save R objects for downstream analysis -------------------------------------

rev_dir("data", "processed")

save(combined_data_list_revision,
     pca_index_list, pca_df, pca_loadings, pca_variance,
     crossfit_df_list, crossfit_weights, crossfit_weight_summary,
     crossfit_specs, crossfit_windows,
     window_pca_df, window_exposure_df,
     composite_correlations, composite_descriptives,
     exp_groups, K_FOLDS,
     file = rev_here("data", "processed",
                     "combined_data_list_revision.RData"))

message("\nComposite construction completed!")
message("Results saved to:")
message("  - ", rev_here("data", "processed",
                         "combined_data_list_revision.RData"))
message("  - ", rev_here("tables", "composites"))
message("  - ", rev_here("figures", "composites"))

#--------------------------------End of the code--------------------------------
