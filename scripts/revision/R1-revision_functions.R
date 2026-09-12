## ---------------------------
##
## Script name: R1-revision_functions.R
## Purpose of script: To create shared paths, configuration and helper
##                    functions for the R1 revision analyses
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
## Notes: This file is sourced by every R*-*_revision.R script AFTER
##        scripts/1-functions.R (which calls rm(list = ls())).
##
##        It provides:
##        1. rev_here() - every revision output goes under revision_output/,
##           mirroring the layout of the primary pipeline so that
##           revision_output/tables/mwas_results/... is the exact analogue of
##           tables/mwas_results/...
##        2. The exposure groupings, covariate sets and study/population
##           structure shared with scripts 3-7
##        3. make_pca_index()  - unsupervised PC1 of the quartile-scored
##           exposure matrix (the new outcome-free primary exposure)
##        4. crossfit_composite() - K-fold cross-fitted WQS / QGcomp indices,
##           so no participant's exposure score is built from their own
##           dementia/CIND outcome
##
## ---------------------------

# Revision output paths ------------------------------------------------------

## Everything the revision writes lives under revision_output/. Nothing under
## tables/, figures/, data/ or metaboAnalyst/ is touched - those stay as the
## record of the submitted analysis.

REV_ROOT <- Sys.getenv("SALSA_REVISION_ROOT", unset = "revision_output")

rev_here <- function(...) here::here(REV_ROOT, ...)

rev_dir <- function(...) {
  path <- rev_here(...)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

## Reduced replicate counts for a smoke test. Set
## Sys.setenv(SALSA_REVISION_QUICK = "true") to validate a script end to end
## in minutes; nothing produced under QUICK is reportable.
REV_QUICK <- tolower(Sys.getenv("SALSA_REVISION_QUICK", unset = "false")) %in%
  c("true", "1", "yes")

rev_n <- function(full, quick) if (REV_QUICK) quick else full

rev_announce <- function(script_name) {
  message(strrep("=", 78))
  message("  ", script_name)
  message("  output root: ", rev_here())
  if (REV_QUICK) {
    message("  *** QUICK MODE - reduced replicates, NOT reportable ***")
  }
  message(strrep("=", 78))
}


# Shared configuration -------------------------------------------------------

## Covariate sets. demcind is carried in the data but dropped from the model
## covariates by covars_list_new below.
##
## BATCH IS NOT ADJUSTED FOR (changed 2026-08-28, was in scripts 3-7).
## Three reasons, in order of weight:
##
##   1. ComBat already removed batch effects from the feature matrices, so a
##      batch term in the regression is a second adjustment for the same thing
##      -- Reviewer 1 comment 7.
##   2. It was entering as a NUMERIC term across 39 nominal levels, which fits
##      one linear slope on batch index rather than adjusting for batch at all.
##      It added nothing over wave: anova(~wave, ~wave + batch) gives
##      F = 1.81, p = 0.178 for comp_pca_all.
##   3. Coded correctly as a factor it would be over-adjustment: batch is
##      confounded with collection time (factor(batch) explains 77% of calendar
##      year, 78.5% of wave, and 18 of 39 batches hold a single wave), so it
##      absorbs 45% of the comp_pca_all contrast and 60% of butadiene's.
##
## Dropping it does not move the results: refitting C18 and HILIC for
## comp_pca_all in total/all/covar gives coefficient correlations of 0.999 and
## 1.000 against the batch-adjusted fit.
##
## WAVE IS RETAINED, AND AS A NUMERIC TERM. Specimens were collected at
## different visits spanning 1998-2007, and wave carries that collection-time
## structure; without any time adjustment the same model returns 754 (C18) and
## 1,430 (HILIC) FDR-significant features, which is secular confounding rather
## than signal.
covar_list <- list(
  covar = quote_all(age_at_blooddraw, gender, edu_year, mh62,
                    ruca_metro, nses,
                    wave, demcind),

  covar_sen = quote_all(age_at_blooddraw, gender, edu_year, mh62,
                        ruca_metro, nses,
                        alcohol_drinking, pa3_met_if_ca,
                        bmi_at_blooddraw, diab_at_blooddraw,
                        wave, demcind)
)

covars_list_new <- covar_list |>
  purrr::map(function(covars){
    covars |>
      purrr::discard(~ str_detect(.x, regex("demcind", ignore_case = TRUE)))
  })

## Covariates used when deriving the mixture weights. As in 3-clean_data.R the
## weight models drop demcind (it is the outcome), and wave/batch (they are
## metabolomics run variables, not exposure-model variables).
covars_weight_list <- covar_list |>
  purrr::map(function(covars){
    covars |>
      purrr::discard(~ str_detect(.x, "demcind|dcst|wave|batch"))
  })

## Mixture groupings, matching scripts/3-clean_data.R exactly. Ozone and zinc
## were dropped upstream, so the eight exp_ columns are the mixture components.
make_exp_groups <- function(data) {
  exp_vars <- data |>
    dplyr::select(dplyr::starts_with("exp_")) |>
    dplyr::select(-dplyr::matches("o3|zinc|_iqr$")) |>
    names()

  ## `all` only. The `traffic` and `metal` subgroupings were built and
  ## diagnosed in earlier revision runs and are not carried anywhere now: the
  ## metal QGcomp weights were not identified (weights changed sign between
  ## cross-fitting folds), and traffic was sound but redundant. Collapsing to a
  ## single mixture definition also cuts the analytical configurations behind
  ## Reviewer 1 major comment 5. Their diagnostics remain in
  ## revision_output/archive/ as the record.
  list(all = exp_vars)
}

## Exposures carried into the revision MWAS, by study.
##
##   total : the pooled repeated-measures analysis. Unsupervised PCA is the
##           PRIMARY exposure; the cross-fitted WQS and QGcomp logistic indices
##           are the outcome-informed comparators.
##   cox   : the incident cohort, carrying only the cross-fitted QGcomp Cox
##           index. PCA is deliberately NOT repeated here -- it is a single
##           participant-level score, so running it in the incident cohort as
##           well would only re-estimate the same exposure in a subsample and
##           add three more models to every multiplicity count. This mirrors
##           scripts/4-mwas_analysis.R, which also excludes comp_pca from the
##           cox arm.
##
## MIXTURE GROUPINGS: `all` only. The `traffic` and `metal` subgroupings have
## been removed from the revision entirely (see make_exp_groups above) -- they
## are no longer built, diagnosed or carried, so nothing downstream needs to
## filter them out. What remains is a single mixture definition weighted four
## ways: unsupervised PCA, feature-wise QGcomp, cross-fitted WQS and
## cross-fitted QGcomp logistic. Their earlier diagnostics are preserved under
## revision_output/archive/.
## SINGLE POLLUTANTS: the eight mixture components are also carried as
## exposures in their own right -- Reviewer 1 major comment 1, which asks for
## single-pollutant MWAS results alongside the mixture. They were run in the
## submitted analysis (`tables/mwas_results/`), but under the pre-revision
## covariate set, so pairing those numbers with the revision composites would
## compare two different models. Running them here puts both sides of the
## attribution argument on the same covariates, the same imputation and the
## same duplicateCorrelation, which is what lets the manuscript say whether the
## taurine association and the alanine/aspartate enrichment belong to any one
## pollutant or only to the mixture.
##
## They are IQR-scaled, as in the submitted analysis, so a coefficient is the
## difference in log2 abundance per interquartile-range increase in that
## pollutant. Note this is NOT the per-SD scale the composites were put on in
## R4 -- SD/IQR ranges from 0.73 (benzene) to 1.53 (butadiene) across the
## eight, so the two scales are not interchangeable and any figure placing a
## pollutant next to a composite has to say which it is using.

## Pollutants entering the feature-wise quantile g-computation contrast, in the
## order R3 uses for the `all` grouping. The same eight are the single-pollutant
## exposures.
qgcomp_pollutants <- c("exp_benzene", "exp_butadiene", "exp_chromium",
                       "exp_lead", "exp_no2", "exp_nickel", "exp_pm2.5",
                       "exp_nox")

single_pollutant_exposures <- paste0(qgcomp_pollutants, "_iqr")

## EXPOSURE WINDOWS: Reviewer 1 major comment 13. The primary exposure averages
## each pollutant over the five years preceding the draw. R2 rebuilds those
## averages over a 1-, 3- and 10-year window.
##
## comp_pca_all_w1 and comp_pca_all_w3 are SEVEN-pollutant indices. CALINE4
## NOx exists for model years 1998-2002 only, so 1,442 of 1,546 specimens (93%)
## have no NOx year inside a one-year window and 701 (45%) none inside a
## three-year window -- they would enter the mixture as the pipeline's NA -> 0
## fill. NOx is dropped from both, and both must be reported as
## seven-component indices. comp_pca_all_w10 carries all eight. R2 decides this
## from a coverage threshold rather than hardcoding it.
window_pca_exposures <- c("comp_pca_all_w1", "comp_pca_all_w3",
                          "comp_pca_all_w10")

## Backwards-compatible alias. R7 reads `window_exposures` when it derives the
## pollutant composition of a window index from R2's loadings table.
window_exposures <- window_pca_exposures


## THE FULL EXPOSURE SET AT THE 3- AND 10-YEAR WINDOWS (added 2026-09-03)
## ---------------------------------------------------------------------------
##
## Rebuilding only the PCA index answered "does the window matter for the
## primary exposure", but it could not answer the question the window analysis
## is actually asked in the same breath as major comment 1 -- whether the
## attribution to the mixture, and to any one pollutant inside it, is a
## property of the five-year averaging. That needs the OTHER exposures at the
## same windows, so the 3- and 10-year windows now carry the whole set:
##
##   comp_wqs_cf_all_{w}      cross-fitted WQS, weights re-derived on that
##                            window's pollutant matrix
##   comp_qgcomp_cf_all_{w}   cross-fitted QGcomp logistic, likewise
##   comp_qgcomp_fw_all_{w}   feature-wise QGcomp contrast over that window's
##                            quartile-scored pollutants
##   exp_*_{w}_iqr            that window's single pollutants
##
## The cross-fitted weights are RE-DERIVED per window rather than transported
## from the five-year fit. A weight is a property of the pollutant matrix it
## was estimated on, and carrying the five-year weights onto a ten-year matrix
## would confound the window with the weighting -- the opposite of what the
## PCA arm was designed to avoid.
##
## The 1-year window is deliberately not extended. It reaches NOx for 6.7% of
## specimens, and at that coverage the cross-fitted weight models and the
## single-pollutant NOx model are estimated on a near-constant column.
MWAS_WINDOWS <- c("w3", "w10")

## The 5-year primary carries no window tag. "" is that window everywhere a
## window has to be named.
ALL_MWAS_WINDOWS <- c("", MWAS_WINDOWS)

WINDOW_YEARS <- c(w1 = 1, w3 = 3, w5 = 5, w10 = 10)

## COLUMN-NAME CONVENTION. The window tag is a suffix on the RAW pollutant
## column, so every derived name is a plain append and one rule covers all of
## them:
##
##   5-year   exp_benzene       exp_benzene_iqr       exp_benzene_q
##   3-year   exp_benzene_w3    exp_benzene_w3_iqr    exp_benzene_w3_q
##
## Composites tag at the end: comp_wqs_cf_all_w3, comp_qgcomp_fw_all_w10.

## Which pollutants each window carries. Read back from R2's own coverage
## decision so the rule lives in ONE place -- R2 decides it from
## WINDOW_MIN_COVERAGE against the data, and a change there has to reach the
## exposure list rather than being restated here. The fallback matches the
## measured coverage documented in R2 and is used only before R2 has run.
window_pollutant_map <- function() {
  fallback <- list(w3  = setdiff(qgcomp_pollutants, "exp_nox"),
                   w10 = qgcomp_pollutants)[MWAS_WINDOWS]

  path <- rev_here("data", "processed", "exposure_windows_revision.RData")
  if (!file.exists(path)) return(fallback)

  env <- new.env()
  ok <- try(load(path, envir = env), silent = TRUE)
  if (inherits(ok, "try-error") || is.null(env$window_pca_loadings)) {
    return(fallback)
  }

  m <- env$window_pca_loadings |>
    dplyr::filter(window %in% MWAS_WINDOWS) |>
    dplyr::distinct(window, pollutant) |>
    (\(d) split(paste0("exp_", d$pollutant), d$window))()

  if (!all(MWAS_WINDOWS %in% names(m))) return(fallback)
  ## R2's table is alphabetical; keep qgcomp_pollutants' order so the contrast
  ## terms and the panel C columns read the same at every window.
  m[MWAS_WINDOWS] |> purrr::map(~ intersect(qgcomp_pollutants, .x))
}

WINDOW_POLLUTANTS <- window_pollutant_map()

## The raw pollutant columns of one window. window = "" is the 5-year primary.
pollutants_for_window <- function(window = "") {
  if (identical(window, "")) return(qgcomp_pollutants)
  paste0(WINDOW_POLLUTANTS[[window]], "_", window)
}

single_pollutant_exposures_for <- function(window = "") {
  paste0(pollutants_for_window(window), "_iqr")
}

qgcomp_q_names_for <- function(window = "") {
  paste0(pollutants_for_window(window), "_q")
}

qgcomp_fw_exposure_for <- function(window = "") {
  if (identical(window, "")) "comp_qgcomp_fw_all"
  else paste0("comp_qgcomp_fw_all_", window)
}

crossfit_exposure_for <- function(prefix, window = "") {
  if (identical(window, "")) paste0(prefix, "all")
  else paste0(prefix, "all_", window)
}

## Every window's set, flattened, for the membership tests below.
all_single_pollutant_exposures <- ALL_MWAS_WINDOWS |>
  purrr::map(single_pollutant_exposures_for) |>
  unlist(use.names = FALSE)

all_qgcomp_fw_exposures <- ALL_MWAS_WINDOWS |>
  purrr::map_chr(qgcomp_fw_exposure_for)

all_raw_pollutants <- ALL_MWAS_WINDOWS |>
  purrr::map(pollutants_for_window) |>
  unlist(use.names = FALSE)

## The window an exposure name belongs to, "" for the 5-year primary.
##
## Only MWAS_WINDOWS are matched, so comp_pca_all_w1 reports "" -- there is no
## 1-year single-pollutant or cross-fitted arm for it to be matched against,
## and reporting "w1" here would send the combined panel looking for models
## that were never fitted.
exposure_window <- function(exposure_var) {
  x <- as.character(exposure_var)
  pat <- paste0("_(", paste(MWAS_WINDOWS, collapse = "|"), ")(_iqr|_q)?$")
  tag <- stringr::str_match(x, pat)[, 2]
  ifelse(is.na(tag), "", tag)
}

## The averaging window of an exposure IN YEARS, whatever its name shape.
##
## Distinct from exposure_window(), which reports only the windows that carry
## their own fitted exposure set. comp_pca_all_w1 is a one-year index but has
## no one-year single-pollutant arm to be matched against, so exposure_window()
## calls it "" while this reports 1 -- which is what a figure has to say.
exposure_window_years <- function(exposure_var) {
  tag <- stringr::str_match(as.character(exposure_var),
                            "_(w[0-9]+)(_iqr|_q)?$")[, 2]
  ifelse(is.na(tag), 5, unname(WINDOW_YEARS[tag]))
}

## The same exposure at the 5-year primary window: the name with its tag
## removed. Used for labels and for the window-versus-window comparisons.
untagged_exposure <- function(exposure_var) {
  x <- as.character(exposure_var)
  w <- exposure_window(x)
  out <- purrr::map2_chr(x, w, function(nm, tag){
    if (identical(tag, "")) return(nm)
    stringr::str_remove(nm, paste0("_", tag, "(?=(_iqr|_q)?$)"))
  })
  out
}

## The 3- and 10-year arms, in the same order the 5-year exposures appear in
## so a table or a facet strip reads window-by-window.
window_mwas_exposures <- MWAS_WINDOWS |>
  purrr::map(function(w){
    c(qgcomp_fw_exposure_for(w),
      crossfit_exposure_for("comp_wqs_cf_", w),
      crossfit_exposure_for("comp_qgcomp_cf_", w),
      single_pollutant_exposures_for(w))
  }) |>
  unlist(use.names = FALSE)

rev_exposure_vars_list <- list(
  total = c("comp_pca_all", "comp_qgcomp_fw_all",
            "comp_wqs_cf_all", "comp_qgcomp_cf_all",
            single_pollutant_exposures, window_pca_exposures,
            window_mwas_exposures),
  cox   = "comp_qgcomp_cox_cf_all"
)

## Populations each exposure is fitted in.
##
## The three cognitive strata -- `all`, `no demcind`, `demcind` -- carry
## everything. `all predx` (the post-diagnosis sensitivity population, R1
## comment 2) carries the 5-year composites and the PCA window indices only.
##
## THE SINGLE POLLUTANTS ARE NO LONGER `all`-ONLY (changed 2026-09-03). They
## were, on the argument that the attribution question is about the primary
## stratum. But the single-pollutant half of the combined panel then showed
## full-cohort estimates beside a stratum-specific index, and -- once the 3-
## and 10-year windows carried their own pollutants -- five-year estimates
## beside a ten-year index. Both mismatches are invisible in the figure. They
## are fitted in all three strata at all three windows now, so panel C is
## matched to the stratum AND the window of the index it sits beside.
##
## `all predx` is deliberately left out of the widening. It is a row filter on
## `all` whose only purpose is the post-diagnosis sensitivity check, and every
## exposure added there is another cell in the pathway grid and another model
## in the multiplicity count Reviewer 1 comment 4 is about.
STRATA_POPULATIONS <- c("all", "no demcind", "demcind")

## ONE EXCEPTION, added 2026-09-09: `all predx` additionally carries the
## FIVE-YEAR single pollutants, for the MWAS only.
##
## The response to Reviewer 1 comment 1 argues that benzene and 1,3-butadiene
## are the components surviving the retention-time check, and the response to
## comment 2 argues that post-diagnosis specimens do not drive the findings.
## The second claim was demonstrated for the composites only, so the letter was
## asserting for the pollutants what it had shown for the indices. All eight
## are added rather than only the two the letter highlights: running the check
## just where a finding is expected would be the outcome-informed selection
## comment 1 objects to, in a smaller form.
##
## Deliberately NOT a widening of STRATA_POPULATIONS. That constant gates
## `restricted_exposures`, which is the single pollutants at ALL THREE windows
## PLUS the windowed composites, so flipping it would add 24 pollutant arms and
## the whole window grid rather than 8 cells.
##
## Deliberately MWAS-only. These cells are limma-only anyway (pls_eligible()
## discards single pollutants), and pathway_exposures_for() below keeps them
## out of the Mummichog and Metapone grids, so the multiplicity argument above
## still holds for the pathway results.
PREDX_POPULATION      <- "all predx"
PREDX_EXTRA_EXPOSURES <- single_pollutant_exposures

## Kept as an alias: R4 and R7 refer to the old name in comments and a stale
## reference should not be a silent NULL.
RESTRICTED_POPULATIONS <- STRATA_POPULATIONS

restricted_exposures <- c(single_pollutant_exposures, window_mwas_exposures)

is_restricted_exposure <- function(exposure_var) {
  as.character(exposure_var) %in% restricted_exposures
}

is_single_pollutant <- function(exposure_var) {
  as.character(exposure_var) %in% all_single_pollutant_exposures
}

## The exposure set for one study x population cell.
##
## Every stage of the revision -- design matrices, PLS, Mummichog input,
## figures -- iterates over exposures inside a population loop, so this is the
## one place that decides which exposures a cell carries. `available` lets a
## caller intersect with what a particular object actually holds.
exposures_for <- function(study, population, available = NULL) {
  ## R4 narrows rev_exposure_vars_list to the columns actually present in the
  ## analysis frames and saves the result as exposure_vars_list, which R6 and
  ## R7 load back. Prefer it when it exists so a missing column is never
  ## silently promised downstream.
  exps <- if (exists("exposure_vars_list", inherits = TRUE)) {
    get("exposure_vars_list")[[study]]
  } else {
    rev_exposure_vars_list[[study]]
  }
  if (!is.null(available)) exps <- intersect(exps, available)
  if (!population %in% STRATA_POPULATIONS) {
    ## `all predx` keeps the 5-year single pollutants; every other
    ## non-stratum population drops the restricted set entirely.
    allowed <- if (identical(population, PREDX_POPULATION)) {
      PREDX_EXTRA_EXPOSURES
    } else {
      character(0)
    }
    exps <- exps[!is_restricted_exposure(exps) | exps %in% allowed]
  }
  exps
}

## The exposure set for one pathway cell.
##
## Same as exposures_for() everywhere except `all predx`, where the single
## pollutants added on 2026-09-09 are MWAS-only: eight pollutants x two
## covariate sets x two tools would be 32 new pathway cells for a sensitivity
## population, and the post-diagnosis question is answered by the MWAS
## coefficient agreement rather than by a second pathway table. Nickel is the
## specific hazard -- its HILIC hits are 69.7% void-region artifact (R13), and
## a pre-diagnosis nickel pathway table would read as corroboration of it.
pathway_exposures_for <- function(study, population, available = NULL) {
  exps <- exposures_for(study, population, available = available)
  if (identical(population, PREDX_POPULATION)) {
    exps <- exps[!exps %in% PREDX_EXTRA_EXPOSURES]
  }
  exps
}

QGCOMP_FW_EXPOSURE <- "comp_qgcomp_fw_all"
QGCOMP_FW_Q        <- 4
qgcomp_q_names     <- paste0(qgcomp_pollutants, "_q")

## `comp_qgcomp_fw_all` is deliberately NOT a column in the analysis frames: it
## is a contrast. Every other exposure here is a single scored index whose
## weights came from a model fitted to dementia/CIND; this one inverts the
## direction of the weighting entirely, taking each metabolite as the outcome
## and the eight quartile-scored pollutants as the exposures, so no outcome
## model is involved at any point. psi is the sum of the eight coefficients --
## which is exactly what qgcomp::qgcomp.glm.noboot returns for a Gaussian model
## with no product terms, verified to machine precision. Estimating it as a
## limma contrast rather than by calling qgcomp per feature keeps
## duplicateCorrelation for the repeated draws and the eBayes moderation,
## neither of which qgcomp provides.
## One per window: comp_qgcomp_fw_all is the 5-year contrast,
## comp_qgcomp_fw_all_w3 and _w10 the same estimand over their own window's
## quartile-scored pollutants.
is_qgcomp_fw <- function(exposure_var) {
  as.character(exposure_var) %in% all_qgcomp_fw_exposures
}

## Quartile-score each pollutant (0 .. q-1) within the frame it will be modelled
## in, mirroring qgcomp's default breaks. Ties collapse the upper categories
## rather than erroring, which matters for the small demcind strata.
##
## Pollutants absent from the frame are skipped rather than erroring: the cox
## frames carry no windowed exposure columns, and the 3-year window carries no
## NOx.
add_qgcomp_quantiles <- function(data, pollutants = qgcomp_pollutants,
                                 q = QGCOMP_FW_Q) {
  pollutants <- intersect(pollutants, names(data))
  if (length(pollutants) == 0) return(data)

  data |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(pollutants),
      function(x) {
        brk <- stats::quantile(x, probs = seq(0, 1, length.out = q + 1),
                               na.rm = TRUE)
        as.numeric(cut(x, breaks = unique(brk), include.lowest = TRUE,
                       labels = FALSE)) - 1
      },
      .names = "{.col}_q"))
}

## Every window's quartile scores, in one pass over the frame.
add_all_qgcomp_quantiles <- function(data, windows = ALL_MWAS_WINDOWS,
                                     q = QGCOMP_FW_Q) {
  purrr::reduce(windows,
                function(d, w) add_qgcomp_quantiles(d, qgcomp_pollutants_raw(w),
                                                    q = q),
                .init = data)
}

## Raw pollutant columns of a window, guarded so a caller can ask for a window
## the run does not carry.
qgcomp_pollutants_raw <- function(window = "") {
  if (!identical(window, "") && is.null(WINDOW_POLLUTANTS[[window]])) {
    return(character(0))
  }
  pollutants_for_window(window)
}

## Contrast that sums the quantized-pollutant coefficients into psi.
##
## `window` names which set of quantized terms to sum. The design for
## comp_qgcomp_fw_all_w3 carries exp_*_w3_q and no others, so summing the
## 5-year terms there would silently sum nothing -- hence the hard stop.
qgcomp_psi_contrast <- function(design_matrix, window = "") {
  q_names <- qgcomp_q_names_for(window)
  w <- as.numeric(colnames(design_matrix) %in% q_names)
  if (sum(w) != length(q_names)) {
    stop("Design is missing quantized pollutant terms: ",
         paste(setdiff(q_names, colnames(design_matrix)), collapse = ", "))
  }
  matrix(w, ncol = 1,
         dimnames = list(colnames(design_matrix), "psi"))
}

## Human-readable labels used in figures and tables
rev_exposure_labels <- c(
  comp_pca_all             = "PCA all toxicants (unsupervised)",
  comp_qgcomp_fw_all       = "QGcomp all toxicants (feature-wise)",
  comp_wqs_cf_all          = "WQS all toxicants (cross-fitted)",
  comp_qgcomp_cf_all       = "QGcomp all toxicants (cross-fitted)",
  comp_qgcomp_cox_cf_all   = "QGcomp Cox all toxicants (cross-fitted)",
  ## Single pollutants, per IQR. Labels match scripts/7-visualization.R so the
  ## revision figures read the same as the submitted ones.
  exp_benzene_iqr          = "Benzene",
  exp_butadiene_iqr        = "1,3-Butadiene",
  exp_chromium_iqr         = "Chromium",
  exp_lead_iqr             = "Lead",
  exp_no2_iqr              = "NO2",
  exp_nickel_iqr           = "Nickel",
  `exp_pm2.5_iqr`          = "PM2.5",
  exp_nox_iqr              = "NOx",
  ## Exposure-window sensitivity (R1 comment 13). The 5-year window is the
  ## primary comp_pca_all above.
  comp_pca_all_w1          = "PCA, 1-year window (7 pollutants)",
  comp_pca_all_w3          = "PCA, 3-year window (7 pollutants)",
  comp_pca_all_w10         = "PCA, 10-year window"
)

## The 3- and 10-year arms take their 5-year label with the window appended,
## derived rather than restated so a change to a base label reaches all three
## windows. A composite whose window carries fewer than eight pollutants says
## so, exactly as the PCA window labels above do -- comparing a seven-component
## index with an eight-component one as though the window were the only
## difference is the mistake the annotation is there to prevent.
rev_exposure_labels <- c(
  rev_exposure_labels,
  MWAS_WINDOWS |>
    purrr::map(function(w){
      exps <- c(qgcomp_fw_exposure_for(w),
                crossfit_exposure_for("comp_wqs_cf_", w),
                crossfit_exposure_for("comp_qgcomp_cf_", w),
                single_pollutant_exposures_for(w))
      base <- unname(rev_exposure_labels[untagged_exposure(exps)])
      n_p  <- length(WINDOW_POLLUTANTS[[w]])
      note <- ifelse(startsWith(exps, "comp_") & n_p < length(qgcomp_pollutants),
                     paste0(" (", n_p, " pollutants)"), "")
      stats::setNames(
        paste0(base, ", ", WINDOW_YEARS[[w]], "-year window", note),
        exps)
    }) |>
    unlist()
)

## Always character, including for zero-length input. ifelse() returns
## logical(0) when its test is empty, which is not a labelling problem until
## the result is joined against a character column -- dplyr then refuses with
## "Can't join `x$pollutant` with `y$pollutant` due to incompatible types",
## from a figure whose only fault was that nothing passed its filter.
rev_label <- function(exp_name) {
  if (length(exp_name) == 0) return(character(0))
  as.character(ifelse(exp_name %in% names(rev_exposure_labels),
                      rev_exposure_labels[exp_name], exp_name))
}

## Population labels, short enough for a panel subtitle or a facet strip.
rev_population_labels <- c(
  "all"        = "All",
  "no demcind" = "No dem/CIND",
  "demcind"    = "Dem/CIND",
  "all predx"  = "All, pre-diagnosis"
)

rev_population_label <- function(population) {
  if (length(population) == 0) return(character(0))
  as.character(ifelse(population %in% names(rev_population_labels),
                      rev_population_labels[population], population))
}

## Markdown variant, for the figure text ggtext renders --------------------
##
## The subscripts in NO2, NOx and PM2.5 are part of the species name, not
## decoration, and "PM2.5" written flat is wrong. Only the labels that need
## one differ from rev_label(); everything else passes through, so a caller
## can use rev_label_md() everywhere it renders markdown without special-
## casing.
##
## Plain-text callers -- table columns, glue strings, file names -- must keep
## using rev_label(). An <sub> tag in a spreadsheet cell is worse than a flat
## label.
## Built from rev_exposure_labels rather than listed, so the windowed variants
## (exp_no2_w3_iqr, exp_pm2.5_w10_iqr, ...) get their subscript too. Listing
## them by hand is how a windowed label silently reverts to flat "PM2.5".
rev_exposure_labels_md <- (function(){
  keep <- stringr::str_detect(names(rev_exposure_labels),
                              "^exp_(no2|nox|pm2\\.5)(_w[0-9]+)?_iqr$")
  labs <- rev_exposure_labels[keep]
  ## NOx before NO2: "^NO2" cannot match "NOx", but keeping the more specific
  ## species first makes the intent obvious if another NO-something is added.
  stats::setNames(
    labs |>
      stringr::str_replace("^NOx",     "NO<sub>x</sub>") |>
      stringr::str_replace("^NO2",     "NO<sub>2</sub>") |>
      stringr::str_replace("^PM2\\.5", "PM<sub>2.5</sub>"),
    names(labs))
})()

rev_label_md <- function(exp_name) {
  if (length(exp_name) == 0) return(character(0))
  as.character(ifelse(exp_name %in% names(rev_exposure_labels_md),
                      rev_exposure_labels_md[exp_name], rev_label(exp_name)))
}


## Shared significance palette -------------------------------------------------
##
## ONE definition, used by the volcano, the Manhattan and the effect-estimate
## panel alike. They previously carried near-identical palettes that differed
## in the red only (#BE3F42 vs #B73F42). patchwork collects guides by comparing
## them for equality, so that one-character difference was enough to make it
## treat the two Significance legends as distinct and draw both -- which is
## what crowded the bottom of the combined panel. Keep this shared.
##
## THE MIDDLE BAND IS FDR < 0.10, NOT VIP > 2 (changed 2026-09-02).
##
## It used to be "P < 0.05 & VIP > 2". That put a quantity with no error
## control in the middle of a legend whose other two entries are calibrated:
## R14 measured the empirical FDR of VIP > 2 at ~1, which is why
## filter_significant() dropped it on 2026-09-01. Shading a band with it was
## the last place it still carried visual weight, and a reader has no way to
## tell a display class from a selection rule by looking. A second FDR tier
## says something a reader can act on -- these are the features that would
## survive a more permissive error rate -- and it nests properly inside
## P < 0.05, which VIP never did.
##
## VIP has not been discarded: create_vip_manhattan() still plots it, as its
## own figure, with its own axis, where it is labelled for what it is.
SIG_COLORS <- c("FDR < 0.05" = "#B73F42",
                "FDR < 0.10" = "#436C85",
                "P < 0.05"   = "#DE9960")

## Class order, shared so the volcano, Manhattan and forest cannot drift.
SIG_LEVELS <- c(names(SIG_COLORS), "NS")

## The one place the rule itself lives. Nested and calibrated: FDR < 0.05 is a
## subset of FDR < 0.10, which is (essentially) a subset of P < 0.05.
sig_class <- function(p_value, adj_p_value) {
  factor(
    dplyr::case_when(
      adj_p_value < 0.05 ~ "FDR < 0.05",
      adj_p_value < 0.10 ~ "FDR < 0.10",
      p_value     < 0.05 ~ "P < 0.05",
      TRUE               ~ "NS"),
    levels = SIG_LEVELS)
}


# Identification confidence --------------------------------------------------

## Schymanski et al., Environ Sci Technol 2014;48:2097-2098 -- the scheme
## reviewers in this journal family expect:
##
##   1  confirmed structure (reference standard / MS-MS)
##   2  probable structure (MS-MS spectral match to a library)
##   3  tentative candidate(s) -- a structure is proposed but the evidence
##      does not resolve one exact structure
##   4  unequivocal molecular formula, no structure proposed
##   5  exact mass only, formula not assignable
##
## What the two annotation routes in this study can support:
##
##   In-house library -- accurate mass (<=10 ppm) AND retention time (<=30 s)
##     against authentic reference standards run on the same C18-neg /
##     HILIC-pos methods (scripts/5-annotation.R:184). Reported as LEVEL 1.
##     The Methods must state that MS/MS was not acquired and that
##     identification rests on accurate mass and retention time against
##     authentic standards -- Schymanski's Level 1 lists MS, MS/MS and RT, so
##     the claim needs that sentence to stand on its own.
##
##   xMSannotator (HMDB / KEGG / LIPID MAPS) -- accurate mass plus adduct and
##     isotope consistency, no authentic standard and no MS/MS. LEVEL 3: a
##     structure is proposed, but mass alone does not resolve it.
##
## LEVEL 2 IS DELIBERATELY EMPTY. It requires an MS/MS spectral match and this
## study acquired MS1 full scan only. An empty level is informative -- it says
## the scheme was applied rather than the features distributed across all five
## bins.
##
## LEVELS 4 AND 5 ARE NOT ASSIGNED HERE. They describe features for which no
## structure is proposed at all, which is the ~14,500 unannotated features of
## the ~20,125 in the feature space, not the annotated ones. Grading an
## annotated feature down to 4 or 5 because the database returned several
## names would say it has weaker FORMULA confidence than an unannotated
## feature, which is backwards.
##
## Multiplicity is reported separately, as n_candidates, rather than folded
## into the level: it says exactly how ambiguous an annotation is, which a
## level number cannot.
##
## NOTE: the `confidence` column on the annotation frames is xMSannotator's own
## 0-3 quality score, NOT one of these levels. Do not conflate them.

## An in-house hit is not automatically one compound: 48 of the 153 in-house
## annotations list several candidates sharing the matched mass and retention
## time. Some are the same substance written twice or differing only in
## stereochemistry ("Ornithine; D-Ornithine"), which accurate mass and RT
## cannot separate anyway; others are genuinely different structures
## ("L-Leucine; L-Isoleucine; L-Norleucine"). The first collapses to one
## compound and keeps Level 1; the second does not, and is Level 3 -- several
## tentative candidates, which is what Level 3 is for.
STEREO_PREFIX <- "^(d|l|dl|r|s|rs|\\(r\\)|\\(s\\)|\\(rs\\)|\\(\\+\\)|\\(-\\)|cis|trans|allo)-"

candidate_list <- function(compound) {
  compound |>
    stringr::str_split(";") |>
    purrr::map(~ stringr::str_squish(.x)) |>
    purrr::map(~ .x[.x != "" & !is.na(.x)])
}

## Distinct compounds after folding away duplicates and stereodescriptors.
n_candidate_compounds <- function(compound) {
  candidate_list(compound) |>
    purrr::map_int(function(v){
      if (length(v) == 0) return(NA_integer_)
      length(unique(stringr::str_remove_all(
        stringr::str_remove(stringr::str_to_lower(v),
                            stringr::regex(STEREO_PREFIX)),
        "[[:space:]\\-]")))
    })
}

## What to print on a figure. The shortest candidate with its stereodescriptor
## removed: mass and retention time cannot assign D from L, so "Tyrosine" is a
## more honest label than "L-Tyrosine".
candidate_display_name <- function(compound) {
  candidate_list(compound) |>
    purrr::map_chr(function(v){
      if (length(v) == 0) return(NA_character_)
      short <- v[which.min(nchar(v))]
      out <- stringr::str_remove(short,
                                 stringr::regex(STEREO_PREFIX,
                                                ignore_case = TRUE))
      paste0(toupper(substr(out, 1, 1)), substr(out, 2, nchar(out)))
    })
}

annotation_confidence_level <- function(reference, compound) {
  n_cand <- n_candidate_compounds(compound)
  dplyr::case_when(
    is.na(reference)                              ~ NA_integer_,
    reference == "In House Library" & n_cand == 1 ~ 1L,
    reference == "In House Library"               ~ 3L,
    TRUE                                          ~ 3L
  )
}


# Small helpers --------------------------------------------------------------

## demcind is a labelled factor. Passing it to glm() or survival::Surv()
## silently does the wrong thing, so every model outcome goes through as01().
as01 <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  dplyr::case_when(
    as.character(x) == "Dementia/CIND"    ~ 1,
    as.character(x) == "No Dementia/CIND" ~ 0,
    TRUE ~ NA_real_
  )
}

## Quartile score (0-3) using supplied breakpoints. Values outside the training
## range are clamped into the end bins rather than returned as NA.
quantile_score <- function(x, breaks) {
  breaks <- unique(breaks)
  if (length(breaks) < 2) return(rep(0L, length(x)))
  q <- as.integer(cut(x, breaks = breaks,
                      include.lowest = TRUE, labels = FALSE)) - 1L
  q[!is.na(x) & x <= min(breaks, na.rm = TRUE)] <- 0L
  q[!is.na(x) & x >= max(breaks, na.rm = TRUE)] <- length(breaks) - 2L
  q
}

## Quartile-score a whole exposure matrix on its own quantiles
quantile_score_matrix <- function(data, vars) {
  data |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(vars),
                    ~ quantile_score(.x, stats::quantile(.x, probs = 0:4 / 4,
                                                         na.rm = TRUE)))
    )
}

rev_save_table <- function(tbl, name, topic) {
  path <- file.path(rev_dir("tables", topic), paste0(name, ".xlsx"))
  writexl::write_xlsx(tbl, path)
  message("  table -> ", path)
  invisible(path)
}

## Reviewer 1 minor comment 14: "Figure text is difficult to read at the
## current size. Larger labels and clearer legends ... would improve
## accessibility." Rather than chase the size arguments scattered through the
## individual theme() blocks, every figure gets this override appended at save
## time, so the sizes are set in one place and stay consistent across the
## volcano, Manhattan, scatter, heatmap and pathway figures.
##
## Sizes are chosen for the panel context: the combined panels are saved at
## 8 inches per column, so a 20-pt axis title here renders at roughly the same
## visual weight as a 12-pt label in a single-column journal figure.
REV_TEXT <- theme(
  plot.title    = element_text(face = "bold", size = 18, hjust = 0.5),
  plot.subtitle = element_text(size = 14),
  plot.caption  = element_text(size = 12),
  axis.title    = element_text(face = "bold", size = 16),
  axis.text     = element_text(size = 14),
  strip.text    = element_text(face = "bold", size = 16),
  legend.title  = element_text(face = "bold", size = 14),
  legend.text   = element_text(size = 13),
  plot.tag      = element_text(face = "bold", size = 26)
)

## Percent difference, for Reviewer 1 minor comment 6 -------------------------
##
## "consider reporting percent differences, which would also permit direct
## comparison with Hu et al. and Qi et al., both cited in percent terms."
##
## Feature intensities are log2-transformed, so the limma coefficient is the
## difference in log2 abundance per one-SD increase in the exposure index
## (see the standardization step in R4). Back-transforming,
##
##     percent difference = (2^beta - 1) * 100
##
## The interval is built from the coefficient's own standard error, recovered
## as beta / t, and back-transformed on the log2 scale before conversion --
## the transform is monotone, so the bounds map directly and the interval is
## asymmetric on the percent scale, as it should be.
add_percent_difference <- function(df, conf = 0.95) {
  if (!all(c("logFC", "t") %in% names(df))) return(df)
  z  <- stats::qnorm(1 - (1 - conf) / 2)
  se <- ifelse(is.finite(df$t) & df$t != 0, df$logFC / df$t, NA_real_)
  df |>
    dplyr::mutate(
      pct_diff    = (2^.data$logFC - 1) * 100,
      pct_diff_lo = (2^(.data$logFC - z * se) - 1) * 100,
      pct_diff_hi = (2^(.data$logFC + z * se) - 1) * 100
    ) |>
    dplyr::relocate(pct_diff, pct_diff_lo, pct_diff_hi, .after = "logFC")
}

## The combined panels are 30 inches wide, so they carry larger tick numbers
## and legend text than the standalone figures without crowding. Applied on
## top of REV_TEXT to the individual panels, not to the assembled patchwork,
## so the standalone volcano and Manhattan keep their own sizes.
REV_TEXT_PANEL <- theme(
  axis.title   = element_text(face = "bold", size = 19),
  axis.text    = element_text(size = 18),
  strip.text   = element_text(face = "bold", size = 19),
  legend.title = element_text(face = "bold", size = 17),
  legend.text  = element_text(size = 16)
)

## Long axis titles at these sizes run into neighbouring panels -- panel C's
## y-axis title collided with the figure title before this. Wrap them instead
## of shrinking the text back down.
##
## Plotmath labels must be left alone. The volcano y-axis is
## expression(-log[10](P - value)) and the Manhattan x-axis is
## expression(bold("Mass-to-charge ratio (" * italic(m/z) * ")")); these are
## language objects, and str_wrap() coerces them to their deparsed SOURCE,
## so the axis then reads "-log[10](P - value)" literally. Wrap only plain
## character labels and pass everything else through untouched.
rev_wrap_axis_titles <- function(p, width = 38) {
  wrap_one <- function(lab) {
    if (is.character(lab) && length(lab) == 1 && !is.na(lab)) {
      stringr::str_wrap(lab, width)
    } else {
      lab
    }
  }
  if (!is.null(p$labels$y)) p <- p + ggplot2::labs(y = wrap_one(p$labels$y))
  if (!is.null(p$labels$x)) p <- p + ggplot2::labs(x = wrap_one(p$labels$x))
  p
}

## `post` is a theme applied AFTER REV_TEXT, for the few elements REV_TEXT
## would otherwise overwrite. The case that needs it is ggtext: REV_TEXT sets
## axis.text as a plain element_text, and ggplot2 refuses to merge that over
## an element_markdown ("Only elements of the same class can be merged"), so a
## figure wanting markdown axis labels has to re-assert them last rather than
## first.
rev_save_plot <- function(plot, name, topic, width = 10, height = 7,
                          post = NULL) {
  path <- file.path(rev_dir("figures", topic), paste0(name, ".png"))
  p <- plot & REV_TEXT
  if (!is.null(post)) p <- p & post
  ggplot2::ggsave(path, p, width = width, height = height,
                  dpi = 300, bg = "white")
  message("  figure -> ", path)
  invisible(path)
}


# A. Unsupervised PCA exposure index -----------------------------------------

## PC1 of the quartile-scored exposure matrix.
##
## PCA is applied to the SAME quartile-scored (0-3) matrix that WQS and QGcomp
## use, so the three indices differ only in how the pollutants are weighted,
## not in how they are scaled. PCA finds the direction of maximum variance in
## the exposure data and never touches dementia/CIND, which is exactly what
## Reviewer 1 (major comment 1) and Reviewer 2 (comment 1) asked for.
##
## The sign of PC1 is arbitrary, so it is flipped when needed so that the
## majority of loadings are positive - higher score means higher overall
## exposure burden, matching the WQS convention.
make_pca_index <- function(data, vars, name,
                           id_cols = c("rand_id", "blood_date")) {

  dat <- data |>
    dplyr::select(dplyr::all_of(c(id_cols, vars))) |>
    tidyr::drop_na()

  dat_q <- quantile_score_matrix(dat, vars)

  pca_fit <- stats::prcomp(dat_q[vars], center = TRUE, scale. = TRUE)

  loadings_pc1 <- pca_fit$rotation[, 1]
  scores_pc1   <- pca_fit$x[, 1]

  if (sum(loadings_pc1 < 0) > sum(loadings_pc1 > 0)) {
    loadings_pc1 <- -loadings_pc1
    scores_pc1   <- -scores_pc1
    message("    PC1 sign flipped so most loadings are positive")
  }

  pve <- summary(pca_fit)$importance[2, ]

  message("    ", name, ": PC1 explains ", round(100 * pve[1], 1),
          "% of the variance in ", length(vars), " quartile-scored pollutants")
  message("    PC1 loadings: ",
          paste(names(loadings_pc1), round(loadings_pc1, 3),
                sep = " = ", collapse = ", "))

  composites <- dat |>
    dplyr::select(dplyr::all_of(id_cols)) |>
    dplyr::mutate(!!name := as.numeric(scores_pc1))

  list(pca_fit = pca_fit, composites = composites,
       loadings = loadings_pc1, pve = pve, name = name, vars = vars)
}


# B. Cross-fitted outcome-informed indices ------------------------------------

## Folds are assigned to PARTICIPANTS, not specimens, so every repeated blood
## draw from one person stays in the same fold and a participant's exposure
## score never depends on their own outcome.
make_participant_folds <- function(ids, k, seed = 42) {
  set.seed(seed)
  tibble::tibble(rand_id = ids,
                 fold = sample(rep_len(seq_len(k), length(ids))))
}

## One row per participant for weight derivation. Exposure is taken at the
## participant's first blood draw; the outcome is participant level anyway.
## The submitted analysis fit the weight models at the specimen level, which
## double-counts participants with repeated draws.
participant_frame <- function(data, vars, covars_weight) {
  data |>
    dplyr::arrange(rand_id, blood_date) |>
    dplyr::distinct(rand_id, .keep_all = TRUE) |>
    dplyr::mutate(event = as01(demcind)) |>
    dplyr::select(dplyr::any_of(c("rand_id", "blood_date", "event", "dcst",
                                  vars, covars_weight))) |>
    tidyr::drop_na()
}

## Score a held-out fold with the TRAINING fold's quartile breakpoints and
## weights. Using the test fold's own quantiles would leak information back in
## and defeat the purpose of cross-fitting.
score_holdout <- function(train, test, vars, weights) {
  scored <- purrr::map(vars, function(v) {
    brks <- stats::quantile(train[[v]], probs = 0:4 / 4, na.rm = TRUE)
    quantile_score(test[[v]], brks) * unname(weights[v])
  })
  Reduce(`+`, scored)
}

## Fit the mixture model in K-1 folds and score the held-out fold.
##
## method:
##   wqs        - gWQS logistic WQS on dementia/CIND
##   qgcomp_glm - qgcomp logistic on dementia/CIND
##   qgcomp_cox - qgcomp Cox on time to dementia/CIND
##
## Weight conventions match the naive composites in 3-clean_data.R:
##   WQS    weights come from final_weights$Estimate (they sum to 1)
##   QGcomp weights are coefficient / sum(|coefficient|), the definition used
##          by get_weights() in scripts/qgcomp_modified.R
crossfit_composite <- function(data, vars, covars_weight,
                               method = c("wqs", "qgcomp_glm", "qgcomp_cox"),
                               k = 5, b_wqs = 200, rh_wqs = 5, seed = 42) {

  method <- match.arg(method)

  pdat  <- participant_frame(data, vars, covars_weight)
  folds <- make_participant_folds(unique(pdat$rand_id), k, seed = seed)

  pdat <- pdat |> dplyr::left_join(folds, by = "rand_id")
  sdat <- data |>
    dplyr::left_join(folds, by = "rand_id") |>
    dplyr::filter(!is.na(fold))

  rhs <- paste(c(vars, covars_weight), collapse = " + ")

  fold_fits <- seq_len(k) |>
    purrr::map(function(f) {
      message("    ", method, " fold ", f, "/", k,
              " (train n = ", sum(pdat$fold != f), " participants)")

      train <- pdat |> dplyr::filter(fold != f)
      test  <- sdat |> dplyr::filter(fold == f)

      weights <- switch(
        method,
        wqs = {
          fit <- run_wqs(
            data = train |> dplyr::mutate(event = factor(event)),
            outcome = "event", mix_name = vars, covariates = covars_weight,
            id_cols = "rand_id", q = 4, validation = 0.6, b = b_wqs,
            b1_pos = TRUE, b_constr = FALSE, rh = rh_wqs,
            family = "binomial", seed = seed
          )
          w <- fit$model$final_weights
          stats::setNames(w$Estimate, as.character(w$mix_name))[vars]
        },
        qgcomp_glm = {
          fit <- qgcomp::qgcomp.glm.noboot(
            stats::as.formula(paste("event ~", rhs)),
            expnms = vars, data = as.data.frame(train), q = 4,
            family = stats::binomial()
          )
          b <- fit$fit$coefficients[vars]
          b / sum(abs(b))
        },
        qgcomp_cox = {
          fit <- qgcomp::qgcomp.cox.noboot(
            stats::as.formula(paste("survival::Surv(dcst, event) ~", rhs)),
            expnms = vars, data = as.data.frame(train), q = 4
          )
          b <- fit$fit$coefficients[vars]
          b / sum(abs(b))
        }
      )

      scores <- test |>
        dplyr::select(rand_id, blood_date) |>
        dplyr::mutate(score = score_holdout(train, test, vars, weights),
                      fold = f)

      weight_tbl <- tibble::tibble(
        fold      = f,
        pollutant = vars,
        weight    = unname(weights[vars])
      )

      list(scores = scores, weights = weight_tbl)
    })

  list(
    scores  = fold_fits |> purrr::map("scores") |> purrr::list_rbind(),
    weights = fold_fits |> purrr::map("weights") |> purrr::list_rbind()
  )
}

#--------------------------------End of the code--------------------------------
