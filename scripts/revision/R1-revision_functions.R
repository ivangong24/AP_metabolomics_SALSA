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
rev_exposure_vars_list <- list(
  total = c("comp_pca_all", "comp_qgcomp_fw_all",
            "comp_wqs_cf_all", "comp_qgcomp_cf_all"),
  cox   = "comp_qgcomp_cox_cf_all"
)

## Pollutants entering the feature-wise quantile g-computation contrast, in the
## order R3 uses for the `all` grouping.
qgcomp_pollutants <- c("exp_benzene", "exp_butadiene", "exp_chromium",
                       "exp_lead", "exp_no2", "exp_nickel", "exp_pm2.5",
                       "exp_nox")

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
is_qgcomp_fw <- function(exposure_var) {
  identical(as.character(exposure_var), QGCOMP_FW_EXPOSURE)
}

## Quartile-score each pollutant (0 .. q-1) within the frame it will be modelled
## in, mirroring qgcomp's default breaks. Ties collapse the upper categories
## rather than erroring, which matters for the small demcind strata.
add_qgcomp_quantiles <- function(data, pollutants = qgcomp_pollutants,
                                 q = QGCOMP_FW_Q) {
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

## Contrast that sums the quantized-pollutant coefficients into psi.
qgcomp_psi_contrast <- function(design_matrix) {
  w <- as.numeric(colnames(design_matrix) %in% qgcomp_q_names)
  if (sum(w) != length(qgcomp_q_names)) {
    stop("Design is missing quantized pollutant terms: ",
         paste(setdiff(qgcomp_q_names, colnames(design_matrix)),
               collapse = ", "))
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
  comp_qgcomp_cox_cf_all   = "QGcomp Cox all toxicants (cross-fitted)"
)

rev_label <- function(exp_name) {
  ifelse(exp_name %in% names(rev_exposure_labels),
         rev_exposure_labels[exp_name], exp_name)
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

## Long axis titles at these sizes run into neighbouring panels -- panel C's
## y-axis title collided with the figure title before this. Wrap them instead
## of shrinking the text back down.
rev_wrap_axis_titles <- function(p, width = 38) {
  if (!is.null(p$labels$y)) {
    p <- p + ggplot2::labs(y = stringr::str_wrap(p$labels$y, width))
  }
  if (!is.null(p$labels$x)) {
    p <- p + ggplot2::labs(x = stringr::str_wrap(p$labels$x, width))
  }
  p
}

rev_save_plot <- function(plot, name, topic, width = 10, height = 7) {
  path <- file.path(rev_dir("figures", topic), paste0(name, ".png"))
  ggplot2::ggsave(path, plot & REV_TEXT, width = width, height = height,
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
