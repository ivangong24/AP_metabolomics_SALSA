## ---------------------------
##
## Script name: R4-mwas_revision.R
## Purpose of script: To perform the R1 revision Metabolome-Wide Association
##                    Study using the unsupervised PCA index and the
##                    cross-fitted WQS / QGcomp indices
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
## Notes: This is the revision analogue of scripts/4-mwas_analysis.R. Same
##        models, same significance criteria, different exposures:
##        1. limma for linear model fitting with empirical Bayes
##        2. PLS with VIP scores for feature selection
##
##        Exposures come from R3-composites_revision.R:
##          total : comp_pca_*        (unsupervised, PRIMARY)
##                  comp_qgcomp_fw_all (feature-wise QGcomp, SECONDARY)
##                  comp_wqs_cf_*     (cross-fitted)
##                  comp_qgcomp_cf_*  (cross-fitted)
##          cox   : comp_pca_*        (unsupervised, incident cohort)
##                  comp_qgcomp_cox_cf_* (cross-fitted)
##
##        FEATURE-WISE QGCOMP. comp_qgcomp_fw_all is not a column produced by
##        R3; it is built here and is a CONTRAST rather than a scored index.
##        Each metabolite is the outcome and the eight quartile-scored
##        pollutants are the exposures, so no outcome model enters at any
##        point -- this is the estimand of Yu et al., Environ Sci Technol
##        2025;59(1):212-223, who used qgcomp.glm.noboot per metabolite.
##        psi is the sum of the eight pollutant coefficients, which for a
##        Gaussian model with no product terms is exactly what
##        qgcomp.glm.noboot returns (verified to machine precision), so it is
##        estimated with contrasts.fit after lmFit. Doing it this way keeps
##        duplicateCorrelation for the repeated draws and the eBayes
##        moderation, neither of which qgcomp offers, and costs one model
##        rather than ~20,000 per stratum.
##
##        It carries limma results only: PLS needs a single exposure vector
##        as Y and a contrast does not provide one, so its VIP columns are NA.
##
##        DUPLICATE CORRELATION. 4-mwas_analysis.R estimates the consensus
##        within-participant correlation separately for every exposure, which
##        is the single most expensive step in the pipeline. The consensus
##        correlation is a property of the blocking structure and the residual
##        variance, and changes only in the fourth decimal place when one
##        exposure column of the design is swapped for another. This script
##        therefore estimates it ONCE per platform x study x population x
##        covariate set and reuses it across the exposures in that stratum.
##        Set DUPCOR_PER_EXPOSURE <- TRUE to reproduce the original behaviour;
##        either way the consensus correlations used are written to
##        revision_output/tables/mwas/duplicate_correlations.xlsx.
##
##        EXPOSURE-RESPONSE EXTRACT. The script also writes
##        revision_output/data/processed/exposure_response_extract_revision.RData
##        -- the in-house-library feature abundances plus the standardized
##        analysis frames -- so that R7 can draw the exposure-response panels
##        Reviewer 1 asked for without reloading the 2.8 GB feature tables.
##
##        Outputs -> revision_output/{data,tables}/...
##        Downstream: R5-annotation_revision.R, R6-pathway_revision.R,
##                    R7-visualization_revision.R (exposure-response extract)
##
##        RUNTIME: hours to days, dominated by duplicateCorrelation and PLS.
##        Validate first with Sys.setenv(SALSA_REVISION_QUICK = "true").
##
## ---------------------------

# Load required packages -----------------------------------------------------

source(here::here("scripts", "1-functions.R"))

# Load exposure and metabolomics data ----------------------------------------

source(here::here("scripts", "2-load_data.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R4-mwas_revision.R")

load(here::here("data", "processed", "salsa_clean.RData"))
load(here::here("data", "processed", "air_toxicants_exposure.RData"))
load(rev_here("data", "processed", "combined_data_list_revision.RData"))

DUPCOR_PER_EXPOSURE <- FALSE

## Workers for the stratum-level duplicateCorrelation stage. Each multisession
## worker receives one feature matrix by value, so this is the memory knob:
##   base R + limma + packages         ~400 MB  (measured in a worker)
##   one feature matrix                ~120 MB  (C18 9,992 x 1,546 = 118 MB;
##                                               HILIC 10,588 x 1,546 = 125 MB)
##   duplicateCorrelation working copy  ~120 MB
##   => ~650 MB per worker, ~10 GB for 16, on top of ~7 GB in the parent.
## That is ~18 GB peak against 48 GB of RAM. There are only 24 strata, so
## going past 16 buys nothing: 16 workers finish in two waves, 8 in three.
DUPCOR_WORKERS      <- 16

## RESUME = TRUE reuses any completed stage already saved under
## revision_output/data/metabolomics/results/, so a crash in PLS does not cost
## you the day of duplicateCorrelation that came before it. Each checkpoint is
## validated against the CURRENT exposure set before it is trusted (see
## checkpoint_ok below) -- a checkpoint built from a different exposure list is
## ignored and recomputed rather than silently mixed in. Set to FALSE to force
## a clean recomputation of everything.
RESUME              <- TRUE

## ADOPT_UNKEYED_CHECKPOINTS is an ESCAPE HATCH, and it should normally be
## FALSE.
##
## Model-key tracking (see partial_checkpoint) was added on 2026-09-03, at the
## same time as the 3- and 10-year exposure arms. Checkpoints written before it
## carry no keys, so there is no way to show that a stored fit belongs to the
## exposure column of the same name in the current run -- and R3 re-estimates
## the cross-fitted indices every time it runs, so "same name" is not "same
## exposure". They are therefore refused.
##
## SALSA_ADOPT_UNKEYED=true accepts them on NAME ALONE. Do that only after
## checking that the exposure columns really did not move -- for example by
## comparing the composites in combined_data_list_revision.RData against a copy
## taken before R3 was re-run. It is an environment variable rather than an
## edit here precisely so that the default in the file stays safe and nobody
## inherits the switch by forgetting to change it back. Every run from then on
## writes keys, so this is a one-time bridge.
ADOPT_UNKEYED_CHECKPOINTS <-
  tolower(Sys.getenv("SALSA_ADOPT_UNKEYED", unset = "false")) %in%
  c("true", "1", "yes")

if (ADOPT_UNKEYED_CHECKPOINTS) {
  message("SALSA_ADOPT_UNKEYED is set: keyless limma / PLS checkpoints will ",
          "be reused on exposure name alone.")
}

PLS_NCOMP           <- 3
FDR_THRESHOLD       <- 0.05
VIP_THRESHOLD       <- 2


# Checkpointing --------------------------------------------------------------

## All stage objects share the nesting [study][population][covar_set][exposure],
## so the exposure names of a saved object can be read back and compared with
## what this run expects.
## Exposures that get a PLS model. comp_qgcomp_fw_all is a contrast, not a
## scored column, and the single pollutants are limma-only by design (see
## SINGLE_POLLUTANT_POPULATIONS in R1), so neither appears in the VIP lists.
## The VIP checkpoint has to be validated against this set rather than against
## the full exposure list, or it is rejected on every resume and PLS is rerun
## for nothing.
pls_eligible <- function(exposures) {
  exposures |>
    purrr::discard(is_qgcomp_fw) |>
    purrr::discard(is_single_pollutant)
}

checkpoint_ok <- function(path, object_names, expected_fn = identity) {
  if (!RESUME || !file.exists(path)) return(FALSE)

  env <- new.env()
  ok <- try(load(path, envir = env), silent = TRUE)
  if (inherits(ok, "try-error") || !all(object_names %in% ls(env))) {
    warning("Checkpoint ", basename(path),
            " is unreadable or incomplete; recomputing.", call. = FALSE)
    return(FALSE)
  }

  if (!signature_ok(env, path)) return(FALSE)

  ## EVERY population, not just the first.
  ##
  ## This used to read names(pop_list[[1]][[1]]) -- the first population of
  ## each study, which is `all` -- and compare it against exposures_for() for
  ## that same population. The exposure set is NOT uniform across populations
  ## (see RESTRICTED_POPULATIONS in R1), so a change that touched only the
  ## strata passed validation and the stale checkpoint was loaded: on
  ## 2026-09-02 the exposure windows were widened to the cognitive strata and
  ## R4 died in "Extracting MWAS results for total_no demcind" with a
  ## vctrs size error, because the loaded fit carried four exposures there and
  ## the design carried seven. Comparing the whole grid is the only check that
  ## catches a change confined to one population.
  exposure_grid <- function(obj, fn) {
    obj |>
      purrr::imap(function(pop_list, study) {
        pop_list |>
          purrr::imap(function(covar_ls, population) {
            sort(as.character(fn(study, population, covar_ls)))
          })
      })
  }

  for (nm in object_names) {
    obj <- env[[nm]]
    stored <- try(
      exposure_grid(obj, function(study, population, covar_ls)
        names(covar_ls[[1]])),
      silent = TRUE)
    if (inherits(stored, "try-error") || is.null(names(obj))) {
      warning("Checkpoint ", basename(path), " has an unexpected structure; ",
              "recomputing.", call. = FALSE)
      return(FALSE)
    }
    expected <- exposure_grid(obj, function(study, population, covar_ls)
      expected_fn(exposures_for(study, population)))

    if (!identical(stored, expected)) {
      warning("Checkpoint ", basename(path), " was built with a different ",
              "exposure set; recomputing.", call. = FALSE)
      return(FALSE)
    }
  }
  TRUE
}

## Every checkpoint carries the covariate specification it was built under.
## The exposure set is validated separately (checkpoint_ok), so the signature
## deliberately covers COVARIATES only: that is what lets a stage be reused
## when an exposure is added but rejected when a covariate is dropped -- which
## is exactly what happened when `batch` was removed on 2026-08-28.
checkpoint_covar_signature <- covars_list_new

signature_ok <- function(env, path) {
  stored <- env[["checkpoint_covar_signature"]]
  if (is.null(stored)) {
    warning("Checkpoint ", basename(path), " predates covariate-signature ",
            "tracking and may have been built under different covariates; ",
            "recomputing.", call. = FALSE)
    return(FALSE)
  }
  if (!identical(stored, checkpoint_covar_signature)) {
    warning("Checkpoint ", basename(path), " was built under a different ",
            "covariate set; recomputing.", call. = FALSE)
    return(FALSE)
  }
  TRUE
}

## Looser variant for stage objects whose values do not depend on the exposure
## set: it checks that the study / population / covariate-set grid matches and
## ignores the innermost exposure names.
checkpoint_strata_ok <- function(path, object_names) {
  if (!RESUME || !file.exists(path)) return(FALSE)

  env <- new.env()
  ok <- try(load(path, envir = env), silent = TRUE)
  if (inherits(ok, "try-error") || !all(object_names %in% ls(env))) {
    warning("Checkpoint ", basename(path),
            " is unreadable or incomplete; recomputing.", call. = FALSE)
    return(FALSE)
  }

  strata_of <- function(obj) {
    obj |>
      purrr::imap(function(pop_list, study) {
        pop_list |>
          purrr::imap(function(covar_ls, population) {
            paste(study, population, names(covar_ls), sep = "|")
          }) |>
          unlist(use.names = FALSE)
      }) |>
      unlist(use.names = FALSE) |>
      sort()
  }

  if (!signature_ok(env, path)) return(FALSE)

  expected <- strata_of(design_c18_list)
  for (nm in object_names) {
    stored <- try(strata_of(env[[nm]]), silent = TRUE)
    if (inherits(stored, "try-error") || !identical(stored, expected)) {
      warning("Checkpoint ", basename(path), " covers a different stratum ",
              "grid; recomputing.", call. = FALSE)
      return(FALSE)
    }
  }
  TRUE
}

resume_load <- function(path, object_names) {
  message("RESUME: loading ", paste(object_names, collapse = ", "),
          " from ", basename(path))
  load(path, envir = globalenv())
}


## Partial resume -------------------------------------------------------------
##
## checkpoint_ok() is all-or-nothing: it compares the stored exposure grid with
## the expected one and rejects on any difference. That is the right rule when
## the stored VALUES would change, and the wrong one when the exposure list
## simply grew. Widening it to the 3- and 10-year windows on 2026-09-03 would
## otherwise have discarded every model that had not changed and recomputed the
## whole grid to add the new cells -- days of PLS for no difference in the
## numbers.
##
## A limma fit and a PLS/VIP result are per exposure and do not depend on which
## OTHER exposures the run carries, so the stored ones can be kept and only the
## missing cells computed. What must still match exactly is everything the
## values DO depend on:
##
##   covariates    signature_ok(), as before
##   feature space hence the REV_QUICK guard -- quick mode subsamples features
##                 and its fits must never be mixed into a full run
##   the exposure  a cell is reused only when the CURRENT run carries an
##                 exposure of that name in that stratum, and the columns
##                 behind a name are set once, above, from the study's `all`
##                 frame
##
## Set RESUME <- FALSE to force everything to be recomputed.
## MODEL KEYS. Matching names is not enough to reuse a fit. `comp_wqs_cf_all`
## is re-estimated by R3 every time it runs, and a cross-fitted index that came
## back even slightly different would silently pair an old fit with a new
## exposure -- a wrong number that no check downstream could catch. So each
## cell carries a hash of everything its result is a function of, and a stored
## result is reused only when that hash still matches.
##
## For limma: the design matrix (which contains the exposure column and every
## covariate), the consensus correlation, and the identity of the feature
## matrix. For PLS: the exposure vector, the covariate frame and the sample
## order.
##
## The feature matrix enters as its dimnames rather than its 118 MB of values.
## Those name the features and the samples, so a changed feature set or sample
## set is caught; a changed VALUE in an unchanged matrix is not, and cannot be
## without hashing gigabytes on every run. The matrices come from the same
## saved .RData on every run, and QUICK -- the one thing that does subsample
## them -- is refused outright above.
feature_key <- function(metabo) {
  rlang::hash(list(rownames(metabo), colnames(metabo)))
}

partial_checkpoint <- function(path, object_name, key_name) {
  if (!RESUME || REV_QUICK || !file.exists(path)) return(NULL)

  env <- new.env()
  ok <- try(load(path, envir = env), silent = TRUE)
  if (inherits(ok, "try-error") || !object_name %in% ls(env)) {
    warning("Checkpoint ", basename(path),
            " is unreadable or incomplete; recomputing.", call. = FALSE)
    return(NULL)
  }
  if (!signature_ok(env, path)) return(NULL)

  if (!key_name %in% ls(env)) {
    if (!ADOPT_UNKEYED_CHECKPOINTS) {
      warning("Checkpoint ", basename(path), " predates model-key tracking, ",
              "so its fits cannot be shown to belong to the current ",
              "exposures; recomputing. Set SALSA_ADOPT_UNKEYED=true to accept ",
              "it on name alone, but only after checking the exposure ",
              "columns did not move.", call. = FALSE)
      return(NULL)
    }
    message("ADOPTING keyless checkpoint ", basename(path),
            " on exposure NAME alone -- this run asserts the exposure ",
            "columns are unchanged.")
    return(list(obj = env[[object_name]], keys = NULL, adopted = TRUE))
  }

  list(obj = env[[object_name]], keys = env[[key_name]], adopted = FALSE)
}

## What a checkpoint already holds for one study x population x covariate cell.
prior_cell <- function(prior, study, population, covar_name) {
  if (is.null(prior)) return(list())
  out <- tryCatch(prior$obj[[study]][[population]][[covar_name]],
                  error = function(e) NULL)
  if (is.null(out) || !is.list(out) || is.null(names(out))) list() else out
}

prior_keys <- function(prior, study, population, covar_name) {
  if (is.null(prior)) return(character(0))
  out <- tryCatch(prior$keys[[study]][[population]][[covar_name]],
                  error = function(e) NULL)
  if (is.null(out)) character(0) else unlist(out)
}

## `x[[name]]` raises "subscript out of bounds" for a name that is not there --
## for a named list as well as for an atomic vector -- so every key lookup goes
## through this. A missing key is a reason to recompute the cell, not to abort
## the run.
key_of <- function(x, nm) {
  if (is.null(x) || !nm %in% names(x)) return(NA_character_)
  v <- x[[nm]]
  if (length(v) != 1 || is.na(v)) NA_character_ else as.character(v)
}

## The exposures of one cell that still have to be computed, and the stored
## results for the rest, in the order the current run expects them.
##
## `keys_now` is the current model key per exposure; an exposure whose key has
## moved is recomputed even though its name is unchanged.
split_cell <- function(prior, study, population, covar_name, wanted,
                       keys_now) {
  keep <- prior_cell(prior, study, population, covar_name)
  old  <- prior_keys(prior, study, population, covar_name)

  reusable <- if (isTRUE(prior$adopted)) {
    intersect(wanted, names(keep))
  } else {
    intersect(wanted, names(keep)) |>
      purrr::keep(function(e){
        a <- key_of(old, e)
        b <- key_of(keys_now, e)
        !is.na(a) && !is.na(b) && identical(a, b)
      })
  }

  stale <- setdiff(intersect(wanted, names(keep)), reusable)
  if (length(stale) > 0) {
    message("    ", length(stale), " stored result(s) no longer match their ",
            "model and will be refitted: ", paste(stale, collapse = ", "))
  }

  list(keep = keep[reusable], todo = setdiff(wanted, reusable))
}

## Say what a partial resume is reusing before it starts, so a run that quietly
## recomputes everything is visible in the log rather than only in the clock.
report_partial <- function(prior, label, design_data_list, keys) {
  if (is.null(prior)) {
    message(label, ": no reusable checkpoint - computing every cell")
    return(invisible(NULL))
  }
  have <- 0L
  want <- 0L
  design_data_list |> purrr::iwalk(function(pop_list, study){
    pop_list |> purrr::iwalk(function(covar_ls, population){
      covar_ls |> purrr::iwalk(function(designls, covar_name){
        want <<- want + length(designls)
        have <<- have + length(split_cell(
          prior, study, population, covar_name, names(designls),
          tryCatch(keys[[study]][[population]][[covar_name]],
                   error = function(e) NULL))$keep)
      })
    })
  })
  message(label, ": reusing ", have, " of ", want,
          " exposure x stratum cells from the checkpoint, computing ",
          want - have)
  invisible(NULL)
}


# Prepare exposure data for MWAS ---------------------------------------------

## R3-composites_revision.R has already restricted each study to the exposures
## it should carry, so no further column pruning is needed here.
combined_data_list_new <- combined_data_list_revision

## Post-diagnosis specimens, as a sensitivity population ----------------------
##
## 138 of the 1,546 total-cohort specimens (8.9%), from 102 participants, were
## drawn AFTER a dementia/CIND diagnosis. In those the metabolome may reflect
## the disease, its treatment, or altered diet and activity rather than the
## exposure -- Reviewer 1 comment 2. They are KEPT in the primary analysis and
## removed in an additional `all predx` population, which leaves 1,408
## specimens from 873 participants.
##
## Diagnosis date is baseline + dcst years for incident cases. Participants who
## are dementia/CIND in the total cohort but absent from the incident cohort
## are prevalent cases, diagnosed at or before enrolment, so every specimen of
## theirs is post-diagnosis.
##
## Added to `total` only: the `cox` arm is already restricted to the incident
## cohort and carries one exposure, and each extra population costs four
## duplicateCorrelation strata. Set PREDX_STUDIES to add it elsewhere.
PREDX_STUDIES <- "total"

dx_dates <- salsa_clean_cox |>
  dplyr::select(rand_id, bl_date, dcst, demcind_incident = demcind) |>
  dplyr::mutate(dx_date = dplyr::if_else(
    demcind_incident == "Dementia/CIND",
    as.Date(bl_date) + dcst * 365.25, as.Date(NA)))

drop_post_diagnosis <- function(data) {
  data |>
    dplyr::left_join(dplyr::select(dx_dates, rand_id, dx_date,
                                   demcind_incident),
                     by = "rand_id") |>
    dplyr::filter(!(
      ## prevalent case: in the total cohort as demcind, absent from the
      ## incident cohort, so diagnosed at or before baseline
      (is.na(demcind_incident) & demcind == "Dementia/CIND") |
      (!is.na(dx_date) & as.Date(blood_date) > dx_date)
    )) |>
    dplyr::select(-dx_date, -demcind_incident)
}

combined_data_list_new[PREDX_STUDIES] <-
  combined_data_list_new[PREDX_STUDIES] |>
  purrr::map(function(datalist){
    predx <- datalist[["all"]] |> purrr::map(drop_post_diagnosis)
    message("  post-diagnosis exclusion: ", nrow(datalist[["all"]][["covar"]]),
            " -> ", nrow(predx[["covar"]]), " specimens (",
            dplyr::n_distinct(predx[["covar"]]$rand_id), " participants)")
    c(datalist, list(`all predx` = predx))
  })

## Standardize the scored composites to SD = 1 -------------------------------
##
## Reviewer 1 minor comments 5 and 6. The composites are built on unrelated
## scales -- comp_pca_all has SD 1.92 (prcomp scales the inputs, not the
## resulting scores), comp_wqs_cf_all 0.79, comp_qgcomp_cf_all 0.42 -- so
## their coefficients and discovery counts are not comparable, and a
## coefficient "per unit" means something different for each. Dividing by the
## SD makes every coefficient the difference in log2 abundance per one-SD
## increase in that index, which is what the response letter states.
##
## This is a pure linear rescaling of the exposure: the coefficient and its
## standard error are divided by the same constant, so t statistics, p-values
## and FDR are unchanged to machine precision. Only the effect-size scale
## moves.
##
## The SD is taken ONCE per study, from that study's `all` frame, and applied
## to every population within it. Standardizing each stratum by its own SD
## would make "one SD" a different quantity in each, which would silently
## distort the subgroup-versus-full-cohort comparison in the combined panels.
##
## comp_qgcomp_fw_all is deliberately excluded: it is a contrast, not a scored
## column, and psi is already in an interpretable unit (a one-quartile
## increase in every pollutant simultaneously). Its scale is therefore not
## comparable with the standardized indices, and the Methods should say so.
composite_sds <- combined_data_list_new |>
  purrr::map(function(datalist){
    ref <- datalist[["all"]][["covar"]]
    ref |>
      dplyr::select(dplyr::starts_with("comp_")) |>
      purrr::map_dbl(~ stats::sd(.x, na.rm = TRUE))
  })

combined_data_list_new <- list(combined_data_list_new, composite_sds) |>
  purrr::pmap(function(datalist, sds){
    datalist |>
      purrr::map(function(data_list){
        data_list |>
          purrr::map(function(d){
            for (nm in names(sds)) {
              if (nm %in% names(d) && is.finite(sds[[nm]]) && sds[[nm]] > 0) {
                d[[nm]] <- d[[nm]] / sds[[nm]]
              }
            }
            d
          })
      })
  })

purrr::iwalk(composite_sds, function(sds, study){
  message("  ", study, " composite SDs used for standardization: ",
          paste(names(sds), round(sds, 4), sep = " = ", collapse = ", "))
})

rev_dir("tables", "composites")
writexl::write_xlsx(
  composite_sds |>
    purrr::imap(function(sds, study){
      tibble::tibble(study = study, composite = names(sds),
                     sd_used = unname(sds))
    }) |>
    purrr::list_rbind(),
  rev_here("tables", "composites", "composite_standardization_sd.xlsx"))


## Re-derive the single-pollutant IQR scaling, once per study ------------------
##
## 3-clean_data.R builds exp_*_iqr as `.x / IQR(.x)` INSIDE each study x
## population frame, so the divisor is the stratum's own IQR and
## IQR(exp_*_iqr) == 1 in every stratum. "One IQR" is therefore a different
## physical quantity in a subgroup than in the full cohort -- the raw NOx IQR
## is 2.208 in total/all against 2.750 in total/demcind, a 25% gap -- which
## would silently distort any subgroup-versus-full-cohort comparison of a
## pollutant coefficient.
##
## The composites avoid this immediately above by taking the SD once from the
## study's `all` frame and applying it everywhere. This does the same for the
## eight pollutants, so the two conventions match.
##
## For total/all this is an arithmetic no-op -- the divisor IS that frame's
## IQR -- so it does not change the current results, in which the single
## pollutants are fitted in `all` only. It is what makes widening
## SINGLE_POLLUTANT_POPULATIONS safe.
##
## Every window, not only the 5-year primary: the 3- and 10-year windows carry
## their own single-pollutant exposures, and they get the same convention --
## the divisor is that window's IQR in the study's `all` frame, applied
## unchanged to every population within it.
pollutant_iqrs <- combined_data_list_new |>
  purrr::map(function(datalist){
    ref <- datalist[["all"]][["covar"]]
    intersect(all_raw_pollutants, names(ref)) |>
      purrr::set_names() |>
      purrr::map_dbl(~ stats::IQR(ref[[.x]], na.rm = TRUE))
  })

combined_data_list_new <- list(combined_data_list_new, pollutant_iqrs) |>
  purrr::pmap(function(datalist, iqrs){
    datalist |>
      purrr::map(function(data_list){
        data_list |>
          purrr::map(function(d){
            for (nm in names(iqrs)) {
              col <- paste0(nm, "_iqr")
              if (nm %in% names(d) && is.finite(iqrs[[nm]]) && iqrs[[nm]] > 0) {
                d[[col]] <- d[[nm]] / iqrs[[nm]]
              }
            }
            d
          })
      })
  })

purrr::iwalk(pollutant_iqrs, function(iqrs, study){
  message("  ", study, " pollutant IQRs used for scaling (from the `all` ",
          "frame): ", paste(names(iqrs), signif(iqrs, 4), sep = " = ",
                            collapse = ", "))
})

writexl::write_xlsx(
  pollutant_iqrs |>
    purrr::imap(function(iqrs, study){
      tibble::tibble(study = study, pollutant = names(iqrs),
                     iqr_used = unname(iqrs))
    }) |>
    purrr::list_rbind(),
  rev_here("tables", "composites", "pollutant_standardization_iqr.xlsx"))


## Add the quartile-scored pollutant columns the feature-wise QGcomp contrast
## needs. Breaks are taken within each analysis frame, so a stratum is scored
## against its own exposure distribution rather than the pooled one.
## One set per window, so comp_qgcomp_fw_all sums exp_*_q and
## comp_qgcomp_fw_all_w10 sums exp_*_w10_q. Windows whose columns are absent
## (the cox frames carry none) are skipped rather than erroring.
combined_data_list_new <- combined_data_list_new |>
  purrr::map(function(datalist){
    datalist |>
      purrr::map(function(data_list){
        data_list |> purrr::map(add_all_qgcomp_quantiles)
      })
  })

exposure_vars_list <- combined_data_list_new |>
  purrr::imap(function(datalist, study){
    frame   <- datalist[["all"]][["covar"]]
    ## Composites, plus the IQR-scaled single pollutants. The latter are not
    ## built by R3 -- they come through from 3-clean_data.R -- so they are
    ## matched by name against the frame rather than by prefix.
    present <- frame |> dplyr::select(starts_with("comp_")) |> colnames()
    present <- c(present,
                 intersect(all_single_pollutant_exposures, colnames(frame)))
    ## The feature-wise QGcomp contrasts are not columns, so one is "present"
    ## whenever the quantized pollutants of its window are. Checked per window:
    ## the 3-year window carries seven quantized columns and the 10-year eight,
    ## and a window whose columns never arrived must not be promised
    ## downstream.
    for (w in ALL_MWAS_WINDOWS) {
      q_names <- qgcomp_q_names_for(w)
      if (length(q_names) > 0 && all(q_names %in% colnames(frame))) {
        present <- c(present, qgcomp_fw_exposure_for(w))
      }
    }
    intersect(rev_exposure_vars_list[[study]], present)
  })


## finalize the metabolite matrix

list(
  list(salsa_blood_date_c18, salsa_blood_date_hilic),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |>
  purrr::pmap(function(data1, data2){
    data1 |>
      select(rand_id, blood_date, wave, batch, file.name_new) |>
      filter(file.name_new %in% colnames(data2)) |>
      filter(!is.na(rand_id)) |>
      distinct()
  }) |>
  set_names("sample_link_c18", "sample_link_hilic") |>
  list2env(envir = .GlobalEnv)


list(
  list(sample_link_c18, sample_link_hilic),
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed)
) |>
  purrr::pmap(function(sample_link, metabo){
    combined_data_list_new |>
      purrr::map(function(datalist){
        datalist |>
          purrr::map(function(data_list){
            link <- data_list[["covar"]] |>
              dplyr::select(rand_id, blood_date) |>
              dplyr::inner_join(sample_link, by = c("rand_id", "blood_date")) |>
              dplyr::arrange(match(file.name_new, colnames(metabo)))

            metabo |>
              dplyr::select(all_of(link$file.name_new))
          })
      })
  }) |>
  purrr::set_names("metabo_list_c18_final", "metabo_list_hilic_final") |>
  list2env(.GlobalEnv)


## In QUICK mode keep a random subset of features so the whole script can be
## validated end to end in minutes. The same features are kept in every
## study x population slice. Nothing produced this way is reportable.
if (REV_QUICK) {
  QUICK_N_FEATURES <- 50
  set.seed(42)

  list(
    list(metabo_list_c18_final, metabo_list_hilic_final),
    list("C18", "HILIC")
  ) |>
    purrr::pmap(function(metabo_data_list, mode){
      keep <- rownames(metabo_data_list[["total"]][["all"]]) |>
        sample(min(QUICK_N_FEATURES,
                   nrow(metabo_data_list[["total"]][["all"]])))
      message("QUICK: keeping ", length(keep), " ", mode, " features")

      metabo_data_list |>
        purrr::map(function(metabo_list){
          metabo_list |> purrr::map(~ .x[keep, , drop = FALSE])
        })
    }) |>
    purrr::set_names("metabo_list_c18_final", "metabo_list_hilic_final") |>
    list2env(.GlobalEnv)
}


## Merge air toxicants exposure with sample link files

list(
  list(sample_link_c18, sample_link_hilic),
  list(metabo_list_c18_final, metabo_list_hilic_final)
) |>
  purrr::pmap(function(sample_link, metabo_list){
    list(combined_data_list_new,
         metabo_list) |>
      purrr::pmap(function(datalist, metabolist){
        list(datalist, metabolist) |>
          purrr::pmap(function(data_list, metabo){
            data_list |>
              purrr::map(function(data){
                data |>
                  dplyr::inner_join(sample_link, by = c("rand_id", "blood_date",
                                                        "wave", "batch")) |>
                  dplyr::arrange(match(file.name_new, colnames(metabo)))
              })
          })
      })
  }) |>
  purrr::set_names("combined_data_list_c18", "combined_data_list_hilic") |>
  list2env(.GlobalEnv)


# Verify sample ordering -----------------------------------------------------

list(
  list("C18", "HILIC"),
  list(combined_data_list_c18, combined_data_list_hilic),
  list(metabo_list_c18_final, metabo_list_hilic_final)
) |>
  purrr::pmap(function(mode, combined_data_list, metabo_list) {
    message(paste0(mode, " sample ordering check:"))
    list(combined_data_list, metabo_list) |>
      purrr::pmap(function(datalist, metabolist){
        list(datalist, metabolist) |>
          purrr::pmap(function(combined_dflist, metabo) {
            combined_dflist |>
              purrr::map(function(combined_data){
                print(table(combined_data$file.name_new == colnames(metabo)))
              })
          })
      })
  }) |>
  invisible()

# Exposure-response extract for the identified metabolites ---------------------
#
## Reviewer 1 comment 14 asks that the freed space in the main figure carry
## biological content, one option being "exposure-response plots for the Level 1
## metabolites". Drawing those needs three things that exist only here: the
## feature abundances, the sample-to-participant link, and the analysis frames
## AFTER the SD / IQR standardization applied above. R7 loads neither the
## feature tables (2.8 GB) nor the pre-standardization frames, so this section
## writes a compact extract for it.
##
## Only the features carrying an Emory in-house-library annotation are kept.
## That is the superset from which Level 1 is drawn -- an accurate-mass-only
## database hit can never reach Level 1 (see annotation_confidence_level() in
## R1) -- so the extract stays small (~150 features) and does not have to be
## rebuilt if the confidence rule is refined in R5.
##
## Abundances are on the same log2 scale limma models, and the exposure columns
## are the standardized ones the design matrices were built from, so a slope
## drawn from this extract is in the units the coefficients are reported in.

message("\n=== Exposure-response extract for the annotated features ===")

load(here::here("data", "metabolomics", "annotation",
                "annotation_cleaned_wide.RData"))

inhouse_feature_ids <- list(annotation_c18_wide, annotation_hilic_wide) |>
  purrr::map(function(annot){
    annot |>
      dplyr::filter(reference == "In House Library") |>
      dplyr::pull(id) |>
      unique()
  }) |>
  purrr::set_names("c18", "hilic")

purrr::iwalk(inhouse_feature_ids, function(ids, platform){
  message("  ", platform, ": ", length(ids),
          " in-house-library features requested for the extract")
})

## Wide frame: one row per sample, one column per in-house feature, keyed by
## the same rand_id / blood_date pair every analysis frame carries.
build_er_abundance <- function(metabo, sample_link, feature_ids) {
  keep <- intersect(feature_ids, rownames(metabo))
  cols <- intersect(sample_link$file.name_new, colnames(metabo))
  if (length(keep) == 0 || length(cols) == 0) return(NULL)

  metabo[keep, cols, drop = FALSE] |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("file.name_new") |>
    dplyr::inner_join(
      sample_link |> dplyr::select(rand_id, blood_date, file.name_new),
      by = "file.name_new") |>
    dplyr::relocate(rand_id, blood_date, file.name_new)
}

list(
  list(med_c18_raw_combat_processed, med_hil_raw_combat_processed),
  list(sample_link_c18, sample_link_hilic),
  inhouse_feature_ids
) |>
  purrr::pmap(build_er_abundance) |>
  purrr::set_names("er_abundance_c18", "er_abundance_hilic") |>
  list2env(.GlobalEnv)

purrr::iwalk(list(c18 = er_abundance_c18, hilic = er_abundance_hilic),
             function(df, platform){
               message("  ", platform, " extract: ", nrow(df), " samples x ",
                       ncol(df) - 3, " features")
             })

## The analysis frames go with it. combined_data_list_* is what the design
## matrices were built from, so it carries the standardized composites, the
## re-derived pollutant IQR columns, the quartile scores and every covariate --
## which is what lets R7 residualize on exactly the model's adjustment set.
er_data_list_c18   <- combined_data_list_c18
er_data_list_hilic <- combined_data_list_hilic
er_covars_list     <- covars_list_new

er_extract_path <- rev_here("data", "processed",
                            "exposure_response_extract_revision.RData")

save(er_abundance_c18, er_abundance_hilic,
     er_data_list_c18, er_data_list_hilic, er_covars_list,
     file = er_extract_path)

message("  extract -> ", er_extract_path)


message("Exposure variables for the revision MWAS:")
print(exposure_vars_list)

message("Covariates for MWAS:")
print(covars_list_new)


# =============================================================================
# SECTION 1: LIMMA-BASED MWAS
# =============================================================================

# Create design matrices for each exposure -----------------------------------

## Function to create design matrix for a single exposure
##
## For every scored index this is one exposure column plus covariates. For
## comp_qgcomp_fw_all the eight quantized pollutants enter as separate terms
## and are collapsed to psi by a contrast after lmFit (see fit_limma).
create_design_matrix <- function(combined_data, exposure_var, covars) {
  exposure_terms <- if (is_qgcomp_fw(exposure_var)) {
                      qgcomp_q_names_for(exposure_window(exposure_var))
                    } else exposure_var
  formula_matrix <- as.formula(str_c("~ ",
                                     paste(exposure_terms, collapse = " + "),
                                     " + ", paste(covars, collapse = " + ")))
  model.matrix(formula_matrix, data = combined_data)
}


## Create design matrices for all exposures
list(
  list(combined_data_list_c18, combined_data_list_hilic),
  list("C18", "Hilic")
) |>
  purrr::pmap(function(combined_data_list, mode){
    list(combined_data_list, names(combined_data_list)) |>
      purrr::pmap(function(combined_df_list, study){
        combined_df_list |>
          purrr::imap(function(combined_dflist, population){
            list(combined_dflist, covars_list_new, names(covars_list_new)) |>
              purrr::pmap(function(combined_data, covars, covar_name){
                message(paste0("Creating design matrices for ", study,
                               "_", population, " in ", mode,
                               " with covariates set: ", covar_name, " ..."))

                ## exposures_for(), not exposure_vars_list[[study]]: the single
                ## pollutants are fitted in `all` only, so the exposure set now
                ## varies by population and this is the one place that decides
                ## it. Everything downstream reads names(designls).
                exposures_for(study, population) |>
                  purrr::set_names() |>
                  purrr::map(~ create_design_matrix(combined_data, .x, covars))
              })
          })
      })
  }) |>
  purrr::set_names("design_c18_list", "design_hilic_list") |>
  list2env(.GlobalEnv)


# Estimate within-subject correlation for duplicate measures -------------------

## Set up parallel backend
n_workers <- max(1, future::availableCores() - 1)
message(paste0("Setting up parallel plan with ", n_workers, " workers..."))

## The whole stage is checkpointed: this is the expensive one, and a crash
## anywhere downstream should not cost it. See RESUME at the top.
dupcor_path <- rev_here("data", "metabolomics", "results",
                        "duplicate_correlation_revision.RData")

## When the consensus correlation is estimated per stratum rather than per
## exposure it is a property of the blocking structure and the residual
## variance, not of which exposure column the design carries. Adding an
## exposure therefore does not invalidate it -- so accept a checkpoint whose
## stratum grid matches and simply re-broadcast it over the current exposure
## set, rather than discarding a day of duplicateCorrelation over a name.
## Per-exposure estimates are exposure-specific and get the strict check.
dupcor_reusable <- if (DUPCOR_PER_EXPOSURE) {
  checkpoint_ok(dupcor_path, c("dupcor_c18_list", "dupcor_hilic_list"))
} else {
  checkpoint_strata_ok(dupcor_path, c("dupcor_c18_list", "dupcor_hilic_list"))
}

if (dupcor_reusable) {

  resume_load(dupcor_path, c("dupcor_c18_list", "dupcor_hilic_list"))

  if (!DUPCOR_PER_EXPOSURE) {
    rebroadcast <- function(dupcor_data_list, design_data_list) {
      list(dupcor_data_list, design_data_list) |>
        purrr::pmap(function(dupcor_list, design_list){
          list(dupcor_list, design_list) |>
            purrr::pmap(function(dupcor_ls, design_ls){
              list(dupcor_ls, design_ls) |>
                purrr::pmap(function(dupcorls, designls){
                  designls |> purrr::map(~ dupcorls[[1]])
                })
            })
        })
    }
    dupcor_c18_list   <- rebroadcast(dupcor_c18_list,   design_c18_list)
    dupcor_hilic_list <- rebroadcast(dupcor_hilic_list, design_hilic_list)
    message("RESUME: stratum-level consensus correlations re-broadcast over ",
            length(unlist(purrr::map_depth(design_c18_list, 3, names))),
            " exposure x stratum cells")
  }

} else {

  future::plan(future::multisession, workers = n_workers)

  if (DUPCOR_PER_EXPOSURE) {

    system.time({
      list(
        list("C18", "HILIC"),
        list(metabo_list_c18_final, metabo_list_hilic_final),
        list(design_c18_list, design_hilic_list),
        list(combined_data_list_c18, combined_data_list_hilic)
      ) |>
        purrr::pmap(function(mode, metabo_data_list, design_data_list,
                             combined_data_list) {
          list(design_data_list, metabo_data_list,
               combined_data_list, names(combined_data_list)) |>
            purrr::pmap(function(design_list, metabo_list,
                                 combined_df_list, study){
              list(design_list, metabo_list,
                   combined_df_list, names(combined_df_list)) |>
                purrr::pmap(function(design_ls, metabo,
                                     combined_dflist, population){
                  list(design_ls, combined_dflist, names(combined_dflist)) |>
                    purrr::pmap(function(designls, combined_data, covar_name){
                      message(paste0("Calculating duplicate correlation for ",
                                     mode, " in ", study, "_", population,
                                     " with covariates set: ", covar_name, " ..."))

                      block <- combined_data$rand_id
                      stopifnot(
                        length(block) == ncol(metabo),
                        all(combined_data$file.name_new == colnames(metabo)))

                      designls |>
                        furrr::future_map(function(design) {
                          limma::duplicateCorrelation(metabo, design,
                                                      block = block)
                        }, .options = furrr_options(seed = TRUE),
                        .progress = TRUE)
                    })
                })
            })
        }) |>
        purrr::set_names("dupcor_c18_list", "dupcor_hilic_list") |>
        list2env(.GlobalEnv)
    })

  } else {

    ## One job per platform x study x population x covariate set. The design of
    ## the first exposure in the stratum carries the estimate.
    dupcor_jobs <- list(
      list(mode = "C18",   metabo = metabo_list_c18_final,
           design = design_c18_list,   combined = combined_data_list_c18),
      list(mode = "HILIC", metabo = metabo_list_hilic_final,
           design = design_hilic_list, combined = combined_data_list_hilic)
    ) |>
      purrr::map(function(platform){
        platform$design |>
          purrr::imap(function(design_list, study){
            design_list |>
              purrr::imap(function(design_ls, population){
                design_ls |>
                  purrr::imap(function(designls, covar_name){
                    metabo <- platform$metabo[[study]][[population]]
                    combined_data <-
                      platform$combined[[study]][[population]][[covar_name]]

                    stopifnot(
                      nrow(combined_data) == ncol(metabo),
                      all(combined_data$file.name_new == colnames(metabo)))

                    list(platform   = platform$mode,
                         study      = study,
                         population = population,
                         covar_set  = covar_name,
                         exposures  = names(designls),
                         design     = designls[[1]],
                         metabo     = metabo,
                         block      = combined_data$rand_id)
                  })
              })
          })
      }) |>
      purrr::list_flatten() |>
      purrr::list_flatten() |>
      purrr::list_flatten()

    names(dupcor_jobs) <- dupcor_jobs |>
      purrr::map_chr(~ paste(.x$platform, .x$study, .x$population, .x$covar_set,
                             sep = "|"))

    message("Estimating duplicate correlation for ", length(dupcor_jobs),
            " strata ...")

    ## Each future carries one metabolite matrix, so raise the globals ceiling
    ## and keep the worker count modest to bound peak memory.
    options(future.globals.maxSize = 4 * 1024^3)
    future::plan(future::multisession, workers = min(n_workers, DUPCOR_WORKERS))

    system.time({
      dupcor_flat <- dupcor_jobs |>
        furrr::future_map(function(job) {
          limma::duplicateCorrelation(job$metabo, job$design, block = job$block)
        },
        # one stratum per chunk: they differ in sample size, so static chunking
        # leaves workers idle at the end
        .options = furrr_options(seed = TRUE, chunk_size = 1),
        .progress = TRUE)
    })

    ## duplicateCorrelation returns a NaN consensus when the trimmed mean of the
    ## per-feature atanh correlations hits non-finite values -- it happens in the
    ## small strata (demcind) and whenever few features are analysed. lmFit()
    ## then dies on `if (abs(correlation) >= 1)`, so repair it here: first
    ## recompute the trimmed mean ignoring the non-finite features, then, only if
    ## that still fails, borrow the median consensus from the other strata.
    repair_consensus <- function(dupcor, fallback = NA_real_) {
      if (is.finite(dupcor$consensus.correlation)) return(dupcor)

      ac <- dupcor$atanh.correlations
      ac <- ac[is.finite(ac)]
      if (length(ac) > 0) {
        dupcor$consensus.correlation <- tanh(mean(ac, trim = 0.15))
      }
      if (!is.finite(dupcor$consensus.correlation)) {
        dupcor$consensus.correlation <- fallback
      }
      dupcor$cor <- dupcor$consensus.correlation
      dupcor
    }

    dupcor_flat <- dupcor_flat |> purrr::map(repair_consensus)

    consensus_median <- dupcor_flat |>
      purrr::map_dbl(~ .x$consensus.correlation) |>
      stats::median(na.rm = TRUE)

    dupcor_flat <- dupcor_flat |>
      purrr::imap(function(dupcor, stratum){
        if (is.finite(dupcor$consensus.correlation)) return(dupcor)
        warning("Non-finite consensus correlation for stratum ", stratum,
                "; substituting the median across strata (",
                round(consensus_median, 5), ").", call. = FALSE)
        repair_consensus(dupcor, fallback = consensus_median)
      })

    stopifnot(all(purrr::map_lgl(dupcor_flat,
                                 ~ is.finite(.x$consensus.correlation))))

    ## Fan the stratum-level estimate back out over that stratum's exposures, so
    ## the nested structure downstream is identical to the per-exposure version.
    expand_dupcor <- function(design_data_list, mode) {
      design_data_list |>
        purrr::imap(function(design_list, study){
          design_list |>
            purrr::imap(function(design_ls, population){
              design_ls |>
                purrr::imap(function(designls, covar_name){
                  dupcor <- dupcor_flat[[paste(mode, study, population,
                                               covar_name, sep = "|")]]
                  message(mode, " ", study, "_", population, " [", covar_name,
                          "]: consensus correlation = ",
                          round(dupcor$consensus.correlation, 5),
                          " (estimated on ", names(designls)[1],
                          ", reused for ", length(designls), " exposures)")
                  designls |> purrr::map(~ dupcor)
                })
            })
        })
    }

    dupcor_c18_list   <- expand_dupcor(design_c18_list,   "C18")
    dupcor_hilic_list <- expand_dupcor(design_hilic_list, "HILIC")

    rm(dupcor_jobs, dupcor_flat)
    gc()
  }

  rev_dir("data", "metabolomics", "results")
  save(dupcor_c18_list, dupcor_hilic_list, checkpoint_covar_signature,
       file = dupcor_path)
  message("Checkpoint written: ", dupcor_path)
}


## Record the consensus correlations actually used
dupcor_table <- list(
  list("C18", "HILIC"),
  list(dupcor_c18_list, dupcor_hilic_list)
) |>
  purrr::pmap(function(mode, dupcor_data_list){
    dupcor_data_list |>
      purrr::imap(function(dupcor_list, study){
        dupcor_list |>
          purrr::imap(function(dupcor_ls, population){
            dupcor_ls |>
              purrr::imap(function(dupcorls, covar_name){
                tibble::tibble(
                  platform    = mode,
                  study       = study,
                  population  = population,
                  covar_set   = covar_name,
                  exposure    = names(dupcorls),
                  consensus_correlation = purrr::map_dbl(
                    dupcorls, ~ .x$consensus.correlation)
                )
              }) |>
              purrr::list_rbind()
          }) |>
          purrr::list_rbind()
      }) |>
      purrr::list_rbind()
  }) |>
  purrr::list_rbind()

rev_save_table(dupcor_table, "duplicate_correlations", "mwas")


# Fit limma models -----------------------------------------------------------

## Function to fit limma model with duplicate correlation
##
## For comp_qgcomp_fw_all the eight quantized-pollutant coefficients are summed
## into psi by contrasts.fit BEFORE eBayes, so the moderated t and its FDR
## refer to the mixture effect rather than to any single pollutant. The
## resulting fit carries one coefficient, named "psi".
fit_limma <- function(metabolome_matrix, design_matrix, block, correlation,
                      exposure_var = NULL) {
  fit <- limma::lmFit(metabolome_matrix, design_matrix,
                      block = block, correlation = correlation)
  if (is_qgcomp_fw(exposure_var)) {
    fit <- limma::contrasts.fit(
      fit, qgcomp_psi_contrast(design_matrix, exposure_window(exposure_var)))
  }
  fit <- limma::eBayes(fit)
  return(fit)
}

## Checkpointed: lmFit + eBayes over 144 models is cheap next to
## duplicateCorrelation but still worth not repeating.
fit_c18_path   <- rev_here("data", "metabolomics", "results",
                           "limma_fit_c18_revision.RData")
fit_hilic_path <- rev_here("data", "metabolomics", "results",
                           "limma_fit_hilic_revision.RData")

## The model key of every limma cell this run would fit: the design matrix, the
## consensus correlation it is fitted with, and the identity of the feature
## matrix. See feature_key() above for why the matrix enters by its dimnames.
build_limma_keys <- function(design_data_list, dupcor_data_list,
                             metabo_data_list) {
  list(design_data_list, dupcor_data_list, metabo_data_list,
       names(design_data_list)) |>
    purrr::pmap(function(design_list, dupcor_list, metabo_list, study){
      list(design_list, dupcor_list, metabo_list, names(design_list)) |>
        purrr::pmap(function(design_ls, dupcor_ls, metabo, population){
          fkey <- feature_key(metabo)
          list(design_ls, dupcor_ls, names(design_ls)) |>
            purrr::pmap(function(designls, dupcorls, covar_name){
              list(designls, dupcorls, names(designls)) |>
                purrr::pmap(function(design, dupcor, exposure_var){
                  rlang::hash(list(design, dupcor$consensus.correlation,
                                   fkey, exposure_var))
                })
            }) |>
            purrr::set_names(names(design_ls))
        }) |>
        purrr::set_names(names(design_list))
    }) |>
    purrr::set_names(names(design_data_list))
}

limma_keys_c18   <- build_limma_keys(design_c18_list, dupcor_c18_list,
                                     metabo_list_c18_final)
limma_keys_hilic <- build_limma_keys(design_hilic_list, dupcor_hilic_list,
                                     metabo_list_hilic_final)

## Partial resume: keep every fit the checkpoint already holds whose model key
## still matches, and fit only the rest.
prior_limma <- list(
  C18   = partial_checkpoint(fit_c18_path,   "limma_fit_c18",
                             "limma_keys_c18"),
  HILIC = partial_checkpoint(fit_hilic_path, "limma_fit_hilic",
                             "limma_keys_hilic")
)

limma_keys <- list(C18 = limma_keys_c18, HILIC = limma_keys_hilic)

report_partial(prior_limma$C18,   "limma C18",   design_c18_list,
               limma_keys_c18)
report_partial(prior_limma$HILIC, "limma HILIC", design_hilic_list,
               limma_keys_hilic)

future::plan(future::multicore, workers = n_workers)

system.time({
  list(
    list("C18", "HILIC"),
    list(design_c18_list, design_hilic_list),
    list(metabo_list_c18_final, metabo_list_hilic_final),
    list(dupcor_c18_list, dupcor_hilic_list),
    list(combined_data_list_c18, combined_data_list_hilic)
  ) |>
    purrr::pmap(function(mode, design_data_list, metabo_data_list,
                         dupcor_data_list, combined_data_list) {
      prior <- prior_limma[[mode]]

      list(design_data_list, metabo_data_list,
           dupcor_data_list, combined_data_list, names(combined_data_list)) |>
        purrr::pmap(function(design_list, metabo_list, dupcor_list,
                             combined_df_list, study){
          list(design_list, metabo_list, dupcor_list,
               combined_df_list, names(combined_df_list)) |>
            purrr::pmap(function(design_ls, metabo, dupcor_ls,
                                 combined_dflist, population){
              list(design_ls, dupcor_ls, combined_dflist,
                   names(combined_dflist)) |>
                purrr::pmap(function(designls, dupcorls,
                                     combined_data, covar_name){

                  parts <- split_cell(
                    prior, study, population, covar_name, names(designls),
                    limma_keys[[mode]][[study]][[population]][[covar_name]])

                  message(paste0("Fitting limma models for ", mode,
                                 " in ", study, "_", population,
                                 " with covariates set: ", covar_name,
                                 " (", length(parts$todo), " to fit, ",
                                 length(parts$keep), " reused) ..."))

                  if (length(parts$todo) == 0) return(parts$keep)

                  block <- combined_data$rand_id

                  fits <- list(designls[parts$todo], dupcorls[parts$todo],
                               parts$todo) |>
                    furrr::future_pmap(function(design, dupcor, exposure_var) {
                      fit_limma(metabo, design,
                                block = block,
                                correlation = dupcor$consensus.correlation,
                                exposure_var = exposure_var)
                    }, .options = furrr_options(seed = TRUE), .progress = TRUE)

                  c(parts$keep, fits)[names(designls)]
                })
            })
        })
    }) |>
    purrr::set_names("limma_fit_c18", "limma_fit_hilic") |>
    list2env(.GlobalEnv)
})

rm(prior_limma)
gc()

## Reset to sequential plan
plan(sequential)

save(limma_fit_c18,   limma_keys_c18,   checkpoint_covar_signature,
     file = fit_c18_path)
save(limma_fit_hilic, limma_keys_hilic, checkpoint_covar_signature,
     file = fit_hilic_path)
message("Checkpoint written: ", fit_c18_path)
message("Checkpoint written: ", fit_hilic_path)


# Extract MWAS results -------------------------------------------------------

## Function to extract topTable results
extract_toptable <- function(fit, design, metabolome_matrix) {
  ## Scored indices sit in the second design column (the first is the
  ## intercept). After contrasts.fit the feature-wise QGcomp fit carries psi as
  ## its only coefficient, so resolve by name when it is present -- coef = 2
  ## would be out of bounds there.
  coef_use <- if ("psi" %in% colnames(fit$coefficients)) "psi" else 2

  limma::topTable(
    fit,
    coef = coef_use,
    sort.by = "p",
    number = nrow(metabolome_matrix),
    adjust.method = "BH"  # Benjamini-Hochberg FDR adjustment
  )
}

list(
  list("C18", "HILIC"),
  list(limma_fit_c18, limma_fit_hilic),
  list(design_c18_list, design_hilic_list),
  list(metabo_list_c18_final, metabo_list_hilic_final)
) |>
  purrr::pmap(function(mode, fit_data_list, design_data_list, metabo_data_list){
    list(fit_data_list, design_data_list,
         metabo_data_list, names(metabo_data_list)) |>
      purrr::pmap(function(fit_list, design_list, metabo_list, study){
        list(fit_list, design_list, metabo_list, names(metabo_list)) |>
          purrr::pmap(function(fit_ls, design_ls, metabo, population) {
            list(fit_ls, design_ls, names(design_ls)) |>
              purrr::pmap(function(fitls, designls, covar_name){
                message(paste0("Extracting MWAS results for ", mode,
                               " in ", study, "_", population,
                               " with covariates set: ",
                               covar_name, " ..."))
                list(fitls, designls) |>
                  purrr::pmap(function(fit, design) {
                    extract_toptable(fit, design, metabo)
                  })
              })
          })
      })
  }) |>
  purrr::set_names("mwas_results_list_c18", "mwas_results_list_hilic") |>
  list2env(.GlobalEnv) |>
  invisible()


# Summarize significant results ----------------------------------------------

## Function to summarize significant metabolites
summarize_significant <- function(results_list, fdr_threshold = 0.05) {
  results_list |>
    purrr::imap(function(tbl, exp_name) {
      sig_count <- sum(tbl$adj.P.Val < fdr_threshold, na.rm = TRUE)
      tibble(
        exposure = exp_name,
        n_significant = sig_count,
        n_total = nrow(tbl),
        pct_significant = round(sig_count / nrow(tbl) * 100, 2)
      )
    }) |>
    purrr::list_rbind()
}

list(
  list("C18", "HILIC"),
  list(mwas_results_list_c18, mwas_results_list_hilic)
) |>
  purrr::pmap(function(mode, mwas_data_list){
    message(paste0(mode, " MWAS Summary (FDR < 0.05):"))
    print(
      mwas_data_list |>
        purrr::map(function(mwas_list){
          mwas_list |>
            purrr::map(function(result_list){
              result_list |>
                purrr::map(~ summarize_significant(.x,
                                                   fdr_threshold = FDR_THRESHOLD))
            })
        })
    )
  }) |>
  invisible()


# =============================================================================
# SECTION 2: PLS WITH VIP SCORES
# =============================================================================

# Prepare metabolomics matrices for PLS --------------------------------------

## Transpose the metabolite matrices (samples as rows). Built one platform at a
## time, and only when that platform's PLS actually has to run -- each
## transposed slice is ~125 MB and there are six per platform.
make_metabo_matrix_list <- function(metabo_data_list, sample_link, mode) {
  metabo_data_list |>
    purrr::imap(function(metabo_list, study){
      metabo_list |>
        purrr::imap(function(metabo, population){
          message(paste0("Preparing metabolomics matrix for ", mode,
                         " in ", study, "_", population, "..."))
          metabo |>
            t() |>
            as.data.frame() |>
            tibble::rownames_to_column("file.name_new") |>
            dplyr::inner_join(sample_link, by = "file.name_new") |>
            dplyr::select(-c(rand_id, blood_date, wave, batch)) |>
            tibble::column_to_rownames("file.name_new") |>
            as.matrix()
        })
    })
}

## Confirm the PLS matrix rows line up with the exposure/covariate rows
check_pls_ordering <- function(combined_data_list, metabo_matrix_list, mode) {
  message(paste0(mode, " PLS sample ordering check:"))
  list(combined_data_list, metabo_matrix_list) |>
    purrr::pmap(function(combined_df_list, metabo_list){
      list(combined_df_list, metabo_list) |>
        purrr::pmap(function(combined_dflist, metabo_matrix) {
          combined_dflist |>
            purrr::map(function(combined_data) {
              stopifnot(all(combined_data$file.name_new ==
                              rownames(metabo_matrix)))
            })
        })
    }) |>
    invisible()
}


# Fit PLS models for each exposure -------------------------------------------

## Function to fit PLS and extract VIP (with covariate adjustment)
fit_pls_vip <- function(X, Y, covars_df = NULL, ncomp = 3) {
  # Remove NA values
  valid_idx <- complete.cases(Y, covars_df)
  X_valid <- X[valid_idx, ]
  Y_valid <- Y[valid_idx]

  # Residualize X and Y on covariates to adjust for confounders
  if (!is.null(covars_df)) {
    covars_valid <- covars_df[valid_idx, , drop = FALSE]
    covar_matrix <- model.matrix(~ ., data = covars_valid)

    # Residualize Y
    Y_valid <- residuals(lm.fit(covar_matrix, Y_valid))

    # Residualize each column of X
    X_valid <- apply(X_valid, 2, function(x) {
      residuals(lm.fit(covar_matrix, x))
    })
  }

  # Fit PLS
  pls_fit <- mixOmics::pls(X_valid, Y_valid, ncomp = ncomp)

  # Extract VIP scores
  vip_scores <- mixOmics::vip(pls_fit) |>
    as.data.frame() |>
    arrange(desc(comp1))

  return(list(pls_fit = pls_fit, vip = vip_scores))
}


## Run PLS for one platform, returning [study][population][covar_set][exposure]
##
## `prior_vip` is a VIP checkpoint from an earlier run. Exposures it already
## covers are skipped, so widening the exposure list costs only the new fits
## (see partial_checkpoint above). The returned object therefore holds the NEW
## results only; merge_pls_cells() puts them back together with the stored
## ones.
run_pls_platform <- function(combined_data_list, metabo_matrix_list, mode,
                             prior_vip = NULL, keys = NULL) {
  list(combined_data_list, metabo_matrix_list,
       names(combined_data_list)) |>
    purrr::pmap(function(combined_df_list, metabo_list, study) {
      list(combined_df_list, metabo_list, names(combined_df_list)) |>
        purrr::pmap(function(combined_dflist, metabo_matrix, population) {
          list(combined_dflist, covars_list_new, names(covars_list_new)) |>
            purrr::pmap(function(combined_data, covars, covar_name) {
              message(paste0("Run PLS in ", mode, " for ", study, "_",
                             population, " with covariates set: ",
                             covar_name, " ..."))
              covars_df <- combined_data |>
                dplyr::select(all_of(covars))
              ## PLS needs a single exposure vector as Y. comp_qgcomp_fw_all is
              ## a contrast across eight columns, not a scored index, so there
              ## is nothing to regress the feature matrix on -- it is skipped
              ## here and carries limma results only.
              ##
              ## The single pollutants are skipped for a different reason: they
              ## are limma-only by design. VIP > 2 is a ranking heuristic on the
              ## first PLS component, and it is the FDR-adjusted coefficient,
              ## not the ranking, that carries the attribution argument. PLS is
              ## also the long pole in this script, and eight more exposures
              ## would roughly triple it. combine_mwas_vip fills their VIP
              ## columns with NA, exactly as for comp_qgcomp_fw_all.
              wanted <- exposures_for(study, population) |> pls_eligible()
              todo   <- split_cell(
                prior_vip, study, population, covar_name, wanted,
                keys[[study]][[population]][[covar_name]])$todo

              if (length(todo) < length(wanted)) {
                message("  reusing ", length(wanted) - length(todo),
                        " stored PLS result(s), fitting ", length(todo))
              }

              todo |>
                purrr::set_names() |>
                purrr::map(function(exp_var) {
                  fit_pls_vip(
                    X = metabo_matrix,
                    Y = combined_data[[exp_var]],
                    covars_df = covars_df,
                    ncomp = PLS_NCOMP
                  )
                })
            })
        })
    })
}

## The model key of every PLS cell: the exposure vector, the covariate frame
## and the sample order. fit_pls_vip() is a function of exactly those three
## plus the feature matrix, whose rows check_pls_ordering() asserts are
## file.name_new in this order -- so hashing the frame's own column covers it
## without building the matrix first, which is the point (a fully resumed
## platform must never pay to transpose 125 MB slices).
build_pls_keys <- function(combined_data_list) {
  combined_data_list |>
    purrr::imap(function(combined_df_list, study){
      combined_df_list |>
        purrr::imap(function(combined_dflist, population){
          list(combined_dflist, covars_list_new, names(covars_list_new)) |>
            purrr::pmap(function(combined_data, covars, covar_name){
              covars_df <- combined_data |> dplyr::select(dplyr::all_of(covars))
              exposures_for(study, population) |>
                pls_eligible() |>
                purrr::set_names() |>
                purrr::map(function(exp_var){
                  rlang::hash(list(combined_data[[exp_var]], covars_df,
                                   combined_data$file.name_new, exp_var))
                })
            }) |>
            purrr::set_names(names(covars_list_new))
        })
    })
}

## Stored VIP results plus newly computed ones, in the order this run expects.
merge_pls_cells <- function(prior, new_obj, keys) {
  new_obj |>
    purrr::imap(function(pop_list, study){
      pop_list |>
        purrr::imap(function(covar_ls, population){
          covar_ls |>
            purrr::imap(function(cell, covar_name){
              wanted <- exposures_for(study, population) |> pls_eligible()
              keep   <- split_cell(prior, study, population, covar_name,
                                   wanted,
                                   keys[[study]][[population]][[covar_name]])$keep
              c(keep, cell)[wanted]
            })
        })
    })
}

## The [study][population][covar_set] grid with no exposures in it, so a fully
## resumed platform can go through merge_pls_cells() like any other and come
## out ordered the same way.
empty_pls_grid <- function(combined_data_list) {
  combined_data_list |>
    purrr::map(function(combined_df_list){
      combined_df_list |>
        purrr::map(function(combined_dflist){
          names(covars_list_new) |>
            purrr::set_names() |>
            purrr::map(function(covar_name) list())
        })
    })
}

## How many PLS models a platform still owes, given what the checkpoint holds.
## Zero means the transposed feature matrices -- ~125 MB per study x population
## slice -- never have to be built at all.
pls_todo_count <- function(prior, combined_data_list, keys) {
  total <- 0L
  combined_data_list |> purrr::iwalk(function(combined_df_list, study){
    combined_df_list |> purrr::iwalk(function(combined_dflist, population){
      names(covars_list_new) |> purrr::walk(function(covar_name){
        wanted <- exposures_for(study, population) |> pls_eligible()
        total <<- total + length(split_cell(
          prior, study, population, covar_name, wanted,
          keys[[study]][[population]][[covar_name]])$todo)
      })
    })
  })
  total
}

extract_vip <- function(pls_results) {
  pls_results |>
    purrr::map(function(pls_result_list){
      pls_result_list |>
        purrr::map(function(pls_result_ls){
          pls_result_ls |>
            purrr::map(function(pls_results){
              pls_results |> purrr::map(~ .x$vip)
            })
        })
    })
}

# Set maximum vector size for PLS to avoid memory issues
mem.maxVSize(vsize = 65536)

## Run PLS separately for C18 and HILIC to reduce peak memory usage, and
## checkpoint each platform independently: a crash in HILIC does not cost the
## C18 run. Only the small VIP lists are needed downstream, so those are what
## the resume path loads -- the full PLS objects are saved for reference but
## never read back.

pls_specs <- list(
  list(mode      = "C18",
       vip_name  = "vip_c18_list",
       pls_name  = "pls_results_list_c18",
       metabo    = "metabo_list_c18_final",
       link      = "sample_link_c18",
       combined  = "combined_data_list_c18"),
  list(mode      = "HILIC",
       vip_name  = "vip_hilic_list",
       pls_name  = "pls_results_list_hilic",
       metabo    = "metabo_list_hilic_final",
       link      = "sample_link_hilic",
       combined  = "combined_data_list_hilic")
)

for (spec in pls_specs) {

  vip_path <- rev_here("data", "metabolomics", "results",
                       paste0(spec$vip_name, "_revision.RData"))
  pls_path <- rev_here("data", "metabolomics", "results",
                       paste0(spec$pls_name, "_revision.RData"))

  key_name  <- paste0("pls_keys_", spec$mode)
  pls_keys  <- build_pls_keys(get(spec$combined))
  prior_vip <- partial_checkpoint(vip_path, spec$vip_name, key_name)
  n_todo    <- pls_todo_count(prior_vip, get(spec$combined), pls_keys)

  if (!is.null(prior_vip) && n_todo == 0) {
    message("RESUME: ", spec$vip_name,
            " already covers every exposure this run carries")
    assign(spec$vip_name,
           merge_pls_cells(prior_vip, empty_pls_grid(get(spec$combined)),
                           pls_keys),
           envir = globalenv())
    next
  }

  message("=== Running PLS for ", spec$mode, " (", n_todo, " model",
          if (n_todo == 1) "" else "s", " to fit) ===")

  metabo_matrix_list <- make_metabo_matrix_list(
    get(spec$metabo), get(spec$link), spec$mode)
  check_pls_ordering(get(spec$combined), metabo_matrix_list, spec$mode)

  new_pls <- run_pls_platform(get(spec$combined), metabo_matrix_list,
                              spec$mode, prior_vip = prior_vip,
                              keys = pls_keys)

  ## The full PLS objects are ~9 GB per platform and are never read back --
  ## only the VIP lists are. On a partial resume they are written to their own
  ## file rather than over the existing one, which holds the fits this run
  ## deliberately did not recompute; overwriting it with a subset would destroy
  ## them.
  assign(spec$pls_name, new_pls, envir = globalenv())
  pls_out <- if (is.null(prior_vip)) {
    pls_path
  } else {
    rev_here("data", "metabolomics", "results",
             paste0(spec$pls_name, "_revision_added_",
                    format(Sys.Date(), "%Y%m%d"), ".RData"))
  }
  ## checkpoint_covar_signature has to travel with the object: signature_ok()
  ## rejects any checkpoint that does not carry one, so a PLS run saved without
  ## it could never be resumed and this stage was recomputed on every run.
  save(list = c(spec$pls_name, "checkpoint_covar_signature"),
       file = pls_out, envir = globalenv())
  message("PLS objects written: ", pls_out)

  assign(spec$vip_name,
         merge_pls_cells(prior_vip, extract_vip(get(spec$pls_name)), pls_keys),
         envir = globalenv())
  assign(key_name, pls_keys, envir = globalenv())
  save(list = c(spec$vip_name, key_name, "checkpoint_covar_signature"),
       file = vip_path, envir = globalenv())
  message("Checkpoint written: ", vip_path)

  rm(list = c(spec$pls_name), envir = globalenv())
  rm(metabo_matrix_list, new_pls, prior_vip, pls_keys)
  gc()
}




# Identify high-VIP metabolites ----------------------------------------------

## Function to count VIP > 2 metabolites
count_high_vip <- function(vip_list, threshold = 2) {
  vip_list |>
    purrr::imap(function(vip_df, exp_name) {
      high_vip_count <- sum(vip_df$comp1 > threshold, na.rm = TRUE)
      tibble(
        exposure = exp_name,
        n_vip_gt_2 = high_vip_count,
        n_total = nrow(vip_df)
      )
    }) |>
    purrr::list_rbind()
}


message("\nC18 VIP > 2 Summary:")
print(vip_c18_list |>
        purrr::map(function(vip_data_list){
          vip_data_list |>
            purrr::map(function(vip_list){
              vip_list |>
                purrr::map(~ count_high_vip(.x, threshold = VIP_THRESHOLD))
            })
        })
)


message("\nHILIC VIP > 2 Summary:")
print(vip_hilic_list |>
        purrr::map(function(vip_data_list){
          vip_data_list |>
            purrr::map(function(vip_list){
              vip_list |>
                purrr::map(~ count_high_vip(.x, threshold = VIP_THRESHOLD))
            })
        })
)


# =============================================================================
# SECTION 3: COMBINE LIMMA AND VIP RESULTS
# =============================================================================

# Function to combine MWAS and VIP results -----------------------------------

## Joined by exposure NAME, not by position: comp_qgcomp_fw_all has limma
## results but no PLS model (see the PLS section), so the two lists are no
## longer the same length. Its VIP columns are filled with NA rather than
## dropping the exposure.
combine_mwas_vip <- function(mwas_results, vip_results, annotation_df = NULL) {
  names(mwas_results) |>
    purrr::set_names() |>
    purrr::map(function(exp_name) {
      out <- mwas_results[[exp_name]] |>
        tibble::rownames_to_column("met")

      vip <- vip_results[[exp_name]]
      if (is.null(vip)) {
        out <- out |>
          dplyr::mutate(VIP_comp1 = NA_real_, VIP_comp2 = NA_real_,
                        VIP_comp3 = NA_real_)
      } else {
        out <- out |>
          dplyr::left_join(
            vip |>
              tibble::rownames_to_column("met") |>
              dplyr::select(met, VIP_comp1 = comp1,
                            VIP_comp2 = comp2, VIP_comp3 = comp3),
            by = "met"
          )
      }
      out |> dplyr::arrange(adj.P.Val)
    })
}

list(
  list("C18", "HILIC"),
  list(mwas_results_list_c18, mwas_results_list_hilic),
  list(vip_c18_list, vip_hilic_list)
) |>
  purrr::pmap(function(mode, mwas_data_list, vip_data_list){
    list(mwas_data_list, vip_data_list, names(mwas_data_list)) |>
      purrr::pmap(function(mwas_list, vip_list, study){
        list(mwas_list, vip_list, names(mwas_list)) |>
          purrr::pmap(function(mwas_results_ls, vip_results_ls, population){
            list(mwas_results_ls, vip_results_ls, names(mwas_results_ls)) |>
              purrr::pmap(function(mwas_results, vip_results, covar_name){
                message(paste0("Combining MWAS and VIP results for ", mode,
                               " in ", study, "_", population,
                               " with covariates set: ",
                               covar_name, " ..."))

                combine_mwas_vip(mwas_results, vip_results)
              })
          })
      })
  }) |>
  purrr::set_names("combined_results_list_c18",
                   "combined_results_list_hilic") |>
  list2env(.GlobalEnv)


# Filter significant metabolites ---------------------------------------------

## Function to filter significant metabolites
##
## FDR ONLY. VIP > 2 was removed from the selection criterion on 2026-09-01
## after R14-pls_validation.R showed it does not survive validation
## (Reviewer 2 comment 7):
##
##   - cross-validated Q2 is NEGATIVE at every component count on both
##     platforms (best -0.003 HILIC, -0.092 C18, and worsening with more
##     components), so the PLS has no out-of-sample predictive value;
##   - permuting the exposure across participants 200 times puts a mean of
##     433.5 (C18) and 448.8 (HILIC) features above VIP > 2 under the null,
##     against 456 and 425 observed -- an empirical FDR of 0.95 and 1.06.
##
## VIP > 2 therefore selects noise at the rate it selects signal, and a `_sig`
## table built on it would be ~450 noise features per model wide. The PLS is
## retained as exploratory dimension reduction and its VIP columns are still
## carried for description and for the figure classes, but nothing is SELECTED
## on them. `vip_thresh` is kept in the signature so callers do not break.
filter_significant <- function(combined_results, fdr_thresh = 0.05,
                               vip_thresh = 2) {
  combined_results |>
    purrr::map(function(df) {
      df |>
        dplyr::filter(adj.P.Val < fdr_thresh) |>
        dplyr::arrange(adj.P.Val)
    })
}

list(
  list("C18", "HILIC"),
  list(combined_results_list_c18, combined_results_list_hilic)
) |>
  purrr::pmap(function(mode, combined_results_data_list){
    combined_results_data_list |>
      purrr::imap(function(combined_results_list, study){
        combined_results_list |>
          purrr::imap(function(dflist, population){
            dflist |>
              purrr::imap(function(df, covar_name){
                message(paste0("Filtering significant metabolites for ", mode,
                               " in ", study, "_", population,
                               " with covariates set: ",
                               covar_name, " ..."))
                filter_significant(df, fdr_thresh = FDR_THRESHOLD,
                                   vip_thresh = VIP_THRESHOLD)

              })
          })
      })
  }) |>
  purrr::set_names("significant_results_list_c18",
                   "significant_results_list_hilic") |>
  list2env(.GlobalEnv)


# =============================================================================
# SECTION 4: SAVE RESULTS
# =============================================================================

# Create output directories --------------------------------------------------

combined_data_list_c18 |>
  purrr::imap(function(data, study){
    names(data) |>
      purrr::walk(function(population){
        names(covar_list) |>
          purrr::walk(function(covar_name){
            rev_dir("tables", "mwas_results", study, population, covar_name)
          })
      })
  })

rev_dir("data", "metabolomics", "results")

# Save MWAS results to Excel -------------------------------------------------

## Save combined MWAS + VIP results

list(
  list(combined_results_list_c18, combined_results_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(data_list, mode){
    data_list |>
      purrr::imap(function(datalist, study){
        datalist |>
          purrr::imap(function(dflist, population){
            dflist |>
              purrr::imap(function(df_ls, covar_name){
                message(paste0("Saving combined MWAS + VIP results for ", mode,
                               " in ", study, "_", population,
                               " with covariates set: ",
                               covar_name, " ..."))
                df_ls |>
                  purrr::imap(function(df, exp_name) {
                    writexl::write_xlsx(
                      df,
                      path = rev_here(
                        "tables", "mwas_results", study, population, covar_name,
                        glue::glue(
                          "mwas_{mode}_{exp_name}_{study}_{population}_{covar_name}.xlsx"
                        ))
                    )
                  })
              })
          })
      })
  })


# Save significant metabolites -----------------------------------------------

list(
  list(significant_results_list_c18, significant_results_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(data_list, mode){
    data_list |>
      purrr::imap(function(datalist, study){
        datalist |>
          purrr::imap(function(dflist, population){
            dflist |>
              purrr::imap(function(df_ls, covar_name){
                message(paste0("Saving significant MWAS + VIP results for ",
                               mode, " in ", study, "_", population,
                               " with covariates set: ",
                               covar_name, " ..."))
                df_ls |>
                  purrr::imap(function(df, exp_name) {
                    if (nrow(df) > 0) {
                      writexl::write_xlsx(
                        df,
                        path = rev_here(
                          "tables", "mwas_results", study, population, covar_name,
                          glue::glue(
                            "mwas_{mode}_{exp_name}_{study}_{population}_{covar_name}_sig.xlsx"
                          ))
                      )
                    }
                  })
              })
          })
      })
  })


# Save R objects for downstream analysis -------------------------------------

save(mwas_results_list_c18, mwas_results_list_hilic,
     vip_c18_list, vip_hilic_list,
     combined_results_list_c18, combined_results_list_hilic,
     significant_results_list_c18, significant_results_list_hilic,
     dupcor_c18_list, dupcor_hilic_list,
     limma_fit_c18, limma_fit_hilic,
     exposure_vars_list,
     file = rev_here("data", "metabolomics", "results",
                     "mwas_results_all_revision.RData"))

message("MWAS analysis completed! Results saved to ",
        rev_here("tables", "mwas_results"))


# =============================================================================
# SECTION 5: SUMMARY TABLE
# =============================================================================

# Create summary table for all exposures -------------------------------------

## NOTE: the *_vip_gt2 columns are DESCRIPTIVE only. They count features above
## VIP > 2 but that threshold is not a selection criterion (see
## filter_significant), because R14 shows it has an empirical FDR near 1.
create_summary_table <- function(mwas_c18, mwas_hilic,
                                 vip_c18, vip_hilic, exposure_vars) {
  ## NA rather than 0 where no PLS model exists (comp_qgcomp_fw_all): a zero
  ## would read as "PLS ran and found nothing".
  count_vip <- function(vip_list, exp) {
    if (is.null(vip_list[[exp]])) NA_integer_
    else sum(vip_list[[exp]]$comp1 > 2, na.rm = TRUE)
  }

  exposure_vars |>
    purrr::map(function(exp) {
      tibble(
        exposure = exp,
        label = unname(rev_label(exp)),
        c18_total = nrow(mwas_c18[[exp]]),
        c18_sig_fdr05 = sum(mwas_c18[[exp]]$adj.P.Val < 0.05, na.rm = TRUE),
        c18_sig_fdr10 = sum(mwas_c18[[exp]]$adj.P.Val < 0.10, na.rm = TRUE),
        c18_vip_gt2 = count_vip(vip_c18, exp),
        hilic_total = nrow(mwas_hilic[[exp]]),
        hilic_sig_fdr05 = sum(mwas_hilic[[exp]]$adj.P.Val < 0.05, na.rm = TRUE),
        hilic_sig_fdr10 = sum(mwas_hilic[[exp]]$adj.P.Val < 0.10, na.rm = TRUE),
        hilic_vip_gt2 = count_vip(vip_hilic, exp)
      )
    }) |>
    purrr::list_rbind()
}

summary_table_list <- list(
  mwas_results_list_c18, mwas_results_list_hilic,
  vip_c18_list, vip_hilic_list,
  names(mwas_results_list_c18)
) |>
  purrr::pmap(function(mwas_results_c18_dflist, mwas_results_hilic_dflist,
                       vip_c18_dflist, vip_hilic_dflist, study){
    list(mwas_results_c18_dflist, mwas_results_hilic_dflist,
         vip_c18_dflist, vip_hilic_dflist,
         names(mwas_results_c18_dflist)) |>
      purrr::pmap(function(mwas_results_c18_ls, mwas_results_hilic_ls,
                           vip_c18_ls, vip_hilic_ls, population){
        list(mwas_results_c18_ls, mwas_results_hilic_ls,
             vip_c18_ls, vip_hilic_ls, names(mwas_results_c18_ls)) |>
          purrr::pmap(function(mwas_c18, mwas_hilic, vip_c18, vip_hilic,
                               covar_name){
            message("Creating summary table for ", study, "_", population,
                    " [", covar_name, "] ...")
            create_summary_table(mwas_c18, mwas_hilic,
                                 vip_c18, vip_hilic,
                                 exposures_for(study, population)) |>
              dplyr::mutate(study = study, population = population,
                            covar_set = covar_name, .before = 1)
          })
      })
  })

mwas_summary_table <- summary_table_list |>
  purrr::list_flatten() |>
  purrr::list_flatten() |>
  purrr::list_rbind()

rev_save_table(mwas_summary_table, "mwas_summary_table", "mwas")
print(mwas_summary_table |>
        dplyr::filter(population == "all", covar_set == "covar"), n = 30)

#--------------------------------End of the code--------------------------------
