## ---------------------------
##
## Script name: R9-whicap_concordance.R
## Purpose of script: Cross-study concordance between the SALSA and WHICAP
##                    PM2.5 metabolome-wide association results, on a single
##                    common estimand
##
## Author: Yufan Gong
##
## Date Created: 2026-08-26
##
## Date Modified: 2026-09-13
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: STAGE 4. Addresses Reviewer 1 major comment 10 and Reviewer 2
##        comment 8.
##
##        The submitted analysis pooled SALSA multipollutant-index coefficients
##        with WHICAP PM2.5 coefficients. Those do not estimate the same
##        contrast, so the fixed-effect estimand was undefined. This script
##        compares PM2.5 with PM2.5 and nothing else.
##
##        ------------------------------------------------------------------
##        2026-09-13 REVISION. Four defects in the 2026-08-26 version are
##        fixed here; each changes the numbers, so nothing from the previous
##        run should be quoted.
##
##        (1) HALF THE WHICAP FEATURES WERE SILENTLY DROPPED.
##            In 20260427_PM25_MWAS_results_annotations.txt the `mz` and
##            `time` columns are populated only for features that received an
##            annotation: 4,806 of the 10,134 features carry mz = time = NA.
##            The old join was on (mz, time), so those 4,806 never reached the
##            crosswalk -- including several of WHICAP's strongest signals
##            (e.g. hilpos_320.0038_182.8, p = 4.2e-7). m/z and RT are present
##            for every row inside the `feature` string, so they are parsed
##            from there instead. Aligned pairs go from 1,753 to ~3,300, which
##            is the ~3,343 the Reviewer refers to.
##
##        (2) THE WHICAP COEFFICIENTS WERE RE-STANDARDIZED WHEN THEY WERE
##            ALREADY STANDARDIZED -- the single most consequential error.
##            Kalia et al. log10-transformed, quantile-normalized and
##            AUTO-SCALED every feature before the MWAS, so their beta is
##            already per 1 SD of feature abundance; and the estimates behave
##            numerically as partial correlations, so the exposure is on a
##            1-SD scale too. Section 2 below tests this rather than asserting
##            it: beta / [t / sqrt(t^2 + df)] has median 0.995 (IQR
##            0.96-1.02) across all 10,134 features, and the standard errors
##            are near-constant at 1/sqrt(df). The old script divided those
##            betas by a delta-method SD reconstructed from the untransformed
##            feature summaries, which double-standardized them and injected
##            the CV of every feature as noise. That conversion is deleted.
##            The feature summaries are kept only as the audit in section 2.
##
##        (3) THE SALSA PER-SD CONVERSION WAS INVERTED. SALSA coefficients are
##            per IQR-scaled unit of PM2.5, and sd(exp_pm2.5_iqr) = 0.927, so
##            the per-SD coefficient is logFC * 0.927. The old script divided,
##            inflating every SALSA standardized estimate by 1/0.927^2 = 16%.
##
##        (4) ONE SD IS NOT THE SAME CONTRAST IN THE TWO COHORTS. SALSA and
##            WHICAP have almost the same mean PM2.5 (12.4 vs 12.9 ug/m3) but
##            SALSA's SD is 0.66 ug/m3 against WHICAP's 2.41 -- a 3.6-fold
##            difference in spread. Standardizing per SD makes the two
##            coefficients numerically comparable but leaves them describing
##            different absolute increments. Everything is therefore reported
##            on BOTH scales: per 1 SD of PM2.5 (WHICAP's native scale) and
##            per 1 ug/m3 of PM2.5 (a genuinely common increment). Section 4
##            tabulates the two exposure distributions side by side.
##        ------------------------------------------------------------------
##
##        WHAT THE COMPARISON CAN AND CANNOT BE. Both sides are now PM2.5 and
##        identically standardized, but the two studies still differ in
##        exposure window (SALSA 5-year average before the draw, WHICAP the
##        calendar year before the visit), feature processing, covariate set
##        (WHICAP adjusts for AD status; SALSA stratifies on dementia/CIND
##        rather than adjusting for it in total/all), sample size (952
##        participants / 1,546 specimens vs 107), population, and geography.
##        Both cohorts profiled EDTA plasma at the same laboratory, so the
##        biospecimen matrix is NOT a difference. Section 11 removes the
##        window, transform and covariate differences one at a time.
##
##        TWO ANALYSES ARE REPORTED, one for each route the Reviewer offers.
##        Section 12 is the FIXED-EFFECT META-ANALYSIS on the identically
##        standardized scale (route 1) and is the reported result. Section 8 is
##        the CROSS-STUDY CONCORDANCE on the coefficients as fitted (route 2),
##        the honest counterweight: direction concordance over ALL aligned
##        features, with significance-selected subsets printed only to show how
##        much selection inflates apparent agreement. Neither is described as
##        replication.
##
##        INPUTS (all under data/metabolomics/alignment/):
##          20260427_PM25_MWAS_results_annotations.txt
##              WHICAP PM2.5 MWAS: feature, estimate, std.error, statistic,
##              p_value, fdr_q_value, column_esi ("C18 -" / "HILIC +"), mz,
##              time, name, chemical_formula, confidence, adduct.
##              10,161 rows = 10,134 features + 27 annotation duplicates.
##          {c18neg,hilpos}_feature_level_summary.txt
##              WHICAP per-feature mz / time / mean / sd on UNTRANSFORMED
##              intensities, detection >= 70%, zeros replaced by half the
##              per-feature minimum. Used here only to audit (2).
##          {salsa,whicap}_{c18,hil}_link.RData
##              apLCMS alignment crosswalks from scripts/9-salsa_whicap_alignment.R
##
##        Outputs -> <REV_ROOT>/{tables,figures,data}/whicap/
##
## ---------------------------

# Setup -----------------------------------------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R9-whicap_concordance.R")

# Helpers ----------------------------------------------------------------------
#
# list_mwas_files() defaults to the REVISION output tree, not tables/. The SALSA
# PM2.5 estimates consumed here must be the ones produced under the revision
# covariate set (batch dropped), not the submitted ones.

fs_rel_path <- function(files, base) {
  base_norm <- normalizePath(base, winslash = "/", mustWork = FALSE)
  files |>
    normalizePath(winslash = "/", mustWork = FALSE) |>
    stringr::str_remove(stringr::fixed(paste0(base_norm, "/")))
}

parse_feature_names <- function(met) {
  tibble::tibble(met = as.character(met)) |>
    dplyr::mutate(
      mz = stringr::str_extract(met, "(?<=^mz_rt_)[0-9.]+") |> as.numeric(),
      rt = stringr::str_extract(met, "[0-9.]+$") |> as.numeric()
    )
}

list_mwas_files <- function(base = rev_here("tables", "mwas_results")) {
  files <- list.files(base, pattern = "^mwas_.*_covar(_sen)?\\.xlsx$",
                      full.names = TRUE, recursive = TRUE)
  rel <- fs_rel_path(files, base)
  parts <- stringr::str_split(rel, .Platform$file.sep, simplify = TRUE)

  tibble::tibble(
    path       = files,
    file       = basename(files),
    study      = parts[, 1],
    population = parts[, 2],
    covar_set  = parts[, 3]
  ) |>
    dplyr::mutate(
      platform = stringr::str_extract(file, "(?<=^mwas_)(c18|hilic)"),
      exposure = file |>
        stringr::str_remove("^mwas_(c18|hilic)_") |>
        stringr::str_remove(
          stringr::fixed(paste0("_", study, "_", population, "_",
                                covar_set, ".xlsx"))
        )
    ) |>
    dplyr::select(path, file, platform, exposure, study, population, covar_set)
}

read_mwas_xlsx <- function(path) {
  res <- readxl::read_xlsx(path)
  res |>
    dplyr::bind_cols(parse_feature_names(res$met) |> dplyr::select(mz, rt))
}


PPM_TOL      <- 10
RT_TOL       <- 30
N_PERM_ALIGN <- rev_n(200, 5)

## Tolerance grid for the alignment QC the Reviewer asks for. The middle row is
## the tolerance actually used.
FMR_GRID <- tibble::tribble(
  ~ppm, ~rt,
     5,  15,
    10,  30,
    20,  60
)

## WHICAP model residual degrees of freedom: n = 107 less the intercept,
## PM2.5, age, sex, two race/ethnicity indicators, AD status, year of blood
## draw and education = 107 - 9 = 98. The audit in section 2 is insensitive to
## the exact value -- the median beta / partial-correlation ratio moves from
## 0.986 at df = 98 to 0.995 at df = 100 -- so the nominal df is used.
WHICAP_N  <- 107
WHICAP_DF <- 98

align_dir <- here::here("data", "metabolomics", "alignment")

whicap_mwas_file <- file.path(align_dir,
                              "20260427_PM25_MWAS_results_annotations.txt")

required <- c(
  whicap_mwas_file,
  file.path(align_dir, "c18neg_feature_level_summary.txt"),
  file.path(align_dir, "hilpos_feature_level_summary.txt"),
  file.path(align_dir, "salsa_c18_link.RData"),
  file.path(align_dir, "whicap_c18_link.RData"),
  file.path(align_dir, "salsa_hil_link.RData"),
  file.path(align_dir, "whicap_hil_link.RData")
)

missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  message("\n", strrep("-", 66))
  message("R9 SKIPPED - required inputs are not present:")
  purrr::walk(missing, ~ message("  * ", .x))
  message("Alignment crosswalks come from scripts/9-salsa_whicap_alignment.R.")
  message(strrep("-", 66))
  skip_r9 <- TRUE
} else {
  skip_r9 <- FALSE
}

if (!skip_r9) {

# 1. WHICAP PM2.5 MWAS ---------------------------------------------------------
#
## m/z and RT are taken from the `feature` string, which is populated for every
## row; the mz / time COLUMNS are NA for the 4,806 unannotated features and
## must not be used as the join key (see defect 1 in the header).

whicap_raw <- utils::read.delim(whicap_mwas_file, sep = "\t", header = TRUE,
                                check.names = FALSE,
                                stringsAsFactors = FALSE) |>
  tibble::as_tibble()

message("WHICAP MWAS rows as supplied: ", nrow(whicap_raw))
message("  rows whose mz/time COLUMNS are NA: ",
        sum(is.na(whicap_raw$mz) | whicap_raw$mz %in% c("NA", "")),
        "  <- dropped by the previous version of this script")

whicap_mwas <- whicap_raw |>
  dplyr::mutate(
    platform = dplyr::case_when(
      column_esi == "C18 -"   ~ "c18",
      column_esi == "HILIC +" ~ "hil",
      TRUE ~ NA_character_
    ),
    ## parsed from `feature`, never from the mz / time columns
    whicap_mz   = as.numeric(stringr::str_extract(feature,
                                                  "(?<=_)[0-9.]+(?=_[0-9.]+$)")),
    whicap_time = as.numeric(stringr::str_extract(feature, "[0-9.]+$"))
  ) |>
  dplyr::filter(!is.na(platform), is.finite(whicap_mz), is.finite(whicap_time))

message("  features with a parseable m/z and RT: ",
        dplyr::n_distinct(whicap_mwas$feature))

## Collapse the 27 annotation duplicates to one row per feature, keeping the
## most confident annotation.
whicap_mwas <- whicap_mwas |>
  dplyr::mutate(
    conf_rank = suppressWarnings(
      as.integer(stringr::str_extract(confidence, "[0-9]+"))),
    conf_rank = dplyr::coalesce(conf_rank, 99L)
  ) |>
  dplyr::arrange(platform, feature, conf_rank) |>
  dplyr::distinct(platform, feature, .keep_all = TRUE)

message("WHICAP PM2.5 MWAS features carried forward: ", nrow(whicap_mwas))
print(dplyr::count(whicap_mwas, platform))


# 2. Audit: are the WHICAP coefficients already standardized? -------------------
#
## Kalia et al. state that features were log10-transformed, quantile-normalized
## and auto-scaled. If so, and if the exposure is on a 1-SD scale, then for
## every feature beta must equal the partial correlation t / sqrt(t^2 + df),
## and the standard error must be close to 1/sqrt(df) for all but the largest
## effects. Both are tested here, because the whole comparison rests on it.
## A third check confirms that the estimates carry no trace of the raw
## intensity scale: if they had NOT been auto-scaled, |beta| would track the
## per-feature SD of the untransformed data.

whicap_summary <- list(
  c18 = "c18neg_feature_level_summary.txt",
  hil = "hilpos_feature_level_summary.txt"
) |>
  purrr::imap(function(f, platform) {
    utils::read.delim(file.path(align_dir, f), sep = "\t", header = TRUE,
                      check.names = FALSE) |>
      tibble::as_tibble() |>
      dplyr::transmute(platform, whicap_mz = mz, whicap_time = time,
                       raw_mean = mean, raw_sd = sd, cv = sd / mean)
  }) |>
  purrr::list_rbind()

scale_audit_dat <- whicap_mwas |>
  dplyr::mutate(partial_r = statistic / sqrt(statistic^2 + WHICAP_DF),
                beta_over_r = estimate / partial_r) |>
  dplyr::left_join(whicap_summary, by = c("platform", "whicap_mz",
                                          "whicap_time"))

whicap_scale_audit <- tibble::tibble(
  check = c(
    "beta / partial correlation, median",
    "beta / partial correlation, 25th pct",
    "beta / partial correlation, 75th pct",
    "std.error, median (1/sqrt(df) = expected)",
    "std.error, SD across features",
    "Spearman r( |beta|, untransformed feature SD )",
    "Spearman r( |beta|, untransformed feature mean )"
  ),
  value = c(
    round(stats::median(scale_audit_dat$beta_over_r, na.rm = TRUE), 4),
    round(stats::quantile(scale_audit_dat$beta_over_r, 0.25, na.rm = TRUE), 4),
    round(stats::quantile(scale_audit_dat$beta_over_r, 0.75, na.rm = TRUE), 4),
    round(stats::median(scale_audit_dat$std.error, na.rm = TRUE), 4),
    round(stats::sd(scale_audit_dat$std.error, na.rm = TRUE), 4),
    round(stats::cor(abs(scale_audit_dat$estimate), scale_audit_dat$raw_sd,
                     use = "complete.obs", method = "spearman"), 3),
    round(stats::cor(abs(scale_audit_dat$estimate), scale_audit_dat$raw_mean,
                     use = "complete.obs", method = "spearman"), 3)
  ),
  expected_if_autoscaled = c("1", "~1", "~1",
                             as.character(round(1 / sqrt(WHICAP_DF), 4)),
                             "~0", "~0", "~0")
)

rev_save_table(whicap_scale_audit, "whicap_scale_audit", "whicap")
print(whicap_scale_audit, width = Inf)

message("\nThe WHICAP coefficients are already expressed per 1 SD of PM2.5 and ",
        "per 1 SD of feature abundance. They are used AS SUPPLIED; the ",
        "delta-method rescaling in the previous version of this script has ",
        "been removed.")


# 3. The two SALSA exposures compared with WHICAP --------------------------------
#
## Both are the PRIMARY 5-year models from the revision tree, taken as fitted:
##
##   exp_pm2.5_iqr   PM2.5 alone -- the same pollutant WHICAP estimated, and the
##                   only strictly like-for-like contrast available
##   comp_pca_all    the unsupervised PC1 multipollutant index, which is the
##                   paper's primary exposure
##
## Coefficients are used AS FITTED, with no rescaling on either side. See the
## note at section 8 on why a concordance analysis does not need a common scale.

read_salsa_exposure <- function(exposure_var, label) {
  list_mwas_files() |>
    dplyr::filter(exposure == exposure_var, study == "total",
                  population == "all", covar_set == "covar") |>
    dplyr::mutate(data = purrr::map(path, read_mwas_xlsx)) |>
    dplyr::select(platform, data) |>
    tidyr::unnest(data) |>
    dplyr::mutate(platform = dplyr::if_else(platform == "hilic", "hil", platform),
                  exposure_lab = label)
}

salsa_pm25 <- read_salsa_exposure("exp_pm2.5_iqr", "PM2.5 (5-year)")
salsa_pc1  <- read_salsa_exposure("comp_pca_all",  "PC1 index")

salsa_both <- dplyr::bind_rows(salsa_pm25, salsa_pc1)

message("\nSALSA MWAS features loaded:")
print(dplyr::count(salsa_both, exposure_lab, platform))


# 4. The two exposure distributions --------------------------------------------
#
## Reviewer 1's second bullet asks why the two coefficients are commensurate.
## Both are now PM2.5, but one SD of PM2.5 is a very different increment in the
## two cohorts, so the answer has to be quantitative rather than asserted.

load(here::here("data", "processed", "combined_data_list_new.RData"))
dat <- combined_data_list_new[["total"]][["all"]][["covar"]]

salsa_pm25_ugm3 <- dat$exp_pm2.5
salsa_iqr_ugm3  <- stats::IQR(salsa_pm25_ugm3, na.rm = TRUE)
salsa_sd_ugm3   <- stats::sd(salsa_pm25_ugm3, na.rm = TRUE)

## SALSA coefficients are per one IQR-scaled unit, i.e. per `salsa_iqr_ugm3`
## ug/m3. sd of the scaled exposure is salsa_sd_ugm3 / salsa_iqr_ugm3.
salsa_sd_in_iqr_units <- stats::sd(dat$exp_pm2.5_iqr, na.rm = TRUE)

## The PC1 index needs NO exposure multiplier, and this is easy to get wrong.
## R4 divides every `comp_*` column by its own SD IN MEMORY before fitting
## (`combined_data_list_new <- ... d[[nm]] <- d[[nm]] / sds[[nm]]`), so the
## stored MWAS coefficient is ALREADY per 1 SD of the index. The frame saved to
## disk keeps the UNSCALED column, SD = 1.92, so taking the SD from it and
## multiplying would standardize a second time -- the same error corrected on
## the WHICAP side in section 2, and it inflates every PC1 standardized estimate
## and SE by 1.92x, which distorts the inverse-variance weights.
##
## Verified by refitting: using the stored (unscaled) column reproduces the
## stored MWAS at a ratio of exactly 1.9216, while dividing it by its SD first
## reproduces it to 1.4e-15.
##
## The single pollutants are NOT pre-scaled this way -- R4 rescales them by IQR,
## and `exp_pm2.5_iqr` as stored is what was fitted (confirmed by the step-1
## validation in section 11) -- so PM2.5 does need its multiplier.
pc1_sd_units <- 1

WHICAP_PM25_MEAN <- 12.9   # Kalia et al. 2023, Table 1
WHICAP_PM25_SD   <- 2.41

exposure_comparison <- tibble::tribble(
  ~quantity,                       ~SALSA,                      ~WHICAP,
  "Cohort",                        "SALSA, Sacramento CA",      "WHICAP, New York NY",
  "Participants",                  as.character(dplyr::n_distinct(dat$rand_id)),
                                                                as.character(WHICAP_N),
  "Specimens in the model",        as.character(nrow(dat)),     as.character(WHICAP_N),
  "Matrix",                        "Plasma (EDTA)",             "Plasma",
  "Exposure model",                "CA statewide / CALINE4-linked annual surfaces",
                                                                "Regionalized universal kriging",
  "Averaging window",              "5 years before blood draw", "Calendar year before visit",
  "PM2.5 mean (ug/m3)",            as.character(round(mean(salsa_pm25_ugm3, na.rm = TRUE), 2)),
                                                                as.character(WHICAP_PM25_MEAN),
  "PM2.5 SD (ug/m3)",              as.character(round(salsa_sd_ugm3, 2)),
                                                                as.character(WHICAP_PM25_SD),
  "PM2.5 IQR (ug/m3)",             as.character(round(salsa_iqr_ugm3, 2)), "not reported",
  "Feature transform",             "log2, median-scaled, ComBat",
                                                                "log10, quantile-normalized, auto-scaled",
  "Model",                         "limma, duplicateCorrelation on repeated draws",
                                                                "OLS, one specimen per participant",
  "Covariates",                    "age, sex, education, smoking, RUCA, nSES, wave",
                                                                "age, sex, race/ethnicity, education, year of draw, AD"
)

rev_save_table(exposure_comparison, "exposure_comparison", "whicap")
print(exposure_comparison, width = Inf)

message("\nOne SD of PM2.5 is ", round(salsa_sd_ugm3, 2), " ug/m3 in SALSA and ",
        WHICAP_PM25_SD, " ug/m3 in WHICAP -- a ",
        round(WHICAP_PM25_SD / salsa_sd_ugm3, 1),
        "-fold difference. Results are therefore reported both per 1 SD and ",
        "per 1 ug/m3.")


# 5. Feature SDs, kept for the scale-sensitivity check ---------------------------
#
## The primary analysis uses the coefficients as fitted. These standardized
## versions exist only so section 8 can report what standardizing WOULD do, and
## so the pooling in section 12 -- which genuinely requires a common scale --
## has one.

load(here::here("data", "metabolomics", "processed", "med_c18_raw_combat.Rdata"))
load(here::here("data", "metabolomics", "processed", "med_hil_raw_combat.Rdata"))

salsa_feature_sd <- list(c18 = med_c18_raw_combat, hil = med_hil_raw_combat) |>
  purrr::imap(function(m, platform) {
    tibble::tibble(
      met = rownames(m),
      salsa_sd_log2 = apply(m, 1, stats::sd, na.rm = TRUE),
      platform = platform
    )
  }) |>
  purrr::list_rbind()

## limma's topTable does not emit a standard error, but t = logFC / SE, so it
## is recoverable exactly.
salsa_std <- salsa_both |>
  dplyr::left_join(salsa_feature_sd, by = c("met", "platform")) |>
  dplyr::mutate(
    salsa_se_logfc = abs(logFC / t),
    ## Per 1 SD of exposure and 1 SD of log2 abundance. The exposure constant
    ## is the SD of that exposure in the units its own coefficient is stated in:
    ## IQR-scaled units for PM2.5, index units for PC1. MULTIPLY by it -- logFC
    ## is per one unit of the scaled exposure.
    ## 1 for PC1 (already per SD; see the note in section 4), sd of the
    ## IQR-scaled variable for PM2.5.
    sd_exposure_units = dplyr::if_else(exposure_lab == "PC1 index",
                                       pc1_sd_units, salsa_sd_in_iqr_units),
    salsa_beta_std = logFC * sd_exposure_units / salsa_sd_log2,
    salsa_se_std   = salsa_se_logfc * sd_exposure_units / salsa_sd_log2,
    ## per 1 ug/m3 of PM2.5, per 1 SD of log2 abundance (PM2.5 only)
    salsa_beta_ugm3 = dplyr::if_else(
      exposure_lab == "PC1 index", NA_real_,
      (logFC / salsa_iqr_ugm3) / salsa_sd_log2),
    salsa_se_ugm3 = dplyr::if_else(
      exposure_lab == "PC1 index", NA_real_,
      (salsa_se_logfc / salsa_iqr_ugm3) / salsa_sd_log2)
  )

whicap_std <- whicap_mwas |>
  dplyr::mutate(
    ## used as supplied: already per 1 SD exposure, per 1 SD feature
    whicap_beta_std  = estimate,
    whicap_se_std    = std.error,
    whicap_beta_ugm3 = estimate   / WHICAP_PM25_SD,
    whicap_se_ugm3   = std.error  / WHICAP_PM25_SD
  )


# 6. Map aligned features across cohorts ---------------------------------------

## Use the apLCMS crosswalk produced by scripts/9-salsa_whicap_alignment.R
## rather than re-matching: it is the alignment the manuscript describes.
load_link <- function(f) {
  e <- new.env(); load(file.path(align_dir, f), envir = e)
  get(ls(e)[1], envir = e) |> tibble::as_tibble()
}

align_xwalk <- list(
  c18 = list(salsa = "salsa_c18_link.RData", whicap = "whicap_c18_link.RData"),
  hil = list(salsa = "salsa_hil_link.RData", whicap = "whicap_hil_link.RData")
) |>
  purrr::imap(function(files, platform) {
    salsa_side <- load_link(files$salsa) |>
      dplyr::transmute(alg_match,
                       met = paste0("mz_rt_", mz, "_", time)) |>
      dplyr::distinct(alg_match, .keep_all = TRUE)

    whicap_side <- load_link(files$whicap) |>
      dplyr::transmute(alg_match, whicap_mz = mz, whicap_time = time) |>
      dplyr::distinct(alg_match, .keep_all = TRUE)

    dplyr::inner_join(salsa_side, whicap_side, by = "alg_match") |>
      dplyr::mutate(platform = platform)
  }) |>
  purrr::list_rbind() |>
  dplyr::distinct(platform, met, whicap_mz, whicap_time, .keep_all = TRUE)

message("\nAligned feature pairs in the crosswalk: ", nrow(align_xwalk))
print(dplyr::count(align_xwalk, platform))

## Provenance check: the alignment IS built from the feature-level summaries ----
##
## scripts/9-salsa_whicap_alignment.R feeds apLCMS the two WHICAP feature-level
## summary files, and this asserts that the crosswalk it produced still consists
## only of peaks from those files. The m/z and RT parsed out of the MWAS
## `feature` strings in section 1 are used ONLY to attach each WHICAP RESULT to
## a crosswalk row; they are not an alternative basis for alignment, and this
## check is here so that distinction cannot quietly stop being true.
summary_keys <- whicap_summary |>
  dplyr::mutate(key = paste(sprintf("%.4f", whicap_mz),
                            sprintf("%.1f", whicap_time))) |>
  dplyr::pull(key)

xwalk_from_summary <- align_xwalk |>
  dplyr::mutate(key = paste(sprintf("%.4f", whicap_mz),
                            sprintf("%.1f", whicap_time))) |>
  dplyr::summarise(n = dplyr::n(), in_summary = sum(key %in% summary_keys))

message("Crosswalk WHICAP peaks drawn from the feature-level summaries: ",
        xwalk_from_summary$in_summary, " of ", xwalk_from_summary$n)
if (xwalk_from_summary$in_summary < xwalk_from_summary$n) {
  warning("Some crosswalk peaks are not in {c18neg,hilpos}_feature_level_",
          "summary.txt. The alignment and the summaries have diverged; ",
          "re-run scripts/9-salsa_whicap_alignment.R.", call. = FALSE)
}


## Reconciliation with the 3,343 reported in the submitted manuscript ----------
##
## The submitted Results say "3343 metabolic features were aligned". That was a
## ROW count of the meta-analysis tables, and those tables repeat 18 feature
## pairs: 4 on C18 and 14 on HILIC appear twice, each time against the same
## WHICAP peak. Counting unique SALSA-to-WHICAP pairs gives 3,325. No feature
## was lost between the two numbers -- de-duplicating simply stops 18 pairs from
## contributing twice to a direction-concordance or correlation statistic.
legacy_mixture <- here::here("data", "metabolomics", "alignment",
                             "meta_mixture_salsa_whicap.RData")
if (file.exists(legacy_mixture)) {
  le <- new.env(); load(legacy_mixture, envir = le)
  grab_legacy <- function(obj, platform) {
    d <- try(obj$total$all$covar[[1]], silent = TRUE)
    if (inherits(d, "try-error") || is.null(d)) return(NULL)
    tibble::tibble(
      platform     = platform,
      submitted_rows  = nrow(d),
      unique_pairs    = nrow(dplyr::distinct(d, salsa_met, whicap_mz,
                                             whicap_time)),
      current_pairs   = sum(align_xwalk$platform == platform))
  }
  alignment_reconciliation <- dplyr::bind_rows(
    grab_legacy(le$meta_results_c18,   "c18"),
    grab_legacy(le$meta_results_hilic, "hil"))

  if (!is.null(alignment_reconciliation) && nrow(alignment_reconciliation) > 0) {
    alignment_reconciliation <- alignment_reconciliation |>
      dplyr::mutate(duplicate_rows = submitted_rows - unique_pairs)
    rev_save_table(alignment_reconciliation, "alignment_reconciliation",
                   "whicap")
    message("\n3,343 vs 3,325 -- reconciliation with the submitted tables:")
    print(alignment_reconciliation, width = Inf)
    message("Submitted total: ", sum(alignment_reconciliation$submitted_rows),
            " rows; unique pairs: ",
            sum(alignment_reconciliation$unique_pairs),
            "; current crosswalk: ",
            sum(alignment_reconciliation$current_pairs))
  }
}


# 7. Join ----------------------------------------------------------------------

concordance_all <- align_xwalk |>
  dplyr::inner_join(
    salsa_std |>
      dplyr::select(met, platform, exposure_lab, logFC, P.Value, adj.P.Val,
                    salsa_sd_log2, salsa_se_logfc,
                    salsa_beta_std, salsa_se_std,
                    salsa_beta_ugm3, salsa_se_ugm3),
    by = c("met", "platform"), relationship = "many-to-many"
  ) |>
  dplyr::inner_join(
    whicap_std |>
      dplyr::select(platform, whicap_mz, whicap_time,
                    whicap_feature = feature, whicap_name = name,
                    whicap_confidence = confidence,
                    whicap_estimate = estimate, whicap_se = std.error,
                    whicap_p = p_value, whicap_fdr = fdr_q_value,
                    whicap_beta_std, whicap_se_std,
                    whicap_beta_ugm3, whicap_se_ugm3),
    by = c("platform", "whicap_mz", "whicap_time"),
    relationship = "many-to-one"
  ) |>
  dplyr::distinct(exposure_lab, platform, met, whicap_feature, .keep_all = TRUE) |>
  dplyr::mutate(same_direction = sign(logFC) == sign(whicap_estimate))

message("\nAligned features with an estimate in BOTH cohorts:")
print(dplyr::count(concordance_all, exposure_lab, platform))

rev_save_table(concordance_all, "concordance_all_features", "whicap")


# 8. Cross-study concordance -- the primary analysis ------------------------------
#
## THIS IS A CONCORDANCE ANALYSIS, NOT A META-ANALYSIS, and the coefficients are
## used exactly as each study fitted them. No rescaling is applied to either
## side. That is a deliberate simplification, and it is defensible:
##
##   - Pooling two estimates requires a common scale, because the pooled number
##     is an average and an average of differently-scaled quantities is
##     meaningless. Nothing is pooled here.
##   - Direction concordance -- the statistic this analysis leads with -- is
##     invariant to ANY positive rescaling of either cohort's coefficients,
##     per feature or overall. It cannot be affected by the choice.
##   - The correlations are descriptive. Rescaling changes them, so section 8b
##     reports what standardizing would do rather than leaving it unstated.
##
## Two SALSA exposures are compared with the WHICAP PM2.5 results: PM2.5 alone,
## the like-for-like contrast, and the PC1 index, the paper's primary exposure.
## Both are the 5-year models exactly as reported in the manuscript.

concordance_stats <- function(d, exposure_lab, platform_lab) {
  if (nrow(d) < 3) return(NULL)
  bt <- stats::binom.test(sum(d$same_direction, na.rm = TRUE),
                          sum(!is.na(d$same_direction)), p = 0.5)
  ct <- stats::cor.test(d$logFC, d$whicap_estimate)
  sp <- stats::cor.test(d$logFC, d$whicap_estimate, method = "spearman",
                        exact = FALSE)
  tibble::tibble(
    exposure = exposure_lab, platform = platform_lab, n = nrow(d),
    pct_same_direction = round(100 * mean(d$same_direction, na.rm = TRUE), 1),
    dir_lo = round(100 * bt$conf.int[1], 1),
    dir_hi = round(100 * bt$conf.int[2], 1),
    direction_p = signif(bt$p.value, 3),
    pearson_r = round(unname(ct$estimate), 3),
    pearson_lo = round(ct$conf.int[1], 3),
    pearson_hi = round(ct$conf.int[2], 3),
    pearson_p = signif(ct$p.value, 3),
    spearman_r = round(unname(sp$estimate), 3),
    spearman_p = signif(sp$p.value, 3))
}

by_exposure_platform <- function(data, fn) {
  tidyr::expand_grid(
    exposure_lab = unique(data$exposure_lab),
    platform_lab = c("Both platforms", "C18/neg-", "HILIC/pos+")) |>
    purrr::pmap(function(exposure_lab, platform_lab) {
      d <- dplyr::filter(data, exposure_lab == !!exposure_lab)
      d <- switch(platform_lab,
                  "Both platforms" = d,
                  "C18/neg-"       = dplyr::filter(d, platform == "c18"),
                  "HILIC/pos+"     = dplyr::filter(d, platform == "hil"))
      fn(d, exposure_lab, platform_lab)
    }) |>
    purrr::list_rbind()
}

concordance_summary <- by_exposure_platform(concordance_all, concordance_stats)

rev_save_table(concordance_summary, "concordance_summary", "whicap")
message("\nCross-study concordance, coefficients as fitted:")
print(concordance_summary, width = Inf)


## 8a. Selection inflation ------------------------------------------------------
##
## Reviewer 1's fourth bullet. The submitted correlations were computed on
## features selected at pooled p < 0.05; these rows show how much that inflates
## the apparent agreement. They are reported to quantify the problem, not as
## findings.

selection_check <- list(
  list(lab = "all aligned features (PRIMARY)", f = function(d) d),
  list(lab = "SALSA p < 0.05",   f = function(d) dplyr::filter(d, P.Value < 0.05)),
  list(lab = "WHICAP p < 0.05",  f = function(d) dplyr::filter(d, whicap_p < 0.05)),
  list(lab = "both p < 0.05",    f = function(d) dplyr::filter(d, P.Value < 0.05,
                                                               whicap_p < 0.05))
) |>
  purrr::map(function(x) {
    unique(concordance_all$exposure_lab) |>
      purrr::map(function(e) {
        d <- x$f(dplyr::filter(concordance_all, exposure_lab == e))
        res <- concordance_stats(d, e, "Both platforms")
        if (is.null(res)) return(NULL)
        dplyr::mutate(res, subset = x$lab, .after = exposure)
      }) |>
      purrr::list_rbind()
  }) |>
  purrr::list_rbind() |>
  dplyr::select(exposure, subset, n, pct_same_direction, pearson_r, spearman_r)

rev_save_table(selection_check, "selection_inflation", "whicap")
message("\nWhat significance-based selection does to the same comparison:")
print(selection_check, width = Inf)


## 8b. Does the scale choice matter? ---------------------------------------------
##
## The primary analysis leaves both cohorts' coefficients alone. Standardizing
## each to "per 1 SD of exposure, per 1 SD of feature abundance" is the obvious
## alternative, and it is not neutral: dividing each feature by its own SD
## reorders features, so the correlations move. Direction concordance does not
## move at all, which is the reason it leads.

scale_sensitivity <- unique(concordance_all$exposure_lab) |>
  purrr::map(function(e) {
    d <- concordance_all |>
      dplyr::filter(exposure_lab == e,
                    is.finite(salsa_beta_std), is.finite(whicap_beta_std))
    tibble::tibble(
      exposure = e,
      n = nrow(d),
      pearson_as_fitted   = round(stats::cor(d$logFC, d$whicap_estimate), 3),
      pearson_standardized = round(stats::cor(d$salsa_beta_std,
                                              d$whicap_beta_std), 3),
      spearman_as_fitted  = round(stats::cor(d$logFC, d$whicap_estimate,
                                             method = "spearman"), 3),
      spearman_standardized = round(stats::cor(d$salsa_beta_std,
                                               d$whicap_beta_std,
                                               method = "spearman"), 3),
      pct_same_direction_BOTH_SCALES =
        round(100 * mean(d$same_direction, na.rm = TRUE), 1))
  }) |>
  purrr::list_rbind()

rev_save_table(scale_sensitivity, "scale_sensitivity", "whicap")
message("\nAs fitted versus standardized -- direction concordance is identical ",
        "by construction, the correlations are not:")
print(scale_sensitivity, width = Inf)


# 9. Figures ---------------------------------------------------------------------

pacman::p_load("patchwork")

EXPOSURE_COLS <- c("PM2.5 (5-year)" = "#B73F42",
                   "PC1 index"      = "#436C85")

## Platforms are combined in one panel per exposure. Both carry the same log2
## abundance scale, so there is nothing to separate them on, and pooling them
## makes the shape of the cloud readable. Colour is the project's shared
## sig_class() ladder applied to the SALSA side -- the same rule and the same
## palette as the MWAS volcanoes, so a reader moving between figures is not
## asked to learn a second convention -- and shape carries the platform.

plot_dat <- concordance_all |>
  dplyr::filter(is.finite(logFC), is.finite(whicap_estimate)) |>
  dplyr::mutate(
    exposure = factor(exposure_lab, levels = names(EXPOSURE_COLS)),
    platform_lab = factor(dplyr::if_else(platform == "c18", "C18/neg-",
                                         "HILIC/pos+"),
                          levels = c("C18/neg-", "HILIC/pos+")),
    significant = sig_class(P.Value, adj.P.Val))

message("\nSALSA significance classes in the scatter:")
print(dplyr::count(plot_dat, exposure, significant))

plot_lab <- concordance_summary |>
  dplyr::filter(platform == "Both platforms") |>
  dplyr::mutate(
    exposure = factor(exposure, levels = names(EXPOSURE_COLS)),
    label = paste0("n = ", n,
                   "\nsame direction ", sprintf("%.1f", pct_same_direction),
                   "% (", sprintf("%.1f", dir_lo), "-",
                   sprintf("%.1f", dir_hi), ")",
                   "\nr = ", sprintf("%+.3f", pearson_r),
                   ", rho = ", sprintf("%+.3f", spearman_r)))

p_scatter <- ggplot(plot_dat, aes(x = logFC, y = whicap_estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  ## Drawn in three layers, least to most significant, so the handful of
  ## FDR < 0.05 features -- 17 for PM2.5 and 26 for PC1, the ones the paper
  ## actually reports -- are never buried under 2,800 grey points.
  geom_point(data = ~ dplyr::filter(.x, significant == "NS"),
             aes(colour = significant, shape = platform_lab),
             alpha = 0.35, size = 1.1) +
  geom_point(data = ~ dplyr::filter(.x, significant %in% c("P < 0.05",
                                                           "FDR < 0.10")),
             aes(colour = significant, shape = platform_lab),
             alpha = 0.85, size = 1.8) +
  geom_point(data = ~ dplyr::filter(.x, significant == "FDR < 0.05"),
             aes(colour = significant, shape = platform_lab),
             alpha = 0.95, size = 2.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7,
              colour = "grey25", fill = "grey80") +
  geom_text(data = plot_lab, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.12, size = 3.1, lineheight = 1.05,
            inherit.aes = FALSE) +
  ## NS is given an explicit swatch rather than left as an unexplained grey.
  scale_colour_manual(values = c(SIG_COLORS, "NS" = "grey72"),
                      name = "SALSA significance", drop = FALSE,
                      breaks = SIG_LEVELS) +
  scale_shape_manual(values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
                     name = "Column") +
  facet_wrap(~ exposure, nrow = 1) +
  labs(title = "SALSA against WHICAP, every aligned feature",
       subtitle = paste0("Coefficients as each study fitted them. No rescaling, ",
                         "no significance-based selection.\nColour is SALSA ",
                         "significance; the fitted line uses all features ",
                         "regardless of it."),
       x = "SALSA coefficient (log2 abundance per unit of the exposure shown)",
       y = "WHICAP PM2.5 coefficient") +
  theme_bw(base_size = 11) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        legend.box = "horizontal")

rev_save_plot(p_scatter, "concordance_scatter", "whicap",
              width = 11, height = 6)

## Direction concordance with its binomial CI -- the scale-free statistic.
p_dir <- concordance_summary |>
  dplyr::mutate(exposure = factor(exposure, levels = names(EXPOSURE_COLS)),
                platform = factor(platform,
                                  levels = c("Both platforms", "C18/neg-",
                                             "HILIC/pos+"))) |>
  ggplot(aes(x = pct_same_direction, y = forcats::fct_rev(platform),
             colour = exposure)) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = dir_lo, xmax = dir_hi), orientation = "y",
                width = 0, linewidth = 0.7,
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = EXPOSURE_COLS, name = "SALSA exposure") +
  labs(title = "Direction concordance with WHICAP PM2.5",
       subtitle = paste0("Percent of aligned features agreeing in sign, ",
                         "95% CI. 50% is chance.\nInvariant to how either ",
                         "cohort's coefficients are scaled."),
       x = "Same direction (%)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

rev_save_plot(p_dir, "direction_concordance", "whicap", width = 8, height = 4.5)


# 10. Alignment quality control ------------------------------------------------
#
## Reviewer 1's fifth bullet. Permute the WHICAP retention times within
## platform, leaving the m/z distribution untouched, and re-run the match: any
## surviving hit is a mass coincidence. The ratio of the permuted count to the
## observed count is the expected false-match rate. It is computed on the FULL
## feature list of each cohort, and at three tolerances, so the choice of
## 10 ppm / 30 s can be seen in context.

count_matches <- function(a_mz, a_rt, b_mz, b_rt, ppm, rt_tol) {
  keep_a <- is.finite(a_mz) & is.finite(a_rt)
  a_mz <- a_mz[keep_a]; a_rt <- a_rt[keep_a]
  keep_b <- is.finite(b_mz) & is.finite(b_rt)
  b_mz <- b_mz[keep_b]; b_rt <- b_rt[keep_b]
  if (length(a_mz) == 0 || length(b_mz) == 0) return(0L)

  ord <- order(b_mz)
  b_mz <- b_mz[ord]; b_rt <- b_rt[ord]

  sum(purrr::map_int(seq_along(a_mz), function(i) {
    tol <- a_mz[i] * ppm / 1e6
    lo <- findInterval(a_mz[i] - tol, b_mz) + 1L
    hi <- findInterval(a_mz[i] + tol, b_mz)
    if (hi < lo) return(0L)
    sum(abs(b_rt[lo:hi] - a_rt[i]) <= rt_tol)
  }))
}

## Full feature lists. SALSA from its own MWAS, WHICAP from the feature-level
## summaries, which are the complete post-filtering tables (3,759 C18 and
## 6,375 HILIC, matching Kalia et al.).
salsa_features <- salsa_pm25 |>
  dplyr::distinct(platform, mz, rt)
whicap_features <- whicap_summary |>
  dplyr::distinct(platform, whicap_mz, whicap_time)

## Seeded so the reported false-match rates are reproducible.
set.seed(20260913)

alignment_fmr <- tidyr::expand_grid(
  platform = c("c18", "hil"),
  FMR_GRID
) |>
  purrr::pmap(function(platform, ppm, rt) {
    s <- salsa_features  |> dplyr::filter(platform == !!platform)
    w <- whicap_features |> dplyr::filter(platform == !!platform)

    observed <- count_matches(s$mz, s$rt, w$whicap_mz, w$whicap_time, ppm, rt)
    null <- purrr::map_int(seq_len(N_PERM_ALIGN), function(i) {
      count_matches(s$mz, s$rt, w$whicap_mz, sample(w$whicap_time), ppm, rt)
    })

    tibble::tibble(
      platform         = platform,
      ppm_tol          = ppm,
      rt_tol_s         = rt,
      n_salsa          = nrow(s),
      n_whicap         = nrow(w),
      observed_matches = observed,
      null_mean        = round(mean(null), 1),
      null_p95         = round(unname(stats::quantile(null, 0.95)), 1),
      expected_fmr     = round(mean(null) / max(observed, 1), 4),
      n_permutations   = N_PERM_ALIGN
    )
  }) |>
  purrr::list_rbind()

## What the analysis actually uses is the apLCMS crosswalk, which is stricter
## than naive tolerance matching: it adjusts retention time between the two
## runs and keeps one pair per alignment group.
aplcms_counts <- align_xwalk |>
  dplyr::count(platform, name = "aplcms_aligned_pairs") |>
  dplyr::left_join(
    ## one row per aligned PAIR, not per pair x exposure
    concordance_all |>
      dplyr::distinct(platform, met, whicap_feature) |>
      dplyr::count(platform, name = "pairs_with_both_estimates"),
    by = "platform")

alignment_fmr <- alignment_fmr |>
  dplyr::left_join(aplcms_counts, by = "platform")

rev_save_table(alignment_fmr, "alignment_false_match_rate", "whicap")
print(alignment_fmr, width = Inf)

message("\n`expected_fmr` is the share of naive tolerance matches attributable ",
        "to mass coincidence alone. Quote the ", PPM_TOL, " ppm / ", RT_TOL,
        " s row alongside the aligned-feature count (Reviewer 1 comment 10).")


# 11. Making the SALSA estimate look like WHICAP's ------------------------------
#
## Reviewer 1's second bullet asks why the two coefficients are commensurate.
## Section 4 answers it descriptively; this section answers it by construction,
## rebuilding the SALSA PM2.5 estimate under WHICAP's own choices one at a time
## and re-running the concordance after each change. Four models, cumulative:
##
##   1  SALSA as published   5-year window, SALSA's own feature processing
##                           (log2, median-scaled, ComBat), no dementia term
##   2  + WHICAP's window    1-year window: the calendar year before the draw,
##                           which is exactly Kalia et al.'s definition
##   3  + WHICAP's transform log10 -> quantile-normalized -> auto-scaled, so the
##                           outcome is constructed the same way on both sides
##   4  + WHICAP's covariate dementia/CIND entered as a covariate, as WHICAP
##                           adjusts for AD status
##
## If the null concordance in section 8 were an artefact of the analytical
## choices rather than of the data, it should move along this ladder.
##
## The consensus correlation is reused from the stratum-level estimate in the
## revision checkpoint. That is the pipeline's own convention -- R4 runs with
## DUPCOR_PER_EXPOSURE = FALSE and broadcasts one consensus across every
## exposure in a stratum -- so no duplicateCorrelation is re-run here.
##
## The whole block is cached: it is the only expensive part of this script.

W1_CACHE <- file.path(rev_dir("data", "whicap"), "salsa_pm25_matched_mwas.RData")

if (!file.exists(W1_CACHE)) {

  load(here::here("data", "processed", "salsa_clean.RData"))
  load(rev_here("data", "processed", "combined_data_list_revision.RData"))
  load(rev_here("data", "processed", "exposure_windows_revision.RData"))
  load(rev_here("data", "metabolomics", "results",
                "duplicate_correlation_revision.RData"))

  ## 2-load_data.R lowercases sample names when it reads the feature tables;
  ## the stored matrices do not, and the sample link is on the lowercased form.
  mat_c18 <- med_c18_raw_combat; colnames(mat_c18) <- tolower(colnames(mat_c18))
  mat_hil <- med_hil_raw_combat; colnames(mat_hil) <- tolower(colnames(mat_hil))

  frame <- combined_data_list_revision[["total"]][["all"]][["covar"]] |>
    dplyr::left_join(
      window_wide[["w1"]] |>
        dplyr::select(rand_id, blood_date, exp_pm2.5_w1 = exp_pm2.5),
      by = c("rand_id", "blood_date"))
  stopifnot(!anyNA(frame$exp_pm2.5_w1))

  ## Scaled by its own IQR, the convention R4 applies to every pollutant in
  ## every window.
  w1_iqr <- stats::IQR(frame$exp_pm2.5_w1, na.rm = TRUE)
  frame$exp_pm2.5_w1_iqr <- frame$exp_pm2.5_w1 / w1_iqr

  message("\n1-year PM2.5: mean ", round(mean(frame$exp_pm2.5_w1), 2),
          " ug/m3, SD ", round(stats::sd(frame$exp_pm2.5_w1), 2),
          ", IQR ", round(w1_iqr, 2),
          "; correlation with the 5-year average r = ",
          round(stats::cor(frame$exp_pm2.5_w1, frame$exp_pm2.5), 3))

  ## demcind is the stratifying variable in total/all, so the published model
  ## does NOT carry it; covars_list_new is the set R4 actually fits.
  COVARS_PUBLISHED <- covars_list_new[["covar"]]
  COVARS_WHICAP    <- covar_list[["covar"]]      # the same set plus demcind

  ## WHICAP's feature processing. The stored matrices are log2, so log10 is a
  ## constant rescale (log10(x) = log2(x) * log10(2)); it is applied explicitly
  ## rather than skipped so the step is visible. Quantile normalization is
  ## across samples, auto-scaling is per feature, which leaves every feature at
  ## unit variance -- the reason a model fitted on it needs no further division
  ## by a feature SD.
  whicap_transform <- function(mat) {
    x <- mat * log10(2)
    x <- limma::normalizeQuantiles(x)
    x <- t(scale(t(x)))
    x[!is.finite(x)] <- NA_real_
    x
  }

  fit_one <- function(metabo, sample_link, consensus, exposure, covars) {
    link <- sample_link |>
      dplyr::select(rand_id, blood_date, wave, batch, file.name_new) |>
      dplyr::filter(file.name_new %in% colnames(metabo), !is.na(rand_id)) |>
      dplyr::distinct()
    dat <- frame |>
      dplyr::inner_join(link, by = c("rand_id", "blood_date", "wave", "batch")) |>
      dplyr::arrange(match(file.name_new, colnames(metabo)))
    mat <- metabo[, dat$file.name_new, drop = FALSE]
    stopifnot(ncol(mat) == nrow(dat))

    design <- model.matrix(
      as.formula(paste0("~ ", exposure, " + ", paste(covars, collapse = " + "))),
      data = dat)
    stopifnot(nrow(design) == ncol(mat))

    fit <- limma::lmFit(mat, design, block = dat$rand_id,
                        correlation = consensus)
    fit <- limma::eBayes(fit)
    limma::topTable(fit, coef = 2, sort.by = "p", number = nrow(mat),
                    adjust.method = "BH") |>
      tibble::rownames_to_column("met") |>
      tibble::as_tibble()
  }

  cons <- c(c18 = dupcor_c18_list$total$all$covar[[1]]$consensus.correlation,
            hil = dupcor_hilic_list$total$all$covar[[1]]$consensus.correlation)

  ## The ladder. `feature_sd` says how a coefficient becomes a per-SD-of-feature
  ## coefficient: "native" divides by the feature's own SD, "unit" needs no
  ## division because auto-scaling already set it to 1.
  matched_spec <- tibble::tribble(
    ~step, ~label,                               ~exposure,           ~transform, ~covars,
    1L, "1. SALSA as published (5-year)",        "exp_pm2.5_iqr",     "native",   "published",
    2L, "2. + WHICAP 1-year window",             "exp_pm2.5_w1_iqr",  "native",   "published",
    3L, "3. + WHICAP feature transform",         "exp_pm2.5_w1_iqr",  "whicap",   "published",
    4L, "4. + WHICAP dementia adjustment",       "exp_pm2.5_w1_iqr",  "whicap",   "whicap"
  )

  mats <- list(c18 = mat_c18, hil = mat_hil)
  links <- list(c18 = salsa_blood_date_c18, hil = salsa_blood_date_hilic)
  mats_whicap <- purrr::map(mats, whicap_transform)

  salsa_matched <- matched_spec |>
    purrr::pmap(function(step, label, exposure, transform, covars) {
      covar_set <- if (covars == "whicap") COVARS_WHICAP else COVARS_PUBLISHED
      c("c18", "hil") |>
        purrr::map(function(pf) {
          m <- if (transform == "whicap") mats_whicap[[pf]] else mats[[pf]]
          fit_one(m, links[[pf]], cons[[pf]], exposure, covar_set) |>
            dplyr::mutate(platform = pf)
        }) |>
        purrr::list_rbind() |>
        dplyr::mutate(step = step, label = label, exposure = exposure,
                      transform = transform)
    }) |>
    purrr::list_rbind()

  ## Exposure SD expressed in the units each model's coefficient is stated in.
  exposure_sd_units <- c(
    exp_pm2.5_iqr    = stats::sd(frame$exp_pm2.5_iqr,    na.rm = TRUE),
    exp_pm2.5_w1_iqr = stats::sd(frame$exp_pm2.5_w1_iqr, na.rm = TRUE))

  ## Validation: step 1 must reproduce the stored MWAS exactly, or none of the
  ## other three can be trusted.
  matched_validation <- c("c18", "hil") |>
    purrr::map(function(pf) {
      ref_file <- rev_here("tables", "mwas_results", "total", "all", "covar",
                           paste0("mwas_", if (pf == "c18") "c18" else "hilic",
                                  "_exp_pm2.5_iqr_total_all_covar.xlsx"))
      if (!file.exists(ref_file)) return(NULL)
      ref <- readxl::read_xlsx(ref_file) |> dplyr::select(met, ref_logFC = logFC)
      cmp <- salsa_matched |>
        dplyr::filter(step == 1L, platform == pf) |>
        dplyr::inner_join(ref, by = "met")
      tibble::tibble(platform = pf, n = nrow(cmp),
                     r = signif(stats::cor(cmp$logFC, cmp$ref_logFC), 10),
                     max_abs_diff = signif(max(abs(cmp$logFC - cmp$ref_logFC)), 3))
    }) |>
    purrr::list_rbind()

  ## Calendar-time diagnostic --------------------------------------------------
  ##
  ## The inverse correlation deepens when WHICAP's 1-year window is adopted, so
  ## it needs a cause. Two candidates are separated here.
  ##
  ##   (a) THE CALENDAR-TIME COVARIATE. PM2.5 fell over the years both cohorts
  ##       drew blood, so `wave` and the exposure are correlated (r = -0.51 for
  ##       the 5-year average, -0.17 for the 1-year). Adjusting for a collinear
  ##       term can move a coefficient a long way. Refitting WITHOUT `wave`
  ##       tests whether the adjustment is what produces the disagreement.
  ##
  ##   (b) A SHARED FEATURE-LEVEL TIME AXIS. Take the coefficient for `wave`
  ##       itself out of the SALSA model and correlate it, over the aligned
  ##       features, with each cohort's PM2.5 coefficient. Two coefficient
  ##       vectors that load on a common axis with OPPOSITE signs are inversely
  ##       correlated whatever the biology.
  ##
  ## Because beta_PM2.5 and beta_wave are estimated jointly, part of any
  ## correlation between them is sampling covariance, which for two terms in one
  ## design is approximately -cor(exposure, wave) across features. That
  ## benchmark is reported alongside the observed value so the two can be told
  ## apart.
  COVARS_NO_TIME <- setdiff(COVARS_PUBLISHED, "wave")

  fit_coef_named <- function(metabo, sample_link, consensus, exposure, covars,
                             coef_name) {
    link <- sample_link |>
      dplyr::select(rand_id, blood_date, wave, batch, file.name_new) |>
      dplyr::filter(file.name_new %in% colnames(metabo), !is.na(rand_id)) |>
      dplyr::distinct()
    dat <- frame |>
      dplyr::inner_join(link, by = c("rand_id", "blood_date", "wave", "batch")) |>
      dplyr::arrange(match(file.name_new, colnames(metabo)))
    mat <- metabo[, dat$file.name_new, drop = FALSE]
    design <- model.matrix(
      as.formula(paste0("~ ", exposure, " + ", paste(covars, collapse = " + "))),
      data = dat)
    fit <- limma::eBayes(limma::lmFit(mat, design, block = dat$rand_id,
                                      correlation = consensus))
    limma::topTable(fit, coef = coef_name, sort.by = "p", number = nrow(mat),
                    adjust.method = "BH") |>
      tibble::rownames_to_column("met") |>
      tibble::as_tibble()
  }

  over_platforms <- function(exposure, covars, coef_name) {
    c("c18", "hil") |>
      purrr::map(function(pf) {
        fit_coef_named(mats[[pf]], links[[pf]], cons[[pf]], exposure, covars,
                       coef_name) |>
          dplyr::mutate(platform = pf)
      }) |>
      purrr::list_rbind()
  }

  time_fits <- list(
    pm5_no_time  = over_platforms("exp_pm2.5_iqr",    COVARS_NO_TIME, 2),
    pm1_no_time  = over_platforms("exp_pm2.5_w1_iqr", COVARS_NO_TIME, 2),
    wave_in_pm5  = over_platforms("exp_pm2.5_iqr",    COVARS_PUBLISHED, "wave"),
    wave_in_pm1  = over_platforms("exp_pm2.5_w1_iqr", COVARS_PUBLISHED, "wave")
  )

  ## Collinearity benchmark for the beta-to-beta correlations.
  collinearity_benchmark <- c(
    pm5 = -stats::cor(frame$exp_pm2.5_iqr,    frame$wave),
    pm1 = -stats::cor(frame$exp_pm2.5_w1_iqr, frame$wave))

  save(salsa_matched, matched_spec, matched_validation, exposure_sd_units,
       w1_iqr, time_fits, collinearity_benchmark, file = W1_CACHE)
  rm(mats, mats_whicap, mat_c18, mat_hil)

} else {
  message("\nReusing cached window/transform-matched MWAS: ", W1_CACHE)
}

load(W1_CACHE)

rev_save_table(matched_validation, "matched_step1_validation", "whicap")
message("\nValidation -- step 1 against the stored MWAS:")
print(matched_validation, width = Inf)
if (any(matched_validation$max_abs_diff > 1e-10)) {
  warning("Step 1 does not reproduce the stored MWAS; the matched ladder is ",
          "not trustworthy.", call. = FALSE)
}

## Concordance after each step -------------------------------------------------

matched_concordance <- matched_spec$step |>
  purrr::map(function(st) {
    spec <- matched_spec |> dplyr::filter(step == st)
    res <- salsa_matched |> dplyr::filter(step == st)

    sd_exp <- exposure_sd_units[[spec$exposure]]

    res_std <- res |>
      dplyr::left_join(salsa_feature_sd, by = c("met", "platform")) |>
      dplyr::mutate(
        ## after auto-scaling the feature SD is 1 by construction
        sd_use = if (spec$transform == "whicap") 1 else salsa_sd_log2,
        salsa_beta_std = logFC * sd_exp / sd_use)

    d <- align_xwalk |>
      dplyr::inner_join(res_std |> dplyr::select(met, platform, logFC, P.Value,
                                                 adj.P.Val, salsa_beta_std),
                        by = c("met", "platform"), relationship = "many-to-one") |>
      dplyr::inner_join(whicap_std |>
                          dplyr::select(platform, whicap_mz, whicap_time,
                                        whicap_feature = feature,
                                        whicap_beta_std),
                        by = c("platform", "whicap_mz", "whicap_time"),
                        relationship = "many-to-one") |>
      dplyr::distinct(platform, met, whicap_feature, .keep_all = TRUE) |>
      dplyr::filter(is.finite(salsa_beta_std), is.finite(whicap_beta_std)) |>
      dplyr::mutate(same_direction = sign(salsa_beta_std) == sign(whicap_beta_std))

    one <- function(dd, pf_label) {
      if (nrow(dd) < 3) return(NULL)
      bt <- stats::binom.test(sum(dd$same_direction), nrow(dd), p = 0.5)
      ct <- stats::cor.test(dd$salsa_beta_std, dd$whicap_beta_std)
      tibble::tibble(
        step = st, model = spec$label, platform = pf_label, n = nrow(dd),
        pct_same_direction = round(100 * mean(dd$same_direction), 1),
        direction_p = signif(bt$p.value, 3),
        pearson_r = round(unname(ct$estimate), 3),
        pearson_ci = paste0(round(ct$conf.int[1], 3), " to ",
                            round(ct$conf.int[2], 3)),
        pearson_p = signif(ct$p.value, 3),
        salsa_fdr05 = sum(dd$adj.P.Val < 0.05, na.rm = TRUE))
    }

    dplyr::bind_rows(one(d, "both"),
                     one(dplyr::filter(d, platform == "c18"), "c18"),
                     one(dplyr::filter(d, platform == "hil"), "hil"))
  }) |>
  purrr::list_rbind()

rev_save_table(matched_concordance, "matched_concordance", "whicap")
message("\nConcordance after each step toward WHICAP's analysis:")
print(matched_concordance |> dplyr::filter(platform == "both"), width = Inf)
print(matched_concordance |> dplyr::filter(platform != "both"), width = Inf)

message("\nIf the null concordance were an artefact of how SALSA was analysed ",
        "rather than a property of the data, these rows would move.")


## Where does the inverse correlation come from? -------------------------------

whicap_beta <- whicap_std |>
  dplyr::select(platform, whicap_mz, whicap_time, whicap_feature = feature,
                whicap_beta_std)

## A SALSA coefficient is in log2 abundance units and a WHICAP coefficient is
## per 1 SD of feature, so the SALSA side must be divided by the feature's own
## SD before the two are correlated -- the same standardization the ladder in
## section 11 uses. The exposure-side constant is dropped: it is common to all
## features and a positive constant does not change a correlation.
align_to_whicap <- function(betas) {
  align_xwalk |>
    dplyr::inner_join(betas |> dplyr::select(met, platform, logFC),
                      by = c("met", "platform"), relationship = "many-to-one") |>
    dplyr::left_join(salsa_feature_sd, by = c("met", "platform")) |>
    dplyr::mutate(b = logFC / salsa_sd_log2) |>
    dplyr::inner_join(whicap_beta, by = c("platform", "whicap_mz", "whicap_time"),
                      relationship = "many-to-one") |>
    dplyr::distinct(platform, met, whicap_feature, .keep_all = TRUE) |>
    dplyr::filter(is.finite(b), is.finite(whicap_beta_std))
}

r_vs_whicap <- function(betas) {
  d <- align_to_whicap(betas)
  c(all = stats::cor(d$b, d$whicap_beta_std),
    c18 = stats::cor(d$b[d$platform == "c18"], d$whicap_beta_std[d$platform == "c18"]),
    hil = stats::cor(d$b[d$platform == "hil"], d$whicap_beta_std[d$platform == "hil"]))
}

## (a) Does the calendar-time covariate create the disagreement?
time_covariate_test <- dplyr::bind_rows(
  tibble::tibble(model = "5-year window, wave adjusted",
                 !!!as.list(round(r_vs_whicap(
                   salsa_matched |> dplyr::filter(step == 1L)), 3))),
  tibble::tibble(model = "5-year window, wave REMOVED",
                 !!!as.list(round(r_vs_whicap(time_fits$pm5_no_time), 3))),
  tibble::tibble(model = "1-year window, wave adjusted",
                 !!!as.list(round(r_vs_whicap(
                   salsa_matched |> dplyr::filter(step == 2L)), 3))),
  tibble::tibble(model = "1-year window, wave REMOVED",
                 !!!as.list(round(r_vs_whicap(time_fits$pm1_no_time), 3)))
)
rev_save_table(time_covariate_test, "time_covariate_test", "whicap")
message("\nCorrelation with the WHICAP PM2.5 coefficients, with and without ",
        "the calendar-time covariate:")
print(time_covariate_test, width = Inf)

## (b) Do both cohorts' coefficients load on SALSA's calendar-time axis?
shared_axis_test <- list(
  list(lab = "5-year", wave_fit = time_fits$wave_in_pm5,
       pm_fit = salsa_matched |> dplyr::filter(step == 1L), bench = "pm5"),
  list(lab = "1-year", wave_fit = time_fits$wave_in_pm1,
       pm_fit = salsa_matched |> dplyr::filter(step == 2L), bench = "pm1")
) |>
  purrr::map(function(x) {
    joint <- x$pm_fit |>
      dplyr::select(met, platform, b_pm = logFC) |>
      dplyr::inner_join(x$wave_fit |> dplyr::select(met, platform, b_wave = logFC),
                        by = c("met", "platform"))
    ## SALSA against SALSA: both coefficients come from one model and are
    ## already in the same log2 units, so they are compared unstandardized.
    tibble::tibble(
      window = x$lab,
      r_salsaPM25_vs_salsaTIME = round(stats::cor(joint$b_pm, joint$b_wave), 3),
      collinearity_benchmark   = round(collinearity_benchmark[[x$bench]], 3),
      r_whicapPM25_vs_salsaTIME = round(unname(r_vs_whicap(x$wave_fit)["all"]), 3))
  }) |>
  purrr::list_rbind()

rev_save_table(shared_axis_test, "shared_time_axis_test", "whicap")
message("\nSALSA's coefficient for WAVE against each cohort's PM2.5 ",
        "coefficient, with the collinearity benchmark:")
print(shared_axis_test, width = Inf)

message("\nReading: removing the calendar-time covariate collapses the 5-year ",
        "correlation but leaves the 1-year one intact, so the disagreement ",
        "under WHICAP's window is NOT an adjustment artefact. Both cohorts' ",
        "PM2.5 coefficients are related to SALSA's calendar-time coefficient ",
        "with OPPOSITE signs, and for the 1-year exposure the SALSA side ",
        "exceeds what joint estimation alone predicts. What that shared axis ",
        "physically is -- specimen chronology, storage, run order -- cannot be ",
        "settled from summary statistics alone.")


# 12. Fixed-effect meta-analysis ---------------------------------------------------
#
## Reviewer 1 comment 10 offers two routes. Section 8 reports the concordance
## (route 2). This section takes ROUTE 1 -- "standardize exposures and metabolic
## outcomes identically before pooling".
##
## TWO EXPOSURES ARE POOLED WITH WHICAP'S PM2.5, and the reason for the second
## is a substantive one raised by the co-investigators rather than a hedge:
##
##   PC1 index       PRIMARY. The unsupervised multipollutant index. WHICAP's cohort is
##                   in northern Manhattan, where ambient PM2.5 is itself a
##                   traffic- and combustion-derived mixture rather than a
##                   single agent; SALSA's PC1 is the corresponding
##                   traffic-related axis through eight LUR surfaces (benzene,
##                   1,3-butadiene, NO2, NOx, PM2.5, chromium, nickel, lead).
##                   On that reading the two indices describe the same
##                   underlying source mixture more nearly than two PM2.5 mass
##                   concentrations measured in different airsheds do. Note that
##                   PC1 is UNSUPERVISED -- Reviewer 1's objection was to
##                   "outcome-informed multipollutant indices", and PCA on the
##                   exposure surfaces is not outcome-informed, so the
##                   circularity half of the objection does not apply to it.
##   PM2.5 (5-year)  SECONDARY, reported alongside. Same pollutant on both
##                   sides, so its estimand needs no argument at all; it is the
##                   like-for-like check that the estimand problem is fixed.
##
## THAT ARGUMENT IS QUALITATIVE AND MUST BE STATED AS SUCH. A 1-SD move in PC1
## is not calibrated to a 1-SD move in Manhattan PM2.5; the two are commensurate
## in what they proxy, not in their units. Reviewer 1's objection to the
## SUBMITTED analysis -- that a mixture coefficient and a PM2.5 coefficient were
## pooled without justification -- is answered by making the argument explicitly
## and by reporting the PM2.5-with-PM2.5 pooling alongside it, not by asserting
## equivalence. Both are therefore reported.
##
## WHAT POOLING STILL CANNOT FIX, for either exposure: the alignment false-match
## rate is 40% (C18) and 48% (HILIC) at 10 ppm / 30 s (section 10), and SALSA
## supplies most of the inverse-variance weight, so a pooled estimate is a
## weighted summary rather than independent corroboration. Both travel with the
## results.

pacman::p_load("patchwork", "ggrepel", "ggtext")

## PC1 FIRST: it is the primary comparison (see the header), so it leads every
## table, every facet and the combined panel. PM2.5 follows as the like-for-like
## check.
META_EXPOSURES <- c("PC1 index", "PM2.5 (5-year)")
ANNOT_KEY <- c("PM2.5 (5-year)" = "exp_pm2.5_iqr", "PC1 index" = "comp_pca_all")

## Level 1 identifications, from OUR annotation and per exposure. Never WHICAP's:
## every WHICAP annotation in the aligned set is Schymanski Level 3 or 5, and
## naming a main figure from those while the MWAS figures refuse anything below
## Level 1 is the inconsistency Reviewer 1 comment 8 is about.
annot_path <- rev_here("data", "metabolomics", "results",
                       "mwas_annotation_revision.RData")
level1_names <- function(exposure_var) {
  if (!file.exists(annot_path)) return(tibble::tibble(met = character(),
                                                      compound_display = character()))
  if (!exists(".annot_env", envir = globalenv())) {
    ae <- new.env(); load(annot_path, envir = ae)
    assign(".annot_env", ae, envir = globalenv())
  }
  ae <- get(".annot_env", envir = globalenv())
  grab <- function(lst) {
    d <- lst[["total"]][["all"]][["covar"]][[exposure_var]]
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d |>
      dplyr::filter(confidence_level == 1, !is.na(compound), compound != "") |>
      dplyr::group_by(met) |> dplyr::slice_head(n = 1) |> dplyr::ungroup() |>
      dplyr::transmute(met, compound_display = candidate_display_name(compound))
  }
  dplyr::bind_rows(grab(ae$mwas_full_annotated_list_c18),
                   grab(ae$mwas_full_annotated_list_hilic))
}

pool_one <- function(lab) {
  d <- concordance_all |>
    dplyr::filter(exposure_lab == lab,
                  is.finite(salsa_beta_std), is.finite(whicap_beta_std),
                  is.finite(salsa_se_std),   is.finite(whicap_se_std),
                  salsa_se_std > 0, whicap_se_std > 0)
  d |>
    dplyr::mutate(
      w_salsa  = 1 / salsa_se_std^2,
      w_whicap = 1 / whicap_se_std^2,
      w_total  = w_salsa + w_whicap,
      beta_pooled = (w_salsa * salsa_beta_std + w_whicap * whicap_beta_std) / w_total,
      se_pooled   = sqrt(1 / w_total),
      z_pooled    = beta_pooled / se_pooled,
      p_pooled    = 2 * stats::pnorm(-abs(z_pooled)),
      Q   = w_salsa  * (salsa_beta_std  - beta_pooled)^2 +
            w_whicap * (whicap_beta_std - beta_pooled)^2,
      p_Q = stats::pchisq(Q, df = 1, lower.tail = FALSE),
      I2  = pmax(0, (Q - 1) / Q) * 100,
      pct_weight_salsa = round(100 * w_salsa / w_total, 1)) |>
    dplyr::group_by(platform) |>
    dplyr::mutate(fdr_pooled = stats::p.adjust(p_pooled, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::left_join(level1_names(ANNOT_KEY[[lab]]), by = "met")
}

meta_all <- META_EXPOSURES |> purrr::map(pool_one) |> purrr::list_rbind()
## Kept for the sections that consume the PM2.5 pooling by name.
meta <- dplyr::filter(meta_all, exposure_lab == "PM2.5 (5-year)")

message("\nFeatures entering the meta-analysis:")
print(dplyr::count(meta_all, exposure_lab, platform))
message("Carrying a SALSA Level 1 identification:")
print(meta_all |> dplyr::filter(!is.na(compound_display)) |>
        dplyr::count(exposure_lab))

## Ordering a discrete axis WITHIN each facet. tidytext::reorder_within does
## this, but it is a two-function idiom and not worth a package dependency here,
## so the pair is written out. The suffix makes each level unique across facets
## -- the same compound can appear under both exposures -- and the scale strips
## it again for display.
reorder_within <- function(x, by, within, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) sub(paste0(sep, ".+$"), "", x),
                            ...)
}



## Validation against meta::metagen ---------------------------------------------
##
## The pooling is written out rather than delegated to a package: six vectorised
## lines over the whole table, where metagen() would be one model fit per
## feature. "I implemented it myself" is not a reason for a reviewer to trust
## it, so a sample is refitted with meta::metagen() and the two must agree.
if (requireNamespace("meta", quietly = TRUE)) {
  chk <- meta_all |> dplyr::slice_head(n = 300)
  ref <- purrr::map_dfr(seq_len(nrow(chk)), function(i) {
    m <- meta::metagen(TE   = c(chk$salsa_beta_std[i], chk$whicap_beta_std[i]),
                       seTE = c(chk$salsa_se_std[i],   chk$whicap_se_std[i]),
                       common = TRUE, random = FALSE)
    tibble::tibble(
      te = if (!is.null(m$TE.common))   m$TE.common   else m$TE.fixed,
      se = if (!is.null(m$seTE.common)) m$seTE.common else m$seTE.fixed,
      p  = if (!is.null(m$pval.common)) m$pval.common else m$pval.fixed,
      Q  = m$Q, pQ = m$pval.Q, I2 = m$I2 * 100)
  })
  metagen_check <- tibble::tibble(
    quantity = c("pooled estimate", "pooled SE", "pooled p", "Cochran Q",
                 "p(Q)", "I2"),
    max_abs_diff = c(max(abs(chk$beta_pooled - ref$te), na.rm = TRUE),
                     max(abs(chk$se_pooled   - ref$se), na.rm = TRUE),
                     max(abs(chk$p_pooled    - ref$p),  na.rm = TRUE),
                     max(abs(chk$Q           - ref$Q),  na.rm = TRUE),
                     max(abs(chk$p_Q         - ref$pQ), na.rm = TRUE),
                     max(abs(chk$I2          - ref$I2), na.rm = TRUE)),
    n_features = nrow(chk))
  rev_save_table(metagen_check, "metagen_validation", "whicap")
  message("\nValidation against meta::metagen():")
  print(metagen_check, width = Inf)
  if (any(metagen_check$max_abs_diff > 1e-8)) {
    warning("The hand-coded pooling has diverged from meta::metagen().",
            call. = FALSE)
  }
}


## Tables ------------------------------------------------------------------------

meta_table <- meta_all |>
  dplyr::transmute(
    exposure = exposure_lab, platform, met, compound_display, whicap_feature,
    salsa_beta_std, salsa_se_std, salsa_p = P.Value, salsa_fdr = adj.P.Val,
    whicap_beta_std, whicap_se_std, whicap_p, whicap_fdr,
    beta_pooled, se_pooled,
    ci_lo = beta_pooled - 1.96 * se_pooled,
    ci_hi = beta_pooled + 1.96 * se_pooled,
    z_pooled, p_pooled, fdr_pooled,
    pct_weight_salsa, Q, p_Q, I2, same_direction) |>
  dplyr::arrange(exposure, p_pooled)

rev_save_table(meta_table, "meta_pooled", "whicap")
rev_save_table(dplyr::filter(meta_table, fdr_pooled < 0.05),
               "meta_pooled_significant", "whicap")

meta_summary <- meta_all |>
  dplyr::group_by(exposure_lab, platform) |>
  dplyr::summarise(
    n_features         = dplyr::n(),
    n_salsa_fdr05      = sum(adj.P.Val < 0.05, na.rm = TRUE),
    n_whicap_fdr05     = sum(whicap_fdr < 0.05, na.rm = TRUE),
    n_pooled_fdr05     = sum(fdr_pooled < 0.05),
    n_both_p05_same_dir = sum(P.Value < 0.05 & whicap_p < 0.05 &
                                same_direction, na.rm = TRUE),
    pct_same_direction = round(100 * mean(same_direction), 1),
    median_pct_weight_salsa = round(stats::median(pct_weight_salsa), 1),
    pct_Q_rejects      = round(100 * mean(p_Q < 0.05, na.rm = TRUE), 1),
    median_I2          = round(stats::median(I2, na.rm = TRUE), 1),
    .groups = "drop")

rev_save_table(meta_summary, "meta_pooled_summary", "whicap")
message("\nFixed-effect meta-analysis with WHICAP PM2.5:")
print(meta_summary, width = Inf)

## The named metabolites a reader will look for.
meta_level1 <- meta_table |>
  dplyr::filter(!is.na(compound_display)) |>
  dplyr::arrange(exposure, p_pooled) |>
  dplyr::transmute(exposure, compound = compound_display, platform,
                   salsa_p = signif(salsa_p, 2), whicap_p = signif(whicap_p, 2),
                   pooled_p = signif(p_pooled, 2),
                   pooled_fdr = signif(fdr_pooled, 2),
                   direction = dplyr::if_else(same_direction, "same", "opposite"),
                   pct_weight_salsa)
rev_save_table(meta_level1, "meta_pooled_level1", "whicap")
message("\nLevel 1 identifications, strongest pooled evidence per exposure:")
print(meta_level1 |> dplyr::group_by(exposure) |> dplyr::slice_head(n = 8),
      width = Inf, n = 20)


## 12b. Combined panel --------------------------------------------------------------
##
## Rows are exposures, not panel types: A-B are PM2.5 with WHICAP PM2.5, C-D the
## PC1 index with the same WHICAP estimates. Each row is a volcano of the pooled
## estimate beside the two cohorts' coefficients against each other, so the two
## comparisons are read the same way and can be compared directly. The
## heterogeneity distribution and the Level 1 forests move to the Supporting
## Information, where they are given for both exposures.

PLATFORM_LAB <- c(c18 = "C18/neg-", hil = "HILIC/pos+")

## CROSS-COHORT CONCORDANCE CLASS, which is what the points are coloured by.
##
## Colouring by pooled FDR says only how strong the combined evidence is, and in
## a figure whose subject is whether two cohorts agree that is the less useful
## of the two things a colour can carry -- particularly since SALSA supplies
## ~90% of the pooled weight, so pooled significance largely tracks SALSA alone.
## The classes below say which cohort supports a feature and whether they agree
## in sign, so a reader can see the concordant features directly. The scheme and
## the palette are the ones scripts/9-salsa_whicap_alignment.R used for the
## submitted Figure S9, so the encoding is already familiar to the co-authors.
CONC_LEVELS <- c("Both P < 0.05 (same direction)",
                 "Both P < 0.05 (opposite direction)",
                 "P < 0.05 in WHICAP", "P < 0.05 in SALSA", "NS")
CONC_COLORS <- c("Both P < 0.05 (same direction)"     = "#B73F42",
                 "Both P < 0.05 (opposite direction)" = "#436C85",
                 "P < 0.05 in WHICAP"                 = "#7E9A6C",
                 "P < 0.05 in SALSA"                  = "#DE9960",
                 "NS"                                 = "grey75")

panel_data <- function(lab) {
  meta_all |>
    dplyr::filter(exposure_lab == lab) |>
    dplyr::mutate(
      platform_lab = unname(PLATFORM_LAB[platform]),
      significant  = sig_class(p_pooled, fdr_pooled),
      concordance = factor(dplyr::case_when(
        P.Value < 0.05 & whicap_p < 0.05 & same_direction ~ CONC_LEVELS[1],
        P.Value < 0.05 & whicap_p < 0.05                  ~ CONC_LEVELS[2],
        whicap_p < 0.05                                   ~ CONC_LEVELS[3],
        P.Value  < 0.05                                   ~ CONC_LEVELS[4],
        TRUE                                              ~ "NS"),
        levels = CONC_LEVELS),
      neg_log10_p  = -log10(p_pooled))
}

volcano_panel <- function(d, title_txt, tag_txt) {
  labelable <- dplyr::filter(d, !is.na(compound_display))
  labs_df <- dplyr::bind_rows(
    dplyr::filter(labelable, fdr_pooled < 0.10),
    labelable |> dplyr::filter(p_pooled < 0.05) |>
      dplyr::arrange(p_pooled) |> dplyr::slice_head(n = 8)) |>
    dplyr::distinct(met, .keep_all = TRUE) |>
    dplyr::mutate(label = paste0(stringr::str_trunc(compound_display, 22), " (",
                                 dplyr::if_else(platform == "c18", "C18", "HILIC"),
                                 ")"))

  ggplot(d, aes(x = beta_pooled, y = neg_log10_p)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    ## Layered so the concordant features sit on top of 3,000 grey points.
    ## NS is mapped through the colour aes (rather than given a fixed grey) so
    ## that its legend key is drawn; a fixed colour leaves the key blank.
    geom_point(data = ~ dplyr::filter(.x, concordance == "NS"),
               aes(colour = concordance, shape = platform_lab),
               alpha = 0.40, size = 1.2) +
    geom_point(data = ~ dplyr::filter(.x, concordance %in% CONC_LEVELS[3:4]),
               aes(colour = concordance, shape = platform_lab),
               alpha = 0.80, size = 1.9) +
    geom_point(data = ~ dplyr::filter(.x, concordance == CONC_LEVELS[2]),
               aes(colour = concordance, shape = platform_lab),
               alpha = 0.90, size = 2.4) +
    geom_point(data = ~ dplyr::filter(.x, concordance == CONC_LEVELS[1]),
               aes(colour = concordance, shape = platform_lab),
               alpha = 0.95, size = 2.8) +
    ggrepel::geom_label_repel(
      data = labs_df, aes(label = label), size = 3.3,
      max.overlaps = Inf, box.padding = 0.9, point.padding = 0.4,
      force = 12, max.iter = 50000, min.segment.length = 0,
      segment.colour = "grey50", segment.size = 0.3, seed = 42,
      show.legend = FALSE) +
    ## override.aes goes HERE, on the scale, not in a figure-level `& guides()`:
    ## applied with `&` it overrides the cohort panel's `guide = "none"` and
    ## prints that panel's colour legend a second time under its default name.
    scale_colour_manual(values = CONC_COLORS, name = "Cross-cohort",
                        breaks = CONC_LEVELS, drop = FALSE,
                        guide = guide_legend(
                          override.aes = list(size = 2.8, alpha = 1))) +
    scale_shape_manual(values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
                       name = "Column",
                       guide = guide_legend(override.aes = list(size = 2.8))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(title = title_txt, x = "Pooled standardized coefficient",
         y = expression(-log[10](italic(P))), tag = tag_txt) +
    theme_bw(base_size = 11)
}

cohort_panel <- function(d, title_txt, tag_txt) {
  lab_df <- c("c18", "hil") |>
    purrr::map(function(pf) {
      dd <- dplyr::filter(d, platform == pf)
      tibble::tibble(
        platform_lab = unname(PLATFORM_LAB[pf]),
        label = paste0("n = ", nrow(dd), "\nsame direction ",
                       sprintf("%.1f", 100 * mean(dd$same_direction)), "%",
                       "\nboth P < 0.05, same sign: ",
                       sum(dd$concordance == CONC_LEVELS[1])))
    }) |>
    purrr::list_rbind()

  ggplot(d, aes(x = salsa_beta_std, y = whicap_beta_std)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey55", linewidth = 0.4) +
    geom_point(data = ~ dplyr::filter(.x, concordance == "NS"),
               aes(colour = concordance), alpha = 0.35, size = 1.0) +
    geom_point(data = ~ dplyr::filter(.x, concordance %in% CONC_LEVELS[3:4]),
               aes(colour = concordance), alpha = 0.80, size = 1.8) +
    geom_point(data = ~ dplyr::filter(.x, concordance == CONC_LEVELS[2]),
               aes(colour = concordance), alpha = 0.90, size = 2.2) +
    geom_point(data = ~ dplyr::filter(.x, concordance == CONC_LEVELS[1]),
               aes(colour = concordance), alpha = 0.95, size = 2.6) +
    geom_text(data = lab_df, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.08, vjust = 1.15, size = 3.1, lineheight = 1.05,
              inherit.aes = FALSE) +
    scale_colour_manual(values = CONC_COLORS, breaks = CONC_LEVELS,
                        drop = FALSE, guide = "none") +
    facet_wrap(~ platform_lab, scales = "free") +
    labs(title = title_txt, x = "SALSA standardized coefficient",
         y = "WHICAP standardized coefficient", tag = tag_txt) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"))
}

d_pm <- panel_data("PM2.5 (5-year)")
d_pc <- panel_data("PC1 index")

PANEL_TEXT <- theme(
  panel.grid    = element_blank(),
  plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
  plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey30"),
  axis.title    = element_text(face = "bold", size = 12),
  axis.text     = element_text(size = 11),
  strip.text    = element_text(face = "bold", size = 11),
  legend.title  = element_text(face = "bold", size = 11),
  legend.text   = element_text(size = 10),
  plot.tag      = element_text(face = "bold", size = 22))

## Titles name BOTH sides of the pooling. "Pooled PC1 mixture estimate" left the
## reader to work out what it was pooled with, and since the WHICAP side is
## PM2.5 under both rows that is exactly the thing worth stating.
panels <- list(
  volcano_panel(d_pc, "Pooled analysis: SALSA PC1 mixture & WHICAP PM2.5", "A"),
  cohort_panel(d_pc,  "Cohort-specific estimates: SALSA PC1 mixture", "B"),
  volcano_panel(d_pm, "Pooled analysis: SALSA PM2.5 & WHICAP PM2.5", "C"),
  cohort_panel(d_pm,  "Cohort-specific estimates: SALSA PM2.5", "D")) |>
  purrr::map(~ .x + PANEL_TEXT)

meta_panel <- patchwork::wrap_plots(panels, ncol = 2, widths = c(4.6, 5.4)) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.justification = "center",
        legend.background = element_rect(colour = "grey45", fill = NA,
                                         linewidth = 0.4),
        legend.margin = margin(5, 9, 5, 9),
        legend.spacing.x = unit(10, "pt"),
        legend.box.margin = margin(4, 0, 2, 0))

ggplot2::ggsave(rev_here("figures", "whicap", "combined_panel_meta.png"),
                meta_panel, width = 17, height = 13, dpi = 300, bg = "white")
message("  combined panel -> ",
        rev_here("figures", "whicap", "combined_panel_meta.png"))


## 12c. Supporting figures: heterogeneity and the Level 1 forests -------------------

het_lab <- meta_all |>
  dplyr::mutate(platform_lab = unname(PLATFORM_LAB[platform])) |>
  dplyr::group_by(exposure_lab, platform_lab) |>
  dplyr::summarise(n = dplyr::n(),
                   pct_het = 100 * mean(p_Q < 0.05, na.rm = TRUE),
                   med_I2 = stats::median(I2, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(label = paste0("median *I*<sup>2</sup> = ", round(med_I2), "%<br>",
                               "Q rejects homogeneity: ", round(pct_het),
                               "% of ", n))

p_het <- meta_all |>
  dplyr::mutate(platform_lab = unname(PLATFORM_LAB[platform])) |>
  ggplot(aes(x = I2)) +
  geom_histogram(aes(fill = platform_lab), bins = 30, colour = "white",
                 linewidth = 0.2, show.legend = FALSE) +
  ggtext::geom_richtext(data = het_lab, aes(x = Inf, y = Inf, label = label),
                        hjust = 1.02, vjust = 1.12, size = 3.0,
                        fill = NA, label.colour = NA, inherit.aes = FALSE) +
  scale_fill_manual(values = c("C18/neg-"   = unname(SIG_COLORS[["FDR < 0.10"]]),
                               "HILIC/pos+" = unname(SIG_COLORS[["P < 0.05"]]))) +
  facet_grid(exposure_lab ~ platform_lab, scales = "free_y") +
  labs(title = "Between-cohort heterogeneity",
       subtitle = "Per feature, Cochran Q on 1 df",
       x = expression(italic(I)^2~"(%)"), y = "Features") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"))

rev_save_plot(p_het, "meta_heterogeneity", "whicap", width = 10, height = 7)

## Cohort-specific and pooled estimates side by side, for the Level 1 features.
## Reviewer 1's third bullet made visual. A palette outside SIG_COLORS, because
## red/blue/orange mean significance classes everywhere else in this figure set.
COHORT_COLORS <- c(SALSA = "#6A51A3", WHICAP = "#2C7873", Pooled = "#101010")
COHORT_SHAPES <- c(SALSA = 16, WHICAP = 17, Pooled = 18)

forest_src <- meta_all |>
  dplyr::filter(!is.na(compound_display)) |>
  dplyr::group_by(exposure_lab) |>
  dplyr::arrange(p_pooled, .by_group = TRUE) |>
  dplyr::slice_head(n = 12) |>
  dplyr::ungroup()

if (nrow(forest_src) > 0) {
  fdat <- forest_src |>
    dplyr::mutate(label = paste0(stringr::str_trunc(compound_display, 24), " (",
                                 dplyr::if_else(platform == "c18", "C18", "HILIC"),
                                 ")")) |>
    dplyr::select(exposure_lab, label, SALSA = salsa_beta_std,
                  WHICAP = whicap_beta_std, Pooled = beta_pooled,
                  se_s = salsa_se_std, se_w = whicap_se_std, se_p = se_pooled,
                  p_pooled) |>
    tidyr::pivot_longer(c(SALSA, WHICAP, Pooled), names_to = "cohort",
                        values_to = "beta") |>
    dplyr::mutate(
      se = dplyr::case_when(cohort == "SALSA" ~ se_s,
                            cohort == "WHICAP" ~ se_w, TRUE ~ se_p),
      lo = beta - 1.96 * se, hi = beta + 1.96 * se,
      cohort = factor(cohort, levels = c("SALSA", "WHICAP", "Pooled")),
      label = reorder_within(label, -p_pooled, exposure_lab))

  p_forest_cohorts <- ggplot(fdat, aes(x = beta, y = label, colour = cohort,
                                       shape = cohort)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                  linewidth = 0.6, position = position_dodge(width = 0.7)) +
    geom_point(position = position_dodge(width = 0.7), size = 2.4) +
    scale_colour_manual(values = COHORT_COLORS, name = NULL) +
    scale_shape_manual(values = COHORT_SHAPES, name = NULL) +
    scale_y_reordered() +
    facet_wrap(~ exposure_lab, scales = "free_y", ncol = 1) +
    labs(title = "Cohort-specific and pooled estimates",
         subtitle = "Level 1 identifications, 12 lowest pooled P; 95% CI. WHICAP intervals are wide because n = 107.",
         x = "Standardized coefficient (per 1 SD exposure, per 1 SD abundance)",
         y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  rev_save_plot(p_forest_cohorts, "meta_forest_cohort_specific", "whicap",
                width = 10, height = 9)
}


## 12d. Pathway enrichment on each pooled result ------------------------------------
##
## Mummichog on the pooled statistics, one run per exposure.
##
## BACKGROUND CAVEAT, which must appear in the figure legend because the figure
## carries no caption: mummichog infers enrichment against the feature list it
## is given, and that list is the aligned features, not the ~20,000 of either
## platform. These p-values are NOT comparable with the primary pathway
## analysis; they answer a narrower question -- which pathways are enriched
## among features measured in both cohorts.

source(here::here("scripts", "mummichog_pathway.R"))
MUMMICHOG_PERM <- rev_n(1000, 50)

pooled_pathway_tables <- META_EXPOSURES |>
  purrr::set_names() |>
  purrr::map(function(lab) {
    slug <- if (lab == "PC1 index") "pc1" else "pm25"
    cache <- file.path(rev_dir("data", "whicap"),
                       paste0("meta_pathway_pooled_", slug, ".RData"))
    if (!file.exists(cache)) {
      d <- dplyr::filter(meta_all, exposure_lab == lab)
      inp <- d |>
        dplyr::transmute(
          `m.z` = as.numeric(stringr::str_extract(met, "(?<=^mz_rt_)[0-9.]+")),
          rt    = as.numeric(stringr::str_extract(met, "[0-9.]+$")),
          `p.value` = p_pooled, t.score = z_pooled,
          mode = dplyr::if_else(platform == "c18", "negative", "positive")) |>
        dplyr::filter(is.finite(`m.z`), is.finite(rt), is.finite(`p.value`)) |>
        dplyr::arrange(`p.value`)
      in_dir  <- rev_dir("metaboAnalyst", "Input",  paste0("whicap_meta_", slug))
      out_dir <- rev_dir("metaboAnalyst", "Output", paste0("whicap_meta_", slug))
      in_file <- file.path(in_dir, paste0("mwas_meta_pooled_", slug, ".txt"))
      utils::write.table(inp, file = in_file, row.names = FALSE,
                         col.names = TRUE, quote = FALSE, sep = "\t")
      res <- tryCatch(
        run_mummichog(input_file = in_file, output_dir = out_dir, p_cutoff = 0.1,
                      organism = "hsa_mfn", instrument_ppm = 10.0,
                      ion_mode = "mixed",
                      adducts = c("M-H [1-]", "M-2H [2-]", "M-H2O-H [1-]",
                                  "M [1+]", "M+H [1+]", "M+Na [1+]"),
                      min_hits = 3, num_permutations = MUMMICHOG_PERM),
        error = function(e) { warning("Mummichog failed for ", lab, ": ",
                                      e$message, call. = FALSE); NULL })
      pooled_pathway_table <- if (!is.null(res)) as.data.frame(res$result_table) else NULL
      save(pooled_pathway_table, file = cache)
    }
    e <- new.env(); load(cache, envir = e); e$pooled_pathway_table
  })

path_plot_dat <- META_EXPOSURES |>
  purrr::map(function(lab) {
    tbl <- pooled_pathway_tables[[lab]]
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    pcol <- intersect(c("P(Fisher)", "FET", "P(Gamma)", "Gamma"), names(tbl))[1]
    hcol <- intersect(c("Hits.sig", "Hits_sig"), names(tbl))[1]
    tcol <- intersect(c("Hits.total", "Pathway.total", "Hits.all"), names(tbl))[1]
    all_p <- tbl |> tibble::as_tibble() |>
      dplyr::transmute(pathway, p_value = .data[[pcol]],
                       n_sig = .data[[hcol]], n_total = .data[[tcol]]) |>
      dplyr::filter(is.finite(p_value), n_total > 0)
    ## The background hit rate MUST come from the whole table; computing it after
    ## the p < 0.05 filter rescales every enrichment by the significant subset's
    ## own rate and pushes enriched pathways below 1.
    overall <- sum(all_p$n_sig) / sum(all_p$n_total)
    all_p |>
      dplyr::filter(p_value < 0.05, n_sig >= 2) |>
      dplyr::mutate(exposure = lab,
                    enrichment = (n_sig / n_total) / overall,
                    neg_log10_p = -log10(p_value))
  }) |>
  purrr::list_rbind()

if (!is.null(path_plot_dat) && nrow(path_plot_dat) > 0) {
  rev_save_table(path_plot_dat, "meta_pathway_pooled", "whicap")

  pd <- path_plot_dat |>
    dplyr::mutate(exposure = factor(exposure, levels = META_EXPOSURES),
                  pathway = stringr::str_trunc(pathway, 40),
                  row = reorder_within(pathway, neg_log10_p, exposure))

  ## x is -log10(P), not the enrichment factor: enrichment is a ratio of hit
  ## rates and peaks on pathways with three or four hits, so putting it on x
  ## makes the least reliable rows visually dominant.
  p_pathway <- ggplot(pd, aes(x = neg_log10_p, y = row)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey55") +
    geom_segment(aes(x = 0, xend = neg_log10_p, yend = row),
                 colour = "grey78", linewidth = 0.4) +
    geom_point(aes(size = n_sig, fill = enrichment), shape = 21,
               colour = "grey25", stroke = 0.3) +
    scale_fill_gradient(low = "#F3E4E4", high = "#B73F42",
                        name = "Enrichment\nfactor") +
    scale_size_continuous(range = c(3, 9), name = "Significant\nhits") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    scale_y_reordered() +
    facet_wrap(~ exposure, scales = "free_y", ncol = 1) +
    labs(title = "Pathway enrichment, pooled with WHICAP PM2.5",
         x = expression(bold(-log[10](italic(P)))), y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          legend.position = "right")

  rev_save_plot(p_pathway, "meta_pathway_pooled", "whicap",
                width = 11, height = 8)

  message("\nPooled pathways at P < 0.05:")
  print(pd |> dplyr::select(exposure, pathway, p_value, n_sig, n_total,
                            enrichment) |>
          dplyr::mutate(p_value = signif(p_value, 3),
                        enrichment = round(enrichment, 2)),
        n = 40, width = Inf)
}


# 13. Specification sensitivity ---------------------------------------------------
#
## Everything fitted along the way, on ONE statistic: direction concordance.
## That choice removes the scale question entirely -- sign agreement is
## invariant to any positive rescaling of either cohort's coefficients, per
## feature or overall, so these rows are comparable across specifications whose
## outcomes are on different scales (steps 3 and 4 auto-scale the features,
## the rest do not). Spearman on the coefficients as fitted is carried
## alongside as a descriptive second column.

pacman::p_load("patchwork")

spec_rows <- function(betas, exposure_lab, spec_lab) {
  d <- align_xwalk |>
    dplyr::inner_join(betas |> dplyr::select(met, platform, logFC),
                      by = c("met", "platform"), relationship = "many-to-one") |>
    dplyr::inner_join(whicap_beta_raw,
                      by = c("platform", "whicap_mz", "whicap_time"),
                      relationship = "many-to-one") |>
    dplyr::distinct(platform, met, whicap_feature, .keep_all = TRUE) |>
    dplyr::filter(is.finite(logFC), is.finite(whicap_estimate)) |>
    dplyr::mutate(same_direction = sign(logFC) == sign(whicap_estimate))

  c("Both platforms", "C18/neg-", "HILIC/pos+") |>
    purrr::map(function(pf) {
      dd <- switch(pf, "Both platforms" = d,
                   "C18/neg-"   = dplyr::filter(d, platform == "c18"),
                   "HILIC/pos+" = dplyr::filter(d, platform == "hil"))
      bt <- stats::binom.test(sum(dd$same_direction), nrow(dd), p = 0.5)
      tibble::tibble(
        exposure = exposure_lab, spec = spec_lab, platform = pf, n = nrow(dd),
        pct_dir = 100 * mean(dd$same_direction),
        dir_lo = 100 * bt$conf.int[1], dir_hi = 100 * bt$conf.int[2],
        dir_p = bt$p.value,
        spearman_r = round(stats::cor(dd$logFC, dd$whicap_estimate,
                                      method = "spearman"), 3))
    }) |>
    purrr::list_rbind()
}

whicap_beta_raw <- whicap_std |>
  dplyr::select(platform, whicap_mz, whicap_time, whicap_feature = feature,
                whicap_estimate = estimate)

## The PC1 calendar-time refit is cached by an earlier run; section 12 no longer
## builds it, because the pooling is PM2.5-only. Its rows are included when the
## cache is present and quietly skipped when it is not.
pc1_wave_cache <- file.path(rev_dir("data", "whicap"), "pc1_no_wave_mwas.RData")
have_pc1_no_wave <- file.exists(pc1_wave_cache)
if (have_pc1_no_wave) load(pc1_wave_cache)

spec_correlations <- dplyr::bind_rows(
  spec_rows(salsa_pm25, "PM2.5 (5-year)", "1. as published (PRIMARY)"),
  spec_rows(time_fits$pm5_no_time, "PM2.5 (5-year)", "2. no calendar-time term"),
  spec_rows(dplyr::filter(salsa_matched, step == 2L), "PM2.5 (1-year)",
            "3. WHICAP's 1-year window"),
  spec_rows(time_fits$pm1_no_time, "PM2.5 (1-year)",
            "4. 1-year, no time term"),
  spec_rows(dplyr::filter(salsa_matched, step == 3L), "PM2.5 (1-year)",
            "5. 1-year + WHICAP transform"),
  spec_rows(dplyr::filter(salsa_matched, step == 4L), "PM2.5 (1-year)",
            "6. 1-year + transform + dementia"),
  spec_rows(salsa_pc1, "PC1 index", "7. as published (PRIMARY)"),
  if (have_pc1_no_wave)
    spec_rows(pc1_no_wave, "PC1 index", "8. no calendar-time term")
)

rev_save_table(spec_correlations, "specification_sensitivity", "whicap")
message("\nSpecification sensitivity, direction concordance:")
print(spec_correlations |>
        dplyr::filter(platform == "Both platforms") |>
        dplyr::mutate(dplyr::across(c(pct_dir, dir_lo, dir_hi), ~ round(.x, 1))) |>
        dplyr::select(exposure, spec, n, pct_dir, dir_lo, dir_hi, spearman_r),
      width = Inf)

SPEC_COLS <- c("PM2.5 (5-year)" = "#B73F42",
               "PM2.5 (1-year)" = "#DE9960",
               "PC1 index"      = "#436C85")

p_spec <- spec_correlations |>
  dplyr::mutate(
    exposure = factor(exposure, levels = names(SPEC_COLS)),
    platform = factor(platform, levels = c("Both platforms", "C18/neg-",
                                           "HILIC/pos+")),
    row = forcats::fct_rev(factor(spec, levels = sort(unique(spec))))) |>
  ggplot(aes(x = pct_dir, y = row, colour = exposure)) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = dir_lo, xmax = dir_hi), orientation = "y",
                width = 0, linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = SPEC_COLS, name = "SALSA exposure") +
  facet_wrap(~ platform, nrow = 1) +
  labs(title = "Direction concordance with WHICAP PM2.5, by specification",
       subtitle = paste0("Percent of aligned features agreeing in sign, 95% CI. ",
                         "50% is chance.\nSign agreement is invariant to how ",
                         "either cohort's coefficients are scaled."),
       x = "Same direction (%)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom")

rev_save_plot(p_spec, "specification_sensitivity", "whicap",
              width = 12, height = 5.5)

message("\nFigures written to ", rev_here("figures", "whicap"))


# 14. Pathway-level concordance -------------------------------------------------
#
## THE COMPARISON THAT DOES NOT NEED FEATURE ALIGNMENT.
##
## Everything above compares coefficients feature by feature, and that comparison
## is hobbled before it starts: aligning two independently acquired untargeted
## datasets by mass and retention time carries a 40% (C18) to 48% (HILIC)
## false-match rate (section 10), so a large minority of the pairs are not the
## same metabolite. No amount of care with scales fixes that.
##
## The pathway level sidesteps it completely. Both studies ran MUMMICHOG against
## the SAME human MFN reference library, so the two sets of results are named in
## the same vocabulary and can be compared directly -- no mass matching, no
## retention-time tolerance, no false-match rate. This is the level at which
## both papers actually make their claims ("amino acid metabolism"), and it is
## where the two cohorts do agree.
##
## LIMITATION, stated plainly: Kalia et al. report pathway NAMES in their
## Results and Figure 2C but do not publish per-pathway p-values or enrichment
## factors, and there is no supplementary table of them. So the WHICAP side is
## set membership only, transcribed below from the published text. That allows a
## test of whether WHICAP's pathways sit unusually high in SALSA's ranking, but
## not a correlation of enrichment strengths.

## Transcribed verbatim from Kalia et al. 2023, Results: "The metabolic features
## associated with PM2.5 enriched several metabolic pathways including: alanine
## and aspartate metabolism, the TCA cycle, glutamate metabolism, glycolysis and
## gluconeogenesis, butanoate metabolism, pyruvate metabolism, methionine and
## cysteine metabolism, tyrosine metabolism, fatty acid oxidation, vitamin A
## metabolism, glycerophospholipid metabolism, and aminosugars metabolism."
## (mummichog, human MFN, feature p < 0.01, Fisher p < 0.1.)
##
## Each is given the exact name it carries in the mummichog library, so the join
## is exact rather than fuzzy. All twelve are present in SALSA's own output.
whicap_pathways <- c(
  "Alanine and Aspartate Metabolism",
  "TCA cycle",
  "Glutamate metabolism",
  "Glycolysis and Gluconeogenesis",
  "Butanoate metabolism",
  "Pyruvate Metabolism",
  "Methionine and cysteine metabolism",
  "Tyrosine metabolism",
  "Fatty acid oxidation",
  "Vitamin A (retinol) metabolism",
  "Glycerophospholipid metabolism",
  "Aminosugars metabolism")

pathway_long <- readxl::read_xlsx(
  rev_here("tables", "pathway", "pathway_results_long.xlsx"))

salsa_pathways <- pathway_long |>
  dplyr::filter(algorithm == "mummichog", study == "total",
                population == "all", covar_set == "covar",
                exposure %in% c("comp_pca_all", "exp_pm2.5_iqr")) |>
  dplyr::mutate(exposure_lab = dplyr::if_else(exposure == "comp_pca_all",
                                              "PC1 index", "PM2.5 (5-year)"),
                in_whicap = pathway %in% whicap_pathways)

missing_names <- setdiff(whicap_pathways, unique(salsa_pathways$pathway))
if (length(missing_names) > 0) {
  warning("WHICAP pathway names absent from the SALSA mummichog output -- ",
          "the library versions may differ: ",
          paste(missing_names, collapse = "; "), call. = FALSE)
}
message("\nWHICAP pathways matched into the SALSA library: ",
        length(intersect(whicap_pathways, unique(salsa_pathways$pathway))),
        " of ", length(whicap_pathways))

## Two tests, because the thresholded one throws away most of the information.
##
##   overlap    hypergeometric on SALSA's p < 0.05 set against WHICAP's set
##   ranking    Wilcoxon rank-sum on SALSA's -log10(p) for WHICAP's pathways
##              against every other pathway in the library. This uses all of
##              SALSA's p-values and only needs set membership on the WHICAP
##              side, so it is much better powered than the 2x2.

pathway_concordance <- unique(salsa_pathways$exposure_lab) |>
  purrr::map(function(e) {
    d <- dplyr::filter(salsa_pathways, exposure_lab == e)
    N <- nrow(d)                                   # pathways in the library
    K <- sum(d$in_whicap)                          # WHICAP's, present here
    sig <- dplyr::filter(d, p_value < 0.05)
    n <- nrow(sig); k <- sum(sig$in_whicap)

    wt <- stats::wilcox.test(-log10(d$p_value[d$in_whicap]),
                             -log10(d$p_value[!d$in_whicap]),
                             alternative = "greater")
    tibble::tibble(
      exposure = e,
      pathways_tested = N,
      whicap_pathways_in_library = K,
      salsa_sig_p05 = n,
      overlap = k,
      expected_overlap = round(n * K / N, 2),
      fold_enrichment = round(k / (n * K / N), 2),
      hypergeometric_p = signif(
        stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE), 3),
      median_rank_whicap = stats::median(rank(d$p_value)[d$in_whicap]),
      median_rank_other  = stats::median(rank(d$p_value)[!d$in_whicap]),
      rank_test_p = signif(wt$p.value, 3))
  }) |>
  purrr::list_rbind()

rev_save_table(pathway_concordance, "pathway_concordance", "whicap")
message("\nPathway-level concordance with WHICAP (mummichog, human MFN, both):")
print(pathway_concordance, width = Inf)

## Which pathways, named -- the table a reader will actually want.
pathway_overlap_detail <- salsa_pathways |>
  dplyr::filter(in_whicap) |>
  dplyr::select(exposure = exposure_lab, pathway, salsa_p = p_value,
                salsa_fdr = p_fdr, n_sig, n_total) |>
  dplyr::mutate(salsa_p = signif(salsa_p, 3),
                salsa_fdr = signif(salsa_fdr, 3),
                salsa_sig = salsa_p < 0.05) |>
  dplyr::arrange(exposure, salsa_p)

rev_save_table(pathway_overlap_detail, "pathway_overlap_detail", "whicap")
message("\nWHICAP's twelve pathways, with SALSA's own p-value for each:")
print(pathway_overlap_detail, width = Inf, n = 30)


# 15. What Figure S9 of the submitted manuscript was measuring ---------------------
#
## The submitted Supporting Information reported Pearson correlations of
## 0.38-0.45 between SALSA mixture-composite coefficients and WHICAP PM2.5
## coefficients, and Figure S9 drew them. Those numbers are very different from
## anything the revised analysis produces, and the reason is worth establishing
## rather than asserting, because it also answers whether the unsupervised PC1
## index would reproduce them.
##
## HOW FIGURE S9 WAS BUILT (scripts/9-salsa_whicap_alignment.R):
##   * SALSA contributed `estimate_salsa = logFC` -- the RAW coefficient, in log2
##     abundance per IQR-scaled unit -- and WHICAP contributed `estimate` as
##     supplied, which is already per 1 SD exposure and 1 SD feature. The two
##     sides were pooled by meta::metagen() on those mismatched scales.
##   * The SALSA side was an OUTCOME-INFORMED mixture composite (WQS, QGcomp),
##     not PM2.5, so the pooled contrast was undefined -- Reviewer 1's comment.
##   * The scatter and its printed correlation were computed ONLY on features
##     with pooled p < 0.05 (`dplyr::filter(!is.na(p_fe), p_fe < 0.05)` feeding
##     ggpubr::stat_cor).
##
## The third of these turns out to be the whole story. Selecting on the pooled
## p-value keeps features where the two cohorts' estimates REINFORCE each other,
## because the pooled statistic is a weighted sum of them; the surviving subset
## is concordant by construction. The permutation below makes that concrete: the
## WHICAP features are scrambled so the cohorts are unrelated BY CONSTRUCTION,
## the pooling and the p < 0.05 selection are re-applied, and the correlation is
## recomputed. If selection is what produces 0.38-0.45, the scrambled data will
## produce it too.

N_PERM_SELECT <- rev_n(200, 20)
set.seed(20260914)

fe_pvalue <- function(b1, s1, b2, s2) {
  w1 <- 1 / s1^2; w2 <- 1 / s2^2; wt <- w1 + w2
  2 * stats::pnorm(-abs(((w1 * b1 + w2 * b2) / wt) / sqrt(1 / wt)))
}

selection_artifact_row <- function(d, exposure_lab, scale_lab) {
  if (scale_lab == "as fitted") {
    b1 <- d$logFC; s1 <- d$salsa_se_logfc
    b2 <- d$whicap_estimate; s2 <- d$whicap_se
  } else {
    b1 <- d$salsa_beta_std; s1 <- d$salsa_se_std
    b2 <- d$whicap_beta_std; s2 <- d$whicap_se_std
  }
  ok <- is.finite(b1) & is.finite(s1) & is.finite(b2) & is.finite(s2) &
        s1 > 0 & s2 > 0
  b1 <- b1[ok]; s1 <- s1[ok]; b2 <- b2[ok]; s2 <- s2[ok]

  sel <- fe_pvalue(b1, s1, b2, s2) < 0.05

  ## Null: the same selection applied to cohorts that cannot agree.
  null <- purrr::map_dfr(seq_len(N_PERM_SELECT), function(i) {
    j <- sample(length(b2)); bb <- b2[j]; ss <- s2[j]
    k <- fe_pvalue(b1, s1, bb, ss) < 0.05
    if (sum(k) < 10) return(tibble::tibble(r = NA_real_, dir = NA_real_))
    tibble::tibble(r = stats::cor(b1[k], bb[k]),
                   dir = mean(sign(b1[k]) == sign(bb[k])))
  })

  tibble::tibble(
    exposure = exposure_lab, scale = scale_lab,
    n_all = length(b1),
    r_all = round(stats::cor(b1, b2), 3),
    pct_dir_all = round(100 * mean(sign(b1) == sign(b2)), 1),
    n_selected = sum(sel),
    r_selected = round(stats::cor(b1[sel], b2[sel]), 3),
    pct_dir_selected = round(100 * mean(sign(b1[sel]) == sign(b2[sel])), 1),
    null_r_mean = round(mean(null$r, na.rm = TRUE), 3),
    null_r_sd   = round(stats::sd(null$r, na.rm = TRUE), 3),
    null_pct_dir = round(100 * mean(null$dir, na.rm = TRUE), 1),
    z_vs_null = round((stats::cor(b1[sel], b2[sel]) -
                         mean(null$r, na.rm = TRUE)) /
                        stats::sd(null$r, na.rm = TRUE), 2))
}

selection_artifact <- purrr::pmap(
  tidyr::expand_grid(
    exposure_lab = c("PM2.5 (5-year)", "PC1 index"),
    scale_lab = c("as fitted", "standardized")),
  function(exposure_lab, scale_lab) {
    selection_artifact_row(
      dplyr::filter(concordance_all, exposure_lab == !!exposure_lab),
      exposure_lab, scale_lab)
  }) |>
  purrr::list_rbind()

rev_save_table(selection_artifact, "selection_artifact", "whicap")
message("\nWhat selecting on the pooled p-value does, against a null in which ",
        "the two cohorts cannot agree:")
print(selection_artifact, width = Inf)

## Reproduce the submitted Figure S9's own numbers from its own object, so the
## comparison is against what was actually published rather than a recollection.
legacy_mixture <- here::here("data", "metabolomics", "alignment",
                             "meta_mixture_salsa_whicap.RData")
if (file.exists(legacy_mixture)) {
  le <- new.env(); load(legacy_mixture, envir = le)
  figS9_reproduction <- c("comp_wqs_all", "comp_qgcomp_all") |>
    purrr::map(function(nm) {
      d <- dplyr::bind_rows(le$meta_results_c18$total$all$covar[[nm]],
                            le$meta_results_hilic$total$all$covar[[nm]])
      if (is.null(d) || nrow(d) == 0) return(NULL)
      sel <- d$p_fe < 0.05 & !is.na(d$p_fe)
      tibble::tibble(
        submitted_exposure = nm,
        n_rows = nrow(d),
        r_ALL_features = round(stats::cor(d$estimate_salsa, d$estimate_whicap,
                                          use = "complete.obs"), 3),
        n_pFE_lt_05 = sum(sel),
        r_SELECTED = round(stats::cor(d$estimate_salsa[sel],
                                      d$estimate_whicap[sel],
                                      use = "complete.obs"), 3),
        pct_dir_SELECTED = round(100 * mean(
          sign(d$estimate_salsa[sel]) == sign(d$estimate_whicap[sel]),
          na.rm = TRUE), 1))
    }) |>
    purrr::list_rbind()

  rev_save_table(figS9_reproduction, "figureS9_reproduction", "whicap")
  message("\nThe submitted Figure S9, recomputed from its own saved object:")
  print(figS9_reproduction, width = Inf)
  message("The published correlations (0.38 and 0.45) are the r_SELECTED ",
          "column. Over ALL aligned features the same composites give ",
          paste(figS9_reproduction$r_ALL_features, collapse = " and "), ".")
}


# 16. Save ---------------------------------------------------------------------

save(concordance_all, concordance_summary,
     whicap_scale_audit, exposure_comparison, alignment_fmr,
     matched_concordance, matched_validation, matched_spec,
     time_covariate_test, shared_axis_test,
     pathway_concordance, pathway_overlap_detail, alignment_reconciliation,
     meta_all, meta_table, meta_summary, meta_level1,
     pooled_pathway_tables,
     metagen_check, selection_artifact,
     selection_check,
     scale_sensitivity, spec_correlations,
     whicap_std, salsa_std, align_xwalk,
     file = file.path(rev_dir("data", "whicap"), "whicap_concordance.RData"))

message("\nR9 SALSA-WHICAP PM2.5 analysis complete.")
message("Manuscript: the reported analysis is the PM2.5 fixed-effect ",
        "meta-analysis of section 12 (Reviewer 1 comment 10, route 1: ",
        "standardize identically, then pool). Two statements must travel with ",
        "it -- SALSA carries a median 91% of the inverse-variance weight, and ",
        "feature-level direction concordance is 47.4%, at chance. Do NOT ",
        "describe the comparison as supporting the robustness of the SALSA ",
        "findings.")

}  # end if (!skip_r9)

#--------------------------------End of the code--------------------------------
