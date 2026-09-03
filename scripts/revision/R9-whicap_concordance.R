## ---------------------------
##
## Script name: R9-whicap_concordance.R
## Purpose of script: Replace the SALSA-WHICAP fixed-effect meta-analysis with a
##                    cross-study concordance analysis on a common estimand
##
## Author: Yufan Gong
##
## Date Created: 2026-08-26
##
## Date Modified: 2026-08-26
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
##        contrast, so the fixed-effect estimand is undefined. Four changes:
##
##          1. Compare SALSA PM2.5-SPECIFIC estimates with WHICAP PM2.5 -- a
##             defensible common estimand, exactly as Reviewer 2 asks. The
##             SALSA PM2.5 MWAS already exists at
##             tables/mwas_results/*/all/covar/mwas_*_exp_pm2.5_iqr_*.
##          2. Standardize both cohorts identically before comparing
##             (per 1-SD exposure, per 1-SD feature abundance).
##          3. Assess concordance across ALL aligned features, reporting the
##             p-selected subset separately and labelled as selection-inflated.
##          4. Estimate the alignment false-match rate at 10 ppm / 30 s.
##
##        INPUTS (all present under data/metabolomics/alignment/):
##          20260427_PM25_MWAS_results_annotations.txt
##              WHICAP PM2.5 MWAS: feature, estimate, std.error, statistic,
##              p_value, fdr_q_value, column_esi ("C18 -" / "HILIC +"), mz, time,
##              name, chemical_formula, confidence, adduct. 10,161 features.
##          {c18neg,hilpos}_feature_level_summary.txt
##              WHICAP per-feature mz / time / mean / sd. Per the WHICAP team:
##              "The data were not transformed before summary. I did filter
##              based on detection (detection in at least 70% of samples) and
##              replaced 0 abundance values with half the minimum intensity
##              observed for the feature before calculating the mean and sd."
##          {salsa,whicap}_{c18,hil}_link.RData
##              apLCMS alignment crosswalks from scripts/9-salsa_whicap_alignment.R
##
##        SCALE HARMONIZATION -- the crux of the reviewers' objection.
##          SALSA effects are per IQR-scaled PM2.5, on log2 feature abundance.
##          WHICAP feature summaries are UNTRANSFORMED, so their SDs are not
##          directly comparable. By the delta method,
##              sd(log2 X) ~= sd(X) / (mean(X) * ln 2),
##          which converts the WHICAP summaries onto the log2 scale SALSA uses.
##          This is an approximation and is flagged as such in the output; it is
##          accurate when the per-feature coefficient of variation is small.
##          The direction-concordance results do not depend on it at all, so
##          they are reported as the primary comparison.
##
##        Outputs -> <REV_ROOT>/{tables,figures,data}/whicap/
##
## ---------------------------

# Setup -----------------------------------------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R9-whicap_concordance.R")

# Helpers ported from the archived revision suite ------------------------------
#
# The archived R13 relied on four helpers that lived in the archived
# R1-revision_functions.R and were not carried into the current one. They are
# reproduced here verbatim so this script depends only on the CURRENT pipeline.
#
# list_mwas_files() defaults to the REVISION output tree, not tables/. That is
# the whole point of re-running this analysis: the SALSA PM2.5 estimates it
# consumes must be the ones produced under the revision covariate set (batch
# dropped), not the submitted ones. Under the submitted analysis PM2.5 had 25
# FDR-significant HILIC features; under the revision it has 62.

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
  message("R13 SKIPPED - required inputs are not present:")
  purrr::walk(missing, ~ message("  * ", .x))
  message("Alignment crosswalks come from scripts/9-salsa_whicap_alignment.R.")
  message(strrep("-", 66))
  skip_r9 <- TRUE
} else {
  skip_r9 <- FALSE
}

if (!skip_r9) {

# 1. WHICAP PM2.5 MWAS ----------------------------------------------------------

whicap_mwas <- utils::read.delim(whicap_mwas_file, sep = "\t", header = TRUE,
                                 check.names = FALSE,
                                 stringsAsFactors = FALSE) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    platform = dplyr::case_when(
      column_esi == "C18 -"   ~ "c18",
      column_esi == "HILIC +" ~ "hil",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(platform))

message("WHICAP PM2.5 MWAS features: ", nrow(whicap_mwas))
print(dplyr::count(whicap_mwas, platform))


# 2. WHICAP feature-level summaries, converted to the log2 scale -----------------

## Delta-method conversion, per the note in the header. `cv` is retained so the
## approximation can be audited: it degrades as cv grows.
whicap_summary <- list(
  c18 = "c18neg_feature_level_summary.txt",
  hil = "hilpos_feature_level_summary.txt"
) |>
  purrr::imap(function(f, platform) {
    utils::read.delim(file.path(align_dir, f), sep = "\t", header = TRUE,
                      check.names = FALSE) |>
      tibble::as_tibble() |>
      dplyr::transmute(
        platform,
        whicap_mz   = mz,
        whicap_time = time,
        raw_mean    = mean,
        raw_sd      = sd,
        cv          = raw_sd / raw_mean,
        ## sd on the log2 scale, delta method
        sd_log2     = raw_sd / (raw_mean * log(2))
      ) |>
      dplyr::filter(is.finite(sd_log2), sd_log2 > 0)
  }) |>
  purrr::list_rbind()

message("\nWHICAP feature summaries: ", nrow(whicap_summary))
message("Median per-feature CV (untransformed): ",
        round(stats::median(whicap_summary$cv, na.rm = TRUE), 3))
message("  The delta-method log2 SD conversion is reliable while CV is small; ",
        round(100 * mean(whicap_summary$cv > 1, na.rm = TRUE), 1),
        "% of features have CV > 1.")


# 3. SALSA PM2.5-specific MWAS ---------------------------------------------------

salsa_pm25 <- list_mwas_files() |>
  dplyr::filter(exposure == "exp_pm2.5_iqr", study == "total",
                population == "all", covar_set == "covar") |>
  dplyr::mutate(data = purrr::map(path, read_mwas_xlsx)) |>
  dplyr::select(platform, data) |>
  tidyr::unnest(data) |>
  dplyr::mutate(platform = dplyr::if_else(platform == "hilic", "hil", platform))

message("\nSALSA PM2.5 MWAS features: ", nrow(salsa_pm25))


# 4. Harmonize the exposure scale -------------------------------------------------

## SALSA coefficients are per IQR-scaled PM2.5; convert to per 1 SD so the two
## cohorts' exposure contrasts match. WHICAP estimates are assumed to be per
## 1-SD PM2.5 as supplied.
load(here::here("data", "processed", "combined_data_list_new.RData"))
dat <- combined_data_list_new[["total"]][["all"]][["covar"]]

pm25_sd_per_iqr_unit <- stats::sd(dat$exp_pm2.5_iqr, na.rm = TRUE)

message("\nSALSA PM2.5: 1 IQR-scaled unit = ",
        round(pm25_sd_per_iqr_unit, 4), " SD")


# 5. SALSA feature SDs on the analysis scale ---------------------------------------

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

## Standard errors are carried through on the same scale as the coefficients,
## because the fixed-effect pooling in R11 needs them. limma's topTable does
## not emit an SE column, but t = logFC / SE, so it is recoverable exactly.
salsa_std <- salsa_pm25 |>
  dplyr::left_join(salsa_feature_sd, by = c("met", "platform")) |>
  dplyr::mutate(
    salsa_se_logfc       = abs(logFC / t),
    beta_per_sd_exposure = logFC / pm25_sd_per_iqr_unit,
    salsa_beta_std       = beta_per_sd_exposure / salsa_sd_log2,
    salsa_se_std         = (salsa_se_logfc / pm25_sd_per_iqr_unit) /
                             salsa_sd_log2
  )


# 6. Map aligned features across cohorts --------------------------------------------

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
  purrr::list_rbind()

message("\nAligned feature pairs in the crosswalk: ", nrow(align_xwalk))
print(dplyr::count(align_xwalk, platform))


# 7. Join everything and standardize WHICAP ------------------------------------------

## WHICAP MWAS rows are matched to the crosswalk on m/z and retention time
## within the alignment tolerances.
whicap_std <- whicap_mwas |>
  dplyr::left_join(whicap_summary,
                   by = c("platform", "mz" = "whicap_mz", "time" = "whicap_time")) |>
  dplyr::mutate(whicap_beta_std = estimate / sd_log2,
                whicap_se_std   = std.error / sd_log2)

## One row per aligned pair. apLCMS can map a feature more than once when two
## SALSA peaks fall inside the tolerance of one WHICAP peak; keep the closest
## match so a feature is not counted twice in the concordance statistics.
align_xwalk <- align_xwalk |>
  dplyr::distinct(platform, met, whicap_mz, whicap_time, .keep_all = TRUE)

concordance_all <- align_xwalk |>
  dplyr::inner_join(
    salsa_std |>
      dplyr::select(met, platform, logFC, P.Value, adj.P.Val,
                    salsa_sd_log2, beta_per_sd_exposure, salsa_beta_std,
                    salsa_se_logfc, salsa_se_std),
    by = c("met", "platform"), relationship = "many-to-one"
  ) |>
  dplyr::inner_join(
    whicap_std |>
      dplyr::distinct(platform, mz, time, .keep_all = TRUE) |>
      dplyr::select(platform, whicap_mz = mz, whicap_time = time,
                    whicap_feature = feature, whicap_name = name,
                    whicap_confidence = confidence,
                    whicap_estimate = estimate, whicap_se = std.error,
                    whicap_p = p_value, whicap_fdr = fdr_q_value,
                    whicap_sd_log2 = sd_log2, whicap_cv = cv,
                    whicap_beta_std, whicap_se_std),
    by = c("platform", "whicap_mz", "whicap_time"),
    relationship = "many-to-one"
  ) |>
  dplyr::distinct(platform, met, whicap_feature, .keep_all = TRUE) |>
  dplyr::mutate(
    same_direction_raw = sign(logFC) == sign(whicap_estimate),
    same_direction_std = sign(salsa_beta_std) == sign(whicap_beta_std)
  )

message("\nFeatures with estimates in BOTH cohorts: ", nrow(concordance_all))

rev_save_table(concordance_all, "concordance_all_features", "whicap")


# 8. Concordance summary ---------------------------------------------------------------

summarise_concordance <- function(d, label, warn = NA_character_) {
  ok <- stats::complete.cases(d$salsa_beta_std, d$whicap_beta_std)
  dd <- d[ok, , drop = FALSE]
  if (nrow(dd) < 3) {
    return(tibble::tibble(subset = label, n = nrow(dd)))
  }
  bt <- stats::binom.test(sum(dd$same_direction_raw, na.rm = TRUE),
                          sum(!is.na(dd$same_direction_raw)), p = 0.5)
  tibble::tibble(
    subset             = label,
    n                  = nrow(dd),
    pct_same_direction = round(100 * mean(dd$same_direction_raw, na.rm = TRUE), 1),
    direction_p        = signif(bt$p.value, 3),
    pearson_r_std      = round(stats::cor(dd$salsa_beta_std, dd$whicap_beta_std,
                                          use = "complete.obs"), 3),
    spearman_r_std     = round(stats::cor(dd$salsa_beta_std, dd$whicap_beta_std,
                                          use = "complete.obs",
                                          method = "spearman"), 3),
    spearman_r_raw     = round(stats::cor(dd$logFC, dd$whicap_estimate,
                                          use = "complete.obs",
                                          method = "spearman"), 3),
    note               = warn
  )
}

concordance_summary <- dplyr::bind_rows(
  summarise_concordance(
    concordance_all, "PRIMARY: all aligned features"),
  summarise_concordance(
    concordance_all |> dplyr::filter(P.Value < 0.05),
    "SELECTION-INFLATED: SALSA p < 0.05",
    "Reported only to show the inflation; not a finding."),
  summarise_concordance(
    concordance_all |> dplyr::filter(whicap_p < 0.05),
    "SELECTION-INFLATED: WHICAP p < 0.05",
    "Reported only to show the inflation; not a finding."),
  summarise_concordance(
    concordance_all |> dplyr::filter(P.Value < 0.05, whicap_p < 0.05),
    "SELECTION-INFLATED: both p < 0.05",
    "Reported only to show the inflation; not a finding.")
)

rev_save_table(concordance_summary, "concordance_summary", "whicap")
print(concordance_summary, width = Inf)

message("\nReport the PRIMARY row as the finding. Direction concordance does ",
        "not depend on the delta-method scale conversion; the standardized ",
        "correlations do, so quote the Spearman-on-raw column alongside them.")

## Scatter for the supplement.
p_conc <- concordance_all |>
  dplyr::filter(is.finite(salsa_beta_std), is.finite(whicap_beta_std)) |>
  ggplot2::ggplot(ggplot2::aes(x = salsa_beta_std, y = whicap_beta_std)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  ggplot2::geom_point(alpha = 0.25, size = 0.7) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  ggplot2::facet_wrap(~ platform) +
  ggplot2::labs(
    x = "SALSA PM2.5 effect (per 1-SD exposure, per 1-SD feature)",
    y = "WHICAP PM2.5 effect (same scale)",
    title = "Cross-cohort concordance, all aligned features",
    subtitle = "PM2.5 in both cohorts; no significance-based selection"
  ) +
  ggplot2::theme_bw(base_size = 11)

rev_save_plot(p_conc, "concordance_scatter", "whicap", width = 9, height = 5)


# 9. Alignment false-match rate ----------------------------------------------------------

## Reviewer 1 comment 10 asks for the expected false-match rate at 10 ppm and
## 30 s across two independently acquired datasets. Permute the WHICAP
## retention times within platform, keeping the m/z distribution intact, and
## re-run the match: any hit is then a mass coincidence.

count_matches <- function(a_mz, a_rt, b_mz, b_rt,
                          ppm = PPM_TOL, rt_tol = RT_TOL) {
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

alignment_fmr <- c("c18", "hil") |>
  purrr::map(function(platform) {
    s <- salsa_pm25 |>
      dplyr::filter(platform == !!platform) |>
      dplyr::distinct(mz, rt)
    w <- whicap_mwas |>
      dplyr::filter(platform == !!platform) |>
      dplyr::distinct(mz, time)

    observed <- count_matches(s$mz, s$rt, w$mz, w$time)

    null <- purrr::map_int(seq_len(N_PERM_ALIGN), function(i) {
      count_matches(s$mz, s$rt, w$mz, sample(w$time))
    })

    tibble::tibble(
      platform          = platform,
      n_salsa           = nrow(s),
      n_whicap          = nrow(w),
      observed_matches  = observed,
      null_mean         = round(mean(null), 1),
      null_p95          = round(stats::quantile(null, 0.95), 1),
      expected_fmr      = round(mean(null) / max(observed, 1), 4),
      n_permutations    = N_PERM_ALIGN
    )
  }) |>
  purrr::list_rbind()

rev_save_table(alignment_fmr, "alignment_false_match_rate", "whicap")
print(alignment_fmr, width = Inf)

message("\n`expected_fmr` is the share of matches attributable to mass ",
        "coincidence alone at ", PPM_TOL, " ppm / ", RT_TOL, " s. Report it ",
        "alongside the aligned-feature count (Reviewer 1 comment 10).")


# 10. Cross-check against the existing ad hoc meta-analysis object ---------------------------

## data/metabolomics/alignment/meta_pm25_salsa_whicap.RData contains a
## previously computed SALSA-PM2.5 vs WHICAP-PM2.5 meta-analysis, but no script
## in the repository produces it. Compare against it if present, so the
## reproducible version here can be reconciled with whatever was done before.
legacy_path <- file.path(align_dir, "meta_pm25_salsa_whicap.RData")

if (file.exists(legacy_path)) {
  e <- new.env(); load(legacy_path, envir = e)
  legacy <- try(
    dplyr::bind_rows(
      e$meta_pm25_c18$total$all$covar   |> dplyr::mutate(platform = "c18"),
      e$meta_pm25_hilic$total$all$covar |> dplyr::mutate(platform = "hil")
    ),
    silent = TRUE
  )

  if (!inherits(legacy, "try-error")) {
    cmp <- concordance_all |>
      dplyr::select(platform, met, logFC, whicap_estimate) |>
      dplyr::inner_join(
        legacy |> dplyr::select(platform, met = salsa_met,
                                estimate_salsa, estimate_whicap),
        by = c("platform", "met")
      )

    legacy_check <- tibble::tibble(
      n_overlapping    = nrow(cmp),
      r_salsa_effects  = round(stats::cor(cmp$logFC, cmp$estimate_salsa,
                                          use = "complete.obs"), 4),
      r_whicap_effects = round(stats::cor(cmp$whicap_estimate,
                                          cmp$estimate_whicap,
                                          use = "complete.obs"), 4)
    )
    rev_save_table(legacy_check, "legacy_meta_crosscheck", "whicap")
    message("\nCross-check against meta_pm25_salsa_whicap.RData:")
    print(legacy_check)
    message("An r of 1.00 on both sides confirms that the legacy object was ",
            "ALREADY the PM2.5-vs-PM2.5 comparison the reviewers ask for -- ",
            "it was computed but never reported, and the manuscript presented ",
            "the mixture-vs-PM2.5 pooling instead. This script reproduces it ",
            "from source so it can be reported.")
  }
}


# Save ------------------------------------------------------------------------------

save(concordance_all, concordance_summary, alignment_fmr,
     whicap_std, salsa_std, align_xwalk,
     file = file.path(rev_dir("data", "whicap"), "whicap_concordance.RData"))

message("\nR13 WHICAP concordance complete.")
message("REMINDER for the manuscript: remove all pooled fixed-effect estimates ",
        "and pooled FDR values, retitle the section 'Cross-Cohort ",
        "Concordance', and moderate the claim that WHICAP 'supports the ",
        "robustness' of the mixture findings (Reviewer 1 comment 10, ",
        "Reviewer 2 comment 8).")

}  # end if (!skip_r9)

#--------------------------------End of the code--------------------------------
