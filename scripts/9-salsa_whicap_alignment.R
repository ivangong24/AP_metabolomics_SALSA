## ---------------------------
##
## Script name: 9-salsa_whicap_alignment.R
## Purpose of script: Feature-level alignment between SALSA and WHICAP serum
##                    metabolomics feature tables using apLCMS, as a prerequisite
##                    to feature-level meta-analysis across the two cohorts.
##
## Author: Yufan Gong
##
## Date Created: 2026-05-11
##
## Notes:
##   - Follows the same logic as scripts/archive/S4_only_alignment_JAN26.R
##     (PEG + S4 alignment by Kimberly Paul), generalized to SALSA + WHICAP.
##   - Runs alignment separately for C18-negative and HILIC-positive platforms,
##     iterated via purrr over a single `platform_inputs` config list.
##   - SALSA input: ComBat-corrected feature matrices saved under
##     `data/metabolomics/processed/` as `c18_raw_combat.Rdata` and
##     `hil_raw_combat.Rdata`. These are matrices with rownames
##     `mz_rt_<mz>_<time>` (mz and time rounded to 4 decimals) and per-sample
##     intensity columns. They retain more features than the CV-filtered
##     ComBat text tables under `data/metabolomics/{c18neg,hilicpos}/`. Per
##     Dr. Paul, raw vs. ComBat-corrected mean/SD make little practical
##     difference for alignment; we use the ComBat-corrected matrices here.
##   - WHICAP input: feature-level summary files provided by the WHICAP team
##     (`data/metabolomics/alignment/{c18neg,hilpos}_feature_level_summary.txt`),
##     each with columns mz, time, mean, sd. These summaries were built by
##     filtering to features detected in >=70% of samples and replacing 0
##     abundances with half the per-feature minimum before computing mean/sd
##     on untransformed intensities. This gives apLCMS real per-feature
##     intensity (`area` = mean) and SD weighting for both cohorts.
##
## ---------------------------

# 0. Setup -------------------------------------------------------------

# remotes::install_github("yufree/apLCMS")
pacman::p_load(apLCMS, here, tidyverse,
               matrixStats, doParallel, parallelly,
               writexl, glue, meta)

options(stringsAsFactors = FALSE)

# parallel backend (apLCMS uses foreach internally)
cl <- makeCluster(max(1, parallelly::availableCores() - 1))
registerDoParallel(cl)

# output directory for aligned features
align_dir <- here::here("data", "metabolomics", "alignment")
dir.create(align_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Helpers ----------------------------------------------------------

# Load a SALSA ComBat-corrected feature matrix (`c18_raw_combat` /
# `hil_raw_combat`) from `.Rdata` and return:
#   $df  — a data frame with mz, time, and one column per sample
#   $mat — the apLCMS feature matrix with columns:
#            col 1 = mz, col 2 = pos (retention time),
#            col 3 = SD of peak intensity, col 4 = mean peak area
# Matrix and data frame are both sorted by mz so adjusted RTs from
# adjust.time() line up row-for-row with `df` in the link step.
load_salsa_combat <- function(combat_path) {

  env <- new.env()
  load(combat_path, envir = env)
  obj_names <- ls(env)
  stopifnot(length(obj_names) == 1L)
  combat_mat <- get(obj_names[1], envir = env)

  # rownames format: "mz_rt_<mz>_<time>"  →  parts[3] = mz, parts[4] = time
  parts <- stringr::str_split_fixed(rownames(combat_mat), "_", 4)
  mz  <- as.numeric(parts[, 3])
  pos <- as.numeric(parts[, 4])

  intensities <- as.matrix(combat_mat)
  area <- rowMeans(intensities, na.rm = TRUE)
  sd1  <- matrixStats::rowSds(intensities, na.rm = TRUE)

  keep <- !is.na(mz) & !is.na(pos)
  ord  <- order(mz[keep])

  mat <- cbind(mz = mz[keep][ord],  pos = pos[keep][ord],
               sd1 = sd1[keep][ord], area = area[keep][ord])

  df <- intensities[keep, , drop = FALSE][ord, , drop = FALSE] |>
    as.data.frame(check.names = FALSE) |>
    (\(x) cbind(mz = mat[, "mz"], time = mat[, "pos"], x))()

  list(df = df, mat = mat)
}

# Read a WHICAP feature-level summary file (mz, time, mean, sd) and
# return $df (cleaned, sorted by mz) and $mat (apLCMS feature matrix).
load_whicap_summary <- function(path) {

  df <- read.table(path, sep = "\t", header = TRUE,
                   colClasses = "numeric", check.names = FALSE) |>
    dplyr::filter(!is.na(mz), !is.na(time)) |>
    dplyr::distinct(mz, time, .keep_all = TRUE) |>
    dplyr::arrange(mz)

  mat <- cbind(mz  = df$mz,  pos  = df$time,
               sd1 = df$sd,  area = df$mean)

  list(df = df, mat = mat)
}

# Run apLCMS retention-time correction + feature alignment for one
# platform. SALSA is index 1 (used by apLCMS as template), WHICAP is
# index 2; min.exp = 2 keeps only features present in BOTH cohorts.
# Tolerances: mz.tol is relative (1e-5 = 10 ppm); chr.tol is in the RT
# unit of the input feature tables (seconds here, so 30 = 30 s).
align_platform <- function(salsa_mat, whicap_mat, platform_tag,
                           mz_tol = 1e-5, chr_tol = 30) {

  features  <- list(salsa = salsa_mat, whicap = whicap_mat)
  features2 <- adjust.time(features, mz.tol = mz_tol, chr.tol = chr_tol)
  aligned   <- feature.align(features2, min.exp = 2,
                             mz.tol = mz_tol, chr.tol = chr_tol)

  comb_ftrs <- as.data.frame(aligned$aligned.ftrs)
  comb_rts  <- as.data.frame(aligned$pk.times)

  colnames(comb_ftrs)[colnames(comb_ftrs) == "exp 1"] <- "salsa_area"
  colnames(comb_ftrs)[colnames(comb_ftrs) == "exp 2"] <- "whicap_area"
  colnames(comb_rts)[colnames(comb_rts)   == "exp 1"] <- "salsa_rt"
  colnames(comb_rts)[colnames(comb_rts)   == "exp 2"] <- "whicap_rt"

  comb_ftrs$match <- rownames(comb_ftrs)
  comb_rts$match  <- rownames(comb_rts)

  comb <- merge(comb_ftrs, comb_rts,
                by = c("match", "mz", "chr", "min.mz", "max.mz")) |>
    dplyr::rename(mz.adj     = mz,
                  time.adj   = chr,
                  min.mz.adj = min.mz,
                  max.mz.adj = max.mz) |>
    dplyr::mutate(alg_match = paste0("mz_rt_",
                                     round(mz.adj, 3), "_",
                                     round(time.adj, 3)),
                  platform  = platform_tag)

  list(
    Comb_features = comb,
    features2     = features2,
    mz.tol        = aligned$mz.tol,
    chr.tol       = aligned$chr.tol,
    n_salsa_in    = nrow(salsa_mat),
    n_whicap_in   = nrow(whicap_mat),
    n_aligned     = nrow(comb)
  )
}

# Attach adjusted RT (col 2 of features2 slot) row-for-row onto a
# per-cohort data frame and merge to the alignment table by (rounded
# mz, adjusted RT). Row order in `cohort_df` must already match the
# adjust.time input matrix.
link_cohort <- function(cohort_df, mz_col, features2_slot,
                        comb_features, cohort_rt_col) {

  cohort_df <- cohort_df[order(cohort_df[[mz_col]]), ]
  cohort_df$algRT  <- features2_slot[, 2]
  cohort_df$mz.rnd <- round(cohort_df[[mz_col]], 2)
  comb_features$mz.rnd <- round(comb_features$mz.adj, 2)

  merge(comb_features, cohort_df,
        by.x = c("mz.rnd", cohort_rt_col),
        by.y = c("mz.rnd", "algRT"))
}

# 2. Platform-level inputs --------------------------------------------

platform_inputs <- list(
  c18 = list(
    tag              = "c18neg",
    salsa_combat     = here::here("data", "metabolomics", "processed",
                            "c18_raw_combat.Rdata"),
    whicap_summary   = file.path(align_dir,
                                 "c18neg_feature_level_summary.txt"),
    comb_file        = "Comb_features_c18_salsa_whicap.RData",
    salsa_link_file  = "salsa_c18_link.RData",
    whicap_link_file = "whicap_c18_link.RData"
  ),
  hil = list(
    tag              = "hilicpos",
    salsa_combat     = here::here("data", "metabolomics", "processed",
                            "hil_raw_combat.Rdata"),
    whicap_summary   = file.path(align_dir,
                                 "hilpos_feature_level_summary.txt"),
    comb_file        = "Comb_features_hilic_salsa_whicap.RData",
    salsa_link_file  = "salsa_hil_link.RData",
    whicap_link_file = "whicap_hil_link.RData"
  )
)

# 3. Load SALSA + WHICAP feature tables per platform -------------------

salsa_loaded <- platform_inputs |>
  purrr::imap(function(cfg, platform){
    message(paste0("Loading SALSA ComBat matrix for ", platform, " ..."))
    load_salsa_combat(cfg$salsa_combat)
  })

whicap_loaded <- platform_inputs |>
  purrr::imap(function(cfg, platform){
    message(paste0("Loading WHICAP feature-level summary for ",
                   platform, " ..."))
    load_whicap_summary(cfg$whicap_summary)
  })

# 4. Run apLCMS alignment per platform --------------------------------

alignment_results <- list(
  platform_inputs, salsa_loaded, whicap_loaded, names(platform_inputs)
) |>
  purrr::pmap(function(cfg, salsa, whicap, platform){
    message(paste0("Running apLCMS alignment for ", platform, " ..."))
    aln <- align_platform(salsa$mat, whicap$mat, platform_tag = cfg$tag)

    message(sprintf(
      "[%s] SALSA features: %d | WHICAP features: %d | aligned: %d (mz.tol=%.2g, chr.tol=%.2f)",
      cfg$tag, aln$n_salsa_in, aln$n_whicap_in, aln$n_aligned,
      aln$mz.tol, aln$chr.tol))

    aln
  }) |>
  purrr::set_names(names(platform_inputs))

# Expose Comb_features_{c18,hil} for downstream code + save to disk

alignment_results |>
  purrr::map(~ .x$Comb_features) |>
  purrr::set_names("Comb_features_c18", "Comb_features_hil") |>
  list2env(envir = .GlobalEnv)

list(alignment_results, platform_inputs) |>
  purrr::pwalk(function(aln, cfg){
    Comb_features <- aln$Comb_features
    save(Comb_features, file = file.path(align_dir, cfg$comb_file))
  })

# 5. Link aligned features back to each cohort's records --------------

# SALSA side
list(salsa_loaded, alignment_results, names(platform_inputs)) |>
  purrr::pmap(function(salsa, aln, platform){
    message(paste0("Linking aligned features to SALSA records for ",
                   platform, " ..."))
    link_cohort(salsa$df, mz_col = "mz",
                features2_slot = aln$features2[[1]],
                comb_features  = aln$Comb_features,
                cohort_rt_col  = "salsa_rt")
  }) |>
  purrr::set_names("salsa_c18_link", "salsa_hil_link") |>
  list2env(envir = .GlobalEnv)

# WHICAP side
list(whicap_loaded, alignment_results, names(platform_inputs)) |>
  purrr::pmap(function(whicap, aln, platform){
    message(paste0("Linking aligned features to WHICAP records for ",
                   platform, " ..."))
    link_cohort(whicap$df, mz_col = "mz",
                features2_slot = aln$features2[[2]],
                comb_features  = aln$Comb_features,
                cohort_rt_col  = "whicap_rt")
  }) |>
  purrr::set_names("whicap_c18_link", "whicap_hil_link") |>
  list2env(envir = .GlobalEnv)

# Save link tables to disk
list(
  list("salsa_c18_link",  "whicap_c18_link",
       "salsa_hil_link",  "whicap_hil_link"),
  list(platform_inputs$c18$salsa_link_file,
       platform_inputs$c18$whicap_link_file,
       platform_inputs$hil$salsa_link_file,
       platform_inputs$hil$whicap_link_file)
) |>
  purrr::pwalk(function(obj_name, out_file){
    save(list = obj_name, envir = .GlobalEnv,
         file = file.path(align_dir, out_file))
  })

# 6. Fixed-effect meta-analysis (SALSA mixtures vs WHICAP) -------------

# SALSA-side: per-feature limma logFC + moderated t (+ VIP) across every
# stratum (study × population × covariate set × exposure). The combined
# results were saved by 4-mwas_analysis.R.
load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))

salsa_combined_list <- list(c18 = combined_results_list_c18,
                            hil = combined_results_list_hilic)

# SALSA exposures to meta-analyze, per study. We use only the
# all-toxicant mixture indices: comp_wqs_all and comp_qgcomp_all in the
# cross-sectional analysis (total), and comp_qgcomp_cox_all in the
# time-to-event analysis (cox; WQS is not refit for cox per
# 4-mwas_analysis.R's exposure-exclusion rules).
salsa_exposures_by_study <- list(
  total = c("comp_wqs_all", "comp_qgcomp_all"),
  cox   = c("comp_qgcomp_cox_all")
)

# WHICAP-side: NOTE — the only per-feature WHICAP MWAS file currently
# in hand is a PM2.5 single-exposure analysis. Until WHICAP shares
# mixture-based summary stats, meta-analysis pools SALSA's all-toxicant
# mixture effect against WHICAP's PM2.5 effect at each aligned feature.
whicap_mwas_all <- read.table(
  here::here("data", "metabolomics", "results",
       "20260427_PM25_MWAS_results_annotations.txt"),
  sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE)

whicap_mwas_by_platform <- list(
  c18 = whicap_mwas_all |> dplyr::filter(column_esi == "C18 -"),
  hil = whicap_mwas_all |> dplyr::filter(column_esi == "HILIC +")
)

# Aligned-feature crosswalk per platform: alg_match → SALSA feature ID
# (`salsa_met` = "mz_rt_<mz>_<time>", matching MWAS `met` column) and
# WHICAP feature ID (`whicap_met`, same construction). `whicap_mz` /
# `whicap_time` are also carried through for the pathway-input prep step.
align_xwalk <- list(
  c18 = list(salsa = salsa_c18_link, whicap = whicap_c18_link),
  hil = list(salsa = salsa_hil_link, whicap = whicap_hil_link)
) |>
  purrr::imap(function(links, platform){
    salsa_side <- links$salsa |>
      dplyr::transmute(alg_match,
                       salsa_met = paste0("mz_rt_", mz, "_", time)) |>
      dplyr::distinct(alg_match, .keep_all = TRUE)

    whicap_side <- links$whicap |>
      dplyr::transmute(alg_match,
                       whicap_met  = paste0("mz_rt_", mz, "_", time),
                       whicap_mz   = mz,
                       whicap_time = time) |>
      dplyr::distinct(alg_match, .keep_all = TRUE)

    dplyr::inner_join(salsa_side, whicap_side, by = "alg_match")
  })

# Per-feature meta-analysis via meta::metagen. Returns both common-effect
# (fixed-effect, FE) and random-effect (RE) results in one call; with k = 2
# studies, REML for tau^2 typically collapses to 0, so DerSimonian-Laird
# (`method.tau = "DL"`) is used as the standard 2-study choice. Output:
# one-row tibble with FE + RE point estimates, SEs, z (the metagen test
# statistic), p, 95% CIs, plus heterogeneity (tau^2, Cochran's Q, Q p,
# I^2 as percentage to match the rest of the pipeline).
meta_one_feature <- function(b1, se1, b2, se2) {
  na_row <- tibble::tibble(
    estimate_fe = NA_real_, se_fe = NA_real_,
    z_fe = NA_real_, p_fe = NA_real_,
    ci_lb_fe = NA_real_, ci_ub_fe = NA_real_,
    estimate_re = NA_real_, se_re = NA_real_,
    z_re = NA_real_, p_re = NA_real_,
    ci_lb_re = NA_real_, ci_ub_re = NA_real_,
    tau2  = NA_real_, q_het = NA_real_,
    p_het = NA_real_, i2_pct = NA_real_)

  if (any(!is.finite(c(b1, se1, b2, se2))) || se1 <= 0 || se2 <= 0) {
    return(na_row)
  }

  m <- suppressWarnings(tryCatch(
    meta::metagen(TE       = c(b1, b2),
                  seTE     = c(se1, se2),
                  studlab  = c("SALSA", "WHICAP"),
                  sm       = "MD",
                  common   = TRUE,
                  random   = TRUE,
                  method.tau    = "DL",
                  method.tau.ci = "J",
                  prediction    = FALSE,
                  warn          = FALSE),
    error = function(e) NULL))

  if (is.null(m)) return(na_row)

  # meta returns I2 as a fraction; convert to % to match the rest of the
  # script. CIs already use the 95% default.
  tibble::tibble(
    estimate_fe = m$TE.common,        se_fe = m$seTE.common,
    z_fe        = m$statistic.common, p_fe  = m$pval.common,
    ci_lb_fe    = m$lower.common,     ci_ub_fe = m$upper.common,
    estimate_re = m$TE.random,        se_re = m$seTE.random,
    z_re        = m$statistic.random, p_re  = m$pval.random,
    ci_lb_re    = m$lower.random,     ci_ub_re = m$upper.random,
    tau2        = m$tau2,
    q_het       = m$Q,                p_het = m$pval.Q,
    i2_pct      = m$I2 * 100)
}

# Helpers to reshape one SALSA / WHICAP MWAS table into the columns the
# meta-analysis merge expects.
prep_salsa <- function(mwas_df) {
  mwas_df |>
    dplyr::transmute(salsa_met      = met,
                     estimate_salsa = logFC,
                     se_salsa       = abs(logFC / t),
                     p_salsa        = P.Value,
                     fdr_salsa      = adj.P.Val,
                     vip_salsa      = VIP_comp1)
}

prep_whicap <- function(mwas_df) {
  # `feature` carries the complete mz_rt info even where the separate
  # mz/time/mz_time columns are NA (~47% of the WHICAP MWAS rows). Build
  # whicap_met by stripping the platform prefix so the key matches the
  # `paste0("mz_rt_", mz, "_", time)` form used on the xwalk side.
  mwas_df |>
    dplyr::transmute(whicap_met      = stringr::str_replace(
                                          feature,
                                          "^(hilpos|c18neg)_",
                                          "mz_rt_"),
                     whicap_feature  = feature,
                     estimate_whicap = estimate,
                     se_whicap       = std.error,
                     p_whicap        = p_value,
                     fdr_whicap      = fdr_q_value,
                     whicap_name     = name,
                     whicap_formula  = chemical_formula)
}

# Run meta-analysis: SALSA mixture exposures vs WHICAP, iterated across
# platform → study → population → covariate set → exposure. Strata that
# don't carry the requested exposure (e.g. cross-sectional-only mixtures
# absent from cox results) are silently skipped via purrr::compact().
meta_results <- list(salsa_combined_list, whicap_mwas_by_platform,
                     align_xwalk, names(salsa_combined_list)) |>
  purrr::pmap(function(salsa_results, whicap_mwas, xwalk, platform){
    whicap_df <- prep_whicap(whicap_mwas)

    salsa_results |>
      purrr::imap(function(by_pop, study){
        exposures <- salsa_exposures_by_study[[study]]
        if (is.null(exposures)) return(NULL)

        by_pop |>
          purrr::imap(function(by_covar, population){
            by_covar |>
              purrr::imap(function(exp_list, covar_name){
                exposures |>
                  purrr::set_names() |>
                  purrr::map(function(exposure){
                    if (!exposure %in% names(exp_list)) return(NULL)

                    message(paste0("Meta-analyzing ", exposure,
                                   " for ", platform,
                                   " in ", study, "_", population,
                                   " with covariate set: ",
                                   covar_name, " ..."))

                    salsa_df <- prep_salsa(exp_list[[exposure]])

                    merged <- xwalk |>
                      dplyr::inner_join(salsa_df,  by = "salsa_met") |>
                      dplyr::inner_join(whicap_df, by = "whicap_met")

                    meta <- purrr::pmap(
                      list(merged$estimate_salsa, merged$se_salsa,
                           merged$estimate_whicap, merged$se_whicap),
                      meta_one_feature) |>
                      purrr::list_rbind()

                    dplyr::bind_cols(merged, meta) |>
                      dplyr::mutate(
                        fdr_fe = p.adjust(p_fe, method = "BH"),
                        fdr_re = p.adjust(p_re, method = "BH")) |>
                      dplyr::arrange(p_fe)
                  }) |>
                  purrr::compact()
              })
          })
      }) |>
      purrr::compact()
  }) |>
  purrr::set_names("meta_results_c18", "meta_results_hilic") |>
  list2env(envir = .GlobalEnv)

# Save meta-analysis results -------------------------------------------

save(meta_results_c18, meta_results_hilic,
     file = file.path(align_dir, "meta_mixture_salsa_whicap.RData"))

# Export every (platform × study × population × covar × exposure) df to
# Excel under tables/meta_analysis/.
dir.create(here::here("tables", "meta_analysis"),
           showWarnings = FALSE, recursive = TRUE)

list(c18 = meta_results_c18, hilic = meta_results_hilic) |>
  purrr::iwalk(function(meta_data, platform){
    meta_data |>
      purrr::iwalk(function(by_pop, study){
        by_pop |>
          purrr::iwalk(function(by_covar, population){
            by_covar |>
              purrr::iwalk(function(exp_dfs, covar_name){
                exp_dfs |>
                  purrr::iwalk(function(df, exposure){
                    if (is.null(df) || nrow(df) == 0) return(invisible())
                    writexl::write_xlsx(
                      df,
                      path = here::here("tables", "meta_analysis",
                                  glue::glue("meta_{exposure}_{platform}_",
                                             "{study}_{population}_",
                                             "{covar_name}.xlsx")))
                  })
              })
          })
      })
  })

message("Meta-analysis complete. Aligned-feature counts (all/covar):")
list(c18 = meta_results_c18, hilic = meta_results_hilic) |>
  purrr::iwalk(function(meta_data, platform){
    salsa_exposures_by_study |>
      purrr::iwalk(function(exposures, study){
        exposures |>
          purrr::walk(function(exposure){
            df <- meta_data[[study]][["all"]][["covar"]][[exposure]]
            if (!is.null(df)) {
              message(sprintf("  %s | %s | %s : %d features",
                              platform, study, exposure, nrow(df)))
            }
          })
      })
  })

# Across-feature heterogeneity summary per stratum: median I^2 (IQR),
# fraction of features with p_het < 0.05, n features. Helps gauge how
# noisy the I^2 distribution is for each (platform × study × pop ×
# covar × exposure) cell.
heterogeneity_summary <- list(c18 = meta_results_c18,
                              hilic = meta_results_hilic) |>
  purrr::imap(function(meta_data, platform){
    meta_data |>
      purrr::imap(function(by_pop, study){
        by_pop |>
          purrr::imap(function(by_covar, population){
            by_covar |>
              purrr::imap(function(exp_dfs, covar_name){
                exp_dfs |>
                  purrr::imap(function(df, exposure){
                    if (is.null(df) || nrow(df) == 0) return(NULL)
                    i2  <- df$i2_pct[is.finite(df$i2_pct)]
                    ph  <- df$p_het[is.finite(df$p_het)]
                    tibble::tibble(
                      platform   = platform,
                      study      = study,
                      population = population,
                      covar_set  = covar_name,
                      exposure   = exposure,
                      n_features = nrow(df),
                      i2_median  = if (length(i2)) median(i2) else NA_real_,
                      i2_q25     = if (length(i2)) quantile(i2, 0.25,
                                                            names = FALSE)
                                   else NA_real_,
                      i2_q75     = if (length(i2)) quantile(i2, 0.75,
                                                            names = FALSE)
                                   else NA_real_,
                      i2_zero_pct = if (length(i2))
                                      mean(i2 == 0) * 100 else NA_real_,
                      pct_het    = if (length(ph))
                                     mean(ph < 0.05) * 100 else NA_real_)
                  }) |> purrr::compact() |> purrr::list_rbind()
              }) |> purrr::list_rbind()
          }) |> purrr::list_rbind()
      }) |> purrr::list_rbind()
  }) |>
  purrr::list_rbind()

writexl::write_xlsx(
  heterogeneity_summary,
  path = here::here("tables", "meta_analysis",
                    "heterogeneity_summary.xlsx"))

message("Heterogeneity summary (head):")
print(utils::head(heterogeneity_summary, 12))
message("Full table written to tables/meta_analysis/heterogeneity_summary.xlsx")

# Same summary but restricted to meta-significant features (p_fe < 0.05).
# If `i2_median_sig` is meaningfully higher than the all-feature `i2_median`
# in a given stratum, the top hits in that cell are driven more by one
# cohort than by both — worth flagging in the manuscript.
heterogeneity_summary_sig <- list(c18 = meta_results_c18,
                                  hilic = meta_results_hilic) |>
  purrr::imap(function(meta_data, platform){
    meta_data |>
      purrr::imap(function(by_pop, study){
        by_pop |>
          purrr::imap(function(by_covar, population){
            by_covar |>
              purrr::imap(function(exp_dfs, covar_name){
                exp_dfs |>
                  purrr::imap(function(df, exposure){
                    if (is.null(df) || nrow(df) == 0) return(NULL)
                    sig <- df |>
                      dplyr::filter(is.finite(p_fe), p_fe < 0.05)
                    if (nrow(sig) == 0) return(tibble::tibble(
                      platform   = platform,
                      study      = study,
                      population = population,
                      covar_set  = covar_name,
                      exposure   = exposure,
                      n_sig         = 0L,
                      i2_median_sig = NA_real_,
                      i2_q25_sig    = NA_real_,
                      i2_q75_sig    = NA_real_,
                      i2_zero_pct_sig = NA_real_,
                      pct_het_sig   = NA_real_))
                    i2 <- sig$i2_pct[is.finite(sig$i2_pct)]
                    ph <- sig$p_het[is.finite(sig$p_het)]
                    tibble::tibble(
                      platform   = platform,
                      study      = study,
                      population = population,
                      covar_set  = covar_name,
                      exposure   = exposure,
                      n_sig         = nrow(sig),
                      i2_median_sig = if (length(i2)) median(i2)
                                      else NA_real_,
                      i2_q25_sig    = if (length(i2)) quantile(i2, 0.25,
                                                               names = FALSE)
                                      else NA_real_,
                      i2_q75_sig    = if (length(i2)) quantile(i2, 0.75,
                                                               names = FALSE)
                                      else NA_real_,
                      i2_zero_pct_sig = if (length(i2))
                                          mean(i2 == 0) * 100 else NA_real_,
                      pct_het_sig   = if (length(ph))
                                        mean(ph < 0.05) * 100 else NA_real_)
                  }) |> purrr::compact() |> purrr::list_rbind()
              }) |> purrr::list_rbind()
          }) |> purrr::list_rbind()
      }) |> purrr::list_rbind()
  }) |>
  purrr::list_rbind()

# Side-by-side: all-feature vs meta-significant. The delta column makes the
# "are top hits more heterogeneous than the bulk?" check a one-glance answer.
heterogeneity_summary_combined <- heterogeneity_summary |>
  dplyr::left_join(heterogeneity_summary_sig,
                   by = c("platform", "study", "population",
                          "covar_set", "exposure")) |>
  dplyr::mutate(
    i2_median_delta = i2_median_sig - i2_median,
    pct_het_delta   = pct_het_sig   - pct_het)

writexl::write_xlsx(
  list(all_features = heterogeneity_summary,
       meta_sig     = heterogeneity_summary_sig,
       side_by_side = heterogeneity_summary_combined),
  path = here::here("tables", "meta_analysis",
                    "heterogeneity_summary.xlsx"))

message("Meta-significant (p_fe < 0.05) heterogeneity summary (head):")
print(utils::head(heterogeneity_summary_sig, 12))
message("Combined sheets written to tables/meta_analysis/",
        "heterogeneity_summary.xlsx (all_features | meta_sig | side_by_side)")

# 7. Annotation of meta-analysis features ------------------------------

# Reuse SALSA-side per-feature annotation (xMSannotator stage-5 + Emory
# in-house library, ranked + collapsed by 5-annotation.R) and left-join
# onto every meta-analysis stratum via the SALSA `met` ID.
load(here::here("data", "metabolomics", "annotation",
         "annotation_cleaned_wide.RData"))

load(here::here("data", "metabolomics",
                "alignment", "meta_mixture_salsa_whicap.RData"))

annot_wide_by_platform <- list(c18 = annotation_c18_wide,
                               hil = annotation_hilic_wide)

list(c18 = meta_results_c18, hilic = meta_results_hilic) |>
  purrr::set_names(names(salsa_combined_list)) |>
  purrr::imap(function(meta_data, platform){
    annot_wide <- annot_wide_by_platform[[platform]]
    meta_data |>
      purrr::map(function(by_pop){
        by_pop |>
          purrr::map(function(by_covar){
            by_covar |>
              purrr::map(function(exp_dfs){
                exp_dfs |>
                  purrr::map(function(df){
                    if (is.null(df) || nrow(df) == 0) return(df)
                    df |>
                      dplyr::left_join(annot_wide,
                                       by = c("salsa_met" = "id"))
                  })
              })
          })
      })
  }) |>
  purrr::set_names("meta_annotated_c18", "meta_annotated_hilic") |>
  list2env(envir = .GlobalEnv)

save(meta_annotated_c18, meta_annotated_hilic,
     file = file.path(align_dir, "meta_mixture_annotated.RData"))

# Export annotated meta-analysis tables under tables/meta_analysis/
# {study}/{population}/{covar}/, mirroring tables/mwas_results/ layout.
list(c18 = meta_annotated_c18, hilic = meta_annotated_hilic) |>
  purrr::iwalk(function(meta_data, platform){
    meta_data |>
      purrr::iwalk(function(by_pop, study){
        by_pop |>
          purrr::iwalk(function(by_covar, population){
            by_covar |>
              purrr::iwalk(function(exp_dfs, covar_name){
                exp_dfs |>
                  purrr::iwalk(function(df, exposure){
                    if (is.null(df) || nrow(df) == 0) return(invisible())
                    out_dir <- here::here("tables", "meta_analysis",
                                    study, population, covar_name)
                    dir.create(out_dir, showWarnings = FALSE,
                               recursive = TRUE)
                    writexl::write_xlsx(
                      df,
                      path = file.path(out_dir,
                        glue::glue("meta_{exposure}_{platform}_{study}_",
                                   "{population}_{covar_name}",
                                   "_annotated.xlsx")))
                  })
              })
          })
      })
  })

# 8. Pathway analysis on meta-analysis features ------------------------

source(here::here("scripts", "mummichog_pathway.R"))
source(here::here("scripts", "metapone_pathway.R"))

# Mummichog/Metapone expect: m.z, rt, p.value, t.score, mode. Pull these
# from the meta-analysis table at each stratum; combine C18 + HILIC.
prep_pathway_input <- function(df, mode){
  if (is.null(df) || nrow(df) == 0) return(NULL)
  # Use the fixed-effect pooled p-value + z-score for pathway enrichment
  # (the conventional choice; switch to p_re / z_re if you'd rather drive
  # enrichment off the random-effects estimates).
  df |>
    dplyr::transmute(`m.z`     = whicap_mz,
                     `rt`      = whicap_time,
                     `p.value` = p_fe,
                     `t.score` = z_fe,
                     mode      = mode) |>
    dplyr::filter(!is.na(`m.z`), !is.na(`p.value`)) |>
    dplyr::arrange(`p.value`)
}

# Build one combined (C18 + HILIC) Mummichog input table per stratum ×
# exposure by walking the meta_results structure (study × population ×
# covar × exposure) in lockstep across the two platforms.
pathway_input_list <- meta_results_c18 |>
  purrr::imap(function(by_pop_c18, study){
    by_pop_hil <- meta_results_hilic[[study]]

    by_pop_c18 |>
      purrr::imap(function(by_covar_c18, population){
        by_covar_hil <- by_pop_hil[[population]]
        if (is.null(by_covar_hil)) return(NULL)

        by_covar_c18 |>
          purrr::imap(function(exp_dfs_c18, covar_name){
            exp_dfs_hil <- by_covar_hil[[covar_name]]
            if (is.null(exp_dfs_hil)) return(NULL)

            exposures <- union(names(exp_dfs_c18), names(exp_dfs_hil))
            exposures |>
              purrr::set_names() |>
              purrr::map(function(exposure){
                dplyr::bind_rows(
                  prep_pathway_input(exp_dfs_c18[[exposure]], "negative"),
                  prep_pathway_input(exp_dfs_hil[[exposure]], "positive")
                ) |>
                  dplyr::arrange(`p.value`)
              }) |>
              purrr::compact()
          })
      })
  })

# Write input files under metaboAnalyst/Input/meta/ and Metapone/Input/meta/
pathway_input_list |>
  purrr::iwalk(function(by_pop, study){
    by_pop |>
      purrr::iwalk(function(by_covar, population){
        by_covar |>
          purrr::iwalk(function(exp_inputs, covar_name){
            mum_dir  <- here::here("metaboAnalyst", "Input", "meta",
                             study, population, covar_name)
            mtp_dir  <- here::here("Metapone", "Input", "meta",
                             study, population, covar_name)
            dir.create(mum_dir, showWarnings = FALSE, recursive = TRUE)
            dir.create(mtp_dir, showWarnings = FALSE, recursive = TRUE)

            exp_inputs |>
              purrr::iwalk(function(df, exposure){
                if (is.null(df) || nrow(df) == 0) return(invisible())
                fname <- paste0("meta_", exposure, "_", study, "_",
                                population, "_", covar_name, ".txt")
                purrr::walk(c(file.path(mum_dir, fname),
                              file.path(mtp_dir, fname)),
                            function(out_path){
                              write.table(df, file = out_path,
                                          row.names = FALSE,
                                          col.names = TRUE,
                                          quote = FALSE, sep = "\t")
                            })
              })
          })
      })
  })

message("Pathway input files written to metaboAnalyst/Input/meta/ ",
        "and Metapone/Input/meta/")

# Run Mummichog per stratum × exposure ---------------------------------

message("Running Mummichog pathway analysis on meta-analysis results...")

system.time({
  mummichog_meta_results <- pathway_input_list |>
    purrr::imap(function(by_pop, study){
      by_pop |>
        purrr::imap(function(by_covar, population){
          by_covar |>
            purrr::imap(function(exp_inputs, covar_name){
              exp_inputs |>
                purrr::imap(function(df, exposure){
                  message(paste0("\n--- Mummichog (meta): ", exposure,
                                 " | ", study, "_", population,
                                 " | ", covar_name, " ---"))

                  input_file <- here::here("metaboAnalyst", "Input", "meta",
                                     study, population, covar_name,
                                     paste0("meta_", exposure, "_",
                                            study, "_", population, "_",
                                            covar_name, ".txt"))
                  output_dir <- here::here("metaboAnalyst", "Output", "meta",
                                     study, population, covar_name,
                                     exposure)
                  dir.create(output_dir, showWarnings = FALSE,
                             recursive = TRUE)

                  result <- tryCatch(
                    run_mummichog(
                      input_file = input_file,
                      output_dir = output_dir,
                      p_cutoff = 0.1,
                      organism = "hsa_mfn",
                      instrument_ppm = 10.0,
                      ion_mode = "mixed",
                      adducts = c("M-H [1-]", "M-2H [2-]",
                                  "M-H2O-H [1-]", "M [1+]",
                                  "M+H [1+]", "M+Na [1+]"),
                      min_hits = 3,
                      num_permutations = 100
                    ),
                    error = function(e){
                      warning(paste0("Mummichog (meta) failed for ",
                                     exposure, " | ", study, "_",
                                     population, " | ", covar_name,
                                     ": ", e$message))
                      NULL
                    })
                  mum_rdata <- file.path(output_dir, "mum.RData")
                  if (file.exists(mum_rdata)) file.remove(mum_rdata)
                  result
                })
            })
        })
    })
})

save(mummichog_meta_results,
     file = here::here("data", "metabolomics", "results",
                 "mummichog_meta_results.RData"))

# Run Metapone per stratum × exposure ----------------------------------

message("Running Metapone pathway analysis on meta-analysis results...")

system.time({
  metapone_meta_results <- pathway_input_list |>
    purrr::imap(function(by_pop, study){
      by_pop |>
        purrr::imap(function(by_covar, population){
          by_covar |>
            purrr::imap(function(exp_inputs, covar_name){
              exp_inputs |>
                purrr::imap(function(df, exposure){
                  message(paste0("\n--- Metapone (meta): ", exposure,
                                 " | ", study, "_", population,
                                 " | ", covar_name, " ---"))

                  input_file <- here::here("Metapone", "Input", "meta",
                                     study, population, covar_name,
                                     paste0("meta_", exposure, "_",
                                            study, "_", population, "_",
                                            covar_name, ".txt"))

                  tryCatch(
                    run_metapone(
                      input_file = input_file,
                      p_cutoff = 0.05,
                      num_permutations = 200,
                      match_tol_ppm = 10,
                      pos.adductlist = c("M+H", "M+Na", "M+"),
                      neg.adductlist = c("M-H", "M-2H", "M-H2O-H")
                    ),
                    error = function(e){
                      warning(paste0("Metapone (meta) failed for ",
                                     exposure, " | ", study, "_",
                                     population, " | ", covar_name,
                                     ": ", e$message))
                      NULL
                    })
                })
            })
        })
    })
})

save(metapone_meta_results,
     file = here::here("data", "metabolomics", "results",
                 "metapone_meta_results.RData"))

# Export pathway result tables + bubble plots per stratum × exposure
list(mummichog = mummichog_meta_results,
     metapone  = metapone_meta_results) |>
  purrr::iwalk(function(method_results, method){
    method_results |>
      purrr::iwalk(function(by_pop, study){
        by_pop |>
          purrr::iwalk(function(by_covar, population){
            by_covar |>
              purrr::iwalk(function(exp_results, covar_name){
                tbl_dir <- here::here("tables", paste0(method, "_meta_results"),
                                study, population, covar_name)
                fig_dir <- here::here("figures", paste0(method, "_meta"),
                                study, population, covar_name)
                dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
                dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

                exp_results |>
                  purrr::compact() |>
                  purrr::iwalk(function(res, exposure){
                    if (!is.null(res$result_table)){
                      writexl::write_xlsx(
                        res$result_table,
                        path = file.path(tbl_dir,
                          glue::glue("{method}_meta_{exposure}_{study}_",
                                     "{population}_{covar_name}.xlsx")))
                    }
                    if (!is.null(res$plot)){
                      ggplot2::ggsave(
                        filename = file.path(fig_dir,
                          glue::glue("{method}_meta_{exposure}_{study}_",
                                     "{population}_{covar_name}.png")),
                        plot = res$plot +
                          ggplot2::ggtitle(paste0(method, " (meta): ",
                                                  exposure)),
                        width = 10, height = 8, dpi = 300)
                    }
                  })
              })
          })
      })
  })

message("Pathway analysis on meta-analysis complete.")

# 9. Visualization of meta-analysis features ---------------------------

pacman::p_load(ggrepel, ggpubr, patchwork)

# Volcano plot: pooled estimate vs -log10(meta p), faceted by platform.
# `model` picks which set of meta::metagen columns to plot — "fe" → estimate_fe
# / p_fe / fdr_fe, "re" → estimate_re / p_re / fdr_re.
create_meta_volcano <- function(meta_c18, meta_hil, exposure,
                                model = c("fe", "re"),
                                fdr_threshold = 0.05, n_labels = 6){
  model    <- match.arg(model)
  est_col  <- paste0("estimate_", model)
  p_col    <- paste0("p_",  model)
  fdr_col  <- paste0("fdr_", model)
  model_lab <- if (model == "fe") "fixed-effects" else "random-effects"

  prep <- function(df, platform_label){
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      dplyr::mutate(
        platform = platform_label,
        estimate_plot = .data[[est_col]],
        p_plot        = .data[[p_col]],
        fdr_plot      = .data[[fdr_col]],
        neg_log10_p   = -log10(p_plot),
        significant = dplyr::case_when(
          fdr_plot < fdr_threshold ~ "FDR < 0.05",
          p_plot   < 0.05          ~ "P < 0.05",
          TRUE                     ~ "NS"),
        significant = factor(significant,
                             levels = c("FDR < 0.05", "P < 0.05", "NS")),
        compound = dplyr::if_else(is.na(compound) | compound == "",
                                  NA_character_, compound))
  }

  plot_data <- dplyr::bind_rows(prep(meta_c18,  "C18/neg−"),
                                prep(meta_hil,  "HILIC/pos+"))
  if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)

  top_labs <- plot_data |>
    dplyr::filter(!is.na(compound), fdr_plot < fdr_threshold) |>
    dplyr::arrange(p_plot) |>
    dplyr::group_by(platform) |>
    dplyr::slice_head(n = n_labels) |>
    dplyr::ungroup() |>
    dplyr::mutate(compound_first =
                    stringr::str_squish(stringr::str_split_i(compound, ";", 1)))

  ggplot(plot_data, aes(x = estimate_plot, y = neg_log10_p)) +
    geom_point(data = ~ dplyr::filter(.x, significant == "NS"),
               aes(shape = platform),
               color = "grey70", alpha = 0.5, size = 1.6) +
    geom_point(data = ~ dplyr::filter(.x, significant != "NS"),
               aes(color = significant, shape = platform),
               alpha = 0.7, size = 2.4) +
    ggrepel::geom_label_repel(
      data = top_labs,
      aes(label = compound_first),
      size = 3.6, max.overlaps = Inf, force = 2,
      box.padding = 0.4, segment.color = "grey50",
      inherit.aes = TRUE, show.legend = FALSE) +
    scale_color_manual(values = c("FDR < 0.05" = "#BE3F42",
                                  "P < 0.05"   = "#DE9960"),
                       name = "Significance", drop = FALSE) +
    scale_shape_manual(values = c("C18/neg−" = 16, "HILIC/pos+" = 17),
                       name = "Column") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40") +
    labs(title = gsub("comp_|exp_", "",
                      paste0("Meta-analysis (", model_lab, "): ", exposure)),
         x = paste0("Pooled effect estimate (", toupper(model), ")"),
         y = bquote(-log[10](P[.(toupper(model))]))) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = .5),
          axis.title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold", size = 12),
          legend.position   = "bottom",
          legend.box        = "horizontal",
          legend.background = element_rect(colour = "grey80", fill = "white",
                                           linewidth = 0.5),
          legend.margin     = margin(4, 6, 4, 6),
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10))
}

# Scatter plot: SALSA estimate vs WHICAP estimate at meta-significant
# (p_FE < 0.05) aligned features. Both platforms are overlaid in the same
# panel and distinguished by point shape; the three mixture indices
# (comp_wqs_all, comp_qgcomp_all, comp_qgcomp_cox_all) are arranged as
# side-by-side facets. Points are coloured by cross-cohort agreement
# (same sign + both p < 0.05).
#
# `meta_inputs` is an ordered, named list — one element per mixture —
# each of the form list(c18 = <df>, hil = <df>, label = "<facet title>").
create_cohort_scatter <- function(meta_inputs){
  prep <- function(df, platform_label, exposure_label){
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      dplyr::filter(!is.na(p_fe), p_fe < 0.05) |>
      dplyr::mutate(
        platform = platform_label,
        exposure_label = exposure_label,
        concordant = dplyr::case_when(
          p_salsa < 0.05 & p_whicap < 0.05 &
            sign(estimate_salsa) == sign(estimate_whicap) ~
              "Both P < 0.05 (same direction)",
          p_salsa < 0.05 & p_whicap < 0.05 ~
              "Both P < 0.05 (opposite direction)",
          p_whicap < 0.05                  ~ "P < 0.05 in WHICAP",
          p_salsa  < 0.05                  ~ "P < 0.05 in SALSA",
          TRUE                             ~ "NS"),
        concordant = factor(concordant,
                            levels = c("Both P < 0.05 (same direction)",
                                       "Both P < 0.05 (opposite direction)",
                                       "P < 0.05 in WHICAP",
                                       "P < 0.05 in SALSA",
                                       "NS")))
  }

  plot_data <- meta_inputs |>
    purrr::map(function(mi){
      dplyr::bind_rows(prep(mi$c18, "C18/neg−",  mi$label),
                       prep(mi$hil, "HILIC/pos+", mi$label))
    }) |>
    dplyr::bind_rows()

  if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)

  exposure_levels <- purrr::map_chr(meta_inputs, "label")
  plot_data <- plot_data |>
    dplyr::mutate(exposure_label = factor(exposure_label,
                                          levels = exposure_levels))

  ggplot(plot_data, aes(x = estimate_salsa, y = estimate_whicap)) +
    geom_point(data = ~ dplyr::filter(.x, concordant == "NS"),
               aes(shape = platform),
               color = "grey70", alpha = 0.4, size = 1.5) +
    geom_point(data = ~ dplyr::filter(.x, concordant != "NS"),
               aes(color = concordant, shape = platform),
               alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linewidth = 0.8, alpha = 0.2) +
    ggpubr::stat_cor(method = "pearson",
                     label.x.npc = "left", label.y.npc = "top",
                     size = 5.5, fontface = "italic") +
    scale_color_manual(values = c("Both P < 0.05 (same direction)"     = "#B73F42",
                                  "Both P < 0.05 (opposite direction)" = "#436C85",
                                  "P < 0.05 in WHICAP"                 = "#7E9A6C",
                                  "P < 0.05 in SALSA"                  = "#DE9960"),
                       name = "Cross-cohort", drop = FALSE) +
    scale_shape_manual(values = c("C18/neg−" = 16, "HILIC/pos+" = 17),
                       name = "Column") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    facet_wrap(~ exposure_label, scales = "free") +
    labs(title = "SALSA vs WHICAP: meta-significant features (p_FE < 0.05)",
         x = "SALSA effect estimate (mixture)",
         y = "WHICAP effect estimate (PM2.5)") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = .5),
          axis.title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold", size = 12),
          legend.position   = "bottom",
          legend.box        = "horizontal",
          legend.background = element_rect(colour = "grey80", fill = "white",
                                           linewidth = 0.5),
          legend.margin     = margin(4, 6, 4, 6),
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10))
}

# Manhattan-style plot of meta-analysis results, faceted by platform.
# Mirrors `create_manhattan` in scripts/7-visualization.R: m/z on x,
# -log10(p_fe) on y, point color by FDR_FE / P_FE significance, with a
# top-right feature-count box and a boxed bottom legend. Meta-analysis
# has no VIP equivalent, so the script-7 "P < 0.05 & VIP > 2" tier is
# dropped — just FDR_FE < 0.05, P_FE < 0.05, NS.
create_meta_manhattan <- function(meta_c18, meta_hil, exposure){

  prep <- function(df, platform_label){
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      dplyr::filter(!is.na(p_fe), !is.na(whicap_mz)) |>
      dplyr::mutate(
        platform    = platform_label,
        mz          = whicap_mz,
        neg_log10_p = -log10(p_fe),
        significant = dplyr::case_when(
          fdr_fe < 0.05 ~ "FDR < 0.05",
          p_fe   < 0.05 ~ "P < 0.05",
          TRUE          ~ "NS"),
        significant = factor(significant,
                             levels = c("FDR < 0.05", "P < 0.05", "NS")))
  }

  plot_data <- dplyr::bind_rows(prep(meta_c18, "C18/neg−"),
                                prep(meta_hil, "HILIC/pos+"))
  if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)

  count_labels <- plot_data |>
    dplyr::group_by(platform) |>
    dplyr::summarize(n_fdr = sum(fdr_fe < 0.05, na.rm = TRUE),
                     n_nom = sum(p_fe   < 0.05, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(label = paste0("FDR < 0.05: ", n_fdr,
                                 "\nP < 0.05: ", n_nom),
                  mz          = Inf,
                  neg_log10_p = Inf)

  ggplot(plot_data, aes(x = mz, y = neg_log10_p)) +
    geom_point(data = \(d) dplyr::filter(d, significant == "NS"),
               color = "grey70", alpha = 0.4, size = 1.2) +
    geom_point(data = \(d) dplyr::filter(d, significant != "NS"),
               aes(color = significant), alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c("FDR < 0.05" = "#B73F42",
                                  "P < 0.05"   = "#DE9960"),
                       name = "Significance", drop = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_text(data = count_labels,
              aes(x = mz, y = neg_log10_p, label = label),
              hjust = 1.05, vjust = 1.2,
              size = 3.5, fontface = "bold",
              inherit.aes = FALSE) +
    facet_wrap(~ platform, scales = "free_x") +
    labs(title = gsub("exp_|comp_", "", paste0("Meta-analysis: ", exposure)),
         x = expression(bold("Mass-to-charge ratio (" * italic(m/z) * ")")),
         y = expression(-log[10](P[FE]))) +
    theme_classic() +
    theme(plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.title       = element_text(face = "bold", size = 12),
          axis.text        = element_text(size = 10),
          strip.text       = element_text(face = "bold", size = 12),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.position   = "bottom",
          legend.box        = "horizontal",
          legend.background = element_rect(colour = "grey80", fill = "white",
                                           linewidth = 0.5),
          legend.margin     = margin(4, 6, 4, 6),
          legend.title      = element_text(face = "bold", size = 12),
          legend.text       = element_text(size = 11),
          legend.key.size   = unit(0.5, "cm"))
}

# Per-exposure volcanos (FE and RE) — saved under
# figures/meta_analysis/{study}/{population}/{covar}/{exposure}/.
purrr::iwalk(meta_annotated_c18, function(by_pop, study){
  by_pop_hil <- meta_annotated_hilic[[study]]
  by_pop |>
    purrr::iwalk(function(by_covar, population){
      by_covar_hil <- by_pop_hil[[population]]
      by_covar |>
        purrr::iwalk(function(exp_dfs, covar_name){
          exp_dfs_hil <- by_covar_hil[[covar_name]]
          exp_dfs |>
            purrr::iwalk(function(df_c18, exposure){
              df_hil <- exp_dfs_hil[[exposure]]
              fig_dir <- here::here("figures", "meta_analysis",
                              study, population, covar_name, exposure)
              dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

              p_volcano_fe <- create_meta_volcano(df_c18, df_hil, exposure,
                                                  model = "fe")
              if (!is.null(p_volcano_fe)){
                ggplot2::ggsave(
                  filename = file.path(fig_dir,
                    glue::glue("volcano_meta_fe_{exposure}.png")),
                  plot = p_volcano_fe, width = 12, height = 6, dpi = 300)
              }

              p_volcano_re <- create_meta_volcano(df_c18, df_hil, exposure,
                                                  model = "re")
              if (!is.null(p_volcano_re)){
                ggplot2::ggsave(
                  filename = file.path(fig_dir,
                    glue::glue("volcano_meta_re_{exposure}.png")),
                  plot = p_volcano_re, width = 12, height = 6, dpi = 300)
              }

              p_manhattan <- create_meta_manhattan(df_c18, df_hil, exposure)
              if (!is.null(p_manhattan)){
                ggplot2::ggsave(
                  filename = file.path(fig_dir,
                    glue::glue("manhattan_meta_{exposure}.png")),
                  plot = p_manhattan, width = 14, height = 6, dpi = 300)
              }
            })
        })
    })
})

# Combined-mixture scatter — one figure per (population × covariate set)
# spanning the three mixture indices (total/comp_wqs_all,
# total/comp_qgcomp_all, cox/comp_qgcomp_cox_all) as facets. Both
# platforms appear in the same panel, distinguished by point shape.
mixture_scatter_specs <- list(
  list(study = "total", exposure = "comp_wqs_all",
       label = "WQS (cross-sectional)"),
  list(study = "total", exposure = "comp_qgcomp_all",
       label = "QGcomp (cross-sectional)"),
  list(study = "cox",   exposure = "comp_qgcomp_cox_all",
       label = "QGcomp (Cox)")
)

populations <- meta_annotated_c18 |>
  purrr::map(names) |>
  unlist() |>
  unique()

purrr::walk(populations, function(population){
  covar_names <- meta_annotated_c18 |>
    purrr::map(~ names(.x[[population]])) |>
    unlist() |>
    unique()

  purrr::walk(covar_names, function(covar_name){
    meta_inputs <- mixture_scatter_specs |>
      purrr::map(function(spec){
        df_c18 <- meta_annotated_c18[[spec$study]][[population]][[covar_name]][[spec$exposure]]
        df_hil <- meta_annotated_hilic[[spec$study]][[population]][[covar_name]][[spec$exposure]]
        if (is.null(df_c18) && is.null(df_hil)) return(NULL)
        list(c18 = df_c18, hil = df_hil, label = spec$label)
      }) |>
      purrr::compact()

    if (length(meta_inputs) == 0) return(invisible())

    p_scatter <- create_cohort_scatter(meta_inputs)
    if (is.null(p_scatter)) return(invisible())

    fig_dir <- here::here("figures", "meta_analysis", "combined",
                          population, covar_name)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(
      filename = file.path(fig_dir,
        glue::glue("scatter_salsa_whicap_mixtures_",
                   "{population}_{covar_name}.png")),
      plot = p_scatter, width = 16, height = 7, dpi = 300)

    # Combined-mixture Manhattan panel — 3 stacked Manhattan rows
    # (WQS / QGcomp / QGcomp-Cox) sharing the bottom legend, mirroring
    # the script-7 manhattan_composite_panel structure.
    no_title <- theme(plot.title = element_blank())
    no_xaxis <- theme(axis.title.x = element_blank(),
                      axis.text.x  = element_blank(),
                      axis.ticks.x = element_blank())

    manhattan_rows <- purrr::imap(meta_inputs, function(mi, i){
      p <- create_meta_manhattan(mi$c18, mi$hil, exposure = mi$label)
      if (is.null(p)) return(NULL)
      tag_letter <- LETTERS[i]
      if (i < length(meta_inputs)) {
        p + labs(title = NULL, tag = tag_letter) + no_title + no_xaxis
      } else {
        p + labs(title = NULL, tag = tag_letter) + no_title
      }
    }) |> purrr::compact()

    if (length(manhattan_rows) > 0){
      row_labels <- purrr::map(meta_inputs, function(mi){
        patchwork::wrap_elements(
          grid::textGrob(mi$label, x = 0.02, hjust = 0,
                         gp = grid::gpar(fontface = "bold", fontsize = 12)))
      })

      interleaved <- purrr::map2(row_labels, manhattan_rows,
                                 ~ list(.x, .y)) |>
        purrr::flatten()

      combined_manhattan <- purrr::reduce(interleaved, `/`) +
        patchwork::plot_layout(
          heights = rep(c(0.03, 1), length(manhattan_rows)),
          guides  = "collect") &
        theme(legend.position   = "bottom",
              legend.box        = "horizontal",
              legend.background = element_rect(colour = "grey80",
                                               fill = "white",
                                               linewidth = 0.5),
              legend.margin     = margin(4, 6, 4, 6),
              legend.title      = element_text(face = "bold", size = 12),
              legend.text       = element_text(size = 11),
              legend.key.size   = unit(0.5, "cm"),
              plot.tag          = element_text(face = "bold", size = 14))

      ggplot2::ggsave(
        filename = file.path(fig_dir,
          glue::glue("manhattan_meta_mixtures_",
                     "{population}_{covar_name}.png")),
        plot = combined_manhattan,
        width = 14, height = 6 * length(manhattan_rows), dpi = 300)
    }

    # Combined volcano + scatter panel — 3 rows (one per mixture),
    # each row pairs the FE volcano (left, faceted by platform) with
    # the per-mixture concordance scatter (right, both platforms via
    # shape). Mirrors the layout of `create_combined_panel` in
    # scripts/7-visualization.R.
    combo_rows <- purrr::imap(meta_inputs, function(mi, i){
      p_v <- create_meta_volcano(mi$c18, mi$hil, exposure = mi$label,
                                  model = "fe")
      p_s <- create_cohort_scatter(list(list(c18 = mi$c18, hil = mi$hil,
                                              label = mi$label)))
      if (is.null(p_v) || is.null(p_s)) return(NULL)
      tag_v <- LETTERS[2 * i - 1]
      tag_s <- LETTERS[2 * i]
      p_v <- p_v + labs(title = NULL, tag = tag_v)
      p_s <- p_s + labs(title = NULL, tag = tag_s) +
        guides(shape = "none") +
        theme(strip.text = element_blank(),
              strip.background = element_blank())
      p_v | p_s
    }) |> purrr::compact()

    if (length(combo_rows) > 0){
      combo_labels <- purrr::map(meta_inputs, function(mi){
        patchwork::wrap_elements(
          grid::textGrob(mi$label, x = 0.02, hjust = 0,
                         gp = grid::gpar(fontface = "bold", fontsize = 13)))
      })

      combo_interleaved <- purrr::map2(combo_labels, combo_rows,
                                       ~ list(.x, .y)) |>
        purrr::flatten()

      combined_panel <- purrr::reduce(combo_interleaved, `/`) +
        patchwork::plot_layout(
          heights = rep(c(0.03, 1), length(combo_rows)),
          guides  = "collect") &
        theme(legend.position   = "bottom",
              legend.box        = "horizontal",
              legend.background = element_rect(colour = "grey80",
                                               fill = "white",
                                               linewidth = 0.5),
              legend.margin     = margin(4, 6, 4, 6),
              legend.title      = element_text(face = "bold", size = 11),
              legend.text       = element_text(size = 10),
              legend.key.size   = unit(0.5, "cm"),
              plot.tag          = element_text(face = "bold", size = 13))

      ggplot2::ggsave(
        filename = file.path(fig_dir,
          glue::glue("combined_meta_panel_",
                     "{population}_{covar_name}.png")),
        plot = combined_panel,
        width = 18, height = 6 * length(combo_rows), dpi = 300)
    }
  })
})

message("Meta-analysis visualizations saved under figures/meta_analysis/")

# 10. Bubble heatmap of meta-analysis mummichog pathway results --------
#
# Mirrors `create_pathway_bubble_heatmap` from scripts/pathway_viz_augment.R,
# but reads the meta-analysis mummichog outputs under
# `metaboAnalyst/Output/meta/{study}/{population}/{covar}/{exposure}/
# mummichog_pathway_enrichment_mummichog.csv` and faceting categories ×
# population with the three mixture indices on the x-axis.

pacman::p_load(ggthemes, ggh4x, legendry, ggnewscale, ggtext)

# Pathway category lookup tables (copied verbatim from pathway_viz_augment.R
# so this section is self-contained and doesn't trigger its side effects).
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
                             "Chondroitin", "Heparan", "Nucleotide sugar")
lipid_metabolism <- c("Bile", "Fatty", "lipid", "Phytanic",
                      "Cholesterol", "neuroprostanes", "steroid",
                      "Sphingolipid", "Glycerophospholipid", "Phospholipid",
                      "Triglyceride", "Ceramide")
energy_metabolism <- c("Butanoate", "Carnitine", "Glycolysis",
                       "Pyruvate", "Pentose", "octadecatrienoate",
                       "TCA", "Citrate", "Oxidative", "Glyoxylate",
                       "Propanoate", "Carbon fixation")
inflammation_metabolism <- c("Arachidonic", "Leukotriene", "Prostaglandin",
                             "linoleic", "Linoleate", "Eicosanoid")
vitamin_cofactor_metabolism <- c("Vitamin", "Biopterin", "Lipoate",
                                 "Porphyrin", "Catabolism",
                                 "Caffeine", "Folate", "Riboflavin",
                                 "Thiamine", "Biotin", "Pantothenate")
nucleotide_metabolism <- c("Pyrimidine", "Purine")
signaling_metabolism  <- c("Dynorphin", "Dopamine", "Serotonin",
                           "Catecholamine", "Neurotransmitter")
secondary_metabolite_metabolism <- c("Alkaloid")
xenobiotic_metabolism <- c("Xenobiotic", "Drug", "Benzoate")

categorize_pathway <- function(pathway_name) {
  dplyr::case_when(
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(amino_acid_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "amino acid",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(carbohydrate_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "carbohydrate",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(lipid_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "lipid",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(energy_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "energy",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(inflammation_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "inflammation",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(vitamin_cofactor_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "vitamin/cofactor",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(nucleotide_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "nucleotide",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(signaling_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "signaling",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(secondary_metabolite_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "secondary",
    stringr::str_detect(pathway_name, stringr::regex(
      stringr::str_c(xenobiotic_metabolism, collapse = "|"),
      ignore_case = TRUE)) ~ "xenobiotic",
    TRUE ~ "other"
  )
}

# Read one mummichog meta CSV, extracting (study × population × covar ×
# exposure) from the directory chain `.../meta/{study}/{pop}/{covar}/{exp}/`.
read_meta_mummichog <- function(file_path) {
  parts <- strsplit(file_path, "/", fixed = TRUE)[[1]]
  meta_idx <- which(parts == "meta")
  if (length(meta_idx) == 0L) return(NULL)
  study      <- parts[meta_idx + 1]
  population <- parts[meta_idx + 2]
  covar_set  <- parts[meta_idx + 3]
  exposure   <- parts[meta_idx + 4]

  df <- readr::read_csv(file_path, show_col_types = FALSE) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(
      pathway_name = any_of(c("...1", "pathway", "pathway_name", "name")),
      pathway_size = any_of(c("pathway total", "pathway_total", "total")),
      hits_total   = any_of(c("hits.total", "hits_total", "total_hits")),
      hits_sig     = any_of(c("hits.sig", "hits_sig", "sig_hits")),
      expected     = any_of(c("expected")),
      p_value      = any_of(c("p(fisher)", "p.value", "pvalue", "p_value"))
    ) |>
    dplyr::mutate(
      study      = study,
      population = population,
      covar_set  = covar_set,
      exposure   = exposure,
      category   = categorize_pathway(pathway_name)
    )
}

mum_meta_files <- list.files(
  here::here("metaboAnalyst", "Output", "meta"),
  pattern = "mummichog_pathway_enrichment.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) |>
  purrr::keep(~ grepl("comp_qgcomp_all|comp_wqs_all|comp_qgcomp_cox_all", .x))

mum_meta_all <- mum_meta_files |>
  purrr::map(read_meta_mummichog) |>
  purrr::compact() |>
  purrr::list_rbind()

mum_meta_sig <- mum_meta_all |>
  dplyr::filter(!is.na(p_value), p_value < 0.1) |>
  dplyr::arrange(p_value) |>
  dplyr::mutate(
    enrichment_factor = hits_sig / expected,
    exposure_short = dplyr::case_when(
      exposure == "comp_wqs_all"        & study == "total" ~ "WQS",
      exposure == "comp_qgcomp_all"     & study == "total" ~ "QGcomp",
      exposure == "comp_qgcomp_cox_all" & study == "cox"   ~ "QGcomp-Cox",
      TRUE ~ exposure
    ),
    population = factor(population,
                        levels = c("all", "no demcind", "demcind"))
  ) |>
  dplyr::filter(is.finite(enrichment_factor))

# Bubble heatmap: tile shading by row stripe, point size = enrichment
# factor, point fill = -log10(p). Faceted category × population, x-axis
# is the three mixture indices grouped by Cross-sectional vs Time-to-event.
create_meta_pathway_bubble_heatmap <- function(data) {

  if (nrow(data) == 0) return(NULL)

  pathway_summary <- data |>
    dplyr::group_by(pathway_name, category) |>
    dplyr::summarise(
      min_p       = min(p_value, na.rm = TRUE),
      total_score = sum(-log10(p_value), na.rm = TRUE),
      .groups = "drop")

  category_order <- pathway_summary |>
    dplyr::group_by(category) |>
    dplyr::summarise(cat_min_p = min(min_p),
                     cat_total = sum(total_score),
                     .groups = "drop") |>
    dplyr::arrange(cat_min_p, dplyr::desc(cat_total)) |>
    dplyr::pull(category)

  pathway_sorted <- pathway_summary |>
    dplyr::mutate(category = factor(category, levels = category_order)) |>
    dplyr::arrange(category, dplyr::desc(min_p), total_score) |>
    dplyr::mutate(pathway_name = forcats::fct_inorder(pathway_name))

  exposure_levels <- c("WQS", "QGcomp", "QGcomp-Cox")
  exposure_display_labels <- c("WQS" = "WQS",
                               "QGcomp" = "QGcomp",
                               "QGcomp-Cox" = "QGcomp")

  plot_data <- data |>
    dplyr::mutate(
      category       = factor(category, levels = category_order),
      pathway_name   = factor(pathway_name,
                              levels = levels(pathway_sorted$pathway_name)),
      neg_log10_p    = -log10(p_value),
      exposure_short = factor(exposure_short, levels = exposure_levels))

  cats <- category_order
  pal  <- ggthemes::tableau_color_pal("Tableau 10")(length(cats))
  names(pal) <- cats
  strip_colors <- lapply(pal[cats], function(col) {
    element_rect(fill = scales::alpha(col, 0.3), colour = "grey80")
  })

  all_combos <- tidyr::expand_grid(
    pathway_name   = levels(pathway_sorted$pathway_name),
    exposure_short = factor(exposure_levels, levels = exposure_levels)) |>
    dplyr::mutate(pathway_name = factor(
      pathway_name, levels = levels(pathway_sorted$pathway_name))) |>
    dplyr::left_join(
      pathway_sorted |> dplyr::select(pathway_name, category),
      by = "pathway_name")

  bubble_data <- all_combos |>
    dplyr::left_join(
      plot_data |>
        dplyr::select(pathway_name, exposure_short,
                      p_value, neg_log10_p, enrichment_factor),
      by = c("pathway_name", "exposure_short"))

  shade_df <- pathway_sorted |>
    dplyr::arrange(category, dplyr::desc(as.integer(pathway_name))) |>
    dplyr::mutate(visual_idx = dplyr::row_number(),
                  shade = ifelse(visual_idx %% 2 == 0, "even", "odd")) |>
    dplyr::select(pathway_name, shade)

  shade_data <- bubble_data |>
    dplyr::left_join(shade_df, by = "pathway_name")

  ggplot(shade_data, aes(x = exposure_short, y = pathway_name)) +
    geom_tile(aes(fill = shade), width = 1, height = 1, alpha = 0.4,
              show.legend = FALSE) +
    scale_fill_manual(values = c("even" = "grey93", "odd" = "white"),
                      guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_point(data = ~ dplyr::filter(.x, !is.na(p_value)),
               aes(size = enrichment_factor, fill = neg_log10_p),
               shape = 21, color = "grey20", stroke = 0.4, alpha = 0.9) +
    scale_fill_gradientn(
      colours = c("#FFF7BC", "#FEC44F", "#F46D43", "#D73027", "#A50026"),
      name    = expression(-log[10](italic(p))),
      limits  = c(-log10(0.1), NA),
      na.value = "grey85") +
    scale_size_continuous(range = c(2, 9), name = "Enrichment factor") +
    scale_x_discrete(
      labels = exposure_display_labels,
      guide  = legendry::guide_axis_nested(
        key = legendry::key_range_manual(
          start = c("WQS", "QGcomp-Cox"),
          end   = c("QGcomp", "QGcomp-Cox"),
          name  = c("Cross-sectional", "Time-to-event")))) +
    guides(
      fill = guide_colorbar(barwidth = 10, barheight = 0.8,
                            title.position = "top", title.hjust = 0.5,
                            frame.colour = "grey40",
                            ticks.colour = "grey40"),
      size = guide_legend(title.position = "top", title.hjust = 0.5,
                          override.aes = list(fill = "#F46D43",
                                              alpha = 0.8))) +
    ggh4x::facet_grid2(
      category ~ .,
      scales = "free_y", space = "free_y",
      strip  = ggh4x::strip_themed(
        background_y = strip_colors)) +
    labs(y = NULL, x = NULL, title = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.margin   = margin(t = 8),
      legend.spacing.x = unit(1, "cm"),
      legend.title    = element_text(face = "bold", size = 13),
      legend.text     = element_text(size = 12),
      axis.line       = element_blank(),
      panel.border    = element_rect(colour = "grey80", fill = NA,
                                     linewidth = 0.4),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      axis.text.y     = element_text(size = 13, color = "grey20"),
      axis.text.x     = ggtext::element_textbox_simple(
        face = "bold", size = 13, color = "grey20",
        halign = 0.5,
        padding = margin(3, 6, 3, 6),
        margin  = margin(t = 3, b = 3),
        box.color = "grey40",
        linewidth = 0.5,
        linetype  = 1,
        fill = "grey95",
        r = unit(2, "pt")),
      strip.text.y    = element_text(face = "bold.italic", size = 13,
                                     angle = 0, color = "grey20"),
      panel.spacing   = unit(0.3, "lines"),
      plot.margin     = margin(10, 15, 10, 10))
}

# Render + save one bubble heatmap per covariate set (primary + sensitivity).
mum_meta_fig_dir <- here::here("figures", "meta_analysis", "pathway_bubble")
dir.create(mum_meta_fig_dir, showWarnings = FALSE, recursive = TRUE)

c("covar", "covar_sen") |>
  purrr::walk(function(cs){
    # Restrict the meta-analysis pathway plot to the "all" population.
    # Subgroup meta-pools become WHICAP-dominated as SALSA SEs grow with
    # smaller n (e.g., demcind n=141 → SALSA SE ~2.6× larger), so the
    # cognitive-subgroup pathway columns reflect WHICAP-PM2.5 biology
    # rather than SALSA-specific signal and should not be reported.
    d <- mum_meta_sig |>
      dplyr::filter(covar_set == cs, population == "all")
    p <- create_meta_pathway_bubble_heatmap(d)
    if (is.null(p)) return(invisible())
    ggsave(
      filename = file.path(mum_meta_fig_dir,
                           glue::glue("mummichog_meta_bubble_heatmap_{cs}.png")),
      plot   = p,
      width  = 9,
      height = max(8, dplyr::n_distinct(d$pathway_name) * 0.35),
      dpi    = 300)
  })

message("Mummichog meta bubble heatmaps saved under ",
        "figures/meta_analysis/pathway_bubble/")

stopCluster(cl)

# ----------------------------------------------------------------------
# DATA REQUEST — what is still needed from the WHICAP team
# ----------------------------------------------------------------------
# The feature-level summary files now in hand (mz, time, mean, sd per
# feature, per platform) are sufficient to run the apLCMS alignment above
# with real intensity/SD weighting on both cohorts. They are NOT, however,
# enough for a per-sample (individual-participant-data) meta-analysis. To
# go beyond summary-statistic meta-analysis we would still need from the
# WHICAP team:
#
#   1. ComBat-corrected (or otherwise normalized) feature tables, one per
#      platform (C18-neg, HILIC-pos), in the same layout as the SALSA files
#      under data/metabolomics/{c18neg,hilicpos}/ComBat_*.txt — i.e. columns
#      for mz, time, mz.min, mz.max, QC scores, and ONE COLUMN PER SAMPLE
#      with peak-area intensities.
#
#   2. The corresponding sample manifest / ID-mapping file linking each
#      LC-MS file name to a participant ID, plus a phenotype file with at
#      least the PM2.5 exposure, dementia status, and the covariates used
#      in the original WHICAP MWAS (analogous to SALSA's S4 manifest +
#      S_4_all_*.csv pair). This is required to run the per-sample
#      analyses needed for inverse-variance meta-analysis or for
#      individual-participant-data pooling.
#
#   3. (Nice to have) The xMSannotator or in-house library annotation
#      table for WHICAP so the merged annotation set is consistent across
#      cohorts.
#
# With only per-feature summary statistics (estimate + std.error), a
# fixed-effects inverse-variance meta-analysis can still be run on the
# aligned features produced by this script — but per-sample harmonization
# (batch effects, matching covariate parameterization, sensitivity
# analyses) will not be possible without items 1 and 2.
# ----------------------------------------------------------------------

#--------------------------------End of the code--------------------------------
