## ---------------------------
##
## Script name: R10-imputation_revision.R
##
## Purpose of script:
##        Reviewer 1 major comment 6 and Reviewer 2 comment 5 -- proper multiple
##        imputation with Rubin's rules for the primary MWAS, plus a
##        complete-case sensitivity analysis.
##
##        scripts/3-clean_data.R:62 generates m = 5 imputations and keeps
##        `complete(1)`. That is single imputation: it ignores between-
##        imputation uncertainty and understates standard errors, which is
##        exactly what the Reviewer objected to. This script refits the primary
##        MWAS in EVERY completed dataset and combines with Rubin's rules --
##        the defect was never the number of imputations, it was discarding
##        four of the five that were already made.
##
## Author: Yufan Gong
##
## Notes:
##        WHICH VARIABLES ARE ACTUALLY IMPUTED. Among the 952 participants in
##        the analysis frame the missingness is: education 6 (0.6%), smoking 8
##        (0.8%), physical activity 41 (4.3%), alcohol 8 (0.8%), baseline age 2,
##        neighbourhood SES 0. Only education and smoking are in the PRIMARY
##        covariate set; physical activity and alcohol belong to the sensitivity
##        set. `age_at_blooddraw` is not imputed by the pipeline at all -- it is
##        selected from the visit-specific age variables
##        (3-clean_data.R:307), not derived from the imputed `blage` -- so the
##        two participants missing baseline age do not propagate into the
##        primary model.
##
##        WHY m = 10. Two constraints, and 10 satisfies both.
##
##        The statistical one is the m >= 100 x FMI rule of thumb, which is
##        stated in the fraction of missing INFORMATION, not the fraction of
##        missing rows; the two are two orders of magnitude apart here.
##        R20-fmi_pilot.R measures the FMI directly and finds 0.0016 for the
##        exposure coefficient this script reports, and 0.037 for the physical
##        activity coefficient, the most affected parameter anywhere in the
##        model. The rule therefore asks for m >= 4 at its most conservative.
##
##        The second is that Reviewer 1 asked specifically for more than the
##        five imputations the submitted pipeline generated. Five would have
##        cleared the statistical bar; ten clears the Reviewer's as well, at a
##        cost of a few minutes, and there is no reason to spend the exchange
##        arguing about it.
##
##        The reason the FMI is so much smaller than the missingness rate is
##        that the incomplete variables are ADJUSTMENT COVARIATES, while the
##        exposure and the outcome are complete for every specimen. Missingness
##        in a covariate propagates to the exposure coefficient only through the
##        covariate's partial association with it, which is weak here.
##
##        The m = 41 outputs are kept under tables/imputation/m41_check/ so the
##        response letter can state that a much larger m was run and changed
##        nothing: m = 5, 10 and 41 agree on every FDR-significant count.
##
##        IMPUTATION IS AT THE PARTICIPANT LEVEL. The covariates are baseline
##        participant attributes; imputing them once per specimen would treat
##        543 participants' repeated draws as independent observations and
##        understate the between-imputation variance. Completed values are
##        broadcast back to that participant's specimens.
##
##        duplicateCorrelation IS ESTIMATED ONCE, not per imputation. It is a
##        property of the blocking structure and the residual variance, and
##        re-estimating it 41 times per platform would cost days for a quantity
##        that moves in the fourth decimal. This is the same argument R4 already
##        makes for reusing it across exposures, and the value used is the one
##        R4 stored for total/all/covar.
##
##        The composites do not depend on the imputation: comp_pca_all is built
##        from the pollutant surfaces, not from the covariates, so it is
##        identical across the 41 completed datasets.
##
##        EXPOSURE WINDOWS (Reviewer 1 comment 13). The exposure is read from
##        the environment variable R10_EXPOSURE (default comp_pca_all, the
##        primary analysis). Setting it to comp_pca_all_w3 or comp_pca_all_w10
##        refits the same Rubin-pooled MWAS for the 3- and 10-year indices, so
##        the window comparison in the letter is pooled-vs-pooled:
##          R10_EXPOSURE=comp_pca_all_w3  Rscript scripts/revision/R10-imputation_revision.R
##          R10_EXPOSURE=comp_pca_all_w10 Rscript scripts/revision/R10-imputation_revision.R
##        Window runs write mwas_pooled_<window>.xlsx and
##        pooled_vs_single_<window>.xlsx only; the m = 5 check, complete-case
##        analysis and .RData are produced by the primary run alone. The
##        window index enters the imputation model in place of comp_pca_all,
##        so each run imputes with the exposure it analyses.
##
## ---------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R10-imputation_revision.R")

M_IMPUTATIONS <- rev_n(10, 3)
MAXIT         <- rev_n(50, 5)
SEED          <- 42
PRIMARY_EXPOSURE <- Sys.getenv("R10_EXPOSURE", "comp_pca_all")
IS_PRIMARY       <- PRIMARY_EXPOSURE == "comp_pca_all"
OUT_SUFFIX       <- if (IS_PRIMARY) "" else
  paste0("_", sub("^comp_pca_all_", "", PRIMARY_EXPOSURE))
message("Exposure: ", PRIMARY_EXPOSURE,
        if (IS_PRIMARY) " (primary)" else " (window sensitivity)")

lc <- function(d) dplyr::rename_all(d, stringr::str_to_lower)


# 1. Frames -------------------------------------------------------------------

load(rev_here("data", "processed", "combined_data_list_revision.RData"))
analysis <- combined_data_list_revision[["total"]][["all"]][["covar"]] |>
  dplyr::mutate(rand_id = as.character(rand_id))

message("Analysis frame: ", nrow(analysis), " specimens from ",
        dplyr::n_distinct(analysis$rand_id), " participants")

## Standardize the composite to SD = 1, exactly as R4 does before fitting
## (Reviewer 1 minor comments 5 and 6). Without this the coefficients here sit
## on the raw PC1 scale and every comparison against the published single-
## imputation result is off by a constant factor of 1 / SD -- which is what the
## first run of this script showed as a uniform SE ratio of 0.520 against
## comp_pca_all's SD of 1.9216. p-values and FDR are unaffected either way,
## but the agreement table is only interpretable on a common scale.
composite_sd <- stats::sd(analysis[[PRIMARY_EXPOSURE]], na.rm = TRUE)
analysis[[PRIMARY_EXPOSURE]] <- analysis[[PRIMARY_EXPOSURE]] / composite_sd
message("  ", PRIMARY_EXPOSURE, " standardized by SD = ",
        round(composite_sd, 4))

## Re-attach the covariate values as they were BEFORE 3-clean_data.R imputed
## them. The analysis frame carries the single-imputation values, so the
## missingness has to be restored from source to impute it properly.
raw <- lc(haven::read_sas(
  here::here("data", "salsa", "original_salsa_data",
             "salsa_data_04212016.sas7bdat")))
pa <- lc(haven::read_sas(
  here::here("data", "salsa", "original_salsa_data",
             "pa_nses_08042023.sas7bdat")))

raw_covars <- raw |>
  dplyr::transmute(
    rand_id  = as.character(rand_id),
    blage_raw = blage,
    edu_year_raw = ses3,
    mh62_raw = mh62,
    mh65, mh66, mh67
  ) |>
  dplyr::mutate(alcohol_raw = dplyr::case_when(
    mh65 == 1 | mh66 == 1 | mh67 == 1 ~ 1,
    mh65 == 0 & mh66 == 0 & mh67 == 0 ~ 0,
    TRUE ~ NA_real_)) |>
  dplyr::select(-mh65, -mh66, -mh67) |>
  dplyr::left_join(
    pa |> dplyr::transmute(rand_id = as.character(rand_id),
                           pa_raw = pa3_met_if_ca),
    by = "rand_id")

## One row per participant, for the imputation model.
impute_frame <- analysis |>
  dplyr::distinct(rand_id, .keep_all = TRUE) |>
  dplyr::select(rand_id, gender, ruca_metro, nses, wave, demcind,
                dplyr::all_of(PRIMARY_EXPOSURE)) |>
  dplyr::left_join(raw_covars, by = "rand_id")

missingness <- impute_frame |>
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~ sum(is.na(.x)))) |>
  tidyr::pivot_longer(dplyr::everything(), names_to = "variable",
                      values_to = "n_missing") |>
  dplyr::mutate(n_total = nrow(impute_frame),
                pct_missing = round(100 * n_missing / n_total, 2)) |>
  dplyr::arrange(dplyr::desc(n_missing))

message("\nPre-imputation missingness among the ", nrow(impute_frame),
        " participants:")
print(missingness |> dplyr::filter(n_missing > 0))
if (IS_PRIMARY) rev_save_table(missingness, "missingness_summary", "imputation")


# 2. Multiple imputation ------------------------------------------------------

message("\nRunning mice: m = ", M_IMPUTATIONS, ", maxit = ", MAXIT,
        ", method = pmm")

mice_input <- impute_frame |>
  dplyr::select(-rand_id) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

mids <- mice::mice(mice_input, m = M_IMPUTATIONS, maxit = MAXIT,
                   method = "pmm", seed = SEED, printFlag = FALSE)

if (IS_PRIMARY) rev_save_table(
  tibble::tibble(variable = names(mice_input),
                 method = as.character(mids$method),
                 in_primary_model = names(mice_input) %in%
                   c("edu_year_raw", "mh62_raw", "gender", "ruca_metro",
                     "nses", "wave")),
  "imputation_model", "imputation")


# 3. Metabolite matrices and the stored consensus correlation -----------------

load(here::here("data", "metabolomics", "processed", "med_c18_raw_combat.Rdata"))
load(here::here("data", "metabolomics", "processed", "med_hil_raw_combat.Rdata"))
metabo <- list(c18 = med_c18_raw_combat, hilic = med_hil_raw_combat) |>
  purrr::map(function(m) { colnames(m) <- tolower(colnames(m)); m })

load(rev_here("data", "metabolomics", "results",
              "duplicate_correlation_revision.RData"))
consensus <- list(
  c18   = dupcor_c18_list[["total"]][["all"]][["covar"]][[1]]$consensus.correlation,
  hilic = dupcor_hilic_list[["total"]][["all"]][["covar"]][[1]]$consensus.correlation
)
message("\nConsensus within-participant correlation reused from R4: C18 = ",
        round(consensus$c18, 5), ", HILIC = ", round(consensus$hilic, 5))

load(here::here("data", "processed", "salsa_clean.RData"))
sample_link <- list(c18 = salsa_blood_date_c18, hilic = salsa_blood_date_hilic) |>
  purrr::map(~ .x |>
               dplyr::transmute(rand_id = as.character(rand_id), blood_date,
                                file.name_new = tolower(file.name_new)) |>
               dplyr::distinct())

COVARS <- c("age_at_blooddraw", "gender", "edu_year", "mh62",
            "ruca_metro", "nses", "wave")


# 4. Fit the primary MWAS in one completed dataset ----------------------------

fit_one <- function(imp_index, platform) {
  completed <- mice::complete(mids, imp_index) |>
    dplyr::mutate(rand_id = impute_frame$rand_id) |>
    dplyr::select(rand_id, edu_year_imp = edu_year_raw, mh62_imp = mh62_raw)

  df <- analysis |>
    dplyr::left_join(completed, by = "rand_id") |>
    ## the imputed values replace the single-imputation ones
    dplyr::mutate(edu_year = edu_year_imp, mh62 = mh62_imp) |>
    dplyr::left_join(sample_link[[platform]], by = c("rand_id", "blood_date"))

  mat <- metabo[[platform]]
  df <- df |> dplyr::filter(file.name_new %in% colnames(mat))
  mat <- mat[, df$file.name_new, drop = FALSE]

  design <- stats::model.matrix(
    stats::as.formula(paste("~", PRIMARY_EXPOSURE, "+",
                            paste(COVARS, collapse = " + "))),
    data = df)
  keep <- rownames(design)
  mat <- mat[, as.integer(keep), drop = FALSE]

  fit <- limma::lmFit(mat, design, block = df$rand_id[as.integer(keep)],
                      correlation = consensus[[platform]]) |>
    limma::eBayes()

  tibble::tibble(
    met      = rownames(fit$coefficients),
    estimate = fit$coefficients[, 2],
    se       = fit$stdev.unscaled[, 2] * sqrt(fit$s2.post),
    df       = fit$df.total,
    imp      = imp_index,
    platform = platform
  )
}


# 5. Run every imputation and pool with Rubin's rules -------------------------

rubin_pool <- function(d) {
  d |>
    dplyr::group_by(met, platform) |>
    dplyr::summarise(
      m          = dplyr::n(),
      qbar       = mean(estimate),                       # pooled estimate
      ubar       = mean(se^2),                           # within-imputation var
      b          = stats::var(estimate),                 # between-imputation var
      total_var  = ubar + (1 + 1 / m) * b,
      se_pooled  = sqrt(total_var),
      lambda     = ((1 + 1 / m) * b) / total_var,
      ## Barnard-Rubin (1999) degrees of freedom. nu_obs is the OBSERVED-data
      ## df -- the complete-data df shrunk by the information lost to
      ## missingness -- not the complete-data df itself. With lambda of order
      ## 1e-4 the two agree to three decimal places here, but the response
      ## letter now prints this formula, so the code has to be the formula.
      nu_com    = mean(df),
      nu_obs    = ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - lambda),
      nu_old    = ifelse(lambda > 0, (m - 1) / lambda^2, Inf),
      df_br     = ifelse(lambda > 0,
                         (nu_old * nu_obs) / (nu_old + nu_obs), nu_com),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ## r is the relative increase in variance; fmi the fraction of missing
      ## information. lambda omits the df term and is NOT the FMI -- reporting
      ## lambda where the m >= 100 x FMI rule expects gamma understates it.
      r_increase = ((1 + 1 / m) * b) / ubar,
      fmi        = (r_increase + 2 / (df_br + 3)) / (r_increase + 1),
      t_pooled = qbar / se_pooled,
      p_pooled = 2 * stats::pt(-abs(t_pooled), df = df_br)
    ) |>
    dplyr::group_by(platform) |>
    dplyr::mutate(fdr_pooled = stats::p.adjust(p_pooled, method = "BH")) |>
    dplyr::ungroup()
}

message("\nFitting the primary MWAS in each of ", M_IMPUTATIONS,
        " completed datasets, per platform ...")
system.time({
  per_imp <- tidyr::expand_grid(imp = seq_len(M_IMPUTATIONS),
                                platform = c("c18", "hilic")) |>
    purrr::pmap(function(imp, platform) fit_one(imp, platform)) |>
    purrr::list_rbind()
})

mwas_pooled <- rubin_pool(per_imp)
rev_save_table(mwas_pooled, paste0("mwas_pooled", OUT_SUFFIX), "imputation")
save_pooled_annotated(mwas_pooled, paste0("mwas_pooled", OUT_SUFFIX))

message("\nRubin-pooled MWAS (", PRIMARY_EXPOSURE, "), FDR < 0.05 by platform:")
print(mwas_pooled |> dplyr::group_by(platform) |>
        dplyr::summarise(n_features = dplyr::n(),
                         n_fdr05 = sum(fdr_pooled < 0.05),
                         n_fdr10 = sum(fdr_pooled < 0.10),
                         .groups = "drop"))

if (IS_PRIMARY) {
## The same pooling over the FIRST FIVE of the same completed datasets, so the
## effect of m can be read off directly. Using a nested subset rather than a
## separate mice run is deliberate: it holds the imputations themselves fixed,
## so any difference is attributable to m alone and not to a different draw.
mwas_pooled_m5 <- per_imp |>
  dplyr::filter(.data$imp <= 5) |>
  rubin_pool()
rev_save_table(mwas_pooled_m5, "mwas_pooled_m5", "imputation")
save_pooled_annotated(mwas_pooled_m5, "mwas_pooled_m5")

m_concordance <- mwas_pooled |>
  dplyr::select(met, platform, est10 = .data$qbar, se10 = .data$se_pooled,
                fdr10 = .data$fdr_pooled) |>
  dplyr::inner_join(
    mwas_pooled_m5 |>
      dplyr::select(met, platform, est5 = .data$qbar, se5 = .data$se_pooled,
                    fdr5 = .data$fdr_pooled),
    by = c("met", "platform"))

rev_save_table(m_concordance, "m_concordance_5_vs_10", "imputation")

message("\nm = 5 against m = 10, over the same imputations:")
print(m_concordance |> dplyr::group_by(platform) |>
        dplyr::summarise(r_estimates = round(stats::cor(est5, est10), 5),
                         median_se_ratio = round(stats::median(se5 / se10), 5),
                         n_fdr05_m5 = sum(fdr5 < 0.05),
                         n_fdr05_m10 = sum(fdr10 < 0.05),
                         n_disagree = sum((fdr5 < 0.05) != (fdr10 < 0.05)),
                         .groups = "drop"))
}


# 6. Agreement with the single-imputation result ------------------------------

single <- c("c18", "hilic") |>
  purrr::set_names() |>
  purrr::map(function(pl) {
    readxl::read_xlsx(rev_here(
      "tables", "mwas_results", "total", "all", "covar",
      paste0("mwas_", pl, "_", PRIMARY_EXPOSURE, "_total_all_covar.xlsx"))) |>
      dplyr::transmute(met, single_estimate = logFC,
                       single_p = P.Value, single_fdr = adj.P.Val,
                       single_se = abs(logFC / t))
  }) |>
  purrr::list_rbind(names_to = "platform")

pooled_vs_single <- mwas_pooled |>
  dplyr::inner_join(single, by = c("met", "platform")) |>
  dplyr::mutate(se_ratio = se_pooled / single_se,
                est_ratio = qbar / single_estimate)

agreement <- pooled_vs_single |>
  dplyr::group_by(platform) |>
  dplyr::summarise(
    n                = dplyr::n(),
    r_estimates      = round(stats::cor(qbar, single_estimate), 4),
    median_se_ratio  = round(stats::median(se_ratio, na.rm = TRUE), 4),
    p95_se_ratio     = round(stats::quantile(se_ratio, 0.95, na.rm = TRUE), 4),
    max_se_ratio     = round(max(se_ratio, na.rm = TRUE), 4),
    n_fdr05_pooled   = sum(fdr_pooled < 0.05),
    n_fdr05_single   = sum(single_fdr < 0.05),
    n_sign_disagree  = sum(sign(qbar) != sign(single_estimate)),
    .groups = "drop"
  )

message("\nPooled versus single imputation:")
print(agreement)
rev_save_table(pooled_vs_single, paste0("pooled_vs_single", OUT_SUFFIX),
               "imputation")
rev_save_table(agreement, paste0("pooled_vs_single_agreement", OUT_SUFFIX),
               "imputation")

if (!IS_PRIMARY) {
  message("\nWindow run complete; m = 5 check, complete-case analysis and ",
          ".RData are produced by the primary run only.")
  quit(save = "no", status = 0)
}


# 7. Complete-case sensitivity analysis ---------------------------------------

## Participants with any missing PRIMARY covariate are dropped entirely, so a
## participant contributes all of their specimens or none.
incomplete_ids <- impute_frame |>
  dplyr::filter(is.na(edu_year_raw) | is.na(mh62_raw)) |>
  dplyr::pull(rand_id)

message("\nComplete-case: dropping ", length(incomplete_ids),
        " participants with a missing primary covariate")

fit_complete_case <- function(platform) {
  df <- analysis |>
    dplyr::filter(!rand_id %in% incomplete_ids) |>
    dplyr::left_join(sample_link[[platform]], by = c("rand_id", "blood_date"))
  mat <- metabo[[platform]]
  df <- df |> dplyr::filter(file.name_new %in% colnames(mat))
  mat <- mat[, df$file.name_new, drop = FALSE]

  design <- stats::model.matrix(
    stats::as.formula(paste("~", PRIMARY_EXPOSURE, "+",
                            paste(COVARS, collapse = " + "))), data = df)
  keep <- as.integer(rownames(design))
  fit <- limma::lmFit(mat[, keep, drop = FALSE], design,
                      block = df$rand_id[keep],
                      correlation = consensus[[platform]]) |>
    limma::eBayes()
  limma::topTable(fit, coef = 2, number = Inf, sort.by = "none",
                  adjust.method = "BH") |>
    tibble::rownames_to_column("met") |>
    dplyr::mutate(platform = platform,
                  n_specimens = length(keep),
                  n_participants = dplyr::n_distinct(df$rand_id[keep]))
}

mwas_cc <- c("c18", "hilic") |> purrr::map(fit_complete_case) |>
  purrr::list_rbind()
rev_save_table(mwas_cc, "mwas_completecase", "imputation")

message("\nComplete-case primary MWAS:")
print(mwas_cc |> dplyr::group_by(platform) |>
        dplyr::summarise(n_specimens = dplyr::first(n_specimens),
                         n_participants = dplyr::first(n_participants),
                         n_fdr05 = sum(adj.P.Val < 0.05), .groups = "drop"))

rev_dir("data", "processed")
save(mids, mwas_pooled, pooled_vs_single, agreement, mwas_cc, missingness,
     file = rev_here("data", "processed", "imputation_revision.RData"))

message("\nMultiple imputation completed.")
message("  tables -> ", rev_here("tables", "imputation"))

#--------------------------------End of the code--------------------------------
