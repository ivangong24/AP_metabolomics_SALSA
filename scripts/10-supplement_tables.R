## ---------------------------
##
## Script name: 10-supplement_tables.R
## Purpose of script: Build the consolidated supplement-table workbook for
##                    SALSA composite MWAS results and the cross-cohort
##                    (SALSA × WHICAP) fixed-effect meta-analysis.
##
## Author: Yufan Gong
##
## Date Created: 2026-06-05
##
## Notes:
##   - Output: tables/Supplement tables_composite_meta.xlsx
##       Sheets:
##         Table S1–S3: SALSA full / no demcind / demcind, primary covariates
##         Table S4–S6: SALSA full / no demcind / demcind, sensitivity covariates
##         Table S7–S8: Meta full cohort, primary / sensitivity covariates
##   - SALSA filter: composite exposures only
##                   (comp_wqs_all, comp_qgcomp_all for cross-sectional;
##                    comp_qgcomp_cox_all for time-to-event)
##                   AND (P.Value < 0.05 OR VIP_comp1 >= 2)
##                   Unannotated features ARE retained.
##   - Meta filter: same composite exposures, "all" population only,
##                  p_fe < 0.05.
##   - Visual format mirrors tables/Supplement tables_composite.xlsx:
##     row 1 holds a merged title; row 2 is the column header.
## ---------------------------

pacman::p_load(tidyverse, here, openxlsx)

load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))
load(here::here("data", "metabolomics", "annotation",
                "annotation_cleaned_wide.RData"))
load(here::here("data", "metabolomics", "alignment",
                "meta_mixture_annotated.RData"))

# 1. Lookup tables ----------------------------------------------------------

composite_exposures <- list(
  total = c("comp_wqs_all", "comp_qgcomp_all"),
  cox   = c("comp_qgcomp_cox_all")
)

analysis_label_lookup <- c(
  comp_wqs_all        = "Cross-sectional",
  comp_qgcomp_all     = "Cross-sectional",
  comp_qgcomp_cox_all = "Time-to-event"
)

# Parse m/z and RT from the "mz_rt_<m/z>_<rt>" feature ID.
parse_met <- function(met){
  parts <- stringr::str_split_fixed(met, "_", 4)
  tibble::tibble(mz = as.numeric(parts[, 3]),
                 rt = as.numeric(parts[, 4]))
}

# 2. SALSA table builder ----------------------------------------------------

build_salsa_tab <- function(population, covar) {
  pull_one <- function(combined_list, esi, annot_wide){
    purrr::imap(composite_exposures, function(exposures, study){
      dfl <- combined_list[[study]][[population]][[covar]]
      if (is.null(dfl)) return(NULL)
      exposures |>
        purrr::set_names() |>
        purrr::map(function(exp){
          d <- dfl[[exp]]
          if (is.null(d)) return(NULL)
          d |>
            dplyr::filter(P.Value < 0.05 | VIP_comp1 >= 2) |>
            dplyr::mutate(exposure = exp,
                          analysis = analysis_label_lookup[[exp]],
                          esi = esi)
        }) |>
        purrr::compact() |>
        purrr::list_rbind()
    }) |>
      purrr::list_rbind() |>
      dplyr::left_join(annot_wide, by = c("met" = "id"))
  }

  c18 <- pull_one(combined_results_list_c18,   "C18/neg-",   annotation_c18_wide)
  hil <- pull_one(combined_results_list_hilic, "HILIC/pos+", annotation_hilic_wide)

  dplyr::bind_rows(c18, hil) |>
    dplyr::mutate(parse_met(met)) |>
    dplyr::transmute(
      `m/z`                = mz,
      `RT(sec)`            = rt,
      annotation           = compound,
      `chemical ID`        = chemical_id,
      `column/ESI`         = esi,
      `confidence level`   = confidence,
      `reference database` = reference,
      logFC                = logFC,
      p                    = P.Value,
      FDR                  = adj.P.Val,
      VIP                  = VIP_comp1,
      exposure             = exposure,
      analysis             = analysis
    ) |>
    dplyr::arrange(p)
}

# 3. Meta-analysis table builder --------------------------------------------

build_meta_tab <- function(covar) {
  pull_meta <- function(meta_list, esi){
    purrr::imap(composite_exposures, function(exposures, study){
      dfl <- meta_list[[study]][["all"]][[covar]]
      if (is.null(dfl)) return(NULL)
      exposures |>
        purrr::set_names() |>
        purrr::map(function(exp){
          d <- dfl[[exp]]
          if (is.null(d)) return(NULL)
          d |>
            dplyr::filter(!is.na(p_fe), p_fe < 0.05) |>
            dplyr::mutate(exposure = exp,
                          analysis = analysis_label_lookup[[exp]],
                          esi = esi)
        }) |>
        purrr::compact() |>
        purrr::list_rbind()
    }) |>
      purrr::list_rbind()
  }

  c18 <- pull_meta(meta_annotated_c18,   "C18/neg-")
  hil <- pull_meta(meta_annotated_hilic, "HILIC/pos+")

  dplyr::bind_rows(c18, hil) |>
    dplyr::transmute(
      `m/z`                = whicap_mz,
      `RT(sec)`            = whicap_time,
      annotation           = compound,
      `chemical ID`        = chemical_id,
      `column/ESI`         = esi,
      `confidence level`   = confidence,
      `reference database` = reference,
      `SALSA logFC`        = estimate_salsa,
      `SALSA p`            = p_salsa,
      `SALSA VIP`          = vip_salsa,
      `WHICAP logFC`       = estimate_whicap,
      `WHICAP p`           = p_whicap,
      `FE estimate`        = estimate_fe,
      `FE 95% CI`          = dplyr::if_else(
                               is.na(ci_lb_fe) | is.na(ci_ub_fe),
                               NA_character_,
                               sprintf("(%.4f, %.4f)", ci_lb_fe, ci_ub_fe)),
      `FE p`               = p_fe,
      `FE FDR`             = fdr_fe,
      `I^2 (%)`            = i2_pct,
      `Cochran Q p`        = p_het,
      exposure             = exposure,
      analysis             = analysis
    ) |>
    dplyr::arrange(`FE p`)
}

# 4. Workbook assembly ------------------------------------------------------

title_style  <- openxlsx::createStyle(textDecoration = "bold",
                                      fontSize = 12,
                                      halign = "left")
header_style <- openxlsx::createStyle(textDecoration = "bold",
                                      border = "bottom",
                                      halign = "left")

write_sheet <- function(wb, sheet_name, title, data) {
  openxlsx::addWorksheet(wb, sheet_name)
  ncols <- ncol(data)
  openxlsx::writeData(wb, sheet_name, x = title,
                      startRow = 1, startCol = 1)
  openxlsx::mergeCells(wb, sheet_name, cols = 1:ncols, rows = 1)
  openxlsx::addStyle(wb, sheet_name, style = title_style,
                     rows = 1, cols = 1)
  openxlsx::writeData(wb, sheet_name, x = data,
                      startRow = 2, startCol = 1)
  openxlsx::addStyle(wb, sheet_name, style = header_style,
                     rows = 2, cols = 1:ncols, gridExpand = TRUE)
  openxlsx::freezePane(wb, sheet_name, firstActiveRow = 3,
                       firstActiveCol = 1)
}

salsa_specs <- list(
  list(sheet = "Table S1", pop = "all",        covar = "covar",
       title = "Table S1. Summary of MWAS-identified features (P < 0.05 or VIP ≥ 2) for air-toxicant mixture composites in the full cohort using the primary covariates"),
  list(sheet = "Table S2", pop = "no demcind", covar = "covar",
       title = "Table S2. Summary of MWAS-identified features (P < 0.05 or VIP ≥ 2) for air-toxicant mixture composites in participants without dementia/CIND using the primary covariates"),
  list(sheet = "Table S3", pop = "demcind",    covar = "covar",
       title = "Table S3. Summary of MWAS-identified features (P < 0.05 or VIP ≥ 2) for air-toxicant mixture composites in participants with dementia/CIND using the primary covariates"),
  list(sheet = "Table S4", pop = "all",        covar = "covar_sen",
       title = "Table S4. Summary of MWAS-identified features (P < 0.05 or VIP ≥ 2) for air-toxicant mixture composites in the full cohort using the sensitivity covariates"),
  list(sheet = "Table S5", pop = "no demcind", covar = "covar_sen",
       title = "Table S5. Summary of MWAS-identified features (P < 0.05 or VIP ≥ 2) for air-toxicant mixture composites in participants without dementia/CIND using the sensitivity covariates"),
  list(sheet = "Table S6", pop = "demcind",    covar = "covar_sen",
       title = "Table S6. Summary of MWAS-identified features (P < 0.05 or VIP ≥ 2) for air-toxicant mixture composites in participants with dementia/CIND using the sensitivity covariates")
)

meta_specs <- list(
  list(sheet = "Table S7", covar = "covar",
       title = "Table S7. Summary of cross-cohort (SALSA × WHICAP) fixed-effect meta-analysis features at P_FE < 0.05 in the full cohort using the primary covariates"),
  list(sheet = "Table S8", covar = "covar_sen",
       title = "Table S8. Summary of cross-cohort (SALSA × WHICAP) fixed-effect meta-analysis features at P_FE < 0.05 in the full cohort using the sensitivity covariates")
)

wb <- openxlsx::createWorkbook()

purrr::walk(salsa_specs, function(spec) {
  message("Building SALSA ", spec$sheet, " (", spec$pop, " / ",
          spec$covar, ") ...")
  d <- build_salsa_tab(spec$pop, spec$covar)
  message("  ", nrow(d), " rows")
  write_sheet(wb, spec$sheet, spec$title, d)
})

purrr::walk(meta_specs, function(spec) {
  message("Building Meta ", spec$sheet, " (all / ", spec$covar, ") ...")
  d <- build_meta_tab(spec$covar)
  message("  ", nrow(d), " rows")
  write_sheet(wb, spec$sheet, spec$title, d)
})

out_path <- here::here("tables", "Supplement tables_composite_meta.xlsx")
openxlsx::saveWorkbook(wb, file = out_path, overwrite = TRUE)
message("Supplement workbook saved to: ", out_path)

#--------------------------------End of the code--------------------------------
