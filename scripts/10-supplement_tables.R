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
##         Table S1–S3:   SALSA full / no demcind / demcind, primary covariates
##         Table S4–S6:   SALSA full / no demcind / demcind, sensitivity covariates
##         Table S7–S8:   Meta full cohort, primary / sensitivity covariates
##         Table S9–S10:  Mummichog enriched pathways (P(Fisher) < 0.05),
##                        main analysis, primary / sensitivity covariates
##         Table S11–S12: Mummichog enriched pathways (P(Fisher) < 0.05),
##                        meta-analysis, primary / sensitivity covariates
##         Table S13–S14: Feature-to-pathway mapping for enriched pathways,
##                        primary / sensitivity covariates
##   - SALSA filter: composite exposures only
##                   (comp_wqs_all, comp_qgcomp_all for cross-sectional;
##                    comp_qgcomp_cox_all for time-to-event)
##                   AND (P.Value < 0.05 OR VIP_comp1 >= 2)
##                   Unannotated features ARE retained.
##   - Meta filter: same composite exposures, "all" population only,
##                  p_fe < 0.05.
##   - S1–S8 exposure column merges the previous (exposure, analysis) pair
##     into a single label: WQS / QGcomp / QGcomp (COX).
##   - Visual format mirrors tables/Supplement tables_composite.xlsx:
##     row 1 holds a merged title; row 2 is the column header.
## ---------------------------

pacman::p_load(tidyverse, here, openxlsx, KEGGREST)

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

# Single-label combining the exposure + analysis columns. Used both for the
# merged "exposure" column in S1-S8 and as the row label in the pathway sheets.
exposure_label_lookup <- c(
  comp_wqs_all        = "WQS",
  comp_qgcomp_all     = "QGcomp",
  comp_qgcomp_cox_all = "QGcomp (COX)"
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
            dplyr::mutate(exposure = exposure_label_lookup[[exp]],
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
      exposure             = exposure
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
            dplyr::mutate(exposure = exposure_label_lookup[[exp]],
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
      exposure             = exposure
    ) |>
    dplyr::arrange(`FE p`)
}

# 4. Mummichog pathway builders --------------------------------------------

# Where the MetaboAnalyst Mummichog outputs live:
#   main analysis:  metaboAnalyst/Output/{total|cox}/{population}/{covar}/{exp}
#   meta analysis:  metaboAnalyst/Output/meta/{total|cox}/all/{covar}/{exp}
mum_root <- here::here("metaboAnalyst", "Output")

# Map of which exposures live under "total" (cross-sectional) vs "cox".
mummi_exposures <- list(
  total = c("comp_wqs_all", "comp_qgcomp_all"),
  cox   = c("comp_qgcomp_cox_all")
)
mummi_populations <- c("all", "no demcind", "demcind")

read_pathway_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw <- utils::read.csv(path, check.names = FALSE,
                         stringsAsFactors = FALSE)
  # MetaboAnalyst writes the pathway name as an unnamed first column.
  names(raw)[1] <- "pathway_name"
  raw
}

read_matched_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Build the "List of enriched pathways" sheet for either main or meta scope.
#   scope: "main" -> iterates total/cox × populations, P(Fisher) < 0.05
#          "meta" -> iterates meta/{total|cox} × "all" only, P(Fisher) < 0.1
build_pathway_tab <- function(covar, scope = c("main", "meta")) {
  scope <- match.arg(scope)
  pops    <- if (scope == "main") mummi_populations else "all"
  base    <- if (scope == "main") mum_root         else file.path(mum_root, "meta")
  p_cut   <- if (scope == "main") 0.05             else 0.1

  purrr::imap(mummi_exposures, function(exposures, study) {
    purrr::map(pops, function(pop) {
      purrr::map(exposures, function(exp) {
        path <- file.path(base, study, pop, covar, exp,
                          "mummichog_pathway_enrichment_mummichog.csv")
        d <- read_pathway_csv(path)
        if (is.null(d)) return(NULL)
        d |>
          dplyr::transmute(
            `pathway name`   = pathway_name,
            `pathway size`   = `Pathway total`,
            `# of peaks`     = Hits.total,
            `# of sig peaks` = Hits.sig,
            p                = `P(Fisher)`,
            empirical        = Empirical,
            adj.p            = AdjP.Fisher,
            exposure         = exposure_label_lookup[[exp]],
            population       = pop
          ) |>
          dplyr::filter(!is.na(p), p < p_cut)
      }) |> purrr::compact() |> purrr::list_rbind()
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind() |>
    dplyr::arrange(p)
}

# KEGG compound-name lookup. Queries the REST API in batches of 10 and caches
# results in-memory; ids that aren't KEGG-formatted (e.g. CE1556, BiGG abbrevs)
# or that 404 on KEGG return NA so callers can drop them.
.kegg_name_cache <- new.env(parent = emptyenv())
get_kegg_names <- function(ids) {
  ids <- unique(ids)
  to_query <- ids[grepl("^C[0-9]+$", ids) & !vapply(ids,
                  function(x) exists(x, envir = .kegg_name_cache),
                  logical(1))]
  if (length(to_query) > 0) {
    message("  KEGG: looking up ", length(to_query), " compound names ...")
    chunks <- split(to_query, ceiling(seq_along(to_query) / 10))
    for (ch in chunks) {
      res <- tryCatch(KEGGREST::keggGet(ch), error = function(e) NULL)
      hit_ids <- character(0)
      if (!is.null(res)) {
        for (e in res) {
          nm <- if (!is.null(e$NAME)) trimws(gsub(";$", "", e$NAME[1])) else NA_character_
          if (is.null(nm) || !nzchar(nm)) nm <- NA_character_
          assign(e$ENTRY, nm, envir = .kegg_name_cache)
          hit_ids <- c(hit_ids, e$ENTRY)
        }
      }
      # Cache misses so we don't retry them.
      for (miss in setdiff(ch, hit_ids)) {
        assign(miss, NA_character_, envir = .kegg_name_cache)
      }
    }
  }
  out <- vapply(ids, function(x) {
    if (exists(x, envir = .kegg_name_cache))
      get(x, envir = .kegg_name_cache)
    else NA_character_
  }, character(1))
  setNames(out, ids)
}

# Build the feature-to-pathway mapping sheet for the requested covariate set.
# Only mappings for pathways that pass P(Fisher) < 0.05 are kept, matching the
# enriched-pathway filter in build_pathway_tab().
build_feature_pathway_tab <- function(covar) {
  pops <- mummi_populations
  purrr::imap(mummi_exposures, function(exposures, study) {
    purrr::map(pops, function(pop) {
      purrr::map(exposures, function(exp) {
        dir_path <- file.path(mum_root, study, pop, covar, exp)
        path_csv <- read_pathway_csv(file.path(dir_path,
                                  "mummichog_pathway_enrichment_mummichog.csv"))
        comp_csv <- read_matched_csv(file.path(dir_path,
                                  "mummichog_matched_compound_all.csv"))
        if (is.null(path_csv) || is.null(comp_csv)) return(NULL)

        sig <- path_csv |>
          dplyr::filter(!is.na(`P(Fisher)`), `P(Fisher)` < 0.05) |>
          dplyr::select(`pathway name` = pathway_name,
                        cpd.hits)
        if (nrow(sig) == 0) return(NULL)

        # Explode the semicolon-separated empirical-compound list and join
        # to the per-feature matched table.
        long <- sig |>
          dplyr::mutate(
            Empirical.Compound = stringr::str_split(cpd.hits, ";")
          ) |>
          tidyr::unnest(Empirical.Compound) |>
          dplyr::select(-cpd.hits) |>
          dplyr::inner_join(comp_csv, by = "Empirical.Compound",
                            relationship = "many-to-many")

        long |>
          dplyr::transmute(
            `pathway name`     = `pathway name`,
            `matched compound` = Matched.Compound,
            met = paste0("mz_rt_",
                         round(as.numeric(Query.Mass), 4), "_",
                         round(as.numeric(Retention.Time), 4)),
            adduct           = Matched.Form,
            exposure         = exposure_label_lookup[[exp]],
            population       = pop
          ) |>
          dplyr::distinct()
      }) |> purrr::compact() |> purrr::list_rbind()
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind() |>
    # Attach KEGG compound names; drop rows where the matched compound is not
    # KEGG-formatted or the ID can't be resolved against KEGG.
    {\(d) {
      if (nrow(d) == 0) return(dplyr::mutate(d, `compound name` = character(0)))
      lookup <- get_kegg_names(d$`matched compound`)
      d |>
        dplyr::mutate(`compound name` = unname(lookup[d$`matched compound`])) |>
        dplyr::filter(!is.na(`compound name`))
    }}() |>
    dplyr::select(`pathway name`, `matched compound`, met, adduct,
                  `compound name`, exposure, population)
}

# 5. Workbook assembly ------------------------------------------------------

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

pathway_specs <- list(
  list(sheet = "Table S9",  scope = "main", covar = "covar",
       title = "Table S9. List of enriched pathways (P < 0.05) from the Mummichog main analysis using the primary covariates"),
  list(sheet = "Table S10", scope = "main", covar = "covar_sen",
       title = "Table S10. List of enriched pathways (P < 0.05) from the Mummichog main analysis using the sensitivity covariates"),
  list(sheet = "Table S11", scope = "meta", covar = "covar",
       title = "Table S11. List of enriched pathways (P < 0.05) from the cross-cohort (SALSA × WHICAP) Mummichog meta-analysis using the primary covariates"),
  list(sheet = "Table S12", scope = "meta", covar = "covar_sen",
       title = "Table S12. List of enriched pathways (P < 0.05) from the cross-cohort (SALSA × WHICAP) Mummichog meta-analysis using the sensitivity covariates")
)

feature_pathway_specs <- list(
  list(sheet = "Table S13", covar = "covar",
       title = "Table S13. Feature-to-pathway mapping table for enriched pathways (P < 0.05) using the primary covariates"),
  list(sheet = "Table S14", covar = "covar_sen",
       title = "Table S14. Feature-to-pathway mapping table for enriched pathways (P < 0.05) using the sensitivity covariates")
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

purrr::walk(pathway_specs, function(spec) {
  message("Building Pathway ", spec$sheet, " (", spec$scope, " / ",
          spec$covar, ") ...")
  d <- build_pathway_tab(spec$covar, scope = spec$scope)
  message("  ", nrow(d), " rows")
  write_sheet(wb, spec$sheet, spec$title, d)
})

purrr::walk(feature_pathway_specs, function(spec) {
  message("Building Feature-to-pathway ", spec$sheet, " (", spec$covar, ") ...")
  d <- build_feature_pathway_tab(spec$covar)
  message("  ", nrow(d), " rows")
  write_sheet(wb, spec$sheet, spec$title, d)
})

out_path <- here::here("tables", "Supplement tables_composite_meta.xlsx")
openxlsx::saveWorkbook(wb, file = out_path, overwrite = TRUE)
message("Supplement workbook saved to: ", out_path)

#--------------------------------End of the code--------------------------------
