## ---------------------------
##
## Script name: 10b-supplement_tables_concise.R
## Purpose of script: Build a *concise* companion to the main supplement
##                    workbook produced by 10-supplement_tables.R. The
##                    concise version reports only FDR < 0.05 features and
##                    pools all platforms, exposures, and populations into
##                    a single sheet per covar set, with platform / exposure /
##                    population recorded as added columns so each row can
##                    still be traced back to its slice.
##
## Author: Yufan Gong
##
## Date Created: 2026-06-10
##
## Notes:
##   - Output: Supplement/Supplement tables_concise.xlsx
##         Table S1 / S1a: SALSA, primary covar, FDR < 0.05
##         Table S2 / S2a: SALSA, sensitivity covar, FDR < 0.05
##         Table S3 / S3a: Meta (SALSA × WHICAP), primary covar, FE FDR < 0.05
##         Table S4 / S4a: Meta, sensitivity covar, FE FDR < 0.05
##         Table S5-S8:    Mummichog pathway enrichment sheets
##                         (main primary, main sensitivity, meta primary,
##                         meta sensitivity); carried over from the main
##                         workbook with renumbered IDs.
##         Table S9-S10:   Feature-to-pathway mapping sheets (primary,
##                         sensitivity); carried over with renumbered IDs.
##   - SALSA stats columns: met, logFC, 95% CI, p, FDR, VIP, platform,
##                          exposure, population.
##   - Meta stats columns:  met, FE estimate, FE 95% CI, FE p, FE FDR,
##                          I^2 (%), Cochran Q p, platform, exposure.
##   - Annotated sheets (*a) carry only rows with a non-empty annotation
##     and only the annotation columns + platform / exposure / population.
##   - This script intentionally duplicates helpers from
##     10-supplement_tables.R rather than sourcing them, so it can be run
##     standalone.
## ---------------------------

pacman::p_load(tidyverse, here, openxlsx, KEGGREST, limma)

load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))
load(here::here("data", "metabolomics", "annotation",
                "annotation_cleaned_wide.RData"))
load(here::here("data", "metabolomics", "alignment",
                "meta_mixture_annotated.RData"))
load(here::here("data", "metabolomics", "results", "limma_fit_c18.RData"))
load(here::here("data", "metabolomics", "results", "limma_fit_hilic.RData"))

# 1. Lookup tables ----------------------------------------------------------

composite_exposures <- list(
  total = c("comp_wqs_all", "comp_qgcomp_all"),
  cox   = c("comp_qgcomp_cox_all")
)

exposure_label_lookup <- c(
  comp_wqs_all        = "WQS",
  comp_qgcomp_all     = "QGcomp",
  comp_qgcomp_cox_all = "QGcomp (COX)"
)

slice_plan <- tibble::tribble(
  ~idx, ~esi_key, ~esi_label, ~exp_key,              ~study,
  1,    "c18",    "C18",      "comp_wqs_all",        "total",
  2,    "c18",    "C18",      "comp_qgcomp_all",     "total",
  3,    "c18",    "C18",      "comp_qgcomp_cox_all", "cox",
  4,    "hilic",  "HILIC",    "comp_wqs_all",        "total",
  5,    "hilic",  "HILIC",    "comp_qgcomp_all",     "total",
  6,    "hilic",  "HILIC",    "comp_qgcomp_cox_all", "cox"
)
slice_plan$exp_label <- exposure_label_lookup[slice_plan$exp_key]

salsa_populations <- c("all", "no demcind", "demcind")
population_label_lookup <- c(
  "all"        = "full cohort",
  "no demcind" = "no dementia/CIND",
  "demcind"    = "dementia/CIND"
)

parse_met <- function(met){
  parts <- stringr::str_split_fixed(met, "_", 4)
  tibble::tibble(mz = as.numeric(parts[, 3]),
                 rt = as.numeric(parts[, 4]))
}

format_ci <- function(lo, hi) {
  dplyr::if_else(is.na(lo) | is.na(hi),
                 NA_character_,
                 sprintf("(%.2f, %.2f)", lo, hi))
}

# 2. logFC CI extraction ----------------------------------------------------

extract_logfc_ci <- function(fit, coef) {
  tt <- limma::topTable(fit,
                        coef     = coef,
                        confint  = TRUE,
                        number   = Inf,
                        sort.by  = "none")
  tibble::tibble(
    met   = rownames(tt),
    CI.L  = tt$CI.L,
    CI.R  = tt$CI.R
  )
}

# 3. SALSA pooled builder ---------------------------------------------------
# Walks every slice × population for a given covar set and returns a single
# long tibble of FDR < 0.05 features. Adds platform / exposure / population
# columns so rows from different slices are distinguishable in the pooled
# sheet.

build_salsa_pooled_wide <- function(covar) {
  purrr::map(salsa_populations, function(pop) {
    purrr::map(seq_len(nrow(slice_plan)), function(i) {
      slice <- as.list(slice_plan[i, ])

      combined_list <- if (slice$esi_key == "c18") {
        combined_results_list_c18
      } else {
        combined_results_list_hilic
      }
      annot_wide <- if (slice$esi_key == "c18") {
        annotation_c18_wide
      } else {
        annotation_hilic_wide
      }
      fit_list <- if (slice$esi_key == "c18") {
        limma_fit_c18
      } else {
        limma_fit_hilic
      }

      d <- combined_list[[slice$study]][[pop]][[covar]][[slice$exp_key]]
      if (is.null(d)) return(NULL)

      fit <- fit_list[[slice$study]][[pop]][[covar]][[slice$exp_key]]
      ci  <- extract_logfc_ci(fit, slice$exp_key)

      d |>
        dplyr::filter(adj.P.Val < 0.05) |>
        dplyr::left_join(ci,         by = "met") |>
        dplyr::left_join(annot_wide, by = c("met" = "id")) |>
        dplyr::mutate(parse_met(met)) |>
        dplyr::transmute(
          met         = paste0("mz_rt_", round(mz, 4), "_", round(rt, 4)),
          annotation  = compound,
          `chemical ID`        = chemical_id,
          `confidence level`   = confidence,
          `reference database` = reference,
          logFC       = logFC,
          `95% CI`    = format_ci(CI.L, CI.R),
          p           = P.Value,
          FDR         = adj.P.Val,
          VIP         = VIP_comp1,
          platform    = slice$esi_label,
          exposure    = slice$exp_label,
          population  = population_label_lookup[[pop]]
        )
    }) |> purrr::compact() |> purrr::list_rbind()
  }) |> purrr::list_rbind() |>
    dplyr::arrange(population, exposure, platform, FDR)
}

salsa_stats_cols <- c("met", "logFC", "95% CI", "p", "FDR", "VIP",
                      "platform", "exposure", "population")
salsa_anno_cols  <- c("met", "annotation", "chemical ID",
                      "confidence level", "reference database",
                      "platform", "exposure", "population")

trim_salsa_stats <- function(d) dplyr::select(d, dplyr::all_of(salsa_stats_cols))
trim_salsa_anno  <- function(d) {
  d |>
    dplyr::filter(!is.na(annotation), nzchar(annotation)) |>
    dplyr::select(dplyr::all_of(salsa_anno_cols)) |>
    dplyr::mutate(
      annotation    = gsub(";\\s*", "\n", annotation),
      `chemical ID` = gsub(";\\s*", "\n", `chemical ID`)
    )
}

# 4. Meta pooled builder ----------------------------------------------------

build_meta_pooled_wide <- function(covar) {
  purrr::map(seq_len(nrow(slice_plan)), function(i) {
    slice <- as.list(slice_plan[i, ])
    meta_list <- if (slice$esi_key == "c18") meta_annotated_c18 else meta_annotated_hilic
    d <- meta_list[[slice$study]][["all"]][[covar]][[slice$exp_key]]
    if (is.null(d)) return(NULL)

    d |>
      dplyr::filter(!is.na(fdr_fe), fdr_fe < 0.05) |>
      dplyr::transmute(
        met                  = paste0("mz_rt_",
                                      round(whicap_mz, 4), "_",
                                      round(whicap_time, 4)),
        annotation           = compound,
        `chemical ID`        = chemical_id,
        `confidence level`   = confidence,
        `reference database` = reference,
        `FE estimate`        = estimate_fe,
        `FE 95% CI`          = format_ci(ci_lb_fe, ci_ub_fe),
        `FE p`               = p_fe,
        `FE FDR`             = fdr_fe,
        `I^2 (%)`            = i2_pct,
        `Cochran Q p`        = p_het,
        platform             = slice$esi_label,
        exposure             = slice$exp_label
      )
  }) |> purrr::compact() |> purrr::list_rbind() |>
    dplyr::arrange(exposure, platform, `FE FDR`)
}

meta_stats_cols <- c("met",
                     "FE estimate", "FE 95% CI", "FE p", "FE FDR",
                     "I^2 (%)", "Cochran Q p",
                     "platform", "exposure")
meta_anno_cols  <- c("met", "annotation", "chemical ID",
                     "confidence level", "reference database",
                     "platform", "exposure")

trim_meta_stats <- function(d) dplyr::select(d, dplyr::all_of(meta_stats_cols))
trim_meta_anno  <- function(d) {
  d |>
    dplyr::filter(!is.na(annotation), nzchar(annotation)) |>
    dplyr::select(dplyr::all_of(meta_anno_cols)) |>
    dplyr::mutate(
      annotation    = gsub(";\\s*", "\n", annotation),
      `chemical ID` = gsub(";\\s*", "\n", `chemical ID`)
    )
}

# 5. Mummichog pathway builders --------------------------------------------

mum_root <- here::here("metaboAnalyst", "Output")

mummi_exposures <- list(
  total = c("comp_wqs_all", "comp_qgcomp_all"),
  cox   = c("comp_qgcomp_cox_all")
)
mummi_populations <- c("all", "no demcind", "demcind")

read_pathway_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw <- utils::read.csv(path, check.names = FALSE,
                         stringsAsFactors = FALSE)
  names(raw)[1] <- "pathway_name"
  raw
}

read_matched_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

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
            population       = population_label_lookup[[pop]]
          ) |>
          dplyr::filter(!is.na(p), p < p_cut)
      }) |> purrr::compact() |> purrr::list_rbind()
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind() |>
    dplyr::arrange(p)
}

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
            population       = population_label_lookup[[pop]]
          ) |>
          dplyr::distinct()
      }) |> purrr::compact() |> purrr::list_rbind()
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind() |>
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

# 6. Workbook styles --------------------------------------------------------

title_style  <- openxlsx::createStyle(fontName = "Times New Roman",
                                      fontSize = 12,
                                      textDecoration = "bold",
                                      halign = "left",
                                      valign = "top",
                                      wrapText = TRUE)
header_style <- openxlsx::createStyle(fontName = "Times New Roman",
                                      fontSize = 11,
                                      textDecoration = "bold",
                                      border = c("top", "bottom"),
                                      borderStyle = c("thin", "thin"),
                                      halign = "left")
body_general_style <- openxlsx::createStyle(fontName = "Times New Roman",
                                            fontSize = 11)
body_number_style  <- openxlsx::createStyle(fontName = "Times New Roman",
                                            fontSize = 11,
                                            numFmt = "0.00")
body_sci_style     <- openxlsx::createStyle(fontName = "Times New Roman",
                                            fontSize = 11,
                                            numFmt = "0.00E+00")
last_row_border    <- openxlsx::createStyle(border = "bottom",
                                            borderStyle = "thin")
body_wrap_style    <- openxlsx::createStyle(fontName = "Times New Roman",
                                            fontSize = 11,
                                            wrapText = TRUE,
                                            valign = "top")
body_vcenter_style <- openxlsx::createStyle(fontName = "Times New Roman",
                                            fontSize = 11,
                                            valign = "center")

number_cols <- c("logFC", "VIP",
                 "SALSA logFC", "SALSA VIP",
                 "WHICAP logFC",
                 "FE estimate", "I^2 (%)")
sci_cols    <- c("p", "FDR", "empirical", "adj.p",
                 "SALSA p", "WHICAP p",
                 "FE p", "FE FDR", "Cochran Q p")

write_sheet <- function(wb, sheet_name, title, data) {
  openxlsx::addWorksheet(wb, sheet_name)
  ncols <- max(1L, ncol(data))
  nrows <- nrow(data)

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

  if (nrows > 0) {
    data_rows <- 3:(2 + nrows)
    openxlsx::addStyle(wb, sheet_name, style = body_general_style,
                       rows = data_rows, cols = 1:ncols, gridExpand = TRUE)
    for (ci in seq_along(colnames(data))) {
      cn <- colnames(data)[ci]
      if (cn %in% number_cols) {
        openxlsx::addStyle(wb, sheet_name, style = body_number_style,
                           rows = data_rows, cols = ci, gridExpand = TRUE)
      } else if (cn %in% sci_cols) {
        openxlsx::addStyle(wb, sheet_name, style = body_sci_style,
                           rows = data_rows, cols = ci, gridExpand = TRUE)
      }
    }
    last_row <- 2 + nrows
    openxlsx::addStyle(wb, sheet_name, style = last_row_border,
                       rows = last_row, cols = 1:ncols,
                       gridExpand = TRUE, stack = TRUE)

    wrap_cols <- c("annotation", "chemical ID",
                   "pathway name", "matched compound", "compound name")
    wrap_idx <- which(colnames(data) %in% wrap_cols)
    if (length(wrap_idx)) {
      openxlsx::addStyle(wb, sheet_name, style = body_wrap_style,
                         rows = data_rows, cols = wrap_idx,
                         gridExpand = TRUE, stack = TRUE)
      vcenter_cols <- which(colnames(data) %in%
                              c("met", "confidence level", "reference database",
                                "platform", "exposure", "population"))
      if (length(vcenter_cols)) {
        openxlsx::addStyle(wb, sheet_name, style = body_vcenter_style,
                           rows = data_rows, cols = vcenter_cols,
                           gridExpand = TRUE, stack = TRUE)
      }
    }
  }

  WRAP_CAP        <- 30
  WRAP_CAP_NARROW <- 22
  widths <- vapply(seq_len(ncols), function(ci) {
    cn <- colnames(data)[ci]
    is_wrap        <- cn %in% c("annotation", "chemical ID",
                                "pathway name", "matched compound",
                                "compound name")
    is_wrap_narrow <- cn %in% c("pathway name", "matched compound",
                                "compound name")
    is_num  <- cn %in% number_cols
    is_sci  <- cn %in% sci_cols
    raw     <- data[[ci]]

    display <- if (is_num) {
      formatC(suppressWarnings(as.numeric(raw)), format = "f", digits = 2)
    } else if (is_sci) {
      formatC(suppressWarnings(as.numeric(raw)), format = "e", digits = 2)
    } else {
      as.character(raw)
    }
    display[is.na(display) | display == "NA"] <- ""

    val_len <- if (length(display) == 0) 0L else max(vapply(display, function(v) {
      lines <- strsplit(v, "\n", fixed = TRUE)[[1]]
      if (!length(lines)) 0L else max(nchar(lines))
    }, integer(1)))
    if (is_wrap_narrow) {
      val_len <- min(val_len, WRAP_CAP_NARROW)
    } else if (is_wrap) {
      val_len <- min(val_len, WRAP_CAP)
    }
    as.numeric(max(nchar(cn), val_len) + 2)
  }, numeric(1))
  openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncols), widths = widths)

  chars_per_line <- max(1, sum(widths) * 0.92)
  n_lines        <- max(1, ceiling(nchar(title) / chars_per_line))
  openxlsx::setRowHeights(wb, sheet_name, rows = 1,
                          heights = n_lines * 13)
}

# 7. Spec lists -------------------------------------------------------------

salsa_concise_specs <- list(
  list(id = "S1", covar = "covar",
       covar_label = "primary covariates"),
  list(id = "S2", covar = "covar_sen",
       covar_label = "sensitivity covariates")
)

meta_concise_specs <- list(
  list(id = "S3", covar = "covar",
       covar_label = "primary covariates"),
  list(id = "S4", covar = "covar_sen",
       covar_label = "sensitivity covariates")
)

pathway_specs <- list(
  list(sheet = "Table S5", scope = "main", covar = "covar",
       title = "Table S5. List of enriched pathways (P < 0.05) from the Mummichog main analysis using the primary covariates"),
  list(sheet = "Table S6", scope = "main", covar = "covar_sen",
       title = "Table S6. List of enriched pathways (P < 0.05) from the Mummichog main analysis using the sensitivity covariates"),
  list(sheet = "Table S7", scope = "meta", covar = "covar",
       title = "Table S7. List of enriched pathways (P < 0.1) from the cross-cohort (SALSA × WHICAP) Mummichog meta-analysis using the primary covariates"),
  list(sheet = "Table S8", scope = "meta", covar = "covar_sen",
       title = "Table S8. List of enriched pathways (P < 0.1) from the cross-cohort (SALSA × WHICAP) Mummichog meta-analysis using the sensitivity covariates")
)

feature_pathway_specs <- list(
  list(sheet = "Table S9", covar = "covar",
       title = "Table S9. Feature-to-pathway mapping table for enriched pathways (P < 0.05) using the primary covariates"),
  list(sheet = "Table S10", covar = "covar_sen",
       title = "Table S10. Feature-to-pathway mapping table for enriched pathways (P < 0.05) using the sensitivity covariates")
)

# 8. Workbook assembly ------------------------------------------------------

wb <- openxlsx::createWorkbook()

# 8a. SALSA pooled FDR < 0.05
purrr::walk(salsa_concise_specs, function(spec) {
  message("Building SALSA concise ", spec$id, " (", spec$covar, ") ...")
  wide <- build_salsa_pooled_wide(spec$covar)

  sheet_s <- sprintf("Table %s", spec$id)
  title_s <- sprintf(
    "%s. MWAS-identified features (FDR < 0.05) across all platforms, composite exposures, and populations in SALSA using the %s",
    sheet_s, spec$covar_label)
  d_s <- if (is.null(wide) || nrow(wide) == 0)
           trim_salsa_stats(tibble::tibble(met = character(),
                                           logFC = numeric(),
                                           `95% CI` = character(),
                                           p = numeric(), FDR = numeric(),
                                           VIP = numeric(),
                                           platform = character(),
                                           exposure = character(),
                                           population = character()))
         else trim_salsa_stats(wide)
  message("  ", sheet_s, ": ", nrow(d_s), " rows")
  write_sheet(wb, sheet_s, title_s, d_s)

  sheet_a <- sprintf("Table %sa", spec$id)
  title_a <- sprintf(
    "%s. Annotated features from %s",
    sheet_a, sheet_s)
  d_a <- if (is.null(wide) || nrow(wide) == 0)
           tibble::tibble(met = character(), annotation = character(),
                          `chemical ID` = character(),
                          `confidence level` = character(),
                          `reference database` = character(),
                          platform = character(),
                          exposure = character(),
                          population = character())
         else trim_salsa_anno(wide)
  message("  ", sheet_a, ": ", nrow(d_a), " rows")
  write_sheet(wb, sheet_a, title_a, d_a)
})

# 8b. Meta pooled FE FDR < 0.05
purrr::walk(meta_concise_specs, function(spec) {
  message("Building Meta concise ", spec$id, " (", spec$covar, ") ...")
  wide <- build_meta_pooled_wide(spec$covar)

  sheet_s <- sprintf("Table %s", spec$id)
  title_s <- sprintf(
    "%s. Cross-cohort (SALSA × WHICAP) fixed-effect meta-analysis features (FE FDR < 0.05) across all platforms and composite exposures in the full cohort using the %s",
    sheet_s, spec$covar_label)
  d_s <- if (is.null(wide) || nrow(wide) == 0)
           trim_meta_stats(tibble::tibble(
             met = character(),
             `FE estimate` = numeric(), `FE 95% CI` = character(),
             `FE p` = numeric(), `FE FDR` = numeric(),
             `I^2 (%)` = numeric(), `Cochran Q p` = numeric(),
             platform = character(), exposure = character()))
         else trim_meta_stats(wide)
  message("  ", sheet_s, ": ", nrow(d_s), " rows")
  write_sheet(wb, sheet_s, title_s, d_s)

  sheet_a <- sprintf("Table %sa", spec$id)
  title_a <- sprintf(
    "%s. Annotated features from %s",
    sheet_a, sheet_s)
  d_a <- if (is.null(wide) || nrow(wide) == 0)
           tibble::tibble(met = character(), annotation = character(),
                          `chemical ID` = character(),
                          `confidence level` = character(),
                          `reference database` = character(),
                          platform = character(), exposure = character())
         else trim_meta_anno(wide)
  message("  ", sheet_a, ": ", nrow(d_a), " rows")
  write_sheet(wb, sheet_a, title_a, d_a)
})

# 8c. Pathway and feature-to-pathway sheets (carried over verbatim)
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

out_path <- here::here("Supplement", "Supplement tables_concise.xlsx")
openxlsx::saveWorkbook(wb, file = out_path, overwrite = TRUE)
message("Concise supplement workbook saved to: ", out_path)

#--------------------------------End of the code--------------------------------
