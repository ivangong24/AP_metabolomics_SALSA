## ---------------------------
##
## Script name: R2-exposure_windows_revision.R
##
## Purpose of script:
##        Reviewer 1 major comment 13 -- exposure window sensitivity.
##
##        The submitted analysis averages each pollutant over the FIVE years
##        preceding the blood draw (scripts/3-clean_data.R:434 is the single
##        line that sets it). Plasma is a comparatively short-term readout, and
##        the Reviewer asks whether the findings depend on that choice. This
##        script rebuilds the exposure averages over a ONE-year and a TEN-year
##        window and derives an unsupervised PC1 index for each, so the primary
##        MWAS can be re-run against them as a sensitivity analysis.
##
##        Only the PCA index is rebuilt. It is the primary exposure of the
##        revision, and it is deterministic and unsupervised, so a window
##        sensitivity built on it isolates the window rather than confounding
##        it with re-drawn cross-fitting folds or a re-fitted outcome model.
##
## Author: Yufan Gong
##
## Notes:
##
##        DATA COVERAGE. The LUR surfaces (benzene, 1,3-butadiene, chromium,
##        lead, nickel, NO2, PM2.5) cover 1989-2007 and support every window.
##        CALINE4 NOx exists for FIVE model years only, 1998-2002, so its
##        coverage depends sharply on the window:
##
##          window   specimens with no NOx model year in range
##            1 yr   1,442 / 1,546  (93%)
##            3 yr     701 / 1,546  (45%)
##            5 yr        16 / 1,546 ( 1%)   <- the submitted analysis
##           10 yr        16 / 1,546 ( 1%)
##
##        The 1-year and 3-year windows therefore CANNOT carry NOx: those
##        specimens would take the NA -> 0 fill of 3-clean_data.R:447 and enter
##        the mixture as a near-constant zero, which would then propagate into
##        the PCA loadings. Both are built from the SEVEN remaining pollutants
##        and must be reported as seven-component indices, not as the
##        eight-component one under a different window. The 3-year shortfall is
##        concentrated in the 2006 and 2007 draws (685 specimens), whose
##        3-year windows fall entirely after the last CALINE4 model year.
##
##        Which pollutants a window carries is decided from the data by a
##        coverage threshold, not hardcoded -- see WINDOW_MIN_COVERAGE below.
##
##        The 10-year window carries all eight, and is in fact BETTER
##        conditioned for NOx than the 5-year primary: every draw from 2004 on
##        averages the full 1998-2002 series, where the 5-year window degrades
##        to two model years for a 2006 draw and one for a 2007 draw.
##
##        MISSING-DATA RULE. The NA -> 0 fill of the primary pipeline is kept
##        deliberately, so that a difference between windows is a difference in
##        the window and not in the missing-data rule. It bites only the 16
##        specimens drawn in 1998 (NOx, 10-year window). The coverage table
##        this script writes makes every such cell visible.
##
## ---------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "2-load_data.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R2-exposure_windows_revision.R")

load(here::here("data", "processed", "air_toxicants_exposure.RData"))


# The specimen list ----------------------------------------------------------

## One row per (participant, blood draw). Taken from the primary pipeline's
## own output so the window analyses cover exactly the specimens the 5-year
## analysis covers, with no re-derivation of the cohort.
specimens <- air_toxicants_avg_list[["total"]][["all"]] |>
  dplyr::select(rand_id, blood_date, wave, batch) |>
  dplyr::distinct()

message("Specimens carried into the window analysis: ", nrow(specimens),
        " from ", dplyr::n_distinct(specimens$rand_id), " participants")


# Rebuild the year-level exposure table --------------------------------------

## Mirrors scripts/3-clean_data.R lines 378-421 exactly: the same objects, the
## same winsorization, the same long shape. Kept as a literal mirror rather
## than refactored into a shared helper, because the point of this script is
## that the ONLY thing differing from the primary analysis is the window.
col_common <- quote_all(rand_id, date, .source_dir)

exp_data_names_new <- c(exp_data_names, "nox")

exp_pattern  <- stringr::regex(stringr::str_c(exp_data_names_new,
                                              collapse = "|"),
                               ignore_case = TRUE)
drop_pattern <- stringr::regex("caline|final", ignore_case = TRUE)

exp_obj_names <- ls(envir = .GlobalEnv, all.names = TRUE) |>
  purrr::keep(~ stringr::str_detect(.x, exp_pattern)) |>
  purrr::discard(~ stringr::str_detect(.x, drop_pattern)) |>
  purrr::keep(~ exists(.x, envir = .GlobalEnv, inherits = TRUE)) |>
  sort()

exp_long <- mget(exp_obj_names, envir = .GlobalEnv, inherits = TRUE) |>
  purrr::imap(function(data, nm) {
    data |>
      dplyr::rename(
        value    = dplyr::all_of(setdiff(names(data), col_common)),
        toxicant = .source_dir
      ) |>
      dplyr::mutate(value = extreme_remove_percentile_win(value))
  }) |>
  purrr::list_rbind()

## CALINE4 NOx, long, as scripts/3-clean_data.R:193 builds it.
##
## The winsorization matters and is easy to miss. In 3-clean_data.R the `nox`
## object is created before exp_long is assembled, and "nox" is one of the
## names in exp_pattern while `drop_pattern` only removes "caline|final" -- so
## `nox` is picked up by mget() along with the LUR surfaces and goes through
## extreme_remove_percentile_win() with them. Building it separately here, as
## this script must, means applying the same winsorization by hand. Without
## it the recomputed 5-year NOx reproduces the stored column at r = 0.998
## rather than 1.000, which the check at the end of this script flags.
nox_long <- caline1789_nox_1998_2002 |>
  dplyr::select(-unique_id) |>
  tidyr::pivot_longer(cols = dplyr::starts_with("nox_"), names_to = "date",
                      names_prefix = "nox_", values_to = "value") |>
  dplyr::mutate(toxicant = "nox",
                date  = lubridate::ymd(stringr::str_c(date, "-01-01")),
                value = extreme_remove_percentile_win(value)) |>
  dplyr::select(rand_id, date, toxicant, value)

exp_long_all <- dplyr::bind_rows(
  exp_long |> dplyr::select(rand_id, date, toxicant, value),
  nox_long
) |>
  ## sdir_merge() sets .source_dir to the raw directory basename, so the
  ## toxicant values arrive capitalised ("Benzene", "PM2.5", "NO2").
  ## 3-clean_data.R only lowercases them after pivot_wider, via
  ## rename_all(str_to_lower); this script matches on the value rather than on
  ## a column name, so it has to lowercase here instead.
  dplyr::mutate(toxicant = stringr::str_to_lower(toxicant))

message("Daily exposure rows: ", nrow(exp_long_all), " across ",
        dplyr::n_distinct(exp_long_all$toxicant), " toxicants, ",
        min(lubridate::year(exp_long_all$date)), "-",
        max(lubridate::year(exp_long_all$date)))

## Collapse to one value per participant-toxicant-year BEFORE joining the
## specimens. The LUR surfaces are supplied per day ("AllDays" files), so the
## raw table is ~37M rows and a many-to-many join against 1,546 specimens is
## needlessly enormous. This is exactly the first stage of the primary
## pipeline's two-stage mean, and it does not depend on the blood draw -- the
## rows for a given (rand_id, toxicant, year) are the same set whichever
## specimen they are later matched to -- so pre-aggregating is arithmetically
## identical to grouping by blood_date as well, and vastly cheaper.
exp_year <- exp_long_all |>
  dplyr::mutate(year = lubridate::year(date)) |>
  dplyr::group_by(rand_id, toxicant, year) |>
  dplyr::summarize(value = mean(value, na.rm = TRUE), .groups = "drop")

message("Participant-toxicant-year rows: ", nrow(exp_year))


# Windowed averages ----------------------------------------------------------

WINDOWS <- c(w1 = 1, w3 = 3, w10 = 10)

## The primary analysis's own window, recomputed here so the diagnostics below
## compare like with like rather than against the stored 5-year columns.
WINDOW_REFERENCE <- c(w5 = 5)

ALL_POLLUTANTS <- stringr::str_remove(qgcomp_pollutants, "^exp_")


## Which pollutants each window can carry ------------------------------------
##
## Decided from the data rather than hardcoded, because the answer depends on
## the window and getting it wrong is silent: a pollutant with no model year in
## range takes the pipeline's NA -> 0 fill and enters the mixture as a
## near-constant zero, which then propagates into the PC1 loadings.
##
## A pollutant is carried only if at least WINDOW_MIN_COVERAGE percent of
## specimens have at least one exposure year inside the window. Measured on
## these data the shares are far from the threshold in both directions, so the
## rule is not sensitive to where exactly it is set -- NOx is at 7% for the
## 1-year window and 55% for the 3-year, against 99% for the 5- and 10-year,
## and every LUR surface is at 100% throughout. Any threshold between 56 and
## 98 gives the same answer.
WINDOW_MIN_COVERAGE <- 80

specimen_coverage <- function(w) {
  exp_year |>
    dplyr::right_join(specimens, by = "rand_id",
                      relationship = "many-to-many") |>
    dplyr::mutate(blood_year = lubridate::year(blood_date)) |>
    dplyr::filter(year >= blood_year - w, year < blood_year) |>
    dplyr::distinct(rand_id, blood_date, toxicant) |>
    dplyr::count(toxicant, name = "n_with_data") |>
    dplyr::mutate(pct_with_data = round(100 * n_with_data / nrow(specimens), 1))
}

window_pollutant_coverage <- c(WINDOWS, WINDOW_REFERENCE) |>
  purrr::imap(function(w, tag){
    specimen_coverage(w) |>
      dplyr::filter(toxicant %in% ALL_POLLUTANTS) |>
      dplyr::right_join(tibble::tibble(toxicant = ALL_POLLUTANTS),
                        by = "toxicant") |>
      dplyr::mutate(
        n_with_data   = dplyr::coalesce(n_with_data, 0L),
        pct_with_data = dplyr::coalesce(pct_with_data, 0),
        window        = tag,
        window_years  = w,
        carried       = pct_with_data >= WINDOW_MIN_COVERAGE
      )
  }) |>
  purrr::list_rbind() |>
  dplyr::arrange(window_years, toxicant)

message("\nPollutant coverage by window (carried if >= ",
        WINDOW_MIN_COVERAGE, "% of specimens have an exposure year):")
print(window_pollutant_coverage |>
        dplyr::select(window, window_years, toxicant, pct_with_data, carried),
      n = 40)
rev_save_table(window_pollutant_coverage, "window_pollutant_coverage",
               "windows")

dropped <- window_pollutant_coverage |> dplyr::filter(!carried)
if (nrow(dropped) > 0) {
  message("\nDROPPED from their window's index for insufficient coverage:")
  dropped |>
    purrr::pwalk(function(toxicant, window, window_years, pct_with_data, ...){
      message("  ", window, " (", window_years, "-year): ", toxicant,
              " -- only ", pct_with_data, "% of specimens have any data")
    })
  message("Indices built on a reduced pollutant set MUST be reported as such ",
          "and not compared component-for-component with the full index.")
}

window_pollutants <- function(w) {
  window_pollutant_coverage |>
    dplyr::filter(window_years == w, carried) |>
    dplyr::pull(toxicant)
}

average_window <- function(w) {
  message("  window: ", w, " year(s); pollutants: ",
          paste(window_pollutants(w), collapse = ", "))

  exp_year |>
    dplyr::filter(toxicant %in% window_pollutants(w)) |>
    dplyr::right_join(specimens, by = "rand_id",
                      relationship = "many-to-many") |>
    dplyr::mutate(blood_year = lubridate::year(blood_date)) |>
    dplyr::filter(year >= blood_year - w, year < blood_year) |>
    ## second stage of the two-stage mean: across the years of the window, so
    ## a year with more daily records does not dominate the average.
    dplyr::group_by(rand_id, blood_date, wave, toxicant) |>
    dplyr::summarize(avg_exp  = mean(value, na.rm = TRUE),
                     n_years  = dplyr::n(),
                     .groups  = "drop")
}

window_long <- c(WINDOWS, WINDOW_REFERENCE) |>
  purrr::imap(function(w, tag){
    average_window(w) |> dplyr::mutate(window = tag, window_years = w)
  }) |>
  purrr::list_rbind()


# Coverage -------------------------------------------------------------------

## How many of the window's years each specimen actually had data for. This is
## what makes the NOx problem, and the 1998-draw shortfall, visible rather than
## hidden behind the zero fill.
window_coverage <- window_long |>
  dplyr::mutate(blood_year = lubridate::year(blood_date)) |>
  dplyr::group_by(window, window_years, toxicant, blood_year) |>
  dplyr::summarize(n_specimens   = dplyr::n(),
                   mean_n_years  = round(mean(n_years), 2),
                   min_n_years   = min(n_years),
                   .groups = "drop") |>
  dplyr::arrange(window, toxicant, blood_year)

rev_save_table(window_coverage, "window_coverage", "windows")

## Specimens that got NO year at all for a toxicant never appear in
## window_long, so they have to be counted against the full specimen grid.
window_missing <- c(WINDOWS, WINDOW_REFERENCE) |>
  purrr::imap(function(w, tag){
    tidyr::expand_grid(
      specimens |> dplyr::select(rand_id, blood_date),
      toxicant = window_pollutants(w)
    ) |>
      dplyr::anti_join(
        window_long |> dplyr::filter(window == tag) |>
          dplyr::select(rand_id, blood_date, toxicant),
        by = c("rand_id", "blood_date", "toxicant")
      ) |>
      dplyr::count(toxicant, name = "n_specimens_no_data") |>
      dplyr::mutate(window = tag, window_years = w,
                    pct = round(100 * n_specimens_no_data / nrow(specimens), 1))
  }) |>
  purrr::list_rbind()

if (nrow(window_missing) > 0) {
  message("\nSpecimens with NO exposure year in the window (these take the ",
          "NA -> 0 fill):")
  print(window_missing, n = 40)
} else {
  message("\nEvery specimen has at least one exposure year in every window.")
}
rev_save_table(window_missing, "window_missing", "windows")


# Wide exposure frames, one per window ---------------------------------------

window_wide <- c(WINDOWS, WINDOW_REFERENCE) |>
  purrr::imap(function(w, tag){
    wide <- window_long |>
      dplyr::filter(window == tag) |>
      dplyr::select(rand_id, blood_date, wave, toxicant, avg_exp) |>
      tidyr::pivot_wider(names_from = toxicant, values_from = avg_exp,
                         names_prefix = "exp_")

    ## Same NA -> 0 fill as the primary pipeline, deliberately (see header).
    specimens |>
      dplyr::select(rand_id, blood_date, wave) |>
      dplyr::left_join(wide, by = c("rand_id", "blood_date", "wave")) |>
      dplyr::mutate(dplyr::across(dplyr::starts_with("exp_"),
                                  ~ dplyr::if_else(is.na(.x), 0, .x)))
  })


# PC1 index per window -------------------------------------------------------

message("\n=== Unsupervised PC1 index per exposure window ===")

window_pca <- list(WINDOWS, names(WINDOWS)) |>
  purrr::pmap(function(w, tag){
    vars <- paste0("exp_", window_pollutants(w))
    message("  ", tag, " (", w, "-year, ", length(vars), " pollutants)")
    make_pca_index(window_wide[[tag]], vars,
                   name = paste0("comp_pca_all_", tag))
  }) |>
  purrr::set_names(names(WINDOWS))

window_pca_df <- window_pca |>
  purrr::map("composites") |>
  purrr::reduce(dplyr::full_join, by = c("rand_id", "blood_date"))

window_pca_loadings <- window_pca |>
  purrr::imap(function(obj, tag){
    tibble::tibble(
      window       = tag,
      window_years = unname(WINDOWS[[tag]]),
      n_pollutants = length(obj$vars),
      pollutant    = names(obj$loadings) |> stringr::str_remove("^exp_"),
      pc1_loading  = round(unname(obj$loadings), 4),
      pc1_pve      = round(unname(obj$pve[1]), 4)
    )
  }) |>
  purrr::list_rbind()

rev_save_table(window_pca_loadings, "window_pca_loadings", "windows")
print(window_pca_loadings, n = 40)


# Agreement between windows --------------------------------------------------

## Reviewer 1 comment 13 asks for the correlation between windows. Two things
## are reported: the pollutant-level correlation of the averaged exposures, and
## the correlation of the resulting PC1 indices with the 5-year primary.

## Every pairwise comparison, derived from the window list rather than named
## one at a time -- an earlier version hardcoded w1/w5/w10 and silently omitted
## w3 when it was added.
window_tags  <- names(c(WINDOWS, WINDOW_REFERENCE))
window_pairs <- utils::combn(window_tags, 2, simplify = FALSE)

pair_label <- function(pr) {
  yrs <- c(WINDOWS, WINDOW_REFERENCE)
  paste0(yrs[[pr[1]]], "yr vs ", yrs[[pr[2]]], "yr")
}

wide_pollutant <- window_long |>
  dplyr::select(rand_id, blood_date, toxicant, window, avg_exp) |>
  tidyr::pivot_wider(names_from = window, values_from = avg_exp)

pollutant_window_cor <- window_pairs |>
  purrr::map(function(pr){
    wide_pollutant |>
      dplyr::group_by(toxicant) |>
      dplyr::summarize(
        comparison = pair_label(pr),
        n       = sum(!is.na(.data[[pr[1]]]) & !is.na(.data[[pr[2]]])),
        pearson = round(stats::cor(.data[[pr[1]]], .data[[pr[2]]],
                                   use = "pairwise.complete.obs"), 4),
        .groups = "drop")
  }) |>
  purrr::list_rbind() |>
  tidyr::pivot_wider(id_cols = toxicant, names_from = comparison,
                     values_from = pearson)

message("\nPollutant-level correlation between windows (Pearson):")
print(pollutant_window_cor, n = 20, width = Inf)
rev_save_table(pollutant_window_cor, "window_pollutant_correlations", "windows")

## The 5-year PC1 index, as R3 builds it, for the index-level comparison.
primary_pca <- make_pca_index(window_wide[["w5"]],
                              paste0("exp_", window_pollutants(5)),
                              name = "comp_pca_all_w5")

index_wide <- primary_pca$composites |>
  dplyr::inner_join(window_pca_df, by = c("rand_id", "blood_date"))

index_window_cor <- window_pairs |>
  purrr::map(function(pr){
    a <- index_wide[[paste0("comp_pca_all_", pr[1])]]
    b <- index_wide[[paste0("comp_pca_all_", pr[2])]]
    tibble::tibble(
      comparison = pair_label(pr),
      n          = sum(!is.na(a) & !is.na(b)),
      pearson    = round(stats::cor(a, b, use = "pairwise.complete.obs"), 4),
      spearman   = round(stats::cor(a, b, method = "spearman",
                                    use = "pairwise.complete.obs"), 4)
    )
  }) |>
  purrr::list_rbind()

message("\nPC1 index correlation between windows:")
print(index_window_cor, n = 20)
rev_save_table(index_window_cor, "window_index_correlations", "windows")

## Sanity check that the recomputed 5-year index reproduces the stored one.
## They are built from the same data by the same code, so anything other than
## a near-perfect correlation means this script has drifted from
## 3-clean_data.R and the window comparison cannot be trusted.
stored_5yr <- air_toxicants_avg_list[["total"]][["all"]] |>
  dplyr::select(rand_id, blood_date,
                dplyr::all_of(paste0("exp_", window_pollutants(5)))) |>
  dplyr::distinct()

recheck <- window_wide[["w5"]] |>
  dplyr::inner_join(stored_5yr, by = c("rand_id", "blood_date"),
                    suffix = c("_new", "_stored"))

recheck_cor <- window_pollutants(5) |>
  purrr::set_names() |>
  purrr::map_dbl(function(p){
    stats::cor(recheck[[paste0("exp_", p, "_new")]],
               recheck[[paste0("exp_", p, "_stored")]],
               use = "pairwise.complete.obs")
  })

message("\n5-year reproduction check (recomputed vs stored, should be 1.000):")
message("  ", paste(names(recheck_cor), round(recheck_cor, 5),
                    sep = " = ", collapse = ", "))
if (any(recheck_cor < 0.999, na.rm = TRUE)) {
  warning("The recomputed 5-year exposures do not reproduce the stored ones. ",
          "This script has drifted from scripts/3-clean_data.R; the window ",
          "comparison is not interpretable until that is resolved.",
          call. = FALSE)
}

rev_save_table(
  tibble::tibble(pollutant = names(recheck_cor),
                 r_recomputed_vs_stored = round(unname(recheck_cor), 5)),
  "window_5yr_reproduction_check", "windows")


# The windowed pollutant columns, for the rest of the revision ----------------

## The PC1 index was the only thing R3 and R4 needed while the window arm was a
## sensitivity analysis on the primary exposure alone. It now carries the whole
## exposure set (MWAS_WINDOWS in R1) -- cross-fitted WQS and QGcomp weights are
## re-derived on each window's pollutant matrix, the feature-wise QGcomp
## contrast sums that window's quartile scores, and each pollutant is an
## exposure in its own right -- so the averaged pollutant columns themselves
## have to travel downstream, not just the index built from them.
##
## Named with the window tag as a plain suffix on the raw column
## (exp_benzene_w3), which is the convention every derived name appends to:
## exp_benzene_w3_iqr, exp_benzene_w3_q.
window_exposure_df <- MWAS_WINDOWS |>
  purrr::map(function(tag){
    vars <- paste0("exp_", window_pollutants(WINDOWS[[tag]]))
    window_wide[[tag]] |>
      dplyr::select(rand_id, blood_date, dplyr::all_of(vars)) |>
      dplyr::rename_with(~ paste0(.x, "_", tag), dplyr::all_of(vars))
  }) |>
  purrr::reduce(dplyr::full_join, by = c("rand_id", "blood_date"))

message("\nWindowed pollutant columns carried downstream (",
        nrow(window_exposure_df), " specimens):")
message("  ", paste(setdiff(names(window_exposure_df),
                            c("rand_id", "blood_date")), collapse = ", "))

## A window that lost a pollutant to the coverage rule must lose it here too,
## or the MWAS silently fits a near-constant zero column.
purrr::walk(MWAS_WINDOWS, function(tag){
  expected <- paste0(pollutants_for_window(tag), "_iqr")
  present  <- paste0(intersect(names(window_exposure_df),
                               pollutants_for_window(tag)), "_iqr")
  if (!setequal(expected, present)) {
    warning("Window ", tag, " carries ", length(present), " pollutant ",
            "columns but R1 expects ", length(expected),
            ". WINDOW_POLLUTANTS in R1 is read from this script's own ",
            "loadings table, so this means the two ran against different ",
            "coverage decisions.", call. = FALSE)
  }
})


# Save -----------------------------------------------------------------------

rev_dir("data", "processed")

save(window_pca_df, window_pca, window_pca_loadings, window_wide,
     window_exposure_df,
     window_coverage, window_missing, pollutant_window_cor, index_window_cor,
     WINDOWS,
     file = rev_here("data", "processed", "exposure_windows_revision.RData"))

message("\nExposure window construction completed!")
message("  data  -> ", rev_here("data", "processed",
                                "exposure_windows_revision.RData"))
message("  tables-> ", rev_here("tables", "windows"))

#--------------------------------End of the code--------------------------------
