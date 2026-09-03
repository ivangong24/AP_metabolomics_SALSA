## ---------------------------
##
## Script name: R5-annotation_revision.R
## Purpose of script: To attach compound annotations to the R1 revision MWAS
##                    results
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
## Notes: This is the revision analogue of scripts/5-annotation.R.
##
##        The annotation itself does NOT depend on the exposure. Only the
##        mixture index changed between the submitted analysis and the
##        revision, so the same feature tables, the same Emory in-house
##        library and the same xMSannotator Stage 5 output apply. Re-running
##        multilevelannotation() would take ~4 hours and reproduce identical
##        tables, so this script reuses the cleaned annotation from
##        data/metabolomics/annotation/annotation_cleaned_wide.RData and joins
##        it to the revision MWAS results.
##
##        Re-run scripts/5-annotation.R (Sections "xMSannotator" through
##        "Clean annotation") only if the feature tables themselves change.
##
##        Outputs -> revision_output/tables/mwas_results/... (_sig_annotated)
##                   revision_output/data/metabolomics/results/
##        Downstream: R7-visualization_revision.R
##
## ---------------------------

# Load required packages -----------------------------------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R5-annotation_revision.R")


# Load required datasets -----------------------------------------------------

load(here::here("data", "metabolomics",
                "annotation", "annotation_cleaned_long.RData"))

load(here::here("data", "metabolomics",
                "annotation", "annotation_cleaned_wide.RData"))

load(rev_here("data", "metabolomics", "results",
              "mwas_results_all_revision.RData"))

message("Annotated features available: c18 = ", nrow(annotation_c18_wide),
        ", hilic = ", nrow(annotation_hilic_wide))


# Identification confidence ---------------------------------------------------

## Reviewer 1 major comment 8: identification confidence has to be traceable.
## Every annotated feature carries a Schymanski (2014) level and the number of
## candidate compounds behind it, so a reader can see both how the annotation
## was made and how ambiguous it is. The rule lives in
## R1-revision_functions.R so that this table and the volcano labels in R7
## cannot drift apart.
##
## n_candidates counts DISTINCT compounds after folding away duplicate names
## and stereodescriptors -- "Ornithine; D-Ornithine" is one compound, not two,
## because accurate mass and retention time cannot separate them.
add_confidence_columns <- function(annot_wide) {
  annot_wide |>
    dplyr::mutate(
      n_candidates     = n_candidate_compounds(compound),
      confidence_level = annotation_confidence_level(reference, compound),
      .after = "compound"
    )
}

annotation_c18_wide   <- add_confidence_columns(annotation_c18_wide)
annotation_hilic_wide <- add_confidence_columns(annotation_hilic_wide)

message("\nIdentification confidence assigned:")
dplyr::bind_rows(
  annotation_c18_wide   |> dplyr::mutate(platform = "c18"),
  annotation_hilic_wide |> dplyr::mutate(platform = "hilic")
) |>
  dplyr::count(platform, confidence_level, reference,
               resolves_to_one = n_candidates == 1) |>
  print(n = 40)



# Annotation evidence table ---------------------------------------------------

## Reviewer 1 asks that every retained feature report its proposed ion and
## adduct, mass error in ppm, retention time, confidence level, supporting
## database, and whether MS/MS evidence is available; and separately objects
## that annotations such as a pentabrominated diphenyl ether are "not a species
## that would be annotated from a singly charged ESI feature".
##
## Both are answered from the xMSannotator Stage5 output, which carries the
## adduct, the theoretical m/z, the ppm error, the molecular formula and the
## isotope grouping. The cleaned annotation objects drop those columns, so they
## are read back here.

STAGE5_DIRS <- c(HMDB = "HMDB", KEGG = "KEGG", `LIPID MAPS` = "LipidMaps")

read_stage5 <- function(platform) {
  STAGE5_DIRS |>
    purrr::imap(function(dir_name, ref) {
      f <- here::here("annotation", platform, dir_name, "Stage5.csv")
      if (!file.exists(f)) {
        warning("No Stage5 output at ", f, call. = FALSE)
        return(NULL)
      }
      suppressWarnings(
        readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
      ) |>
        dplyr::transmute(
          chemical_id    = chemical_ID,
          mz, rt = time,
          match_category = MatchCategory,
          theoretical_mz = theoretical.mz,
          delta_ppm,
          formula        = Formula,
          adduct         = Adduct,
          isotope_group  = ISgroup,
          xms_confidence = Confidence,
          reference      = ref
        )
    }) |>
    purrr::list_rbind()
}

## Element count from a molecular formula. The negative lookahead stops "C"
## matching the C of "Cl" and "N" the N of "Na".
formula_element_count <- function(formula, element) {
  m <- stringr::str_match(formula, paste0(element, "(\\d*)(?![a-z])"))
  dplyr::case_when(
    is.na(m[, 1]) ~ 0L,
    m[, 2] == ""  ~ 1L,
    TRUE          ~ suppressWarnings(as.integer(m[, 2]))
  )
}

## Is this candidate a plausible assignment for a SINGLY CHARGED ESI feature?
##
## The test is halogen content. Chlorine and bromine have diagnostic isotope
## signatures -- one Cl puts ~32% of the intensity at M+2, one Br ~97% -- so a
## genuinely poly-halogenated compound is identifiable from its isotopologues.
## When a poly-halogenated formula is proposed on accurate mass alone, with no
## isotopic confirmation, the match is a mass coincidence rather than an
## annotation. The candidates this removes are exactly the ones the Reviewer
## names: brominated diphenyl ethers, 1,2-dichloroethane [M+H]+,
## chloroacetyl chloride, methoxyflurane -- volatile, reactive or non-ionising
## species that would not survive to be measured as a protonated serum
## metabolite in the first place.
##
## Single-chlorine compounds are NOT flagged: chlorinated endogenous and
## drug-derived metabolites are real and common, and a study of air toxicants
## has no business discarding halogenated species wholesale.
## Metals and alkaline earths in a molecular formula mean the database entry is
## a SALT, not the molecular species that would be observed. Calcium propionate
## (C6H10CaO4) is the Reviewer's own example: the ion measured from a serum
## extract would be propionate, and the calcium salt is a formulation of the
## compound rather than an analyte. Sodium and potassium are the same story and
## are also confusable with the M+Na / M+K adducts already modelled separately.
SALT_ELEMENTS <- c("Na", "K", "Ca", "Mg", "Fe", "Zn", "Cu", "Mn",
                   "Al", "Ba", "Sr", "Li", "Ag", "Hg", "Pb", "Pt", "Co")

implausible_esi_reason <- function(formula) {
  n_br <- formula_element_count(formula, "Br")
  n_i  <- formula_element_count(formula, "I")
  n_cl <- formula_element_count(formula, "Cl")

  salt_hits <- SALT_ELEMENTS |>
    purrr::map(~ formula_element_count(formula, .x) > 0) |>
    purrr::set_names(SALT_ELEMENTS)
  salt_which <- purrr::pmap_chr(salt_hits, function(...){
    hit <- c(...)
    if (any(hit)) paste(SALT_ELEMENTS[hit], collapse = ", ") else NA_character_
  })

  dplyr::case_when(
    is.na(formula) ~ NA_character_,
    !is.na(salt_which) ~ paste0("salt form (contains ", salt_which,
                                "); the measured ion would be the counter-ion"),
    n_br > 0 ~ paste0("contains Br x", n_br,
                      "; no isotopic confirmation of bromine"),
    n_i  > 0 ~ paste0("contains I x", n_i,
                      "; iodinated species not expected here"),
    n_cl >= 2 ~ paste0("contains Cl x", n_cl,
                       "; no isotopic confirmation of chlorine"),
    TRUE ~ NA_character_
  )
}

## MS/MS was not acquired. The Reviewer also asks for spectra and diagnostic
## fragments for the principal candidates; that cannot be produced from these
## data and the limitation is reported rather than worked around. The column is
## carried explicitly so the absence is visible in every table rather than
## inferred from silence.
MS2_EVIDENCE <- "Not acquired (MS1 full scan only)"

## Monoisotopic masses, for computing a theoretical m/z from a formula --------
##
## Needed because the in-house library stores its own m/z per entry, and those
## entries are not always exact: the taurine row is 126.022470 against a
## theoretical [M+H]+ of 126.021941, an error of 4.20 ppm in the library
## itself. Referencing the observed mass against the library entry would
## therefore report a mixture of instrument accuracy and library curation
## error. The mass error is computed against the THEORETICAL m/z of the
## proposed adduct, which is what xMSannotator does for its own matches and
## what a reader means by "mass error in ppm".
MONOISOTOPIC <- c(
  C = 12.0000000000,  H = 1.0078250321,  N = 14.0030740052,
  O = 15.9949146221,  P = 30.9737615100, S = 31.9720707000,
  F = 18.9984032000,  Cl = 34.9688527100, Br = 78.9183376000,
  I = 126.9044680000, Na = 22.9897692800, K = 38.9637069000,
  Se = 79.9165218000, Si = 27.9769265327, Ca = 39.9625912000,
  Mg = 23.9850419000, Fe = 55.9349421000, Zn = 63.9291466000,
  Cu = 62.9296011000, Mn = 54.9380496000, Co = 58.9332002000,
  Li = 7.0160040000,  Al = 26.9815384000, B = 11.0093055000
)

PROTON_MASS <- 1.0072764669

formula_monoisotopic_mass <- function(formula) {
  purrr::map_dbl(formula, function(f) {
    if (is.na(f) || f == "") return(NA_real_)
    m <- stringr::str_match_all(f, "([A-Z][a-z]?)(\\d*)")[[1]]
    if (nrow(m) == 0) return(NA_real_)
    el <- m[, 2]
    n  <- ifelse(m[, 3] == "", 1L, suppressWarnings(as.integer(m[, 3])))
    if (any(!el %in% names(MONOISOTOPIC))) return(NA_real_)
    sum(MONOISOTOPIC[el] * n)
  })
}

## Only the two ion forms the in-house library uses.
adduct_mz <- function(neutral_mass, adduct) {
  dplyr::case_when(
    adduct == "M+H" ~ neutral_mass + PROTON_MASS,
    adduct == "M-H" ~ neutral_mass - PROTON_MASS,
    TRUE            ~ NA_real_
  )
}

## Stage5 read once per platform, and reused: it supplies both the per-match
## evidence and the chemical_id -> formula lookup the in-house branch needs.
stage5_by_platform <- list(c18 = read_stage5("c18"),
                           hilic = read_stage5("hilic"))

formula_lookup <- stage5_by_platform |>
  purrr::list_rbind() |>
  dplyr::distinct(chemical_id, formula) |>
  dplyr::filter(!is.na(formula))

## Observed m/z at full precision. The feature id string rounds it to three
## decimals -- +/- 0.0005 Da, which is +/- 4 ppm at m/z 126 and therefore the
## same size as the numbers being reported -- so it must not be parsed out of
## the id.
observed_mz_by_platform <- list(c18 = xms_c18_long, hilic = xms_hilic_long) |>
  purrr::map(~ dplyr::distinct(.x, id, observed_mz = mz))

annotation_evidence <- list(
  list("c18", "hilic"),
  list(xms_c18_long, xms_hilic_long),
  list(inhouse_c18_matched, inhouse_hilic_matched)
) |>
  purrr::pmap(function(platform, xms_long, inhouse_matched) {

    xms_ev <- xms_long |>
      dplyr::left_join(stage5_by_platform[[platform]],
                       by = c("chemical_id", "reference"),
                       relationship = "many-to-many") |>
      ## Stage5 holds every candidate for every database; keep the row that is
      ## this feature's own match.
      dplyr::filter(abs(mz.x - mz.y) < 0.001, abs(rt.x - rt.y) < 1) |>
      dplyr::transmute(
        platform, id, mz = mz.x, rt = rt.x,
        candidate = match_chemical, chemical_id, database = reference,
        adduct, theoretical_mz, delta_ppm,
        delta_ppm_library = NA_real_,
        mass_reference = paste0("theoretical ", adduct),
        formula, match_category, isotope_group, xms_confidence
      )

    ## The in-house library stores one ion form per mode: its m/z values are
    ## [M+H]+ on HILIC-pos and [M-H]- on C18-neg. Verified against taurine
    ## (library 126.0225 against 125.0147 + 1.00728 = 126.0219) and
    ## 12-hydroxydodecanoic acid (library 215.1647 against 216.1725 - 1.00728
    ## = 215.1652). The adduct is therefore known, not missing, and is
    ## reported rather than left blank.
    ##
    ## NOTE the column order below: `theoretical_mz` and `delta_ppm` are
    ## computed from `library_mz`, NOT from `mz`. dplyr evaluates transmute()
    ## expressions in sequence, so a later reference to `mz` would resolve to
    ## the observed m/z assigned earlier in the same call and every ppm error
    ## would come out as exactly zero.
    ## Resolved OUTSIDE the transmute: `platform` becomes a column in the same
    ## call, so an `if` on it there sees a length-n vector, not the scalar.
    library_adduct <- if (platform == "hilic") "M+H" else "M-H"

    inhouse_ev <- inhouse_matched |>
      dplyr::filter(match_chemical != "") |>
      dplyr::rename(library_mz = mz) |>
      dplyr::left_join(observed_mz_by_platform[[platform]], by = "id") |>
      dplyr::left_join(formula_lookup, by = c("hmdbid" = "chemical_id")) |>
      dplyr::mutate(
        theo_mz = adduct_mz(formula_monoisotopic_mass(formula), library_adduct),
        ## Reference the observed mass against theory where the formula is
        ## known (95-96% of in-house matches), and fall back to the library
        ## entry otherwise -- recording which was used, so no reader has to
        ## guess what a ppm figure is relative to.
        mass_reference = dplyr::if_else(
          !is.na(theo_mz),
          paste0("theoretical ", library_adduct),
          "in-house library entry"),
        ref_mz = dplyr::coalesce(theo_mz, library_mz)
      ) |>
      dplyr::transmute(
        platform, id, mz = observed_mz, rt,
        candidate = match_chemical, chemical_id = hmdbid,
        database = "In House Library",
        adduct = library_adduct,
        theoretical_mz = ref_mz,
        delta_ppm = round(1e6 * (observed_mz - ref_mz) / ref_mz, 2),
        ## How far the library entry itself sits from theory. Reported because
        ## it is a property of the library, not of the measurement, and a few
        ## entries are several ppm out.
        delta_ppm_library = round(1e6 * (library_mz - theo_mz) / theo_mz, 2),
        mass_reference,
        formula,
        match_category = NA_character_,
        isotope_group = NA_character_,
        xms_confidence = NA_real_
      )

    dplyr::bind_rows(inhouse_ev, xms_ev)
  }) |>
  purrr::list_rbind() |>
  dplyr::mutate(
    implausible_reason = implausible_esi_reason(formula),
    plausible_esi      = is.na(implausible_reason),
    ms2_evidence       = MS2_EVIDENCE
  )

message("\nAnnotation evidence assembled: ", nrow(annotation_evidence),
        " feature-candidate rows across ",
        dplyr::n_distinct(annotation_evidence$id), " features")
message("Candidates implausible for a singly charged ESI feature: ",
        sum(!annotation_evidence$plausible_esi), " (",
        round(100 * mean(!annotation_evidence$plausible_esi), 2), "%), on ",
        dplyr::n_distinct(annotation_evidence$id[!annotation_evidence$plausible_esi]),
        " features")

## Per-feature roll-up, joined onto the annotation so the evidence travels with
## every result table rather than living only in a separate file.
evidence_by_feature <- annotation_evidence |>
  dplyr::group_by(id, database) |>
  dplyr::summarize(
    adduct        = paste(sort(unique(stats::na.omit(adduct))), collapse = "; "),
    delta_ppm     = round(mean(delta_ppm, na.rm = TRUE), 2),
    formula       = paste(sort(unique(stats::na.omit(formula))), collapse = "; "),
    isotope_group = paste(sort(unique(stats::na.omit(isotope_group))),
                          collapse = "; "),
    n_implausible = sum(!plausible_esi),
    n_candidates_ev = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    adduct        = dplyr::na_if(adduct, ""),
    formula       = dplyr::na_if(formula, ""),
    isotope_group = dplyr::na_if(isotope_group, ""),
    ## FALSE only when EVERY candidate for the feature is implausible: a
    ## feature with one credible candidate among several is still usable, it
    ## just has to be read with its candidate list.
    plausible_esi = n_implausible < n_candidates_ev,
    ms2_evidence  = MS2_EVIDENCE
  ) |>
  dplyr::select(-n_implausible, -n_candidates_ev)

attach_evidence <- function(annot_wide) {
  annot_wide |>
    dplyr::left_join(evidence_by_feature,
                     by = c("id" = "id", "reference" = "database"))
}

annotation_c18_wide   <- attach_evidence(annotation_c18_wide)
annotation_hilic_wide <- attach_evidence(annotation_hilic_wide)

rev_save_table(annotation_evidence, "annotation_evidence", "annotation")


# Link to MWAS results -------------------------------------------------------

list(
  list(annotation_c18_wide, annotation_hilic_wide),
  list(significant_results_list_c18, significant_results_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(annot_wide, sig_mwas_df_list, mode) {
    sig_mwas_df_list |>
      purrr::map(function(datalist) {
        datalist |>
          purrr::map(function(dflist) {
            dflist |>
              purrr::map(function(dfls){
                dfls |>
                  purrr::map(function(df){
                    df |>
                      add_percent_difference() |>
                      dplyr::left_join(annot_wide, by = c("met" = "id"))
                  })
              })
          })
      })
  }) |>
  purrr::set_names("mwas_annotated_list_c18", "mwas_annotated_list_hilic") |>
  list2env(.GlobalEnv)


# Annotate the full (unfiltered) result tables as well ------------------------

## The submitted pipeline only annotated the _sig tables. Reviewer 1 major
## comment 8 asks for identification confidence to be traceable, so the full
## tables carry the annotation columns here too.
list(
  list(annotation_c18_wide, annotation_hilic_wide),
  list(combined_results_list_c18, combined_results_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(annot_wide, mwas_df_list, mode) {
    mwas_df_list |>
      purrr::map(function(datalist) {
        datalist |>
          purrr::map(function(dflist) {
            dflist |>
              purrr::map(function(dfls){
                dfls |>
                  purrr::map(function(df){
                    df |>
                      add_percent_difference() |>
                      dplyr::left_join(annot_wide, by = c("met" = "id"))
                  })
              })
          })
      })
  }) |>
  purrr::set_names("mwas_full_annotated_list_c18",
                   "mwas_full_annotated_list_hilic") |>
  list2env(.GlobalEnv)


# Save annotated MWAS results ------------------------------------------------

list(
  list(mwas_annotated_list_c18, mwas_annotated_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(datalist, mode) {
    datalist |>
      purrr::imap(function(data_list, study){
        data_list |>
          purrr::imap(function(dflist, population){
            dflist |>
              purrr::imap(function(dfls, covar_name){
                rev_dir("tables", "mwas_results", study, population, covar_name)
                dfls |>
                  purrr::imap(function(df, exp_name) {
                    df |>
                      writexl::write_xlsx(
                        rev_here(
                          "tables", "mwas_results", study,
                          population, covar_name,
                          glue::glue("mwas_{mode}_{exp_name}_{study}_{population}_{covar_name}_sig_annotated.xlsx"))
                      )
                    message(paste0("MWAS results with annotation for ",
                                   exp_name, " ", mode, " in ",
                                   study, "_", population,
                                   " with covariate ", covar_name,
                                   " saved successfully!"))
                  })
              })
          })
      })
  })


# Summarize annotation coverage of the significant features -------------------

annotation_coverage <- list(
  list(mwas_annotated_list_c18, mwas_annotated_list_hilic),
  list("c18", "hilic")
) |>
  purrr::pmap(function(datalist, mode) {
    datalist |>
      purrr::imap(function(data_list, study){
        data_list |>
          purrr::imap(function(dflist, population){
            dflist |>
              purrr::imap(function(dfls, covar_name){
                dfls |>
                  purrr::imap(function(df, exp_name) {
                    tibble::tibble(
                      platform    = mode,
                      study       = study,
                      population  = population,
                      covar_set   = covar_name,
                      exposure    = exp_name,
                      label       = unname(rev_label(exp_name)),
                      n_features  = nrow(df),
                      n_fdr05     = sum(df$adj.P.Val < 0.05, na.rm = TRUE),
                      n_vip_gt2   = sum(df$VIP_comp1 > 2, na.rm = TRUE),
                      n_annotated = sum(!is.na(df$compound) & df$compound != "",
                                        na.rm = TRUE),
                      n_inhouse   = sum(df$reference == "In House Library",
                                        na.rm = TRUE),
                      n_multiple_match = sum(df$multiple_match, na.rm = TRUE),
                      ## Schymanski levels. Level 2 is empty by construction
                      ## (no MS/MS acquired) and levels 4-5 are not assigned
                      ## to annotated features -- see R1.
                      n_level1    = sum(df$confidence_level == 1, na.rm = TRUE),
                      n_level3    = sum(df$confidence_level == 3, na.rm = TRUE),
                      n_level1_fdr05 = sum(df$confidence_level == 1 &
                                             df$adj.P.Val < 0.05, na.rm = TRUE)
                    )
                  }) |>
                  purrr::list_rbind()
              }) |>
              purrr::list_rbind()
          }) |>
          purrr::list_rbind()
      }) |>
      purrr::list_rbind()
  }) |>
  purrr::list_rbind()

rev_save_table(annotation_coverage, "annotation_coverage", "annotation")
print(annotation_coverage |>
        dplyr::filter(population == "all", covar_set == "covar"), n = 40)


# Save R objects for downstream analysis -------------------------------------

save(mwas_annotated_list_c18, mwas_annotated_list_hilic,
     mwas_full_annotated_list_c18, mwas_full_annotated_list_hilic,
     annotation_coverage,
     file = rev_here("data", "metabolomics", "results",
                     "mwas_annotation_revision.RData"))

message("\nAnnotation completed!")
message("Results saved to:")
message("  - ", rev_here("tables", "mwas_results"), " (_sig_annotated.xlsx)")
message("  - ", rev_here("tables", "annotation"))

#--------------------------------End of the code--------------------------------
