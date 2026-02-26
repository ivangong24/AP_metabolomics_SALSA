## ---------------------------
##
## Script name: 9-create_table1.R
## Purpose of script: To create Table 1 (demographic characteristics) for
##                    participants with metabolomics data
##
## Author: Yufan Gong
##
## Date Created: 2026-02-03
##
## Date Modified: 2026-02-03
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: This script creates a descriptive statistics table (Table 1)
##        for SALSA participants with metabolomics data.
##
##        Dependencies: Run scripts 1-3 before this script.
## ---------------------------

# Load required packages -----------------------------------------------------

library(tidyverse)
library(gtsummary)
library(labelled)
library(here)

# Source the functions file for table1() function
source(here::here("scripts", "1-functions.R"))

# Load cleaned data ----------------------------------------------------------

load(here::here("data", "processed", "salsa_clean.RData"))
load(here::here("data", "processed", "air_toxicants_exposure.RData"))
load(here::here("data", "links", "processed", "Sample_links.RData"))

# Identify participants with metabolomics data -------------------------------

## Get unique participant IDs from both platforms
participants_c18 <- sample_link_c18 |>
  dplyr::pull(rand_id) |>
  unique()

participants_hilic <- sample_link_hilic |>
  dplyr::pull(rand_id) |>
  unique()

## Participants with data from either platform
participants_metabolomics <- union(participants_c18, participants_hilic)

## Participants with data from both platforms
participants_both <- intersect(participants_c18, participants_hilic)

message("Participants with metabolomics data:")
message(paste0("  C18 only: ", length(setdiff(participants_c18, participants_hilic))))
message(paste0("  HILIC only: ", length(setdiff(participants_hilic, participants_c18))))
message(paste0("  Both platforms: ", length(participants_both)))
message(paste0("  Total unique: ", length(participants_metabolomics)))


# Create analysis dataset for Table 1 ----------------------------------------

## Merge SALSA data with exposure data and filter to metabolomics participants
table1_data <- salsa_clean |>
  dplyr::mutate(rand_id = as.character(rand_id)) |>
  dplyr::filter(rand_id %in% as.character(participants_metabolomics)) |>
  dplyr::left_join(
    air_toxicants_avg_ztrans |>
      dplyr::mutate(rand_id = as.character(rand_id)),
    by = "rand_id"
  ) |>
  dplyr::mutate(
    # Create metabolomics platform indicator
    has_c18 = rand_id %in% as.character(participants_c18),
    has_hilic = rand_id %in% as.character(participants_hilic),
    platform = case_when(
      has_c18 & has_hilic ~ "Both",
      has_c18 ~ "C18 only",
      has_hilic ~ "HILIC only",
      TRUE ~ "None"
    )
  )

message(paste0("\nFinal sample size for Table 1: ", nrow(table1_data)))


# Prepare variables for Table 1 ----------------------------------------------

## Set variable labels
table1_labeled <- table1_data |>
  dplyr::select(
    # Demographics
    blage, gender, edu_year, mh62,
    # Health outcomes
    demcind,
    # Exposure variables (select key ones for table)
    # exp_pm2.5, exp_no2, exp_nox, exp_o3,
    # exp_benzene, exp_butadiene,
    # exp_lead, exp_chromium, exp_nickel, exp_zinc,
    # Platform
    # platform
  ) |>
  labelled::set_variable_labels(
    blage = "Age at baseline (years)",
    gender = "Sex",
    edu_year = "Years of education",
    mh62 = "Smoking status",
    demcind = "Dementia/CIND during follow-up",
    # exp_pm2.5 = "PM2.5 (z-score)",
    # exp_no2 = "NO2 (z-score)",
    # exp_nox = "NOx (z-score)",
    # exp_o3 = "O3 (z-score)",
    # exp_benzene = "Benzene (z-score)",
    # exp_butadiene = "1,3-Butadiene (z-score)",
    # exp_lead = "Lead (z-score)",
    # exp_chromium = "Chromium (z-score)",
    # exp_nickel = "Nickel (z-score)",
    # exp_zinc = "Zinc (z-score)",
    # platform = "Metabolomics platform"
  )
  # labelled::set_value_labels(
  #   gender = c("Male" = 1, "Female" = 2),
  #   mh62 = c("Never smoker" = 1, "Former smoker" = 2, "Current smoker" = 3),
  #   demcind = c("No" = 0, "Yes" = 1)
  # ) |>
  # dplyr::mutate(
  #   platform = factor(platform, levels = c("Both", "C18 only", "HILIC only"))
  # ) |>
  # labelled::unlabelled()


# Create Table 1 - Overall ---------------------------------------------------

tbl1_overall <- table1_labeled |>
  # dplyr::select(-platform) |>
  gtsummary::tbl_summary(
    type = list(
      all_continuous() ~ "continuous2",
      all_categorical() ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ c("{mean} ({sd})", "{median} [{p25}, {p75}]"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "ifany",
    missing_text = "Missing"
  ) |>
  gtsummary::modify_header(label = "**Characteristic**") |>
  gtsummary::modify_caption("**Table 1. Characteristics of SALSA participants with metabolomics data (N = {N})**") |>
  gtsummary::bold_labels()


# Create Table 1 - By demcind status -----------------------------------------

tbl1_by_demcind <- table1_labeled |>
  gtsummary::tbl_summary(
    by = demcind,
    type = list(
      all_continuous() ~ "continuous2",
      c(gender, mh62) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ c("{mean} ({sd})", "{median} [{p25}, {p75}]"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "ifany",
    missing_text = "Missing"
  ) |>
  gtsummary::add_overall() |>
  # gtsummary::add_p(
  #   test = list(
  #     all_continuous() ~ "wilcox.test",
  #     all_categorical() ~ "chisq.test"
  #   )
  # ) |>
  gtsummary::modify_header(label = "**Characteristic**") |>
  gtsummary::modify_spanning_header(c("stat_1", "stat_2") ~ "**Dementia/CIND Status**") |>
  gtsummary::modify_caption("**Table 1. Characteristics by dementia/CIND status (N = {N})**") |>
  gtsummary::bold_labels()


# Convert tables using table1() function -------------------------------------

tbl1_overall_output <- table1(tbl1_overall)
tbl1_by_demcind_output <- table1(tbl1_by_demcind)



# Save tables ----------------------------------------------------------------

## Create output directory
dir.create(here::here("tables", "table 1"), showWarnings = FALSE)

## Save as Word documents (using flextable)
flextable::save_as_docx(
  tbl1_overall_output$flex,
  path = here::here("tables", "table 1", "table1_overall.docx")
)

flextable::save_as_docx(
  tbl1_by_demcind_output$flex,
  path = here::here("tables", "table 1", "table1_by_demcind.docx")
)


## Save as HTML for QMD integration
gtsummary::as_gt(tbl1_overall) |>
  gt::gtsave(here::here("tables", "table 1", "table1_overall.html"))

gtsummary::as_gt(tbl1_by_demcind) |>
  gt::gtsave(here::here("tables", "table 1", "table1_by_demcind.html"))


# Print tables to console ----------------------------------------------------

message("\n=== Table 1: Overall ===\n")
print(tbl1_overall)

message("\n=== Table 1: By Dementia/CIND Status ===\n")
print(tbl1_by_demcind)



# Save R objects for QMD -----------------------------------------------------

save(tbl1_overall, tbl1_by_demcind, 
     table1_data, table1_labeled,
     file = here::here("tables", "table 1", "table1_objects.RData"))

message("\nTable 1 created successfully!")
message("  - tables/table1_overall.docx")
message("  - tables/table1_by_demcind.docx")
message("  - tables/table1_objects.RData")

#--------------------------------End of the code--------------------------------
