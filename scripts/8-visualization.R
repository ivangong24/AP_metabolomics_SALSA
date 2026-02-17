## ---------------------------
##
## Script name: 7-visualization.R
## Purpose of script: To create visualizations for MWAS results
##                    including volcano plots, heatmaps, and scatter plots
##
## Author: Yufan Gong
##
## Date Created: 2026-01-29
##
## Date Modified: 2026-01-29
##
## Copyright (c) Yufan Gong, 2026
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: This script creates publication-quality visualizations for
##        the air toxicants - metabolomics association study.
##
##        Dependencies: Run scripts 1-6 before this script.
## ---------------------------

# Load required packages -----------------------------------------------------

library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(ggnewscale)

# Load MWAS results ----------------------------------------------------------

load(here::here("data", "metabolomics", "results", "mwas_results_all.RData"))


load(here::here("data", "metabolomics", 
                "processed", "mwas_annotation.RData"))


# Create output directory
dir.create(here::here("figures", "mwas"), 
           showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# SECTION 1: VOLCANO PLOTS
# =============================================================================

# Function to create volcano plot --------------------------------------------

create_volcano_plot <- function(mwas_result, vip_result, annotation_result, 
                                exposure_name,
                                 column_type = "C18",
                                 fdr_threshold = 0.05,
                                 vip_threshold = 2,
                                 n_labels = 10) {

  # Prepare data for plotting
  plot_data <- mwas_result |>
    tibble::rownames_to_column("met") |>
    dplyr::left_join(
      vip_result |>
        tibble::rownames_to_column("met") |>
        dplyr::select(met, VIP = comp1),
      by = "met"
    ) |>
    dplyr::left_join(
      annotation_result |>
        dplyr::select(met, chemical_id, compound, multiple_match, reference),
      by = "met"
    ) |>
    dplyr::group_by(met) |>
    dplyr::slice_head(n = 1) |>  # In case of multiple matches
    dplyr::ungroup() |>
    dplyr::mutate(
      neg_log10_p = -log10(P.Value),
      significant = case_when(
        adj.P.Val < fdr_threshold & VIP > vip_threshold ~ "FDR < 0.05 & VIP>2",
        adj.P.Val < fdr_threshold & logFC ~ "FDR < 0.05 only",
        VIP > vip_threshold ~ "VIP>2 only",
        TRUE ~ "NS"
      )
    )

  # Get top metabolites for labeling (by VIP or p-value)
  top_mets <- plot_data |>
    dplyr::filter(
      is.finite(VIP), !is.na(compound), 
      multiple_match == FALSE | reference == "In House Library", compound != "",
      adj.P.Val < fdr_threshold | VIP > vip_threshold) |>
    # dplyr::arrange(P.Value) |>
    # dplyr::slice_head(n = n_labels) |>
    dplyr::slice_max(order_by = VIP, n = n_labels, with_ties = FALSE) |>
    dplyr::pull(met)

  plot_data <- plot_data |>
    dplyr::mutate(label = ifelse(met %in% top_mets, compound, ""))

  # Create volcano plot
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_p)) +
    # Non-significant points
    geom_point(
      data = . %>% filter(significant == "NS"),
      color = "grey70",
      alpha = 0.5,
      size = 1.5
    ) +
    # Significant points
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant),
      alpha = 0.7,
      size = 2.5
    ) +
    scale_color_manual(
      values = c(
        "FDR < 0.05 & VIP>2" = "#FB9A99",
        "FDR < 0.05 only" = "#FDB462",
        "VIP>2 only" = "#984EA3"
      ),
      name = "Significance"
    ) +
    # Add labels
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50"
    ) +
    # Reference lines
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey40"
    ) +
    # Labels
    labs(
      title = paste0(column_type, " - ", gsub("exp_", "", exposure_name)),
      x = "Log Fold Change",
      y = expression(-log[10](P-value))
    ) +
    # Theme
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    )

  return(p)
}


# Create volcano plots for all exposures -------------------------------------

exposure_vars <- names(mwas_results_c18)

## C18 volcano plots
volcano_c18 <- exposure_vars |>
  purrr::set_names() |>
  purrr::map(function(exp) {
    create_volcano_plot(
      mwas_results_c18[[exp]],
      vip_c18[[exp]],
      mwas_c18_annotated[[exp]],
      exposure_name = exp,
      column_type = "C18"
    )
  })

## HILIC volcano plots
volcano_hilic <- exposure_vars |>
  purrr::set_names() |>
  purrr::map(function(exp) {
    create_volcano_plot(
      mwas_results_hilic[[exp]],
      vip_hilic[[exp]],
      mwas_hilic_annotated[[exp]],
      exposure_name = exp,
      column_type = "HILIC"
    )
  })


# Save individual volcano plots ----------------------------------------------

purrr::iwalk(volcano_c18, function(p, exp) {
  ggsave(
    filename = here::here("figures", "mwas",
                          glue::glue("volcano_c18_{exp}.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
})

purrr::iwalk(volcano_hilic, function(p, exp) {
  ggsave(
    filename = here::here("figures", "mwas",
                          glue::glue("volcano_hilic_{exp}.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
})


# =============================================================================
# SECTION 2: VIP VS LOGFC SCATTER PLOTS
# =============================================================================

# Function to create VIP vs logFC scatter plot -------------------------------

create_vip_scatter <- function(mwas_result, vip_result, annotation_result,
                               exposure_name,
                                column_type = "C18",
                                vip_threshold = 2,
                                n_labels = 10) {

  # Prepare data
  plot_data <- mwas_result |>
    tibble::rownames_to_column("met") |>
    dplyr::left_join(
      vip_result |>
        tibble::rownames_to_column("met") |>
        dplyr::select(met, VIP = comp1),
      by = "met"
    ) |>
    dplyr::left_join(
      annotation_result |>
        dplyr::select(met, chemical_id, compound, multiple_match, reference),
      by = "met"
    ) |>
    dplyr::group_by(met) |>
    dplyr::slice_head(n = 1) |>  # In case of multiple matches
    dplyr::ungroup() |>
    dplyr::mutate(
      category = case_when(
        VIP > vip_threshold & logFC > 0 ~ "VIP>2 & Positive",
        VIP > vip_threshold & logFC < 0 ~ "VIP>2 & Negative",
        TRUE ~ "VIP<=2"
      )
    )

  # Top metabolites for labeling
  top_mets <- plot_data |>
    dplyr::filter(
      is.finite(VIP), !is.na(compound), 
      multiple_match == FALSE | reference == "In House Library",
      compound != "",
      VIP > vip_threshold) |>
    # dplyr::arrange(P.Value) |>
    # dplyr::slice_head(n = n_labels) |>
    dplyr::slice_max(order_by = VIP, n = n_labels, with_ties = FALSE) |>
    dplyr::pull(met)
  
  plot_data <- plot_data |>
    dplyr::mutate(label = ifelse(met %in% top_mets, compound, ""))
  


  # Create plot
  p <- ggplot(plot_data, aes(x = logFC, y = VIP)) +
    geom_point(
      data = . %>% filter(category == "VIP<=2"),
      color = "grey70",
      alpha = 0.5,
      size = 1.5
    ) +
    geom_point(
      data = . %>% filter(category != "VIP<=2"),
      aes(color = category),
      alpha = 0.7,
      size = 2.5
    ) +
    scale_color_manual(
      values = c(
        "VIP>2 & Positive" = "#E41A1C",
        "VIP>2 & Negative" = "#377EB8"
      ),
      name = "Category"
    ) +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = vip_threshold, linetype = "solid", color = "red") +
    labs(
      title = paste0(column_type, " - ", gsub("exp_", "", exposure_name)),
      x = "Log Fold Change",
      y = "VIP Score (Component 1)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right"
    )

  return(p)
}


# Create VIP scatter plots for all exposures ---------------------------------

vip_scatter_c18 <- exposure_vars |>
  purrr::set_names() |>
  purrr::map(function(exp) {
    create_vip_scatter(
      mwas_results_c18[[exp]],
      vip_c18[[exp]],
      mwas_c18_annotated[[exp]],
      exposure_name = exp,
      column_type = "C18"
    )
  })

vip_scatter_hilic <- exposure_vars |>
  purrr::set_names() |>
  purrr::map(function(exp) {
    create_vip_scatter(
      mwas_results_hilic[[exp]],
      vip_hilic[[exp]],
      mwas_hilic_annotated[[exp]],
      exposure_name = exp,
      column_type = "HILIC"
    )
  })


# Save VIP scatter plots -----------------------------------------------------

purrr::iwalk(vip_scatter_c18, function(p, exp) {
  ggsave(
    filename = here::here("figures", "mwas",
                          glue::glue("vip_scatter_c18_{exp}.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
})

purrr::iwalk(vip_scatter_hilic, function(p, exp) {
  ggsave(
    filename = here::here("figures", "mwas",
                          glue::glue("vip_scatter_hilic_{exp}.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
})


# =============================================================================
# SECTION 3: COMBINED PANEL FIGURES
# =============================================================================

# Create combined figure for a single exposure -------------------------------

create_combined_panel <- function(exposure_name) {

  p1 <- volcano_c18[[exposure_name]] +
    labs(title = "C18 - Volcano") +
    theme(legend.position = "none")

  p2 <- volcano_hilic[[exposure_name]] +
    labs(title = "HILIC - Volcano") +
    theme(legend.position = "none")

  p3 <- vip_scatter_c18[[exposure_name]] +
    labs(title = "C18 - VIP vs logFC") +
    theme(legend.position = "none")

  p4 <- vip_scatter_hilic[[exposure_name]] +
    labs(title = "HILIC - VIP vs logFC") +
    theme(legend.position = "none")

  # Combine with patchwork
  combined <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      title = paste0("MWAS Results: ", gsub("exp_", "", exposure_name)),
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
      )
    )

  return(combined)
}


# Create and save combined panels for all exposures --------------------------

purrr::walk(exposure_vars, function(exp) {
  p <- create_combined_panel(exp)
  ggsave(
    filename = here::here("figures", "mwas",
                          glue::glue("combined_panel_{exp}.png")),
    plot = p,
    width = 14, height = 12, dpi = 300
  )
})


# =============================================================================
# SECTION 4: HEATMAP OF TOP METABOLITES
# =============================================================================

# Function to create heatmap of significant metabolites ----------------------

create_sig_heatmap <- function(combined_results_list, column_type = "C18",
                                top_n = 50) {

  # Get top metabolites across all exposures
  all_sig <- combined_results_list |>
    purrr::imap(function(df, exp) {
      df |>
        dplyr::filter(adj.P.Val < 0.1 | VIP_comp1 > 2) |>
        dplyr::mutate(exposure = exp)
    }) |>
    purrr::list_rbind()

  if (nrow(all_sig) == 0) {
    message("No significant metabolites found for heatmap")
    return(NULL)
  }

  # Get unique top metabolites
  top_mets <- all_sig |>
    dplyr::group_by(met) |>
    dplyr::summarize(
      min_p = min(adj.P.Val),
      max_vip = max(VIP_comp1, na.rm = TRUE)
    ) |>
    dplyr::arrange(min_p) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(met)

  # Create matrix for heatmap
  heatmap_data <- combined_results_list |>
    purrr::imap(function(df, exp) {
      df |>
        dplyr::filter(met %in% top_mets) |>
        dplyr::select(met, logFC) |>
        dplyr::rename(!!exp := logFC)
    }) |>
    purrr::reduce(dplyr::full_join, by = "met") |>
    tibble::column_to_rownames("met") |>
    as.matrix()

  # Clean column names
  colnames(heatmap_data) <- gsub("exp_", "", colnames(heatmap_data))

  # Create heatmap
  pheatmap::pheatmap(
    heatmap_data,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = paste0("Top ", top_n, " Significant Metabolites - ", column_type),
    filename = here::here("figures", "mwas",
                          glue::glue("heatmap_top{top_n}_{column_type}.png")),
    width = 10,
    height = 12
  )
}


# Create heatmaps ------------------------------------------------------------

create_sig_heatmap(combined_results_c18, "C18", top_n = 50)
create_sig_heatmap(combined_results_hilic, "HILIC", top_n = 50)


# =============================================================================
# SECTION 5: CORRELATION BETWEEN EXPOSURES
# =============================================================================

# Create correlation plot of effect sizes across exposures -------------------

create_exposure_correlation <- function(combined_results_list, column_type = "C18") {

  # Get common significant metabolites
  sig_mets <- combined_results_list |>
    purrr::map(function(df) {
      df |>
        dplyr::filter(adj.P.Val < 0.1 | VIP_comp1 > 2) |>
        dplyr::pull(met)
    }) |>
    purrr::reduce(union)

  if (length(sig_mets) < 10) {
    message("Not enough significant metabolites for correlation plot")
    return(NULL)
  }

  # Create logFC matrix
  logfc_matrix <- combined_results_list |>
    purrr::imap(function(df, exp) {
      df |>
        dplyr::filter(met %in% sig_mets) |>
        dplyr::select(met, logFC) |>
        dplyr::rename(!!gsub("exp_", "", exp) := logFC)
    }) |>
    purrr::reduce(dplyr::full_join, by = "met") |>
    tibble::column_to_rownames("met")

  # Calculate correlation matrix
  cor_matrix <- cor(logfc_matrix, use = "pairwise.complete.obs")

  # Create correlation heatmap
  pheatmap::pheatmap(
    cor_matrix,
    color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 10,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste0("Correlation of Effect Sizes - ", column_type),
    filename = here::here("figures", "mwas",
                          glue::glue("exposure_correlation_{column_type}.png")),
    width = 10,
    height = 8
  )
}

create_exposure_correlation(combined_results_c18, "C18")
create_exposure_correlation(combined_results_hilic, "HILIC")


# =============================================================================
# SECTION 6: MANHATTAN-STYLE PLOT
# =============================================================================

# Function to create Manhattan plot ------------------------------------------

create_manhattan <- function(mwas_result, vip_result, exposure_name,
                              column_type = "C18") {

  plot_data <- mwas_result |>
    tibble::rownames_to_column("met") |>
    dplyr::left_join(
      vip_result |>
        tibble::rownames_to_column("met") |>
        dplyr::select(met, VIP = comp1),
      by = "met"
    ) |>
    dplyr::mutate(
      index = row_number(),
      neg_log10_p = -log10(P.Value),
      significant = adj.P.Val < 0.05 | VIP > 2
    )

  # Significance threshold line
  sig_line <- -log10(0.05 / nrow(plot_data))  # Bonferroni

  p <- ggplot(plot_data, aes(x = index, y = neg_log10_p)) +
    geom_point(
      aes(color = significant, size = VIP),
      alpha = 0.6
    ) +
    scale_color_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
      labels = c("TRUE" = "Significant", "FALSE" = "NS"),
      name = "Status"
    ) +
    scale_size_continuous(range = c(1, 4), name = "VIP Score") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = sig_line, linetype = "dashed", color = "red") +
    labs(
      title = paste0("Manhattan Plot: ", column_type, " - ",
                     gsub("exp_", "", exposure_name)),
      x = "Metabolite Index",
      y = expression(-log[10](P-value))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      legend.position = "right"
    )

  return(p)
}


# Create Manhattan plots for key exposures (NOx-related) ---------------------

# Focus on NOx if available
nox_exposures <- exposure_vars[grepl("nox|no2", exposure_vars, ignore.case = TRUE)]

if (length(nox_exposures) > 0) {
  purrr::walk(nox_exposures, function(exp) {
    p <- create_manhattan(
      mwas_results_c18[[exp]],
      vip_c18[[exp]],
      exposure_name = exp,
      column_type = "C18"
    )
    ggsave(
      filename = here::here("figures", "mwas",
                            glue::glue("manhattan_c18_{exp}.png")),
      plot = p,
      width = 12, height = 6, dpi = 300
    )
  })
}


# =============================================================================
# SECTION 7: PATHWAY ENRICHMENT PLOTS
# =============================================================================

# Load additional packages for pathway visualization --------------------------

library(ggthemes)

# Define pathway categories ---------------------------------------------------

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
                             "Chondroitin", "Heparan")

lipid_metabolism <- c("Bile", "Fatty", "lipid", "Phytanic",
                      "Cholesterol", "neuroprostanes", "steroid",
                      "Sphingolipid", "Glycerophospholipid", "Phospholipid",
                      "Triglyceride", "Ceramide")

energy_metabolism <- c("Butanoate", "Carnitine", "Glycolysis",
                       "Pyruvate", "Pentose", "octadecatrienoate",
                       "TCA", "Citrate", "Oxidative")

inflammation_metabolism <- c("Arachidonic", "Leukotriene", "Prostaglandin",
                             "linoleic", "Linoleate", "Eicosanoid")

vitamin_cofactor_metabolism <- c("Vitamin", "Biopterin", "Lipoate",
                                 "Porphyrin", "Catabolism", "Purine",
                                 "Caffeine", "Folate", "Riboflavin",
                                 "Thiamine", "Biotin", "Pantothenate")

signaling_metabolism <- c("Dynorphin", "Dopamine", "Serotonin",
                          "Catecholamine", "Neurotransmitter")

secondary_metabolite_metabolism <- c("Alkaloid")

xenobiotic_metabolism <- c("Xenobiotic", "Drug", "Benzoate")

# Function to categorize pathways ---------------------------------------------

categorize_pathway <- function(pathway_name) {
  case_when(
    str_detect(pathway_name,
               regex(str_c(amino_acid_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "amino acid",
    str_detect(pathway_name,
               regex(str_c(carbohydrate_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "carbohydrate",
    str_detect(pathway_name,
               regex(str_c(lipid_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "lipid",
    str_detect(pathway_name,
               regex(str_c(energy_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "energy",
    str_detect(pathway_name,
               regex(str_c(inflammation_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "inflammation",
    str_detect(pathway_name,
               regex(str_c(vitamin_cofactor_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "vitamin/cofactor",
    str_detect(pathway_name,
               regex(str_c(signaling_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "signaling",
    str_detect(pathway_name,
               regex(str_c(secondary_metabolite_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "secondary",
    str_detect(pathway_name,
               regex(str_c(xenobiotic_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "xenobiotic",
    TRUE ~ "other"
  )
}


# Read pathway enrichment results ---------------------------------------------

# Get all mummichog pathway enrichment files
pathway_files <- list.files(
  here::here("metaboAnalyst"),
  pattern = "mummichog_pathway_enrichment",
  recursive = TRUE,
  full.names = TRUE
)

# Function to read and process pathway file
read_pathway_file <- function(file_path) {
  # Extract exposure name from file path
  dir_name <- basename(dirname(file_path))

  # Read the file (handle both csv and xlsx)
  if (grepl("\\.csv$", file_path)) {
    df <- readr::read_csv(file_path, show_col_types = FALSE)
  } else if (grepl("\\.xlsx$", file_path)) {
    df <- readxl::read_xlsx(file_path)
  } else {
    return(NULL)
  }

  # Standardize column names
  df <- df |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(
      pathway_name = any_of(c("...1", "pathway", "pathway_name", "name")),
      pathway_size = any_of(c("pathway total", "pathway_total", "total")),
      hits_total = any_of(c("hits.total", "hits_total", "total_hits")),
      hits_sig = any_of(c("hits.sig", "hits_sig", "sig_hits")),
      p_value = any_of(c("p(fisher)", "p.value", "pvalue", "p_value"))
    )

  # Add exposure info
  df |>
    dplyr::mutate(
      exposure = dir_name,
      category = categorize_pathway(pathway_name)
    )
}

# Read all pathway files if they exist
if (length(pathway_files) > 0) {

  pathway_all <- pathway_files |>
    purrr::map(read_pathway_file) |>
    purrr::compact() |>
    purrr::list_rbind()

  if (nrow(pathway_all) > 0) {

    # Filter significant pathways (drop zinc and o3)
    pathway_sig <- pathway_all |>
      dplyr::filter(p_value < 0.05) |>
      dplyr::filter(!str_detect(exposure, regex("zinc|o3", ignore_case = TRUE))) |>
      dplyr::arrange(p_value)

    # Save significant pathways to Excel
    if (nrow(pathway_sig) > 0) {
      writexl::write_xlsx(
        pathway_sig,
        here::here("tables", "mwas_results", "pathway_sig_all.xlsx")
      )
    }

    # Create pathway summary plot if we have significant pathways
    if (nrow(pathway_sig) > 0) {

      # Step 1: Get unique pathway/category combinations with min p-value
      pathway_summary <- pathway_sig |>
        dplyr::group_by(pathway_name, category) |>
        dplyr::summarise(
          min_p = min(p_value, na.rm = TRUE),
          .groups = "drop"
        )

      # Step 2: Order categories by their minimum p-value (best category first)
      # Reverse so smallest p-value category appears at TOP of figure
      category_order <- pathway_summary |>
        dplyr::group_by(category) |>
        dplyr::summarise(cat_min_p = min(min_p), .groups = "drop") |>
        dplyr::arrange(desc(cat_min_p)) |>
        dplyr::pull(category)

      # Step 3: Sort pathways - first by category order, then by p-value within category
      # Within each category, sort descending so smallest p-value appears at TOP
      pathway_sorted <- pathway_summary |>
        dplyr::mutate(category = factor(category, levels = category_order)) |>
        dplyr::arrange(category, desc(min_p)) |>
        dplyr::mutate(pathway_name = fct_inorder(pathway_name))

      # Step 4: Prepare data for plotting with proper factor levels
      pathway_plot_data <- pathway_sig |>
        dplyr::mutate(
          category = factor(category, levels = category_order),
          pathway_name = factor(pathway_name, levels = levels(pathway_sorted$pathway_name)),
          exposure_clean = gsub("exp_", "", exposure)
        )

      # Get unique categories for color palette (in order)
      cats <- category_order
      pal <- ggthemes::tableau_color_pal("Tableau 10")(length(cats))
      names(pal) <- cats


      # Create heatmap of pathway p-values across exposures
      pathway_heatmap <- pathway_plot_data |>
        ggplot() +
        geom_tile(
          aes(x = exposure_clean, y = pathway_name, fill = p_value),
          lwd = 1.2,
          linetype = 1,
          color = "white"
        ) +
        scale_fill_gradient(
          low = "#E41A1C",
          high = "#FFFFB2",
          limits = c(0, 0.05),
          breaks = c(0.01, 0.03, 0.05),
          labels = c("0.01", "0.03", "0.05"),
          name = "p-value"
        ) +
        coord_fixed(0.8) +
        labs(
          y = "Pathway",
          x = "Exposure"
        ) +
        theme_classic() +
        theme(
          legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10),
          axis.line = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
        )


      # Create scatter plot of pathway enrichment by exposure
      pathway_scatter <- pathway_plot_data |>
        ggplot(aes(x = -log10(p_value), y = pathway_name)) +
        geom_point(
          aes(size = hits_sig, color = exposure_clean),
          alpha = 0.7
        ) +
        scale_size_continuous(range = c(2, 8), name = "Hits (sig)") +
        scale_color_tableau("Tableau 10") +
        labs(
          y = "",
          x = expression(-log[10](p-value)),
          color = "Exposure"
        ) +
        theme_classic() +
        theme(
          legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10),
          axis.line = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(face = "bold", size = 14),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12)
        )


      # Create color block for pathway categories
      # Count pathways per category in the SAME order as the heatmap
      category_summary <- pathway_sorted |>
        dplyr::group_by(category) |>
        dplyr::summarise(
          category_count = n(),
          .groups = "drop"
        ) |>
        # Keep the same category order
        dplyr::mutate(category = factor(category, levels = category_order)) |>
        dplyr::arrange(category) |>
        # Calculate positions from bottom to top (matching ggplot y-axis)
        dplyr::mutate(
          ymax = cumsum(category_count),
          ymin = lag(ymax, default = 0),
          pos = (ymin + ymax) / 2
        )

      pathway_color_block <- category_summary |>
        ggplot() +
        geom_rect(
          aes(xmin = 0.5, xmax = 1.5,
              ymin = ymin + 0.1, ymax = ymax - 0.1,
              fill = category),
          alpha = 0.5
        ) +
        geom_text(
          aes(x = 1, y = pos, label = category),
          size = 3,
          fontface = "italic"
        ) +
        scale_fill_manual(values = pal, limits = cats, drop = FALSE) +
        theme_void() +
        guides(fill = "none") +
        scale_y_continuous(
          limits = c(0, nrow(pathway_sorted)),
          expand = c(0, 0)
        ) +
        scale_x_continuous(expand = c(0, 0))


      # Combine plots using patchwork
      pathway_combined <- pathway_color_block + pathway_heatmap + pathway_scatter +
        patchwork::plot_layout(
          widths = c(0.2, 0.4, 0.4),
          guides = "collect"
        ) &
        theme(legend.position = "bottom")

      # Save combined pathway plot
      ggsave(
        filename = here::here("figures", "mwas", "pathway_enrichment_summary.png"),
        plot = pathway_combined,
        width = 16,
        height = max(8, nrow(pathway_sorted) * 0.3),
        dpi = 300
      )


      # Create individual pathway heatmap for each category
      pathway_by_category <- pathway_sig |>
        dplyr::group_by(category) |>
        dplyr::group_split()

      purrr::walk(pathway_by_category, function(cat_data) {
        if (nrow(cat_data) < 2) return(NULL)

        cat_name <- unique(cat_data$category)

        p <- cat_data |>
          dplyr::mutate(
            pathway_name = fct_reorder(pathway_name, p_value),
            exposure_clean = gsub("exp_", "", exposure)
          ) |>
          ggplot() +
          geom_tile(
            aes(x = exposure_clean, y = pathway_name, fill = p_value),
            color = "white",
            lwd = 1
          ) +
          scale_fill_gradient(
            low = "#E41A1C",
            high = "#FFFFB2",
            limits = c(0, 0.05),
            name = "p-value"
          ) +
          labs(
            title = paste0("Pathway Enrichment: ", cat_name, " metabolism"),
            x = "Exposure",
            y = "Pathway"
          ) +
          theme_classic() +
          theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 9),
            axis.title = element_text(face = "bold", size = 12)
          )

        ggsave(
          filename = here::here("figures", "mwas",
                                glue::glue("pathway_{gsub('/', '_', cat_name)}.png")),
          plot = p,
          width = 10,
          height = max(6, nrow(cat_data) * 0.25),
          dpi = 300
        )
      })

      message("Pathway enrichment plots created!")

    } else {
      message("No significant pathways found (p < 0.05)")
    }
  } else {
    message("No pathway data could be read from files")
  }
} else {
  message("No pathway enrichment files found in metaboAnalyst directory")
}


message("All visualizations completed! Figures saved to figures/mwas/")

#--------------------------------End of the code--------------------------------
