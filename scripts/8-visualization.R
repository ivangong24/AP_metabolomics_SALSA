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

# Create output directory
dir.create(here::here("figures", "mwas"), 
           showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# SECTION 1: VOLCANO PLOTS
# =============================================================================

# Function to create volcano plot --------------------------------------------

create_volcano_plot <- function(mwas_result, vip_result, exposure_name,
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
    dplyr::mutate(
      neg_log10_p = -log10(P.Value),
      significant = case_when(
        adj.P.Val < fdr_threshold & VIP > vip_threshold & logFC > 0 ~ "Up & VIP>2",
        adj.P.Val < fdr_threshold & VIP > vip_threshold & logFC < 0 ~ "Down & VIP>2",
        adj.P.Val < fdr_threshold & logFC > 0 ~ "Up",
        adj.P.Val < fdr_threshold & logFC < 0 ~ "Down",
        VIP > vip_threshold ~ "VIP>2 only",
        TRUE ~ "NS"
      )
    )

  # Get top metabolites for labeling (by VIP or p-value)
  top_mets <- plot_data |>
    dplyr::filter(adj.P.Val < fdr_threshold | VIP > vip_threshold) |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = n_labels) |>
    dplyr::pull(met)

  plot_data <- plot_data |>
    dplyr::mutate(label = ifelse(met %in% top_mets, met, ""))

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
        "Up & VIP>2" = "#E41A1C",
        "Down & VIP>2" = "#377EB8",
        "Up" = "#FB9A99",
        "Down" = "#A6CEE3",
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

create_vip_scatter <- function(mwas_result, vip_result, exposure_name,
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
    dplyr::mutate(
      category = case_when(
        VIP > vip_threshold & logFC > 0 ~ "VIP>2 & Positive",
        VIP > vip_threshold & logFC < 0 ~ "VIP>2 & Negative",
        TRUE ~ "VIP<=2"
      )
    )

  # Top metabolites for labeling
  top_mets <- plot_data |>
    dplyr::filter(VIP > vip_threshold) |>
    dplyr::arrange(desc(VIP)) |>
    dplyr::slice_head(n = n_labels) |>
    dplyr::pull(met)

  plot_data <- plot_data |>
    dplyr::mutate(label = ifelse(met %in% top_mets, met, ""))

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
    plot_annotation(
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
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
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


message("All visualizations completed! Figures saved to figures/mwas/")

#--------------------------------End of the code--------------------------------
