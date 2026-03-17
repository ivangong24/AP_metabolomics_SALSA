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
source(here::here("scripts", "1-functions.R"))

# Load results ---------------------------------------------------------------

load(here::here("data", "metabolomics", "results", 
                "mwas_results_all.RData"))

load(here::here("data", "metabolomics", "results",
                "mwas_annotation.RData"))

load(here::here("data", "metabolomics", "results",
                "metapone_results_all.RData"))


covar_list <- list(
  covar = quote_all(age_at_blooddraw, gender, edu_year, mh62, 
                    wave, batch, demcind),
  
  covar_sen = quote_all(age_at_blooddraw, gender, edu_year, mh62,
                        alcohol_drinking, pa3_met_if_ca, nses, 
                        bmi_at_blooddraw, diab_at_blooddraw, 
                        wave, batch, demcind)
  
)

exposure_vars <- names(combined_results_list_c18[["all"]][["covar"]])

# Create output directory
names(combined_results_list_c18) |> 
  purrr::map(function(population){
    names(covar_list) |> 
      purrr::map(function(covar_names){
        exposure_vars |> 
          purrr::map(function(exposure){
            dir.create(here::here("figures", "mwas", population, 
                                  covar_names, exposure), 
                       showWarnings = FALSE, recursive = TRUE)
          })
      })
  })


# =============================================================================
# SECTION 1: VOLCANO PLOTS (Combined C18 + HILIC)
# =============================================================================

# Function to create combined volcano plot ------------------------------------

create_volcano_plot <- function(mwas_result_c18, mwas_result_hilic,
                                annotation_result_c18, annotation_result_hilic,
                                exposure_name,
                                fdr_threshold = 0.05,
                                n_labels = 10) {

  # Helper to prepare data for one column type
  prep_data <- function(mwas_result, annotation_result, column_type) {
    mwas_result |>
      tibble::rownames_to_column("met") |>
      dplyr::left_join(
        annotation_result |>
          dplyr::select(met, chemical_id, compound, multiple_match, reference),
        by = "met"
      ) |>
      dplyr::group_by(met) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::mutate(column_type = column_type)
  }

  # Combine C18 and HILIC data
  plot_data <- dplyr::bind_rows(
    prep_data(mwas_result_c18, annotation_result_c18, "C18/neg-"),
    prep_data(mwas_result_hilic, annotation_result_hilic, "HILIC/pos+")
  ) |>
    dplyr::mutate(
      neg_log10_p = -log10(P.Value),
      significant = case_when(
        adj.P.Val < fdr_threshold ~ "FDR < 0.05",
        P.Value < 0.05 ~ "P < 0.05",
        TRUE ~ "NS"
      ),
      significant = factor(significant,
                           levels = c("FDR < 0.05", "P < 0.05", "NS"))
    )

  # Get top metabolites for labeling (by p-value)
  top_mets <- plot_data |>
    dplyr::filter(
      !is.na(compound), compound != "",
      # multiple_match == FALSE, 
      P.Value < 0.05) |>
    dplyr::arrange(P.Value) |>
    dplyr::slice_head(n = n_labels) |>
    # dplyr::slice_max(order_by = -P.Value, n = n_labels, with_ties = FALSE) |>
    dplyr::pull(met)

  plot_data <- plot_data |>
    dplyr::mutate(label = ifelse(met %in% top_mets, compound, ""))

  # Create volcano plot
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_p)) +
    # Non-significant points
    geom_point(
      data = . %>% filter(significant == "NS"),
      aes(shape = column_type),
      color = "grey70",
      alpha = 0.5,
      size = 1.5
    ) +
    # Significant points
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant, shape = column_type),
      alpha = 0.7,
      size = 2.5
    ) +
    scale_color_manual(
      values = c(
        "FDR < 0.05" = "#BE3F42",
        "P < 0.05" = "#DE9960"
      ),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
      name = "Column"
    ) +
    # Add labels
    geom_label_repel(
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
      title = gsub("exp_", "", exposure_name),
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


# Create volcano plots for all populations x covariate sets ------------------

population_names <- names(mwas_results_list_c18)
covar_names <- names(covar_list)

volcano_plots <- population_names |>
  purrr::set_names() |>
  purrr::map(function(pop) {
    covar_names |>
      purrr::set_names() |>
      purrr::map(function(cov) {
        exposure_vars |>
          purrr::set_names() |>
          purrr::map(function(exp) {
            create_volcano_plot(
              mwas_result_c18 = mwas_results_list_c18[[pop]][[cov]][[exp]],
              mwas_result_hilic = mwas_results_list_hilic[[pop]][[cov]][[exp]],
              annotation_result_c18 = mwas_annotated_list_c18[[pop]][[cov]][[exp]],
              annotation_result_hilic = mwas_annotated_list_hilic[[pop]][[cov]][[exp]],
              exposure_name = exp
            )
          })
      })
  })


# Save volcano plots ---------------------------------------------------------

population_names |>
  purrr::walk(function(pop) {
    covar_names |>
      purrr::walk(function(cov) {
        purrr::iwalk(volcano_plots[[pop]][[cov]], function(p, exp) {
          ggsave(
            filename = here::here("figures", "mwas", pop, cov, exp,
                                  glue::glue("volcano_{exp}.png")),
            plot = p,
            width = 10, height = 6, dpi = 300
          )
        })
      })
  })


# =============================================================================
# SECTION 2: LOGFC COMPARISON SCATTER PLOTS (Combined C18 + HILIC)
# =============================================================================

# 2a. Function to compare logFC: demcind vs all & no demcind vs all ----------
# Faceted plot, only showing points where P < 0.05 in "all"
# Combines C18 and HILIC with different shapes

create_population_scatter <- function(mwas_list_c18, mwas_list_hilic,
                                      covar_set, exposure_name) {

  # Helper to prep data for one column type
  prep_pop_data <- function(mwas_list, column_type) {
    df_all <- mwas_list[["all"]][[covar_set]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_all = logFC, P.Value_all = P.Value)

    df_demcind <- mwas_list[["demcind"]][[covar_set]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_demcind = logFC, P.Value_demcind = P.Value)

    df_no_demcind <- mwas_list[["no demcind"]][[covar_set]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_no_demcind = logFC, P.Value_no_demcind = P.Value)

    # Join all three populations
    df_all |>
      dplyr::inner_join(df_demcind, by = "met") |>
      dplyr::inner_join(df_no_demcind, by = "met") |>
      # Classify significance based on demcind and no demcind
      dplyr::mutate(
        significant = case_when(
          P.Value_demcind < 0.05 & P.Value_no_demcind < 0.05 ~ "Both P < 0.05",
          P.Value_demcind < 0.05 ~ "demcind P < 0.05",
          P.Value_no_demcind < 0.05 ~ "no demcind P < 0.05",
          TRUE ~ "NS"
        )
      ) |>
      # Pivot to long for faceting
      tidyr::pivot_longer(
        cols = c(logFC_demcind, logFC_no_demcind),
        names_to = "comparison",
        values_to = "logFC_sub",
        names_prefix = "logFC_"
      ) |>
      dplyr::mutate(
        comparison = dplyr::recode(comparison,
                                   "demcind" = "demcind vs all",
                                   "no_demcind" = "no demcind vs all"),
        column_type = column_type
      )
  }

  plot_data <- dplyr::bind_rows(
    prep_pop_data(mwas_list_c18, "C18/neg-"),
    prep_pop_data(mwas_list_hilic, "HILIC/pos+")
  ) |>
    # Only keep points where P < 0.05 in "all"
    dplyr::filter(P.Value_all < 0.05) |>
    dplyr::mutate(
      significant = factor(significant,
                           levels = c("Both P < 0.05",
                                      "demcind P < 0.05",
                                      "no demcind P < 0.05",
                                      "NS"))
    )

  # Per-facet correlation on P < 0.05 points (all points shown are P<0.05 in all)
  # cor_labels <- plot_data |>
  #   dplyr::group_by(comparison) |>
  #   dplyr::summarise(
  #     cor_test = list(cor.test(logFC_all, logFC_sub, use = "complete.obs")),
  #     .groups = "drop"
  #   ) |>
  #   dplyr::mutate(
  #     r = purrr::map_dbl(cor_test, ~ .x$estimate),
  #     p_val = purrr::map_dbl(cor_test, ~ .x$p.value),
  #     label = paste0("r = ", round(r, 3), "\np = ",
  #                    ifelse(p_val < 0.001,
  #                           formatC(p_val, format = "e", digits = 2),
  #                           round(p_val, 3)))
  #   )

  p <- ggplot(plot_data, aes(x = logFC_all, y = logFC_sub)) +
    # NS points
    geom_point(
      data = . %>% filter(significant == "NS"),
      aes(shape = column_type),
      color = "grey70", alpha = 0.5, size = 1.5
    ) +
    # Significant points
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant, shape = column_type),
      alpha = 0.7, size = 2
    ) +
    # Fit line on all displayed points (all P < 0.05 in "all")
    geom_smooth(
      method = "lm", se = TRUE,
      color = "black", linewidth = 0.8, linetype = "solid",
      alpha = 0.2
    ) +
    scale_color_manual(
      values = c(
        "Both P < 0.05" = "#B73F42",
        "demcind P < 0.05" = "#DE9960",
        "no demcind P < 0.05" = "#436C85"
      ),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
      name = "Column"
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    # Per-facet correlation labels
    ggpubr::stat_cor(
      method = "pearson",
      label.x.npc = "left", label.y.npc = "top",
      size = 4, fontface = "italic"
    ) +
    # geom_text(
    #   data = cor_labels,
    #   aes(x = Inf, y = -Inf, label = label),
    #   hjust = 1.1, vjust = -0.3, size = 4, fontface = "italic",
    #   inherit.aes = FALSE
    # ) +
    facet_wrap(~ comparison) +
    labs(
      title = gsub("exp_", "", exposure_name),
      x = "MWAS beta coefficients among all participants",
      y = "MWAS beta coefficients among the subgroup population"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    )

  return(p)
}


# 2b. Function to compare logFC between two covariate sets -------------------
# Combines C18 and HILIC with different shapes

create_covariate_scatter <- function(mwas_list_c18, mwas_list_hilic,
                                     population, cov1, cov2,
                                     exposure_name) {

  # Helper to prep and join data for one column type
  prep_cov_data <- function(mwas_list, column_type) {
    df1 <- mwas_list[[population]][[cov1]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_cov1 = logFC, P.Value_cov1 = P.Value)

    df2 <- mwas_list[[population]][[cov2]][[exposure_name]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_cov2 = logFC, P.Value_cov2 = P.Value)

    dplyr::inner_join(df1, df2, by = "met") |>
      dplyr::mutate(column_type = column_type)
  }

  plot_data <- dplyr::bind_rows(
    prep_cov_data(mwas_list_c18, "C18/neg-"),
    prep_cov_data(mwas_list_hilic, "HILIC/pos+")
  ) |>
    dplyr::mutate(
      significant = case_when(
        P.Value_cov1 < 0.05 & P.Value_cov2 < 0.05 ~ "Both P < 0.05",
        P.Value_cov1 < 0.05 ~ paste0("Primary analysis -", " P < 0.05"),
        P.Value_cov2 < 0.05 ~ paste0("Sensitivity analysis -", " P < 0.05"),
        TRUE ~ "NS"
      ),
      significant = factor(significant,
                           levels = c("Both P < 0.05",
                                      paste0("Primary analysis -", " P < 0.05"),
                                      paste0("Sensitivity analysis -", " P < 0.05"),
                                      "NS"))
    )

  # Correlation on P < 0.05 points only
  sig_data <- plot_data |>
    dplyr::filter(significant != "NS")

  # cor_test <- cor.test(sig_data$logFC_cov1, sig_data$logFC_cov2,
  #                      use = "complete.obs")
  # r <- cor_test$estimate
  # p_val <- cor_test$p.value
  # cor_label <- paste0("r = ", round(r, 3), "\np = ",
  #                     ifelse(p_val < 0.001, formatC(p_val, format = "e", digits = 2),
  #                            round(p_val, 3)))

  p <- ggplot(plot_data, aes(x = logFC_cov1, y = logFC_cov2)) +
    # NS points
    geom_point(
      data = . %>% filter(significant == "NS"),
      aes(shape = column_type),
      color = "grey70", alpha = 0.4, size = 1.5
    ) +
    # Significant points
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant, shape = column_type),
      alpha = 0.7, size = 2
    ) +
    # Fit line on P < 0.05 points only
    geom_smooth(
      data = sig_data,
      aes(x = logFC_cov1, y = logFC_cov2),
      method = "lm", se = TRUE,
      color = "black", linewidth = 0.8, linetype = "solid",
      alpha = 0.2, inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Both P < 0.05" = "#B73F42",
        setNames("#DE9960", paste0("Primary analysis -", " P < 0.05")),
        setNames("#436C85", paste0("Sensitivity analysis -", " P < 0.05"))
      ),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
      name = "Column"
    ) +
    # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    ggpubr::stat_cor(
      method = "pearson",
      label.x.npc = "left", label.y.npc = "top",
      size = 4, fontface = "italic"
    ) +
    # annotate("text", x = Inf, y = -Inf,
    #          label = cor_label,
    #          hjust = 1.1, vjust = -0.3, size = 4, fontface = "italic") +
    labs(
      title = paste0(gsub("exp_", "", exposure_name), " (", population, ")"),
      x = "MWAS beta coefficients in the primary analysis",
      y = "MWAS beta coefficients in the sensitivity analysis (additional covariates)"
    ) +
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


# Create and save population comparison scatter plots -------------------------
# Compare demcind vs all and no demcind vs all (faceted)

covar_names |>
  purrr::walk(function(cov) {
    exposure_vars |>
      purrr::walk(function(exp) {
        p <- create_population_scatter(
          mwas_results_list_c18, mwas_results_list_hilic,
          covar_set = cov, exposure_name = exp
        )
        ggsave(
          filename = here::here("figures", "mwas", "all", cov, exp,
                                glue::glue("scatter_pop_{exp}.png")),
          plot = p,
          width = 14, height = 7, dpi = 300
        )
      })
  })


# Create and save covariate comparison scatter plots --------------------------
# Compare covar vs covar_sen within each population and exposure

population_names |>
  purrr::walk(function(pop) {
    exposure_vars |>
      purrr::walk(function(exp) {
        p <- create_covariate_scatter(
          mwas_results_list_c18, mwas_results_list_hilic,
          population = pop,
          cov1 = "covar", cov2 = "covar_sen",
          exposure_name = exp
        )
        ggsave(
          filename = here::here("figures", "mwas", pop, "covar", exp,
                                glue::glue("scatter_cov_{exp}.png")),
          plot = p,
          width = 8, height = 7, dpi = 300
        )
      })
  })


# =============================================================================
# SECTION 2c: COMPOSITE EXPOSURE COMPARISON SCATTER PLOTS
# =============================================================================

# Function to compare logFC between all-toxicant and traffic-related composites
# Combines C18 and HILIC with different shapes

create_composite_scatter <- function(mwas_list_c18, mwas_list_hilic,
                                     exp_all, exp_traffic,
                                     population, covar_set,
                                     method_label = "WQS") {

  # Helper to prep data for one column type
  prep_comp_data <- function(mwas_list, column_type) {
    df_all <- mwas_list[[population]][[covar_set]][[exp_all]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_all = logFC, P.Value_all = P.Value)

    df_traffic <- mwas_list[[population]][[covar_set]][[exp_traffic]] |>
      tibble::rownames_to_column("met") |>
      dplyr::select(met, logFC_traffic = logFC, P.Value_traffic = P.Value)

    dplyr::inner_join(df_all, df_traffic, by = "met") |>
      dplyr::mutate(column_type = column_type)
  }

  plot_data <- dplyr::bind_rows(
    prep_comp_data(mwas_list_c18, "C18/neg-"),
    prep_comp_data(mwas_list_hilic, "HILIC/pos+")
  ) |>
    dplyr::mutate(
      significant = case_when(
        P.Value_all < 0.05 & P.Value_traffic < 0.05 ~ "Both P < 0.05",
        P.Value_all < 0.05 ~ "All toxicants P < 0.05",
        P.Value_traffic < 0.05 ~ "Traffic-related P < 0.05",
        TRUE ~ "NS"
      ),
      significant = factor(significant,
                           levels = c("Both P < 0.05",
                                      "All toxicants P < 0.05",
                                      "Traffic-related P < 0.05",
                                      "NS"))
    )

  # Fit line data: points where at least one is P < 0.05
  sig_data <- plot_data |>
    dplyr::filter(significant != "NS")

  p <- ggplot(plot_data, aes(x = logFC_all, y = logFC_traffic)) +
    # NS points
    geom_point(
      data = . %>% filter(significant == "NS"),
      aes(shape = column_type),
      color = "grey70", alpha = 0.5, size = 1.5
    ) +
    # Significant points
    geom_point(
      data = . %>% filter(significant != "NS"),
      aes(color = significant, shape = column_type),
      alpha = 0.7, size = 2
    ) +
    # Fit line on P < 0.05 points only
    geom_smooth(
      data = sig_data,
      aes(x = logFC_all, y = logFC_traffic),
      method = "lm", se = TRUE,
      color = "black", linewidth = 0.8, linetype = "solid",
      alpha = 0.2, inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Both P < 0.05" = "#B73F42",
        "All toxicants P < 0.05" = "#DE9960",
        "Traffic-related P < 0.05" = "#436C85"
      ),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c("C18/neg-" = 16, "HILIC/pos+" = 17),
      name = "Column"
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    ggpubr::stat_cor(
      method = "pearson",
      label.x.npc = "left", label.y.npc = "top",
      size = 4, fontface = "italic"
    ) +
    labs(
      title = paste0(method_label, " (", population, ")"),
      x = paste0("logFC (", gsub("exp_", "", exp_all), ")"),
      y = paste0("logFC (", gsub("exp_", "", exp_traffic), ")")
    ) +
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


# Define composite exposure pairs
composite_pairs <- list(
  WQS = list(all = exposure_vars[1], traffic = exposure_vars[2]),
  qgcomp = list(all = exposure_vars[4], traffic = exposure_vars[5])
)

# Create and save composite comparison scatter plots --------------------------

population_names |>
  purrr::walk(function(pop) {
    covar_names |>
      purrr::walk(function(cov) {
        purrr::iwalk(composite_pairs, function(pair, method) {
          p <- create_composite_scatter(
            mwas_results_list_c18, mwas_results_list_hilic,
            exp_all = pair$all, exp_traffic = pair$traffic,
            population = pop, covar_set = cov,
            method_label = method
          )
          ggsave(
            filename = here::here("figures", "mwas", pop, cov,
                                  glue::glue("scatter_composite_{method}.png")),
            plot = p,
            width = 8, height = 7, dpi = 300
          )
        })
      })
  })


# =============================================================================
# SECTION 3: COMBINED PANEL FIGURES
# =============================================================================

# Create combined figure for a single exposure -------------------------------
# Row 1: volcano
# Row 2: population comparison | covariate comparison
# Row 3: WQS composite comparison | qgcomp composite comparison

create_combined_panel <- function(exposure_name, population, covar_set) {

  # Common legend theme: boxed legends
  legend_theme <- theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.background = element_rect(colour = "grey80", fill = "white",
                                     linewidth = 0.5),
    legend.margin = margin(4, 6, 4, 6),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

  # p1: volcano — show Column legend (ordered first) + Significance
  p1 <- volcano_plots[[population]][[covar_set]][[exposure_name]] +
    labs(title = NULL) +
    guides(
      shape = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    legend_theme

  # p2: population scatter — color legend only, no Column
  p2 <- create_population_scatter(
    mwas_results_list_c18, mwas_results_list_hilic,
    covar_set = covar_set, exposure_name = exposure_name
  ) +
    labs(title = NULL) +
    guides(shape = "none") +
    legend_theme

  # p3: covariate scatter — color legend only, no Column
  p3 <- create_covariate_scatter(
    mwas_results_list_c18, mwas_results_list_hilic,
    population = population,
    cov1 = "covar", cov2 = "covar_sen",
    exposure_name = exposure_name
  ) +
    labs(title = NULL) +
    guides(shape = "none") +
    legend_theme

  # Combine with patchwork: single row, each plot keeps its own legend
  combined <- (p1 | p2 | p3) +
    patchwork::plot_annotation(
      title = paste0("MWAS Results: ", gsub("exp_", "", exposure_name),
                     " (", population, ", ", covar_set, ")"),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.tag = element_text(face = "bold", size = 14)
      )
    )

  return(combined)
}


# Create and save combined panels for all populations x covariates x exposures

population_names |>
  purrr::walk(function(pop) {
    covar_names |>
      purrr::walk(function(cov) {
        exposure_vars |>
          purrr::walk(function(exp) {
            p <- create_combined_panel(exp, pop, cov)
            ggsave(
              filename = here::here("figures", "mwas", pop, cov, exp,
                                    glue::glue("combined_panel_{exp}.png")),
              plot = p,
              width = 24, height = 8, dpi = 300
            )
          })
      })
  })


# =============================================================================
# SECTION 4: HEATMAP OF TOP METABOLITES
# =============================================================================

# Function to create heatmap of significant metabolites ----------------------

create_sig_heatmap <- function(mwas_results_exposure_list, column_type = "C18",
                                population = "all", covar_set = "covar",
                                top_n = 50) {

  # Get top metabolites across all exposures
  all_sig <- mwas_results_exposure_list |>
    purrr::imap(function(df, exp) {
      df |>
        tibble::rownames_to_column("met") |>
        dplyr::filter(adj.P.Val < 0.1) |>
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
      .groups = "drop"
    ) |>
    dplyr::arrange(min_p) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(met)

  # Create matrix for heatmap
  heatmap_data <- mwas_results_exposure_list |>
    purrr::imap(function(df, exp) {
      df |>
        tibble::rownames_to_column("met") |>
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
    main = paste0("Top ", top_n, " Significant Metabolites - ", column_type,
                  " (", population, ", ", covar_set, ")"),
    filename = here::here("figures", "mwas", population, covar_set,
                          glue::glue("heatmap_top{top_n}_{column_type}.png")),
    width = 10,
    height = 12
  )
}


# Create heatmaps ------------------------------------------------------------

population_names |>
  purrr::walk(function(pop) {
    covar_names |>
      purrr::walk(function(cov) {
        create_sig_heatmap(mwas_results_list_c18[[pop]][[cov]], "C18",
                           population = pop, covar_set = cov, top_n = 50)
        create_sig_heatmap(mwas_results_list_hilic[[pop]][[cov]], "HILIC",
                           population = pop, covar_set = cov, top_n = 50)
      })
  })


# =============================================================================
# SECTION 5: CORRELATION BETWEEN EXPOSURES
# =============================================================================

# Create correlation plot of effect sizes across exposures -------------------

create_exposure_correlation <- function(mwas_results_exposure_list,
                                        column_type = "C18",
                                        population = "all",
                                        covar_set = "covar") {

  # Get common significant metabolites
  sig_mets <- mwas_results_exposure_list |>
    purrr::map(function(df) {
      df |>
        tibble::rownames_to_column("met") |>
        dplyr::filter(adj.P.Val < 0.1) |>
        dplyr::pull(met)
    }) |>
    purrr::reduce(union)

  if (length(sig_mets) < 10) {
    message("Not enough significant metabolites for correlation plot")
    return(NULL)
  }

  # Create logFC matrix
  logfc_matrix <- mwas_results_exposure_list |>
    purrr::imap(function(df, exp) {
      df |>
        tibble::rownames_to_column("met") |>
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
    main = paste0("Correlation of Effect Sizes - ", column_type,
                  " (", population, ", ", covar_set, ")"),
    filename = here::here("figures", "mwas", population, covar_set,
                          glue::glue("exposure_correlation_{column_type}.png")),
    width = 10,
    height = 8
  )
}

population_names |>
  purrr::walk(function(pop) {
    covar_names |>
      purrr::walk(function(cov) {
        create_exposure_correlation(mwas_results_list_c18[[pop]][[cov]], "C18",
                                    population = pop, covar_set = cov)
        create_exposure_correlation(mwas_results_list_hilic[[pop]][[cov]], "HILIC",
                                    population = pop, covar_set = cov)
      })
  })


# =============================================================================
# SECTION 6: MANHATTAN-STYLE PLOT
# =============================================================================

# Function to create Manhattan plot ------------------------------------------

create_manhattan <- function(mwas_result, exposure_name,
                              column_type = "C18") {

  plot_data <- mwas_result |>
    tibble::rownames_to_column("met") |>
    dplyr::mutate(
      index = row_number(),
      neg_log10_p = -log10(P.Value),
      significant = case_when(
        adj.P.Val < 0.05 ~ "FDR < 0.05",
        P.Value < 0.05 ~ "P < 0.05",
        TRUE ~ "NS"
      )
    )

  # Significance threshold line
  sig_line <- -log10(0.05 / nrow(plot_data))  # Bonferroni

  p <- ggplot(plot_data, aes(x = index, y = neg_log10_p)) +
    geom_point(
      aes(color = significant),
      alpha = 0.6,
      size = 1.5
    ) +
    scale_color_manual(
      values = c("FDR < 0.05" = "#E41A1C", "P < 0.05" = "#377EB8",
                 "NS" = "grey60"),
      name = "Significance"
    ) +
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
  population_names |>
    purrr::walk(function(pop) {
      covar_names |>
        purrr::walk(function(cov) {
          purrr::walk(nox_exposures, function(exp) {
            c("C18", "HILIC") |>
              purrr::walk(function(col_type) {
                mwas_list <- if (col_type == "C18") {
                  mwas_results_list_c18
                } else {
                  mwas_results_list_hilic
                }
                p <- create_manhattan(
                  mwas_list[[pop]][[cov]][[exp]],
                  exposure_name = exp,
                  column_type = col_type
                )
                ggsave(
                  filename = here::here("figures", "mwas", pop, cov, exp,
                                        glue::glue("manhattan_{tolower(col_type)}_{exp}.png")),
                  plot = p,
                  width = 12, height = 6, dpi = 300
                )
              })
          })
        })
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
