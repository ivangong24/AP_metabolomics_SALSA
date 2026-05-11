## ---------------------------
## Test script: Alternative pathway enrichment visualizations
## Purpose: Explore cleaner alternatives to the pathway scatter plot
##          in Section 7 of 8-visualization.R
## ---------------------------

# Load required packages and data -------------------------------------------
source(here::here("scripts", "1-functions.R"))

library(ggthemes)
library(ggh4x)
library(legendary)
library(ggnewscale)
library(ggtext)

load(here::here("data", "metabolomics", "results",
                "metapone_results_all.RData"))

# --- Reuse pathway categories and helpers from 8-visualization.R ---

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
                             "Chondroitin", "Heparan", "Nucleotide sugar")

lipid_metabolism <- c("Bile", "Fatty", "lipid", "Phytanic",
                      "Cholesterol", "neuroprostanes", "steroid",
                      "Sphingolipid", "Glycerophospholipid", "Phospholipid",
                      "Triglyceride", "Ceramide")

energy_metabolism <- c("Butanoate", "Carnitine", "Glycolysis",
                       "Pyruvate", "Pentose", "octadecatrienoate",
                       "TCA", "Citrate", "Oxidative", "Glyoxylate",
                       "Propanoate", "Carbon fixation")

inflammation_metabolism <- c("Arachidonic", "Leukotriene", "Prostaglandin",
                             "linoleic", "Linoleate", "Eicosanoid")

vitamin_cofactor_metabolism <- c("Vitamin", "Biopterin", "Lipoate",
                                 "Porphyrin", "Catabolism", 
                                 "Caffeine", "Folate", "Riboflavin",
                                 "Thiamine", "Biotin", "Pantothenate")

nucleotide_metabolism <- c("Pyrimidine", "Purine")

signaling_metabolism <- c("Dynorphin", "Dopamine", "Serotonin",
                          "Catecholamine", "Neurotransmitter")

secondary_metabolite_metabolism <- c("Alkaloid")

xenobiotic_metabolism <- c("Xenobiotic", "Drug", "Benzoate")

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
               regex(str_c(nucleotide_metabolism, collapse = "|"),
                     ignore_case = TRUE)) ~ "nucleotide",
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

# --- Read and prepare pathway data ---

pathway_files <- list.files(
  here::here("metaboAnalyst", "Output"),
  pattern = "mummichog_pathway_enrichment",
  recursive = TRUE,
  full.names = TRUE
) |>
  purrr::keep(~ grepl("comp_qgcomp_all|comp_wqs_all|comp_qgcomp_cox_all", .x))

read_pathway_file <- function(file_path) {
  path_parts <- unlist(strsplit(file_path, "/"))
  output_idx <- which(path_parts == "Output")
  study <- path_parts[output_idx + 1]
  population <- path_parts[output_idx + 2]
  covar_set <- path_parts[output_idx + 3]
  exposure <- path_parts[output_idx + 4]

  if (grepl("\\.csv$", file_path)) {
    df <- readr::read_csv(file_path, show_col_types = FALSE)
  } else if (grepl("\\.xlsx$", file_path)) {
    df <- readxl::read_xlsx(file_path)
  } else {
    return(NULL)
  }

  df <- df |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(
      pathway_name = any_of(c("...1", "pathway", "pathway_name", "name")),
      pathway_size = any_of(c("pathway total", "pathway_total", "total")),
      hits_total = any_of(c("hits.total", "hits_total", "total_hits")),
      hits_sig = any_of(c("hits.sig", "hits_sig", "sig_hits")),
      expected = any_of(c("expected")),
      p_value = any_of(c("p(fisher)", "p.value", "pvalue", "p_value"))
    )

  df |>
    dplyr::mutate(
      exposure = exposure,
      population = population,
      covar_set = covar_set,
      category = categorize_pathway(pathway_name)
    )
}

pathway_all <- pathway_files |>
  purrr::map(read_pathway_file) |>
  purrr::compact() |>
  purrr::list_rbind()

pathway_sig <- pathway_all |>
  dplyr::filter(p_value < 0.1) |>
  dplyr::arrange(p_value) |>
  dplyr::mutate(
    enrichment_factor = hits_sig / expected,
    exposure_clean = dplyr::case_when(
      exposure == "comp_wqs_all" ~ "Air toxicant composite (WQS)",
      exposure == "comp_qgcomp_all" ~ "Air toxicant composite (QG-computation)",
      exposure == "comp_qgcomp_cox_all" ~ "Air toxicant composite (QGcomp-Cox)",
      TRUE ~ exposure
    ),
    exposure_short = dplyr::case_when(
      exposure == "comp_wqs_all" ~ "WQS",
      exposure == "comp_qgcomp_all" ~ "QGcomp",
      exposure == "comp_qgcomp_cox_all" ~ "QGcomp-Cox",
      TRUE ~ exposure
    ),
    population = factor(population,
                        levels = c("all", "no demcind", "demcind"))
  )

# Create output directory for test figures
dir.create(here::here("figures", "test_pathway"), showWarnings = FALSE,
           recursive = TRUE)


# =============================================================================
# OPTION A: BUBBLE HEATMAP
# =============================================================================
# Replace the scatter with a second heatmap where:
#   - tile color = -log10(p-value) (continuous)
#   - point size = enrichment factor (overlaid circles)
# Faceted: category (rows) × population (columns), x-axis = exposure_short
# This mirrors the left heatmap structure and is easy to read side-by-side.

create_pathway_bubble_heatmap <- function(data, covar_label) {

  if (nrow(data) == 0) return(NULL)

  # Pathway and category ordering. min_p alone is a poor tiebreaker for
  # metapone since permutation p-values frequently tie at the recoded
  # floor (0.0005); fall back on total -log10(p) so pathways significant
  # across more (population × exposure) cells rise to the top of a panel.
  pathway_summary <- data |>
    dplyr::group_by(pathway_name, category) |>
    dplyr::summarise(
      min_p = min(p_value, na.rm = TRUE),
      total_score = sum(-log10(p_value), na.rm = TRUE),
      .groups = "drop"
    )

  category_order <- pathway_summary |>
    dplyr::group_by(category) |>
    dplyr::summarise(
      cat_min_p = min(min_p),
      cat_total = sum(total_score),
      .groups = "drop"
    ) |>
    dplyr::arrange(cat_min_p, dplyr::desc(cat_total)) |>
    dplyr::pull(category)

  pathway_sorted <- pathway_summary |>
    dplyr::mutate(category = factor(category, levels = category_order)) |>
    dplyr::arrange(category, dplyr::desc(min_p), total_score) |>
    dplyr::mutate(pathway_name = fct_inorder(pathway_name))

  exposure_levels <- c("WQS", "QGcomp", "QGcomp-Cox")
  exposure_display_labels <- c(
    "WQS" = "WQS",
    "QGcomp" = "QGcomp",
    "QGcomp-Cox" = "QGcomp"
  )

  plot_data <- data |>
    dplyr::mutate(
      category = factor(category, levels = category_order),
      pathway_name = factor(pathway_name,
                            levels = levels(pathway_sorted$pathway_name)),
      neg_log10_p = -log10(p_value),
      exposure_short = factor(exposure_short, levels = exposure_levels)
    )

  # Color palettes
  cats <- category_order
  pal <- ggthemes::tableau_color_pal("Tableau 10")(length(cats))
  names(pal) <- cats
  strip_colors <- lapply(pal[cats], function(col) {
    element_rect(fill = alpha(col, 0.3), colour = "grey80")
  })

  pop_levels <- c("all", "no demcind", "demcind")
  pop_pal <- setNames(c("#436C85", "#B73F42", "#DE9960"), pop_levels)
  pop_strip_colors <- lapply(pop_pal[pop_levels], function(col) {
    element_rect(fill = alpha(col, 0.7), colour = "grey80")
  })

  # Complete grid for white background
  all_combos <- tidyr::expand_grid(
    pathway_name = levels(pathway_sorted$pathway_name),
    population = factor(pop_levels, levels = pop_levels),
    exposure_short = factor(exposure_levels, levels = exposure_levels)
  ) |>
    dplyr::mutate(
      pathway_name = factor(pathway_name,
                            levels = levels(pathway_sorted$pathway_name))
    ) |>
    dplyr::left_join(
      pathway_sorted |> dplyr::select(pathway_name, category),
      by = "pathway_name"
    )

  bubble_data <- all_combos |>
    dplyr::left_join(
      plot_data |>
        dplyr::select(pathway_name, population, exposure_short,
                       p_value, neg_log10_p, enrichment_factor),
      by = c("pathway_name", "population", "exposure_short")
    )

  # Alternating row shading by visual position (top-to-bottom across all
  # category panels) so adjacent rows at panel seams never share a shade.
  shade_df <- pathway_sorted |>
    dplyr::arrange(category, dplyr::desc(as.integer(pathway_name))) |>
    dplyr::mutate(
      visual_idx = dplyr::row_number(),
      shade = ifelse(visual_idx %% 2 == 0, "even", "odd")
    ) |>
    dplyr::select(pathway_name, shade)
  shade_data <- bubble_data |>
    dplyr::left_join(shade_df, by = "pathway_name")

  p <- shade_data |>
    ggplot(aes(x = exposure_short, y = pathway_name)) +
    # Alternating row shading
    geom_tile(
      aes(fill = shade),
      width = 1, height = 1, alpha = 0.4, show.legend = FALSE
    ) +
    scale_fill_manual(values = c("even" = "grey93", "odd" = "white"),
                      guide = "none") +
    # New fill scale for bubbles (using ggnewscale)
    ggnewscale::new_scale_fill() +
    # Bubbles: size = enrichment factor, color = -log10(p)
    geom_point(
      data = . %>% filter(!is.na(p_value)),
      aes(size = enrichment_factor, fill = neg_log10_p),
      shape = 21, color = "grey20", stroke = 0.4, alpha = 0.9
    ) +
    scale_fill_gradientn(
      colours = c("#FFF7BC", "#FEC44F", "#F46D43", "#D73027", "#A50026"),
      name = expression(-log[10](italic(p))),
      limits = c(-log10(0.1), NA),
      na.value = "grey85"
    ) +
    scale_size_continuous(
      range = c(2, 9), name = "Enrichment factor"
    ) +
    scale_x_discrete(
      labels = exposure_display_labels,
      guide = legendry::guide_axis_nested(
        key = legendry::key_range_manual(
          start = c("WQS", "QGcomp-Cox"),
          end = c("QGcomp", "QGcomp-Cox"),
          name = c("Cross-sectional", "Time-to-event")
        )
      )
    ) +
    guides(
      fill = guide_colorbar(
        barwidth = 10, barheight = 0.8,
        title.position = "top", title.hjust = 0.5,
        frame.colour = "grey40", ticks.colour = "grey40"
      ),
      size = guide_legend(
        title.position = "top", title.hjust = 0.5,
        override.aes = list(fill = "#F46D43", alpha = 0.8)
      )
    ) +
    ggh4x::facet_grid2(
      category ~ population,
      scales = "free_y", space = "free_y",
      strip = ggh4x::strip_themed(
        background_y = strip_colors,
        background_x = pop_strip_colors,
        text_x = list(
          element_text(face = "bold", size = 14, color = "white"),
          element_text(face = "bold", size = 14, color = "white"),
          element_text(face = "bold", size = 14, color = "white")
        )
      )
    ) +
    labs(y = NULL, x = NULL, title = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 8),
      legend.spacing.x = unit(1, "cm"),
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "grey80", fill = NA,
                                  linewidth = 0.4),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 13, color = "grey20"),
      axis.text.x = ggtext::element_textbox_simple(
        face = "bold", size = 13, color = "grey20",
        halign = 0.5,
        padding = margin(3, 6, 3, 6),
        margin = margin(t = 3, b = 3),
        box.color = "grey40",
        linewidth = 0.5,
        linetype = 1,
        fill = "grey95",
        r = unit(2, "pt")
      ),
      strip.text.y = element_text(face = "bold.italic", size = 13,
                                  angle = 0, color = "grey20"),
      panel.spacing = unit(0.3, "lines"),
      plot.margin = margin(10, 15, 10, 10)
    )

  return(p)
}

# Test with covar data
covar_data <- pathway_sig |> dplyr::filter(covar_set == "covar")

p_bubble <- create_pathway_bubble_heatmap(covar_data, "Primary analysis")
if (!is.null(p_bubble)) {
  ggsave(
    filename = here::here("figures", "test_pathway",
                          "option_A_bubble_heatmap.png"),
    plot = p_bubble,
    width = 16,
    height = max(8, dplyr::n_distinct(covar_data$pathway_name) * 0.35),
    dpi = 300
  )
}

# Sensitivity analysis
covar_sen_data <- pathway_sig |> dplyr::filter(covar_set == "covar_sen")

p_bubble_sen <- create_pathway_bubble_heatmap(covar_sen_data,
                                              "Sensitivity analysis")
if (!is.null(p_bubble_sen)) {
  ggsave(
    filename = here::here("figures", "test_pathway",
                          "option_A_bubble_heatmap_sen.png"),
    plot = p_bubble_sen,
    width = 16,
    height = max(8, dplyr::n_distinct(covar_sen_data$pathway_name) * 0.35),
    dpi = 300
  )
}


# =============================================================================
# OPTION B: GROUPED BAR CHART OF -log10(p)
# =============================================================================
# Horizontal bars grouped by exposure within each population facet column.
# Bar length = -log10(p), color = exposure. Clean and easy to compare magnitudes.

create_pathway_bar <- function(data, covar_label) {

  if (nrow(data) == 0) return(NULL)

  pathway_summary <- data |>
    dplyr::group_by(pathway_name, category) |>
    dplyr::summarise(min_p = min(p_value, na.rm = TRUE), .groups = "drop")

  category_order <- pathway_summary |>
    dplyr::group_by(category) |>
    dplyr::summarise(cat_min_p = min(min_p), .groups = "drop") |>
    dplyr::arrange(cat_min_p) |>
    dplyr::pull(category)

  pathway_sorted <- pathway_summary |>
    dplyr::mutate(category = factor(category, levels = category_order)) |>
    dplyr::arrange(category, desc(min_p)) |>
    dplyr::mutate(pathway_name = fct_inorder(pathway_name))

  plot_data <- data |>
    dplyr::mutate(
      category = factor(category, levels = category_order),
      pathway_name = factor(pathway_name,
                            levels = levels(pathway_sorted$pathway_name)),
      neg_log10_p = -log10(p_value)
    )

  cats <- category_order
  pal <- ggthemes::tableau_color_pal("Tableau 10")(length(cats))
  names(pal) <- cats
  strip_colors <- lapply(pal[cats], function(col) {
    element_rect(fill = alpha(col, 0.3), colour = "grey80")
  })

  pop_levels <- c("all", "no demcind", "demcind")
  pop_pal <- setNames(c("#436C85", "#B73F42", "#DE9960"), pop_levels)
  pop_strip_colors <- lapply(pop_pal[pop_levels], function(col) {
    element_rect(fill = alpha(col, 0.7), colour = "grey80")
  })

  exposure_pal <- c("WQS" = "#4E79A7", "QGcomp" = "#F28E2B",
                     "QGcomp-Cox" = "#59A14F")

  p <- plot_data |>
    ggplot(aes(x = neg_log10_p, y = pathway_name, fill = exposure_short)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.6, alpha = 0.85) +
    scale_fill_manual(values = exposure_pal, name = "Exposure") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    ggh4x::facet_grid2(
      category ~ population,
      scales = "free_y", space = "free_y",
      strip = ggh4x::strip_themed(
        background_y = strip_colors,
        background_x = pop_strip_colors
      )
    ) +
    labs(
      y = "Pathway",
      x = expression(-log[10](p-value)),
      title = paste0("Pathway enrichment — ", covar_label)
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      axis.line = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y = element_text(face = "bold.italic", size = 10, angle = 0),
      panel.spacing.x = unit(0.3, "lines"),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    )

  return(p)
}

p_bar <- create_pathway_bar(covar_data, "Primary analysis")
if (!is.null(p_bar)) {
  ggsave(
    filename = here::here("figures", "test_pathway",
                          "option_B_grouped_bar.png"),
    plot = p_bar,
    width = 16,
    height = max(8, dplyr::n_distinct(covar_data$pathway_name) * 0.35),
    dpi = 300
  )
}


# =============================================================================
# OPTION C: LOLLIPOP CHART
# =============================================================================
# Lollipop (segment + point): cleaner than bars, less ink.
# Point color = exposure, faceted by category (rows) × population (cols).

create_pathway_lollipop <- function(data, covar_label) {

  if (nrow(data) == 0) return(NULL)

  pathway_summary <- data |>
    dplyr::group_by(pathway_name, category) |>
    dplyr::summarise(min_p = min(p_value, na.rm = TRUE), .groups = "drop")

  category_order <- pathway_summary |>
    dplyr::group_by(category) |>
    dplyr::summarise(cat_min_p = min(min_p), .groups = "drop") |>
    dplyr::arrange(cat_min_p) |>
    dplyr::pull(category)

  pathway_sorted <- pathway_summary |>
    dplyr::mutate(category = factor(category, levels = category_order)) |>
    dplyr::arrange(category, desc(min_p)) |>
    dplyr::mutate(pathway_name = fct_inorder(pathway_name))

  plot_data <- data |>
    dplyr::mutate(
      category = factor(category, levels = category_order),
      pathway_name = factor(pathway_name,
                            levels = levels(pathway_sorted$pathway_name)),
      neg_log10_p = -log10(p_value)
    )

  cats <- category_order
  pal <- ggthemes::tableau_color_pal("Tableau 10")(length(cats))
  names(pal) <- cats
  strip_colors <- lapply(pal[cats], function(col) {
    element_rect(fill = alpha(col, 0.3), colour = "grey80")
  })

  pop_levels <- c("all", "no demcind", "demcind")
  pop_pal <- setNames(c("#436C85", "#B73F42", "#DE9960"), pop_levels)
  pop_strip_colors <- lapply(pop_pal[pop_levels], function(col) {
    element_rect(fill = alpha(col, 0.7), colour = "grey80")
  })

  exposure_pal <- c("WQS" = "#4E79A7", "QGcomp" = "#F28E2B",
                     "QGcomp-Cox" = "#59A14F")

  p <- plot_data |>
    ggplot(aes(x = neg_log10_p, y = pathway_name, color = exposure_short)) +
    geom_segment(
      aes(x = 0, xend = neg_log10_p, yend = pathway_name),
      position = position_dodge(width = 0.6),
      linewidth = 0.6, alpha = 0.6
    ) +
    geom_point(
      aes(size = enrichment_factor),
      position = position_dodge(width = 0.6),
      alpha = 0.8
    ) +
    scale_color_manual(values = exposure_pal, name = "Exposure") +
    scale_size_continuous(
      range = c(1.5, 5), name = "Enrichment\nfactor",
      breaks = scales::pretty_breaks(3)
    ) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    ggh4x::facet_grid2(
      category ~ population,
      scales = "free_y", space = "free_y",
      strip = ggh4x::strip_themed(
        background_y = strip_colors,
        background_x = pop_strip_colors
      )
    ) +
    labs(
      y = "Pathway",
      x = expression(-log[10](p-value)),
      title = paste0("Pathway enrichment — ", covar_label)
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      axis.line = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y = element_text(face = "bold.italic", size = 10, angle = 0),
      panel.spacing.x = unit(0.3, "lines"),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    )

  return(p)
}

p_lollipop <- create_pathway_lollipop(covar_data, "Primary analysis")
if (!is.null(p_lollipop)) {
  ggsave(
    filename = here::here("figures", "test_pathway",
                          "option_C_lollipop.png"),
    plot = p_lollipop,
    width = 16,
    height = max(8, dplyr::n_distinct(covar_data$pathway_name) * 0.35),
    dpi = 300
  )
}


message("Test pathway figures saved to figures/test_pathway/")


# =============================================================================
# METAPONE BUBBLE HEATMAP
# =============================================================================
# Same bubble-heatmap layout as Option A, but built from metapone results in
# tables/metapone_results/{study}/{population}/{covar_set}/metapone_*.xlsx.
# Metapone column conventions (see scripts/metapone_pathway.R):
#   - pathway, p_value, n_significant_metabolites, n_mapped_metabolites,
#     n_metabolites, lfdr, adjust.p
#   - enrichment_factor := n_significant_metabolites / n_mapped_metabolites
# Filename pattern: metapone_{exposure}_{study}_{population}_{covar_set}.xlsx

read_metapone_file <- function(file_path) {
  path_parts <- unlist(strsplit(file_path, "/"))
  results_idx <- which(path_parts == "metapone_results")
  study <- path_parts[results_idx + 1]
  population <- path_parts[results_idx + 2]
  covar_set <- path_parts[results_idx + 3]
  fname <- path_parts[results_idx + 4]

  # Strip "metapone_" prefix and "_{study}_{population}_{covar_set}.xlsx" suffix
  suffix <- paste0("_", study, "_", population, "_", covar_set, ".xlsx")
  exposure <- sub("^metapone_", "", fname)
  exposure <- sub(paste0(suffix, "$"), "", exposure)

  df <- readxl::read_xlsx(file_path)

  df |>
    dplyr::rename(pathway_name = any_of("pathway")) |>
    dplyr::mutate(
      n_significant_metabolites = as.numeric(n_significant_metabolites),
      n_mapped_metabolites = as.numeric(n_mapped_metabolites),
      p_value = as.numeric(p_value),
      study = study,
      exposure = exposure,
      population = population,
      covar_set = covar_set,
      category = categorize_pathway(pathway_name)
    )
}

metapone_files <- list.files(
  here::here("tables", "metapone_results"),
  pattern = "^metapone_.*\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
) |>
  purrr::keep(~ grepl(
    "metapone_comp_qgcomp_all_total|metapone_comp_wqs_all_total|metapone_comp_qgcomp_cox_all_cox",
    .x
  ))

metapone_all <- metapone_files |>
  purrr::map(read_metapone_file) |>
  purrr::compact() |>
  purrr::list_rbind()

metapone_sig <- metapone_all |>
  dplyr::filter(!is.na(p_value), p_value < 0.05,
                n_mapped_metabolites > 0,
                category != "other") |>
  dplyr::mutate(
    # Metapone reports permutation p-values that hit the floor as 0; recode
    # to 0.0005 (below the smallest non-zero value) so -log10(p) is finite
    # and these pathways receive color rather than rendering as grey NA.
    p_value = ifelse(p_value == 0, 0.0005, p_value)
  ) |>
  dplyr::arrange(p_value) |>
  dplyr::mutate(
    enrichment_factor = n_significant_metabolites / n_mapped_metabolites,
    exposure_short = dplyr::case_when(
      exposure == "comp_wqs_all" & study == "total" ~ "WQS",
      exposure == "comp_qgcomp_all" & study == "total" ~ "QGcomp",
      exposure == "comp_qgcomp_cox_all" & study == "cox" ~ "QGcomp-Cox",
      TRUE ~ exposure
    ),
    population = factor(population,
                        levels = c("all", "no demcind", "demcind"))
  ) |>
  dplyr::filter(is.finite(enrichment_factor))

# Primary analysis
metapone_covar <- metapone_sig |> dplyr::filter(covar_set == "covar")
p_metapone <- create_pathway_bubble_heatmap(metapone_covar, "Primary analysis")
if (!is.null(p_metapone)) {
  ggsave(
    filename = here::here("figures", "test_pathway",
                          "metapone_bubble_heatmap.png"),
    plot = p_metapone,
    width = 16,
    height = max(8, dplyr::n_distinct(metapone_covar$pathway_name) * 0.35),
    dpi = 300
  )
}

# Sensitivity analysis
metapone_covar_sen <- metapone_sig |> dplyr::filter(covar_set == "covar_sen")
p_metapone_sen <- create_pathway_bubble_heatmap(metapone_covar_sen,
                                                "Sensitivity analysis")
if (!is.null(p_metapone_sen)) {
  ggsave(
    filename = here::here("figures", "test_pathway",
                          "metapone_bubble_heatmap_sen.png"),
    plot = p_metapone_sen,
    width = 16,
    height = max(8, dplyr::n_distinct(metapone_covar_sen$pathway_name) * 0.35),
    dpi = 300
  )
}

message("Metapone bubble heatmaps saved to figures/test_pathway/")
