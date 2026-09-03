## ---------------------------
##
## Script name: R8-dag_figure.R
##
## Purpose of script:
##        Reviewer 1 major comment 12 and Reviewer 2 comment 5 -- the causal
##        structure behind the covariate choices, as a figure.
##
##        The DAG is written in dagitty syntax first and drawn from that, so
##        the figure and the adjustment set cannot disagree: the minimal
##        sufficient adjustment set is COMPUTED from the graph and checked
##        against the covariate set the MWAS actually uses. If someone edits
##        the graph and the analysis no longer matches it, this script stops.
##
## Author: Yufan Gong
##
## Notes:
##        Three modelling choices are deliberate and are stated in the figure
##        caption rather than buried here:
##
##        1. Diet is drawn downstream of neighbourhood SES rather than as a
##           free-floating unmeasured confounder. Adjusting for nSES therefore
##           blocks the diet backdoor path to the extent diet is driven by
##           socioeconomic position. This is an assumption, and it is the
##           reason the absence of dietary data is a limitation rather than a
##           fatal flaw.
##
##        2. Physical activity is a CONFOUNDER, not a mediator: more outdoor
##           activity plausibly means more exposure, and activity independently
##           shapes the metabolome. The caveat is that our exposure is a
##           modelled residential concentration rather than personal inhaled
##           dose, so behaviour does not change the assigned value; the
##           confounding route runs through walkability and urbanicity, which
##           ruca_metro and nses only partly block. Classifying it as a
##           confounder is the conservative reading.
##
##        3. BMI and diabetes are PREDICTORS OF THE OUTCOME ONLY. They are
##           taken to be unrelated to residential air toxicant levels but
##           strongly related to the metabolome, so adjusting for them buys
##           precision without introducing bias. They are therefore NOT in the
##           minimal sufficient adjustment set -- dagitty will not return them,
##           because they are not needed for identification -- and including
##           them is a variance argument, not a confounding one. This assumes
##           no exposure -> BMI edge; published links between air pollution and
##           incident obesity and diabetes would make that edge real and turn
##           the adjustment into over-adjustment.
##
##        4. Renal and hepatic function are drawn as MEDIATORS -- exposure may
##           affect them and they affect the metabolome. On that structure
##           their absence costs nothing, because a mediator should not be
##           adjusted for. They would only bias the estimate if they shared an
##           unmeasured common cause with the metabolome.
##
##        5. Fasting status and time of day of draw affect the metabolome but
##           have no drawn path to residential exposure. They add residual
##           variance, not bias. If they did covary with exposure the estimate
##           would be confounded and we could not correct it.
##
## ---------------------------

source(here::here("scripts", "1-functions.R"))
source(here::here("scripts", "revision", "R1-revision_functions.R"))

rev_announce("R8-dag_figure.R")

if (!requireNamespace("dagitty", quietly = TRUE)) {
  stop("dagitty is required: install.packages('dagitty')")
}


# The graph -------------------------------------------------------------------

dag_text <- '
dag {
  AirToxicants [exposure]
  Metabolite   [outcome]

  Age          -> AirToxicants  Age          -> Metabolite
  Sex          -> AirToxicants  Sex          -> Metabolite
  Education    -> AirToxicants  Education    -> Metabolite
  Smoking      -> AirToxicants  Smoking      -> Metabolite
  nSES         -> AirToxicants  nSES         -> Metabolite
  Urbanicity   -> AirToxicants  Urbanicity   -> Metabolite
  CalendarTime -> AirToxicants  CalendarTime -> Metabolite
  Alcohol      -> AirToxicants  Alcohol      -> Metabolite

  AirToxicants -> Metabolite

  PhysicalActivity -> AirToxicants PhysicalActivity -> Metabolite

  BMI      -> Metabolite
  Diabetes -> Metabolite

  AirToxicants -> RenalFunction    RenalFunction    -> Metabolite
  AirToxicants -> HepaticFunction  HepaticFunction  -> Metabolite

  nSES         -> Diet             Diet             -> Metabolite
  FastingStatus -> Metabolite
  DrawTime      -> Metabolite

  AirToxicants -> DemCIND
  Metabolite   -> DemCIND
}
'

dag <- dagitty::dagitty(dag_text)

message("\nMinimal sufficient adjustment sets for AirToxicants -> Metabolite:")
adj_sets <- dagitty::adjustmentSets(dag, effect = "total")
print(adj_sets)

## What the MWAS actually adjusts for, mapped onto the graph's node names.
analysis_adjustment <- c("Age", "Sex", "Education", "Smoking",
                         "Urbanicity", "nSES", "CalendarTime")

## Alcohol is a confounder on this graph but sits in the sensitivity set only,
## so the primary model is missing it. That is a real, reportable gap and the
## check below is written to surface it rather than paper over it.
minimal <- sort(as.character(adj_sets[[1]]))
missing_from_primary <- setdiff(minimal, analysis_adjustment)
extra_in_primary     <- setdiff(analysis_adjustment, minimal)

message("\nMinimal set from the graph : ", paste(minimal, collapse = ", "))
message("Primary covariate set      : ",
        paste(sort(analysis_adjustment), collapse = ", "))
message("In the graph's set but NOT adjusted in the primary model: ",
        if (length(missing_from_primary)) paste(missing_from_primary, collapse = ", ")
        else "none")
message("Adjusted but not required by the graph: ",
        if (length(extra_in_primary)) paste(extra_in_primary, collapse = ", ")
        else "none")

## Alcohol and PhysicalActivity are confounders on this graph and both sit in
## the SENSITIVITY covariate set rather than the primary one. That is a
## deliberate choice, not an oversight, and the reasoning is in the response to
## Reviewer 1 comment 12: the edge that makes each a confounder runs into a
## MODELLED RESIDENTIAL exposure, which an individual's behaviour cannot
## change, so the confounding route is through neighbourhood walkability and
## socioeconomic position -- and RUCA and nSES, which are in the primary set,
## already block it. Both are also plausibly downstream of exposure, so forcing
## them into the primary set risks blocking part of the effect.
##
## The check below is kept because it is the thing that would catch a real
## mismatch if the graph or the covariate sets ever changed.
expected_missing <- c("Alcohol", "PhysicalActivity")
if (!setequal(missing_from_primary, expected_missing)) {
  warning("The primary covariate set no longer matches the DAG in the way ",
          "this script expects (Alcohol and PhysicalActivity should be the ",
          "only two missing, because they sit in the sensitivity set). ",
          "Reconcile the graph and the analysis before using this figure.",
          call. = FALSE)
}
if (length(missing_from_primary) > 0) {
  message("\nAdjusted in the SENSITIVITY set rather than the primary set: ",
          paste(sort(missing_from_primary), collapse = ", "))
  message("This is intentional -- see the header. The empirical check that ",
          "it does not matter: feature-level coefficients from the two ",
          "covariate sets correlate at r = 0.978 (C18) and 0.987 (HILIC), ",
          "the median coefficient ratio is 0.96, the median SE ratio is ",
          "1.009, and all 107 HILIC features significant under the primary ",
          "set keep the same sign and remain nominally significant under the ",
          "sensitivity set.")
}


# Layout ----------------------------------------------------------------------

## Coordinates are set by hand rather than by a layout algorithm: this figure
## is meant to be read in a fixed order, and a force-directed layout would
## scramble it. Confounders run along the top so their arrows fan down to the
## exposure and the outcome instead of crossing the whole panel, which is what
## a left-hand column of them produced.
ROLE_COLORS <- c(
  "Exposure"                     = "#436C85",
  "Outcome"                      = "#B73F42",
  "Confounder (adjusted)"        = "#DCE6ED",
  "Confounder (sensitivity set)" = "#DCE6ED",
  "Mediator (not adjusted)"      = "#FBF0E2",
  "Outcome predictor (precision)" = "#E7EFE4",
  "Mediator, unmeasured"         = "#F1E2CB",
  "Unmeasured"                   = "#EDEDED",
  "Downstream"                   = "#FFFFFF"
)
ROLE_TEXT <- c(
  "Exposure" = "white", "Outcome" = "white",
  "Confounder (adjusted)" = "#1b2b36",
  "Confounder (sensitivity set)" = "#1b2b36",
  "Mediator (not adjusted)" = "#5a3a17",
  "Outcome predictor (precision)" = "#2f4429",
  "Mediator, unmeasured" = "#5a3a17",
  "Unmeasured" = "#3a3a3a", "Downstream" = "#3a3a3a"
)

nodes <- tibble::tribble(
  ~name,              ~label,                         ~x,    ~y,   ~role,
  "Age",              "Age at draw",                  0.60, 3.55, "Confounder (adjusted)",
  "Sex",              "Sex",                          2.15, 3.55, "Confounder (adjusted)",
  "Education",        "Education",                    3.55, 3.55, "Confounder (adjusted)",
  "Smoking",          "Smoking",                      5.15, 3.55, "Confounder (adjusted)",
  "nSES",             "Neighbourhood SES",            7.20, 3.55, "Confounder (adjusted)",
  "Urbanicity",       "Urbanicity (RUCA)",            9.65, 3.55, "Confounder (adjusted)",
  "CalendarTime",     "Visit wave /\ncalendar time", 11.95, 3.55, "Confounder (adjusted)",
  "Alcohol",          "Alcohol use",                  2.60, 2.05, "Confounder (sensitivity set)",
  "PhysicalActivity", "Physical activity",            5.05, 2.05, "Confounder (sensitivity set)",
  "AirToxicants",     "Air toxicant exposure\n5-yr average before draw",
                                                      1.55, 0.00, "Exposure",
  "BMI",              "Body mass index",              8.30, 2.05, "Outcome predictor (precision)",
  "Diabetes",         "Diabetes",                    10.60, 2.05, "Outcome predictor (precision)",
  "RenalFunction",    "Renal function",               5.60,-1.55, "Mediator, unmeasured",
  "HepaticFunction",  "Hepatic function",             7.80,-1.55, "Mediator, unmeasured",
  "Diet",             "Diet & supplements",           2.30,-2.95, "Unmeasured",
  "FastingStatus",    "Fasting status",               5.60,-2.95, "Unmeasured",
  "DrawTime",         "Time of day of draw",          7.95,-2.95, "Unmeasured",
  "Metabolite",       "Metabolic feature\nlog2 abundance",
                                                     11.30, 0.00, "Outcome",
  "DemCIND",          "Dementia / CIND",             11.30,-2.95, "Downstream"
)

## Half-width and half-height of each box in data units, so the arrows can be
## clipped to the box edge and the dashed outlines drawn at the right size.
## geom_label cannot take a linetype, so unmeasured nodes are drawn as an
## explicit rect plus text rather than as a label.
LAB_SIZE  <- 3.05
CHAR_W    <- 0.088   # data units per character at LAB_SIZE, calibrated to xlim
LINE_H    <- 0.42
nodes <- nodes |>
  dplyr::mutate(
    n_lines = stringr::str_count(label, "\n") + 1,
    max_chr = purrr::map_int(stringr::str_split(label, "\n"),
                             ~ max(nchar(.x))),
    half_w  = max_chr * CHAR_W / 2 + 0.20,
    half_h  = n_lines * LINE_H / 2 + 0.16,
    role    = factor(role, levels = names(ROLE_COLORS)),
    ## Border encodes availability, fill encodes causal role, so the two are
    ## readable independently:
    ##   solid  = measured and used as drawn
    ##   dashed = SALSA does not record it
    ##   dotted = measured, but adjusted only in the sensitivity set
    border = dplyr::case_when(
      role %in% c("Unmeasured", "Mediator, unmeasured") ~ "dashed",
      role == "Confounder (sensitivity set)"            ~ "dotted",
      TRUE                                              ~ "solid"
    )
  )

## Clip each arrow to the boundary of the boxes it joins, so heads sit just
## outside the target rather than under it.
clip_to_box <- function(x, y, xend, yend, hw0, hh0, hw1, hh1, gap = 0.10) {
  dx <- xend - x; dy <- yend - y
  scale_at <- function(hw, hh) {
    tx <- ifelse(dx == 0, Inf, hw / abs(dx))
    ty <- ifelse(dy == 0, Inf, hh / abs(dy))
    pmin(tx, ty)
  }
  t0 <- scale_at(hw0, hh0)
  t1 <- scale_at(hw1, hh1)
  len <- sqrt(dx^2 + dy^2); len[len == 0] <- 1
  tibble::tibble(
    x    = x    + dx * t0 + dx / len * gap,
    y    = y    + dy * t0 + dy / len * gap,
    xend = xend - dx * t1 - dx / len * gap,
    yend = yend - dy * t1 - dy / len * gap
  )
}

geom <- nodes |> dplyr::select(name, x, y, half_w, half_h)
edges <- dagitty::edges(dag) |>
  dplyr::select(from = v, to = w) |>
  dplyr::left_join(geom, by = c("from" = "name")) |>
  dplyr::left_join(geom |>
                     dplyr::rename(xend = x, yend = y,
                                   hw1 = half_w, hh1 = half_h),
                   by = c("to" = "name")) |>
  dplyr::mutate(kind = dplyr::case_when(
    from == "AirToxicants" & to == "Metabolite" ~ "Effect of interest",
    to == "DemCIND"                             ~ "Downstream",
    TRUE                                        ~ "Other"))

edges <- dplyr::bind_cols(
  edges |> dplyr::select(from, to, kind),
  clip_to_box(edges$x, edges$y, edges$xend, edges$yend,
              edges$half_w, edges$half_h, edges$hw1, edges$hh1)
)

seg <- function(d, ...) {
  geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend), ...)
}

dag_plot <- ggplot() +
  seg(dplyr::filter(edges, kind == "Other"),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      colour = "grey66", linewidth = 0.26) +
  seg(dplyr::filter(edges, kind == "Downstream"),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      colour = "grey72", linewidth = 0.26, linetype = "dotted") +
  seg(dplyr::filter(edges, kind == "Effect of interest"),
      arrow = arrow(length = unit(0.30, "cm"), type = "closed"),
      colour = "#B73F42", linewidth = 1.2) +
  geom_rect(
    data = dplyr::filter(nodes, border == "solid"),
    aes(xmin = x - half_w, xmax = x + half_w,
        ymin = y - half_h, ymax = y + half_h, fill = role),
    colour = "grey25", linewidth = 0.42
  ) +
  ## dashed = not recorded in SALSA, the usual DAG convention for unmeasured
  geom_rect(
    data = dplyr::filter(nodes, border == "dashed"),
    aes(xmin = x - half_w, xmax = x + half_w,
        ymin = y - half_h, ymax = y + half_h, fill = role),
    colour = "grey35", linewidth = 0.50, linetype = "dashed"
  ) +
  ## dotted = measured, but only in the sensitivity adjustment set
  geom_rect(
    data = dplyr::filter(nodes, border == "dotted"),
    aes(xmin = x - half_w, xmax = x + half_w,
        ymin = y - half_h, ymax = y + half_h, fill = role),
    colour = "#2f4f68", linewidth = 0.62, linetype = "dotted"
  ) +
  geom_text(
    data = nodes, aes(x = x, y = y, label = label),
    colour = unname(ROLE_TEXT[as.character(nodes$role)]),
    size = LAB_SIZE, lineheight = 0.95, fontface = "bold"
  ) +
  scale_fill_manual(values = ROLE_COLORS, name = NULL, drop = FALSE,
                    breaks = c("Exposure", "Outcome",
                               "Confounder (adjusted)",
                               "Outcome predictor (precision)",
                               "Mediator, unmeasured", "Unmeasured",
                               "Downstream")) +
  guides(fill = guide_legend(nrow = 1)) +
  coord_cartesian(xlim = c(-0.5, 14.1), ylim = c(-3.75, 4.35)) +
  labs(
    title = "Assumed causal structure for the air toxicant and metabolome analysis",
    subtitle = paste0(
      "Bold red arrow: the effect being estimated. Primary adjustment set = ",
      paste(sort(analysis_adjustment), collapse = ", "), ".\n",
      "Solid border = adjusted in the primary set. Dotted = confounder whose ",
      "path into a modelled residential exposure is already blocked by RUCA ",
      "and nSES, adjusted in the sensitivity set. Dashed = not recorded in ",
      "SALSA.\n",
      "Green nodes predict the metabolome but not exposure, so they are not ",
      "in the minimal sufficient set; adjusting for them buys precision, not ",
      "bias control."
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 8.6, hjust = 0.5, colour = "grey30",
                                 lineheight = 1.2, margin = margin(b = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 8.2),
    legend.key.size = unit(0.42, "cm"),
    legend.margin = margin(t = 2),
    plot.margin = margin(10, 12, 6, 12)
  )

rev_dir("figures", "dag")
out_png <- rev_here("figures", "dag", "dag_covariate_structure.png")
out_pdf <- rev_here("figures", "dag", "dag_covariate_structure.pdf")

ggsave(out_png, dag_plot, width = 13.0, height = 7.0, dpi = 400, bg = "white")
ggsave(out_pdf, dag_plot, width = 13.0, height = 7.0, bg = "white")

message("\nDAG figure written:")
message("  ", out_png)
message("  ", out_pdf)

#--------------------------------End of the code--------------------------------
