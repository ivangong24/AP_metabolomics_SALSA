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
##        2. Alcohol use and physical activity reach the exposure ONLY through
##           neighbourhood position, so the graph draws nSES and Urbanicity as
##           their parents rather than giving either a direct edge into
##           AirToxicants. The reason is that the exposure is a MODELLED
##           RESIDENTIAL concentration, not personal inhaled dose: an
##           individual's drinking or exercising cannot change the LUR surface
##           evaluated at their address. What can produce an association is a
##           shared cause -- socioeconomic position and walkability drive both
##           the behaviour and the ambient level -- and that is what the
##           arrows say.
##
##           The assumption this rests on, stated rather than buried: no
##           residential self-selection on drinking or activity beyond what
##           nSES and RUCA capture. If heavy drinkers sort into neighbourhoods
##           for reasons those two variables miss, the direct edge is real and
##           alcohol belongs in the primary adjustment set. Both are
##           additionally adjusted in the sensitivity set, so the manuscript
##           reports the analysis under either reading.
##
##        3. BMI, diabetes and medication use are PREDICTORS OF THE OUTCOME
##           ONLY. Each is strongly related to the metabolome and none has an
##           edge into a modelled residential exposure, so dagitty does not
##           return any of them in the minimal sufficient set: adjusting for
##           them is a variance argument, not a confounding one.
##
##           Medication is drawn downstream of its own indications (diabetes,
##           BMI, age) rather than as a free-floating covariate, which is what
##           makes it a non-confounder: nothing links prescribing to the
##           ambient concentration at an address except the conditions that
##           prompt it, and those are themselves drawn.
##
##           This assumes no exposure -> BMI and no exposure -> diabetes edge.
##           Published links between air pollution and incident obesity and
##           diabetes -- including prior work in THIS cohort treating T2DM as
##           a mediator between traffic-related exposure and cognition -- would
##           make those edges real, and would make BMI, diabetes and
##           glucose-lowering medication descendants of the exposure. Adjusting
##           for a descendant does not open a backdoor path, so the checks
##           below stay clean either way; what it does is block part of the
##           total effect. That is the reading under which the sensitivity set
##           is partially OVER-adjusted rather than merely more conservative,
##           and it is consistent with the ~4% coefficient attenuation seen
##           when BMI and diabetes are added (median ratio 0.96).
##
##        4. Renal and hepatic function sit on a causal path from exposure to
##           the metabolome, so a mediator should not be adjusted for and their
##           absence from SALSA costs nothing for identification. They are not
##           PURE mediators, though: alcohol, BMI, diabetes, diet, age and
##           smoking also cause them, which makes each a COLLIDER on paths of
##           the form AirToxicants -> HepaticFunction <- Alcohol -> Metabolite.
##           Conditioning on them would therefore open biasing paths rather
##           than close any -- 23 of them on this graph, against 0 under the
##           primary adjustment set. Not having measured them is better than
##           having measured and adjusted for them. Both counts are computed
##           and asserted below rather than stated here and left to rot.
##
##           The failure mode, stated because this reasoning is convenient for
##           us: if renal or hepatic function shares an UNMEASURED common cause
##           with the metabolome -- chronic inflammation and frailty being the
##           obvious candidates at this age -- the mediator reading is wrong
##           and the estimate is confounded in a way we cannot detect.
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

  nSES       -> Alcohol          Alcohol          -> Metabolite
  Urbanicity -> Alcohol
  nSES       -> PhysicalActivity PhysicalActivity -> Metabolite
  Urbanicity -> PhysicalActivity

  AirToxicants -> Metabolite

  BMI      -> Metabolite
  Diabetes -> Metabolite

  Diabetes -> Medication  BMI -> Medication  Age -> Medication
  Medication -> Metabolite

  AirToxicants -> RenalFunction    RenalFunction    -> Metabolite
  AirToxicants -> HepaticFunction  HepaticFunction  -> Metabolite

  Alcohol  -> HepaticFunction      Alcohol  -> RenalFunction
  BMI      -> HepaticFunction      BMI      -> RenalFunction
  Diabetes -> HepaticFunction      Diabetes -> RenalFunction
  Medication -> HepaticFunction    Medication -> RenalFunction
  Diet     -> HepaticFunction
  Age      -> RenalFunction        Smoking  -> RenalFunction

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

## Compare the computed minimal set with what the MWAS actually fits.
## Any disagreement is a hard stop, not a warning -- see below.

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

## Now that alcohol and physical activity reach the exposure only through nSES
## and urbanicity (header note 2), the graph's minimal sufficient set and the
## primary covariate set should agree EXACTLY. Anything else means the graph
## and the analysis have drifted apart, so this is a hard stop rather than a
## warning.
if (length(missing_from_primary) > 0 || length(extra_in_primary) > 0) {
  stop("The DAG and the primary covariate set disagree. Missing from the ",
       "primary model: ",
       if (length(missing_from_primary)) paste(missing_from_primary, collapse = ", ") else "none",
       ". Adjusted but not required: ",
       if (length(extra_in_primary)) paste(extra_in_primary, collapse = ", ") else "none",
       ". Reconcile the graph and the analysis before using this figure.",
       call. = FALSE)
}
message("\nThe primary covariate set IS the graph's minimal sufficient set.")

## Count the biasing paths left open under the primary set, and under the
## primary set plus the two organ-function variables, so header note 4 is a
## computed number rather than a claim. A biasing path is one that is open
## given Z and is not wholly directed from exposure to outcome -- i.e. it
## contains at least one arrowhead pointing back.
biasing_paths <- function(Z) {
  p <- dagitty::paths(dag, "AirToxicants", "Metabolite", Z = Z)
  p$paths[p$open & grepl("<-", p$paths, fixed = TRUE)]
}
open_primary <- biasing_paths(analysis_adjustment)
open_organs  <- biasing_paths(c(analysis_adjustment,
                                "RenalFunction", "HepaticFunction"))
message("Biasing paths open under the primary set: ", length(open_primary))
message("Biasing paths that would be OPENED by adjusting for renal and ",
        "hepatic function: ", length(open_organs))
if (length(open_primary) != 0) {
  stop("The primary adjustment set no longer blocks every backdoor path ",
       "on this graph. Open: ", paste(open_primary, collapse = " | "),
       call. = FALSE)
}
if (length(open_organs) <= length(open_primary)) {
  warning("Adjusting for renal and hepatic function no longer opens paths ",
          "on this graph, so header note 4 is out of date.", call. = FALSE)
}

## Empirical counterpart: the sensitivity set adds alcohol, physical activity,
## BMI and diabetes, and barely moves the estimates. Feature-level
## coefficients from the two covariate sets correlate at r = 0.978 (C18) and
## 0.987 (HILIC), the median coefficient ratio is 0.96, the median SE ratio is
## 1.010, and all 106 HILIC features significant under the primary set keep
## the same sign and stay nominally significant under the sensitivity set (68
## remain at FDR < 0.05). Recomputed in the response to Reviewer 1 comment 12
## from revision_output/tables/imputation/mwas_pooled{,_covar_sen}.xlsx.


# Layout ----------------------------------------------------------------------

## The full graph has 18 nodes and ~50 edges. Drawn node-by-node that is
## unreadable: the confounder arrows alone are 14 long lines fanning across the
## whole panel, and a reader cannot tell which arrow belongs to which variable.
##
## So the FIGURE groups nodes that play the same causal role into plates and
## draws one arrow per group pair, while the dagitty graph above stays the
## full-fidelity object that every adjustment set and path count is computed
## from. A group arrow means "at least one member of this group causes at least
## one member of that group" -- it is a summary, not a claim that every member
## participates, and the caption says so. Coordinates are set by hand because
## the figure is meant to be read in a fixed order.

ROLE_COLORS <- c(
  "Exposure"                        = "#436C85",
  "Outcome"                         = "#B73F42",
  "Confounder (adjusted)"           = "#DCE6ED",
  "Behaviour, blocked by nSES and RUCA" = "#E4DCEA",
  "Outcome predictor (precision)"   = "#E7EFE4",
  "Mediator / collider, unmeasured" = "#F1E2CB",
  "Unmeasured"                      = "#EDEDED",
  "Downstream"                      = "#FFFFFF"
)
ROLE_TEXT <- c(
  "Exposure" = "white", "Outcome" = "white",
  "Confounder (adjusted)" = "#1b2b36",
  "Behaviour, blocked by nSES and RUCA" = "#3b2b47",
  "Outcome predictor (precision)" = "#2f4429",
  "Mediator / collider, unmeasured" = "#5a3a17",
  "Unmeasured" = "#3a3a3a", "Downstream" = "#3a3a3a"
)

nodes <- tibble::tribble(
  ~name,              ~label,                  ~x,     ~y,   ~role,                                 ~grp,
  "Age",              "Age at draw",            1.23,  4.35, "Confounder (adjusted)",               "CONF",
  "Sex",              "Sex",                    2.60,  4.35, "Confounder (adjusted)",               "CONF",
  "Education",        "Education",              3.88,  4.35, "Confounder (adjusted)",               "CONF",
  "Smoking",          "Smoking",                5.33,  4.35, "Confounder (adjusted)",               "CONF",
  "nSES",             "Neighbourhood SES",      7.14,  4.35, "Confounder (adjusted)",               "CONF",
  "Urbanicity",       "Urbanicity (RUCA)",      9.39,  4.35, "Confounder (adjusted)",               "CONF",
  "CalendarTime",     "Visit wave / calendar",  11.60, 4.35, "Confounder (adjusted)",               "CONF",

  "Alcohol",          "Alcohol use",            7.60,  2.25, "Behaviour, blocked by nSES and RUCA", "BEHAV",
  "PhysicalActivity", "Physical activity",      9.85,  2.25, "Behaviour, blocked by nSES and RUCA", "BEHAV",

  "AirToxicants",     "Air toxicant exposure\n5-yr average before draw",
                                                1.55,  0.00, "Exposure",                            "EXPO",
  "Metabolite",       "Metabolic feature\nlog2 abundance",
                                               11.30,  0.00, "Outcome",                             "OUT",

  "RenalFunction",    "Renal function",         4.60, -2.15, "Mediator / collider, unmeasured",     "ORGAN",
  "HepaticFunction",  "Hepatic function",       6.67, -2.15, "Mediator / collider, unmeasured",     "ORGAN",
  "DemCIND",          "Dementia / CIND",        9.90, -2.15, "Downstream",                          "DEM",

  "BMI",              "Body mass index",        1.41, -4.35, "Outcome predictor (precision)",       "OUTARM",
  "Diabetes",         "Diabetes",               3.14, -4.35, "Outcome predictor (precision)",       "OUTARM",
  "Medication",       "Medication use",         4.83, -4.35, "Outcome predictor (precision)",       "OUTARM",
  "Diet",             "Diet & supplements",     6.96, -4.35, "Unmeasured",                          "OUTARM",
  "FastingStatus",    "Fasting status",         9.09, -4.35, "Unmeasured",                          "OUTARM",
  "DrawTime",         "Time of day of draw",   11.26, -4.35, "Unmeasured",                          "OUTARM"
)

LAB_SIZE <- 3.05
CHAR_W   <- 0.088
LINE_H   <- 0.42
nodes <- nodes |>
  dplyr::mutate(
    n_lines = stringr::str_count(label, "\n") + 1,
    max_chr = purrr::map_int(stringr::str_split(label, "\n"), ~ max(nchar(.x))),
    half_w  = max_chr * CHAR_W / 2 + 0.20,
    half_h  = n_lines * LINE_H / 2 + 0.16,
    role    = factor(role, levels = names(ROLE_COLORS)),
    ## Border encodes AVAILABILITY, fill encodes causal role, so the two read
    ## independently: solid = measured and used as drawn, dashed = not recorded
    ## in SALSA, dotted = measured but adjusted only in the sensitivity set.
    border = dplyr::case_when(
      role %in% c("Unmeasured", "Mediator / collider, unmeasured")  ~ "dashed",
      role == "Behaviour, blocked by nSES and RUCA"                 ~ "dotted",
      TRUE                                                          ~ "solid"
    )
  )

## Group plates, sized to their members.
PAD <- 0.30
plates <- nodes |>
  dplyr::filter(grp %in% c("CONF", "BEHAV", "ORGAN", "OUTARM")) |>
  dplyr::group_by(grp) |>
  dplyr::summarise(xmin = min(x - half_w) - PAD, xmax = max(x + half_w) + PAD,
                   ymin = min(y - half_h) - PAD, ymax = max(y + half_h) + PAD,
                   .groups = "drop")

## lx overrides the label's x so that it does not sit under an arrow; NA keeps
## it centred on its plate.
plate_labels <- tibble::tribble(
  ~grp,     ~lab,                                                              ~ly,   ~just, ~lx,
  "CONF",   "ADJUSTED IN THE PRIMARY MODEL  =  the graph's minimal sufficient set",
                                                                                0.36,  "above", NA,
  "BEHAV",  "BEHAVIOUR  —  blocked by nSES and RUCA",
                                                                                0.36,  "above", 7.80,
  "ORGAN",  "MEDIATOR AND COLLIDER  —  adjusting would open biasing paths",
                                                                                0.36,  "below", 4.90,
  "OUTARM", "OUTCOME-ARM ONLY  —  no path to the exposure",
                                                                                0.36,  "below", NA
) |>
  dplyr::left_join(plates, by = "grp") |>
  dplyr::mutate(x = dplyr::coalesce(lx, (xmin + xmax) / 2),
                y = dplyr::if_else(just == "above", ymax + ly, ymin - ly))

## Box geometry for every drawable unit: the four plates plus the three
## singleton nodes that are not in a plate.
boxes <- dplyr::bind_rows(
  plates,
  nodes |>
    dplyr::filter(grp %in% c("EXPO", "OUT", "DEM")) |>
    dplyr::transmute(grp, xmin = x - half_w, xmax = x + half_w,
                     ymin = y - half_h, ymax = y + half_h)
)

## Each display edge names its source and target box and, optionally, the x (or
## y) at which it should leave and enter them. Specifying the anchors by hand is
## what keeps the 15 arrows from piling on top of one another.
## curv != 0 draws the edge with geom_curve instead of geom_segment. Used for
## the one arrow that would otherwise be a full-height straight line down the
## left margin (nSES -> diet); bowing it outward reads as routing around the
## panel rather than as a stray rule.
E <- function(from, to, x0 = NA, y0 = NA, x1 = NA, y1 = NA,
              kind = "Other", curv = 0)
  tibble::tibble(from, to, ax0 = x0, ay0 = y0, ax1 = x1, ay1 = y1, kind, curv)

edges <- dplyr::bind_rows(
  E("CONF",   "EXPO",   x0 =  1.55,               x1 =  1.55),
  E("CONF",   "OUT",    x0 = 11.60,               x1 = 11.55),
  E("CONF",   "BEHAV",  x0 =  9.60,               x1 =  9.60, kind = "Blocked route"),
  E("CONF",   "ORGAN",  x0 =  4.55,               x1 =  4.55),
  E("CONF",   "OUTARM", x0 =  0.12,               x1 =  0.12, curv = 0.30),
  E("BEHAV",  "OUT",    x0 = 10.60,               x1 = 10.90),
  E("BEHAV",  "ORGAN",  x0 =  6.95,               x1 =  6.95),
  E("OUTARM", "OUT",    x0 = 11.10,               x1 = 11.10),
  E("OUTARM", "ORGAN",  x0 =  7.40,               x1 =  7.40),
  E("EXPO",   "ORGAN",  y0 = -0.40,               x1 =  3.70),
  E("ORGAN",  "OUT",    x0 =  7.60,               y1 = -0.40),
  E("EXPO",   "OUT",    y0 =  0.00, y1 =  0.00,   kind = "Effect of interest"),
  E("EXPO",   "DEM",    y0 = -0.34,               y1 = -2.00, kind = "Downstream"),
  E("OUT",    "DEM",    x0 = 10.70,               x1 = 10.45, kind = "Downstream")
)

## Clip the segment joining two boxes to their borders, honouring any anchor
## the edge specified. An anchor fixes the coordinate; the other one is solved
## from the box edge the arrow leaves through.
anchor_on <- function(box, ax, ay, towards_x, towards_y, gap) {
  cx <- (box$xmin + box$xmax) / 2; cy <- (box$ymin + box$ymax) / 2
  x <- if (!is.na(ax)) ax else cx
  y <- if (!is.na(ay)) ay else cy
  if (!is.na(ax) && is.na(ay)) {
    y <- if (towards_y > cy) box$ymax + gap else box$ymin - gap
  } else if (!is.na(ay) && is.na(ax)) {
    x <- if (towards_x > cx) box$xmax + gap else box$xmin - gap
  } else if (is.na(ax) && is.na(ay)) {
    dx <- towards_x - cx; dy <- towards_y - cy
    tx <- if (dx == 0) Inf else (box$xmax - cx) / abs(dx)
    ty <- if (dy == 0) Inf else (box$ymax - cy) / abs(dy)
    t  <- min(tx, ty)
    x <- cx + dx * t; y <- cy + dy * t
    len <- sqrt(dx^2 + dy^2); x <- x + dx / len * gap; y <- y + dy / len * gap
  } else {
    ## both fixed: nudge outward along the line
    y <- y + sign(towards_y - cy) * gap
  }
  c(x, y)
}

edges <- edges |>
  dplyr::rowwise() |>
  dplyr::mutate({
    b0 <- boxes[boxes$grp == from, ]; b1 <- boxes[boxes$grp == to, ]
    c0 <- c((b0$xmin + b0$xmax) / 2, (b0$ymin + b0$ymax) / 2)
    c1 <- c((b1$xmin + b1$xmax) / 2, (b1$ymin + b1$ymax) / 2)
    tgt <- c(if (!is.na(ax1)) ax1 else c1[1], if (!is.na(ay1)) ay1 else c1[2])
    src <- c(if (!is.na(ax0)) ax0 else c0[1], if (!is.na(ay0)) ay0 else c0[2])
    p0 <- anchor_on(b0, ax0, ay0, tgt[1], tgt[2], 0.10)
    p1 <- anchor_on(b1, ax1, ay1, src[1], src[2], 0.12)
    tibble::tibble(x = p0[1], y = p0[2], xend = p1[1], yend = p1[2])
  }) |>
  dplyr::ungroup()

seg <- function(d, ...) geom_segment(
  data = dplyr::filter(d, curv == 0),
  aes(x = x, y = y, xend = xend, yend = yend), ...)

crv <- function(d, ...) {
  d <- dplyr::filter(d, curv != 0)
  if (nrow(d) == 0) return(NULL)
  geom_curve(data = d, aes(x = x, y = y, xend = xend, yend = yend),
             curvature = d$curv[1], ncp = 12, ...)
}

dag_plot <- ggplot() +
  ## plates first, so arrows and boxes sit on top of them
  geom_rect(data = plates,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey97", colour = "grey78", linewidth = 0.34,
            linetype = "solid") +
  geom_text(data = plate_labels, aes(x = x, y = y, label = lab),
            size = 2.72, colour = "grey38", fontface = "bold") +
  seg(dplyr::filter(edges, kind == "Other"),
      arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
      colour = "grey52", linewidth = 0.62) +
  crv(dplyr::filter(edges, kind == "Other"),
      arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
      colour = "grey52", linewidth = 0.62) +
  seg(dplyr::filter(edges, kind == "Downstream"),
      arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
      colour = "grey68", linewidth = 0.50, linetype = "dotted") +
  seg(dplyr::filter(edges, kind == "Blocked route"),
      arrow = arrow(length = unit(0.24, "cm"), type = "closed"),
      colour = "#7B5EA7", linewidth = 0.95) +
  seg(dplyr::filter(edges, kind == "Effect of interest"),
      arrow = arrow(length = unit(0.32, "cm"), type = "closed"),
      colour = "#B73F42", linewidth = 1.25) +
  geom_rect(data = dplyr::filter(nodes, border == "solid"),
            aes(xmin = x - half_w, xmax = x + half_w,
                ymin = y - half_h, ymax = y + half_h, fill = role),
            colour = "grey25", linewidth = 0.42) +
  geom_rect(data = dplyr::filter(nodes, border == "dashed"),
            aes(xmin = x - half_w, xmax = x + half_w,
                ymin = y - half_h, ymax = y + half_h, fill = role),
            colour = "grey35", linewidth = 0.50, linetype = "dashed") +
  geom_rect(data = dplyr::filter(nodes, border == "dotted"),
            aes(xmin = x - half_w, xmax = x + half_w,
                ymin = y - half_h, ymax = y + half_h, fill = role),
            colour = "#2f4f68", linewidth = 0.62, linetype = "dotted") +
  geom_text(data = nodes, aes(x = x, y = y, label = label),
            colour = unname(ROLE_TEXT[as.character(nodes$role)]),
            size = LAB_SIZE, lineheight = 0.95, fontface = "bold") +
  scale_fill_manual(values = ROLE_COLORS, name = NULL, drop = FALSE,
                    breaks = c("Exposure", "Outcome",
                               "Confounder (adjusted)",
                               "Behaviour, blocked by nSES and RUCA",
                               "Outcome predictor (precision)",
                               "Mediator / collider, unmeasured", "Unmeasured",
                               "Downstream")) +
  guides(fill = guide_legend(nrow = 1)) +
  coord_cartesian(xlim = c(-1.25, 13.3), ylim = c(-5.85, 5.80)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8.2),
    legend.key.size = unit(0.42, "cm"),
    legend.margin = margin(t = 4),
    plot.margin = margin(8, 12, 6, 12)
  )

rev_dir("figures", "dag")
out_png <- rev_here("figures", "dag", "dag_covariate_structure.png")
out_pdf <- rev_here("figures", "dag", "dag_covariate_structure.pdf")
out_cap <- rev_here("figures", "dag", "dag_covariate_structure_caption.txt")

ggsave(out_png, dag_plot, width = 13.0, height = 9.9, dpi = 400, bg = "white")
ggsave(out_pdf, dag_plot, width = 13.0, height = 9.9, bg = "white")

## The figure carries no subtitle; everything it used to say lives here, so the
## caption travels with the manuscript rather than being baked into the image.
caption <- paste0(
  "Figure S[X]. Assumed causal structure for the air toxicant and metabolome ",
  "analysis. The bold red arrow is the effect being estimated. Nodes are ",
  "grouped into plates by the causal role they play, and one arrow is drawn ",
  "per pair of plates: a group arrow means that at least one member of the ",
  "source group causes at least one member of the target group, not that ",
  "every member does. The complete graph, from which every adjustment set ",
  "and path count reported here is computed, is given in dagitty syntax in ",
  "the Supporting Information.\n\n",
  "Border encodes availability and fill encodes causal role, so the two read ",
  "independently: solid = measured and adjusted as drawn; dashed = not ",
  "recorded in SALSA; dotted = measured, and additionally adjusted in the ",
  "sensitivity set.\n\n",
  "The primary adjustment set is ",
  paste(sort(analysis_adjustment), collapse = ", "),
  ", which is exactly this graph's minimal sufficient set; it leaves ",
  length(open_primary), " biasing paths open. Purple: alcohol use and ",
  "physical activity have no direct edge into the exposure, because an ",
  "individual's behaviour cannot change the land-use regression surface ",
  "evaluated at their address; their only route runs back up through ",
  "neighbourhood SES and urbanicity, both already adjusted, which is what ",
  "blocks it. Renal and hepatic function are caused by the exposure and also ",
  "by alcohol, BMI, diabetes, medication, diet, age and smoking, so each is ",
  "a mediator and a collider: adjusting for them would open ",
  length(open_organs), " biasing paths rather than close any. The ",
  "outcome-arm variables predict the metabolome but have no path to the ",
  "exposure, so they are not in the minimal sufficient set and adjusting for ",
  "them buys precision rather than bias control. Dementia/CIND is a common ",
  "effect of the exposure and the metabolome and is therefore a collider, ",
  "which is why it is used only for stratification and never as a covariate."
)
writeLines(caption, out_cap)

message("\nDAG figure written:")
message("  ", out_png)
message("  ", out_pdf)
message("  ", out_cap, "  (caption text, for the manuscript)")

#--------------------------------End of the code--------------------------------
