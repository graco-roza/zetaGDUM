# Figure 5 | MS-GDUM predictive performance and environmental contributions
# Run from docs/final_empirical_analysis/:
#   Rscript figure5_environmental_redesign.R tables figures
# Uses committed CSVs only. No common.R, raw data, zetaGDUM, or model refits.
# Dependencies: ggplot2, dplyr, patchwork, scales, svglite (SVG export).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else "tables"
out_dir <- if (length(args) >= 2) args[2] else "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("svglite", quietly = TRUE))

# ---- Style: Figures 2-4 -----------------------------------------------------
ink <- "#252B2D"; petrol <- "#246878"; orange <- "#C45C19"
purple <- "#615A83"; grey <- "#848B8D"; grid_col <- "#E5E9E9"

pred_levels <- c("TP", "TN", "Chl-a", "Geographic distance")
pred_cols <- setNames(c(petrol, orange, purple, grey), pred_levels)
model_levels <- c("Spread only", "Spread + linear position",
                  "Spread + nonlinear position", "Training: linear position")
model_cols <- setNames(c(grey, petrol, orange, ink), model_levels)
model_lines <- setNames(c("solid", "solid", "solid", "dotted"), model_levels)
model_shapes <- setNames(c(16, 17, 15), model_levels[1:3])
ord_cols <- c(`2` = petrol, `5` = orange)

FIG_W <- 183
FIG_H <- 216
base <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(text = element_text(colour = ink),
        axis.text = element_text(size = 7.5, colour = ink),
        axis.title = element_text(size = 8.2, colour = ink),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.title.x = element_text(margin = margin(t = 3)),
        axis.ticks = element_line(colour = ink, linewidth = .25),
        axis.ticks.length = grid::unit(1, "mm"),
        axis.line = element_line(colour = ink, linewidth = .3),
        panel.grid.major = element_line(colour = grid_col, linewidth = .25),
        panel.grid.minor = element_blank(),
        panel.spacing.x = grid::unit(3.5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, colour = ink, face = "bold",
                                  margin = margin(b = 4)),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7.2, colour = ink),
        legend.key.width = grid::unit(6, "mm"),
        legend.key.height = grid::unit(2.7, "mm"),
        legend.spacing.x = grid::unit(1.3, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(1, "mm"),
        plot.title = element_text(size = 9, face = "bold", margin = margin(b = 4)),
        plot.tag = element_text(size = 11, face = "bold", colour = ink),
        plot.margin = margin(2, 3, 1, 2))
no_strip <- theme(strip.text = element_blank(), strip.background = element_blank())
zero <- geom_hline(yintercept = 0, colour = "#9DA5A7", linewidth = .35)

# Four fixed-width columns, with an independent x-range for each order sequence.
# The names are unique internally, while the facet labels show representation only.
lev <- c("BM_inc", "BM_rel", "ZP_inc", "ZP_rel")
repr_lab <- c(BM_inc = "Incidence", BM_rel = "Relative abundance",
              ZP_inc = "Incidence", ZP_rel = "Relative abundance")
add_panel <- function(d) {
  d <- mutate(d, panel = factor(case_when(
    community == "EPA macroinvertebrates" & response == "incidence" ~ "BM_inc",
    community == "EPA macroinvertebrates" & response == "relative" ~ "BM_rel",
    community == "EPA zooplankton" & response == "incidence" ~ "ZP_inc",
    community == "EPA zooplankton" & response == "relative" ~ "ZP_rel",
    TRUE ~ NA_character_), levels = lev))
  stopifnot(!anyNA(d$panel))
  d
}
facets <- function(free_x = TRUE) {
  facet_grid(cols = vars(panel), scales = if (free_x) "free_x" else "fixed",
             space = "fixed", drop = FALSE,
             labeller = labeller(panel = as_labeller(repr_lab)))
}
# Fewer labels on narrow panels; all evaluated orders remain plotted.
order_breaks <- function(lim) {
  if (max(lim) <= 6.5) c(2, 3, 4, 5, 6) else c(2, 4, 6, 8, 10, 12, 14, 16)
}
order_scale <- function() scale_x_continuous(breaks = order_breaks,
                                             expand = expansion(mult = c(.045, .045)))

# Grouped taxon headings; the representation labels are supplied by A's strips.
# A separate header is used instead of requiring a nested-facet extension.
head_data <- data.frame(x = c(.25, .75), label = c("EPA macroinvertebrates", "EPA zooplankton"))
header <- ggplot(head_data, aes(x, .5, label = label)) +
  geom_text(fontface = "bold", colour = ink, size = 3.35) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 3, 0, 2))

# ---- Read and aggregate existing results -----------------------------------
read_input <- function(stem, pool) {
  path <- file.path(in_dir, paste0(stem, "_", pool, ".csv"))
  if (!file.exists(path)) stop("Missing committed input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
min_na <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
max_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

prepare <- function(pool) {
  V <- read_input("fig5_validation_folds", pool)
  I <- read_input("fig5_importance_folds", pool)
  Cv <- read_input("fig5_response_curves", pool)
  stopifnot(all(c("community", "order", "response", "model", "fold", "seed",
                  "D2_heldout", "D2_train") %in% names(V)),
            all(c("community", "order", "response", "fold", "seed", "predictor",
                  "component", "absolute", "D2_heldout") %in% names(I)),
            all(c("community", "order", "response", "quantity", "predictor", "axis",
                  "bin", "x", "predicted", "observed") %in% names(Cv)),
            all(V$pool == pool), all(I$pool == pool), all(Cv$pool == pool))

  # Same aggregation as 21_fig5_figure.R: mean over seeds within each fold,
  # then median held-out D2 and min/max across folds; training reference = mean.
  Dv <- V |>
    group_by(community, order, response, model, fold) |>
    summarise(D2 = mean_na(D2_heldout), tr = mean_na(D2_train), .groups = "drop") |>
    group_by(community, order, response, model) |>
    summarise(heldout = median_na(D2), lo = min_na(D2), hi = max_na(D2),
              train = mean_na(tr), n_valid_folds = sum(is.finite(D2)), .groups = "drop") |>
    mutate(model = factor(model, levels = c("spread_only", "linear_const", "nonlinear_const"),
                          labels = model_levels[1:3])) |> add_panel()
  stopifnot(!anyNA(Dv$model))

  # Original importance aggregation: mean across seeds, then across folds.
  # Negative held-out contributions are retained, not truncated or reweighted.
  Ip <- I |>
    group_by(community, order, response, fold, predictor, component) |>
    summarise(absolute = mean_na(absolute), D2 = mean(D2_heldout), .groups = "drop") |>
    group_by(community, order, response, predictor, component) |>
    summarise(abs_mean = mean_na(absolute), abs_lo = min_na(absolute),
              abs_hi = max_na(absolute), heldout_D2 = mean(D2),
              n_valid_folds = sum(!is.na(absolute)), .groups = "drop") |>
    mutate(rel = ifelse(heldout_D2 >= .02, abs_mean / heldout_D2, NA_real_),
           pred = factor(recode(predictor, TP = "TP", TN = "TN", Chla = "Chl-a",
                                distance = "Geographic distance"), levels = pred_levels),
           term = ifelse(component == "position", "Position", "Spread or distance"),
           term = factor(term, levels = c("Spread or distance", "Position"))) |>
    add_panel()
  stopifnot(!anyNA(Ip$pred), all(Ip$component %in% c("spread", "position", "distance")))

  # The two selected orders and the original ten equal-count bins are unchanged.
  # These are marginal summaries over sampled sets, not isolated partial effects.
  Cc <- Cv |> filter(quantity == "curve", predictor == "TP", order %in% c(2, 5)) |>
    mutate(ord = factor(order, levels = c(2, 5))) |> add_panel()
  stopifnot(nrow(Cc) > 0, all(Cc$axis %in% c("spread", "position")),
            all(is.finite(Cc$x)), all(is.finite(Cc$predicted)), all(is.finite(Cc$observed)),
            all(Cc$predicted >= 0 & Cc$predicted <= 1),
            all(Cc$observed >= 0 & Cc$observed <= 1))
  list(validation = Dv, importance = Ip, curves = Cc)
}
shared <- prepare("shared")
full <- prepare("full")

# Use common y limits across the shared/full figures, including all values.
# The limits are calculated solely from the plotted quantities, not new fits.
all_D <- bind_rows(shared$validation, full$validation)
all_I <- bind_rows(shared$importance, full$importance)
all_C <- bind_rows(shared$curves, full$curves)
padded_limits <- function(x, step) {
  x <- x[is.finite(x)]
  stopifnot(length(x) > 0)
  r <- range(c(x, 0)); c(floor(r[1]/step)*step, ceiling(r[2]/step)*step)
}
lim_D <- padded_limits(c(all_D$heldout, all_D$lo, all_D$hi, all_D$train), .1)
lim_I <- padded_limits(all_I$abs_mean, .02)
lim_R <- padded_limits(all_I$rel, .25)
# Consistent dissimilarity scale for both environmental summaries, but allow a
# tighter window than 0-1 because all displayed values are well above zero.
lim_C <- c(max(0, floor(min(c(all_C$predicted, all_C$observed)) * 10)/10), 1)
# Common within-axis x limits for the two pools and all four community columns.
lim_x_sp <- c(0, ceiling(max(all_C$x[all_C$axis == "spread"]) * 4)/4)
lim_x_po <- c(floor(min(all_C$x[all_C$axis == "position"]) * 4)/4,
              ceiling(max(all_C$x[all_C$axis == "position"]) * 4)/4)

# ---- A: model validation ---------------------------------------------------
make_A <- function(Dv) {
  tr <- Dv |> filter(model == "Spread + linear position") |>
    mutate(model = factor("Training: linear position", levels = model_levels))
  Dv$model <- factor(Dv$model, levels = model_levels)
  ggplot(Dv, aes(order, heldout, colour = model, group = model)) + zero +
    # Fold ranges (not confidence intervals) at the true orders.
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .42, alpha = .48,
                   show.legend = FALSE) +
    geom_line(aes(linetype = model), linewidth = .65) +
    geom_point(aes(shape = model), size = 1.25, show.legend = FALSE) +
    geom_line(data = tr, aes(y = train, linetype = model), linewidth = .55) +
    facets() + order_scale() +
    scale_y_continuous(breaks = breaks_extended(4), expand = expansion(mult = .035)) +
    coord_cartesian(ylim = lim_D) +
    scale_colour_manual(values = model_cols, breaks = model_levels, drop = FALSE) +
    scale_linetype_manual(values = model_lines, breaks = model_levels, drop = FALSE) +
    scale_shape_manual(values = model_shapes, guide = "none") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                 override.aes = list(linetype = unname(model_lines))),
           linetype = "none") +
    labs(tag = "A", x = "Zeta order, k", y = expression("Held-out " * D^2),
         colour = NULL, linetype = NULL) + base
}

# ---- B and supplementary relative importance ------------------------------
make_importance <- function(Ip, relative = FALSE, show_header = FALSE) {
  value <- if (relative) "rel" else "abs_mean"
  ylab <- if (relative) expression("Share of held-out " * D^2) else
    expression("Absolute contribution (" * Delta * D^2 * ")")
  p <- ggplot(Ip, aes(order, .data[[value]], colour = pred, linetype = term,
                      group = interaction(pred, term))) + zero +
    geom_line(linewidth = .65, na.rm = TRUE) +
    geom_point(size = 1.3, na.rm = TRUE, show.legend = FALSE) +
    facets() + order_scale() +
    scale_y_continuous(breaks = breaks_extended(5), expand = expansion(mult = .035)) +
    coord_cartesian(ylim = if (relative) lim_R else lim_I) +
    scale_colour_manual(values = pred_cols, breaks = pred_levels) +
    scale_linetype_manual(values = c("Spread or distance" = "solid", "Position" = "84")) +
    guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = "solid")),
           linetype = guide_legend(order = 2, nrow = 1,
                                    override.aes = list(colour = ink, linetype = c("solid", "dashed")))) +
    labs(tag = if (relative) NULL else "B", x = "Zeta order, k", y = ylab,
         colour = NULL, linetype = NULL) + base +
    theme(legend.box = "vertical")
  if (!show_header) p <- p + no_strip
  if (relative) p <- p + labs(subtitle = expression("Shown only where mean held-out " * D^2 >= 0.02)) +
    theme(plot.subtitle = element_text(size = 8, colour = ink, margin = margin(b = 4)))
  p
}

# ---- C: observed vs fitted dissimilarity along the TP summaries -------------
make_curve <- function(Cc, which_axis, top = FALSE) {
  Cx <- Cc |> filter(axis == which_axis)
  is_spread <- which_axis == "spread"
  xlabel <- if (is_spread) {
    expression("TP spread: mean pairwise " * "|" * Delta * log[10](TP) * "|")
  } else {
    expression("TP position: mean " * log[10](TP~"("*mu*g~L^{-1}*")"))
  }
  p <- ggplot(Cx, aes(x, colour = ord, group = ord)) +
    geom_line(aes(y = predicted), linewidth = .72) +
    geom_point(aes(y = observed), size = 1.5, shape = 1, stroke = .55,
               show.legend = FALSE) +
    facets(free_x = FALSE) +
    scale_x_continuous(breaks = breaks_extended(4), expand = expansion(mult = .025)) +
    scale_y_continuous(breaks = seq(0, 1, .1), labels = label_number(accuracy = .1),
                       expand = expansion(mult = .035)) +
    coord_cartesian(xlim = if (is_spread) lim_x_sp else lim_x_po, ylim = lim_C) +
    scale_colour_manual(values = ord_cols, labels = c("k = 2", "k = 5")) +
    guides(colour = guide_legend(override.aes = list(linetype = "solid"))) +
    labs(tag = if (top) "C" else NULL, x = xlabel, y = "Dissimilarity", colour = NULL) +
    base + no_strip
  if (top) p <- p + theme(legend.position = "none")
  p
}

# A compact, unambiguous shared key for C. Order is encoded by colour;
# lines are fitted bin means; hollow circles are observed bin means.
curve_key_data <- data.frame(x = c(.14, .35, .62, .82),
                             label = c("k = 2", "k = 5", "Fitted", "Observed"))
curve_key <- ggplot() +
  geom_segment(data = data.frame(x = c(.085, .295, .565), xend = c(.125, .335, .605),
                                 col = c(petrol, orange, ink)),
               aes(x = x, xend = xend, y = .5, yend = .5, colour = col), linewidth = .72) +
  geom_point(aes(x = .77, y = .5), colour = ink, shape = 1, size = 1.6, stroke = .55) +
  geom_text(data = curve_key_data, aes(x, .5, label = label), hjust = 0,
            colour = ink, size = 2.65) + scale_colour_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 3, 0, 2))

build <- function(d) {
  A <- make_A(d$validation)
  B <- make_importance(d$importance)
  C1 <- make_curve(d$curves, "spread", top = TRUE)
  C2 <- make_curve(d$curves, "position") + theme(legend.position = "none")
  wrap_plots(header, A, B, C1, C2, curve_key, ncol = 1,
             heights = grid::unit(c(6, 1.22, 1.22, 1, 1, 6),
                                  c("mm", "null", "null", "null", "null", "mm")))
}

export <- function(fig, name, h = FIG_H) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), fig, width = FIG_W, height = h,
         units = "mm", device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".svg")), fig, width = FIG_W, height = h,
         units = "mm", device = svglite::svglite, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".png")), fig, width = FIG_W, height = h,
         units = "mm", dpi = 600, bg = "white")
}

for (pool in c("shared", "full")) {
  d <- if (pool == "shared") shared else full
  name <- if (pool == "shared") "Figure5_main_redesign" else "FigureS5_full_pool_redesign"
  export(build(d), name)
  p_rel <- make_importance(d$importance, relative = TRUE, show_header = TRUE)
  rel_fig <- wrap_plots(header, p_rel, ncol = 1,
                        heights = grid::unit(c(6, 1), c("mm", "null")))
  export(rel_fig, paste0("FigureS5_relative_importance_", pool, "_redesign"), h = 88)
  # Export the exact plotted summaries (including all negative contributions,
  # original fold ranges, and the relative-importance eligibility threshold).
  write.csv(d$validation, file.path(out_dir, paste0(name, "_A_validation.csv")), row.names = FALSE)
  write.csv(d$importance, file.path(out_dir, paste0(name, "_B_importance.csv")), row.names = FALSE)
  write.csv(d$curves, file.path(out_dir, paste0(name, "_C_TP_curves.csv")), row.names = FALSE)
  message("Saved ", name, " and relative-importance supplement for ", pool)
}
