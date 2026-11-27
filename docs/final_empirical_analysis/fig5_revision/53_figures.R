# 53_figures.R: Figure 5 revision (candidate main figure and supplements). Run from final/:
#   Rscript fig5_revision/53_figures.R fig5_revision/tables fig5_revision/figures
# Reads only the CSV tables written by 50-52 (no model refits). Style: Figures 2-4 (see figure5_environmental_redesign.R). 183 mm wide; PDF (cairo), SVG, 600 dpi PNG.
# Candidate main figure (shared 1,077-lake pool): A paired held-out delta-D2 vs the spread-only baseline (log10 predictors; fold lines + medians), B absolute contributions (log10 predictors,
# nonlinear position, the full model; chosen a priori as the most general configuration, not for its results), C conditional TP responses on the link scale (log10, nonlinear position).
# Plotted values are exported to <figure>_<panel>.csv next to each figure.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales) })
args <- commandArgs(trailingOnly = TRUE); in_dir <- if (length(args) >= 1) args[1] else "fig5_revision/tables"; out_dir <- if (length(args) >= 2) args[2] else "fig5_revision/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); stopifnot(requireNamespace("svglite", quietly = TRUE))
ink <- "#252B2D"; petrol <- "#246878"; orange <- "#C45C19"; purple <- "#615A83"; grey <- "#848B8D"; grid_col <- "#E5E9E9"
pred_levels <- c("TP", "TN", "Chl-a", "Geographic distance"); pred_cols <- setNames(c(petrol, orange, purple, grey), pred_levels)
CON <- c(delta_linear_vs_spread = "Spread + linear position", delta_nonlinear_vs_spread = "Spread + nonlinear position"); CON_COL <- setNames(c(petrol, orange), CON)
FIG_W <- 183
base <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(text = element_text(colour = ink), axis.text = element_text(size = 7.5, colour = ink), axis.title = element_text(size = 8.2, colour = ink),
        axis.title.y = element_text(margin = margin(r = 4)), axis.title.x = element_text(margin = margin(t = 3)), axis.ticks = element_line(colour = ink, linewidth = .25),
        axis.ticks.length = grid::unit(1, "mm"), axis.line = element_line(colour = ink, linewidth = .3), panel.grid.major = element_line(colour = grid_col, linewidth = .25),
        panel.grid.minor = element_blank(), panel.spacing.x = grid::unit(3.5, "mm"), strip.background = element_blank(),
        strip.text = element_text(size = 8, colour = ink, face = "bold", margin = margin(b = 4)), legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7.2, colour = ink), legend.key.width = grid::unit(6, "mm"), legend.key.height = grid::unit(2.7, "mm"), legend.spacing.x = grid::unit(1.3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "mm"), plot.tag = element_text(size = 11, face = "bold", colour = ink), plot.margin = margin(2, 3, 1, 2))
no_strip <- theme(strip.text = element_blank(), strip.background = element_blank()); zero <- geom_hline(yintercept = 0, colour = "#9DA5A7", linewidth = .35)
lev <- c("BM_inc", "BM_rel", "ZP_inc", "ZP_rel"); repr_lab <- c(BM_inc = "Incidence", BM_rel = "Relative abundance", ZP_inc = "Incidence", ZP_rel = "Relative abundance")
add_panel <- function(d) { d <- mutate(d, panel = factor(case_when(community == "EPA macroinvertebrates" & response == "incidence" ~ "BM_inc", community == "EPA macroinvertebrates" & response == "relative" ~ "BM_rel",
    community == "EPA zooplankton" & response == "incidence" ~ "ZP_inc", community == "EPA zooplankton" & response == "relative" ~ "ZP_rel", TRUE ~ NA_character_), levels = lev)); stopifnot(!anyNA(d$panel)); d }
facets <- function(free_x = TRUE) facet_grid(cols = vars(panel), scales = if (free_x) "free_x" else "fixed", space = "fixed", drop = FALSE, labeller = labeller(panel = as_labeller(repr_lab)))
order_breaks <- function(lim) if (max(lim) <= 6.5) c(2, 3, 4, 5, 6) else c(2, 4, 6, 8, 10, 12, 14, 16)
order_scale <- function() scale_x_continuous(breaks = order_breaks, expand = expansion(mult = c(.045, .045)))
head_data <- data.frame(x = c(.25, .75), label = c("EPA macroinvertebrates", "EPA zooplankton"))
header <- ggplot(head_data, aes(x, .5, label = label)) + geom_text(fontface = "bold", colour = ink, size = 3.35) + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(0, 3, 0, 2))
rd <- function(stem, pool) { p <- file.path(in_dir, sprintf("%s_%s.csv", stem, pool)); if (!file.exists(p)) stop("Missing ", p); read.csv(p, stringsAsFactors = FALSE) }
export <- function(fig, name, h) { ggsave(file.path(out_dir, paste0(name, ".pdf")), fig, width = FIG_W, height = h, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".svg")), fig, width = FIG_W, height = h, units = "mm", device = svglite::svglite, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".png")), fig, width = FIG_W, height = h, units = "mm", dpi = 600, bg = "white") }
dump <- function(d, name, panel) write.csv(d, file.path(out_dir, sprintf("%s_%s.csv", name, panel)), row.names = FALSE)
stack <- function(plots, hs, hdr = TRUE) wrap_plots(c(if (hdr) list(header), plots), ncol = 1, heights = grid::unit(c(if (hdr) 6, hs), c(if (hdr) "mm", rep("null", length(hs)))))
pred_lab <- function(x) factor(recode(x, TP = "TP", TN = "TN", Chla = "Chl-a", distance = "Geographic distance"), levels = pred_levels)
cfg_lab <- function(tf, m) paste0(ifelse(tf == "log10", "Log10", "Raw"), " + ", ifelse(m == "linear_const", "linear", "nonlinear"), " position")

# ---------------- panel builders ----------------
make_delta <- function(F, S, tag = "A", transform = "log10", cons = names(CON)) {
  F <- F |> filter(transform == !!transform, contrast %in% cons) |> mutate(contrast = factor(CON[contrast], levels = CON), grp = interaction(fold, contrast)) |> add_panel()
  S <- S |> filter(transform == !!transform, contrast %in% cons) |> mutate(contrast = factor(CON[contrast], levels = CON)) |> add_panel()
  ggplot() + zero + geom_line(data = F, aes(order, delta_D2, group = grp, colour = contrast, linewidth = "Individual spatial folds"), alpha = .38, na.rm = TRUE) +
    geom_point(data = F, aes(order, delta_D2, colour = contrast), size = .45, alpha = .38, na.rm = TRUE, show.legend = FALSE) +
    geom_line(data = S, aes(order, median, colour = contrast, linewidth = "Median across folds"), alpha = .9) + geom_point(data = S, aes(order, median, colour = contrast), size = 1.2, show.legend = FALSE) +
    facets() + order_scale() + scale_y_continuous(breaks = breaks_extended(5), expand = expansion(mult = .035)) + scale_colour_manual(values = CON_COL) + scale_linewidth_manual(values = c("Individual spatial folds" = .3, "Median across folds" = .9)) +
    guides(colour = guide_legend(order = 1, nrow = 1), linewidth = guide_legend(order = 2, nrow = 1, override.aes = list(colour = ink, alpha = 1))) +
    labs(tag = tag, x = "Zeta order, k", y = expression("Paired held-out " * Delta * D^2 * " vs spread only")) + base }
make_importance <- function(Ia, tag = "B", relative = FALSE, rows = NULL) {
  Ia <- Ia |> mutate(pred = pred_lab(predictor), term = factor(ifelse(component == "position", "Position", "Spread or distance"), levels = c("Spread or distance", "Position"))) |> add_panel()
  if (relative) Ia <- Ia |> filter(!is.na(relative_share)) |> mutate(ok = relative_flag == "ok")
  val <- if (relative) "relative_share" else "abs_mean"
  p <- ggplot(Ia, aes(order, .data[[val]], colour = pred, linetype = term, group = interaction(pred, term, if (!is.null(rows)) .data[[rows]] else 1))) + zero + geom_line(linewidth = .62, na.rm = TRUE)
  p <- if (relative) p + geom_point(aes(shape = ok), size = 1.1, na.rm = TRUE, show.legend = FALSE) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) else p + geom_point(size = 1.05, na.rm = TRUE, show.legend = FALSE)
  p <- p + (if (is.null(rows)) facets() else facet_grid(rows = vars(.data[[rows]]), cols = vars(panel), scales = "free_x", labeller = labeller(panel = as_labeller(repr_lab)))) + order_scale() +
    scale_y_continuous(breaks = breaks_extended(5), expand = expansion(mult = .035)) + scale_colour_manual(values = pred_cols, breaks = pred_levels) +
    scale_linetype_manual(values = c("Spread or distance" = "solid", "Position" = "42")) +
    guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = "solid")), linetype = guide_legend(order = 2, nrow = 1, override.aes = list(colour = ink, linetype = c("solid", "42")))) +
    labs(tag = tag, x = "Zeta order, k", y = if (relative) expression("Share of held-out " * D^2 * " (not clipped)") else expression("Absolute contribution (" * Delta * D^2 * ")")) + base + theme(legend.box = "vertical", strip.text.y = element_text(size = 7, angle = -90))
  p }
ORD_COL <- scale_colour_gradient(low = "#9CC3CC", high = "#0E2C36", trans = "log2", breaks = c(2, 3, 4, 6, 8, 12), name = "Zeta order, k", guide = guide_colourbar(barwidth = grid::unit(60, "mm"), barheight = grid::unit(2.2, "mm"), title.position = "top", title.hjust = .5))
axis_lab <- function(v, ax, tf = "log10") { vv <- c(TP = "TP", TN = "TN", Chla = "Chl-a")[[v]]
  if (tf == "log10") { if (ax == "spread") bquote(.(vv) ~ "spread: mean pairwise" ~ "|" * Delta * log[10](.(vv)) * "|") else bquote(.(vv) ~ "position: mean" ~ log[10](.(vv)) ~ "(" * mu * g ~ L^-1 * ")") }
  else { if (ax == "spread") bquote(.(vv) ~ "spread: mean pairwise" ~ "|" * Delta * .(vv) * "|" ~ "(" * mu * g ~ L^-1 * ")") else bquote(.(vv) ~ "position: mean" ~ .(vv) ~ "(" * mu * g ~ L^-1 * ")") } }
make_curve <- function(Gd, v, ax, scale = c("link", "response", "centred"), tag = NULL, legend = FALSE, yl = NULL, xl = NULL) { scale <- match.arg(scale)
  d <- Gd |> filter(range_type == "common", feasible, predictor == v, axis == ax) |> add_panel(); yv <- if (scale == "response") "dissimilarity" else "link"
  if (scale == "centred") d <- d |> group_by(panel, order) |> arrange(x, .by_group = TRUE) |> mutate(link = link - first(link)) |> ungroup()
  ylab <- switch(scale, link = "Model link scale (logit)", response = "Dissimilarity", centred = "Link change from lowest x")
  p <- ggplot(d, aes(x, .data[[yv]], colour = order, group = order)) + geom_line(linewidth = .65) + facets(free_x = TRUE) + ORD_COL +
    scale_x_continuous(breaks = breaks_extended(4), expand = expansion(mult = .03)) + scale_y_continuous(breaks = breaks_extended(5), expand = expansion(mult = .03)) +
    coord_cartesian(ylim = yl %||% range(d[[yv]])) + labs(tag = tag, x = axis_lab(v, ax), y = ylab) + base + no_strip
  if (!legend) p <- p + theme(legend.position = "none"); p }
`%||%` <- function(a, b) if (is.null(a)) b else a
curve_block <- function(Gd, v, scale, tag, legend_last = TRUE) { Gd <- Gd; d <- Gd |> filter(range_type == "common", feasible, predictor == v)
  lim <- function(ax, col) range(d[d$axis == ax, col]); yl <- function(ax) { dd <- d[d$axis == ax, ]; if (scale == "centred") { dd <- dd |> add_panel() |> group_by(panel, order) |> arrange(x, .by_group = TRUE) |> mutate(link = link - first(link)); range(dd$link) } else range(dd[[if (scale == "response") "dissimilarity" else "link"]]) }
  list(make_curve(Gd, v, "spread", scale, tag = tag, yl = yl("spread"), xl = lim("spread", "x")), make_curve(Gd, v, "position", scale, legend = legend_last, yl = yl("position"), xl = lim("position", "x"))) }

# ---------------- build all figures for one pool ----------------
build_pool <- function(pool) { sfx <- if (pool == "shared") "" else "_full_pool"; pre <- if (pool == "shared") "Figure5_revised_main" else "FigureS5r_full_pool_main"
  Fp <- rd("sens_paired_folds", pool); Sp <- rd("sens_paired_summary", pool); Ia <- rd("sens_importance", pool); Pf <- rd("sens_performance", pool); Gd <- rd("sens_curves_grid", pool) |> filter(transform == "log10", position_model == "nonlinear"); Sm <- rd("sens_curves_summary", pool)
  Ibl <- Ia |> filter(transform == "log10", model == "nonlinear_const")
  # main figure
  A <- make_delta(Fp, Sp, "A"); B <- make_importance(Ibl, "B") + no_strip; cb <- curve_block(Gd |> filter(response %in% c("incidence", "relative")), "TP", "link", "C")
  main <- stack(list(A, B, cb[[1]], cb[[2]]), c(1.15, 1.15, 1, 1.15)); export(main, pre, if (pool == "shared") 232 else 232)
  dump(Fp |> filter(transform == "log10", contrast %in% names(CON)), pre, "A_fold_values"); dump(Sp |> filter(transform == "log10", contrast %in% names(CON)), pre, "A_medians"); dump(Ibl, pre, "B_importance"); dump(Gd |> filter(predictor == "TP"), pre, "C_TP_link_curves")
  # S: performance across the four configurations (+ spread-only baselines), training D2, nonlinear vs linear
  Pf <- Pf |> mutate(series = factor(model_label, levels = c("Spread only", "Spread + linear position", "Spread + nonlinear position")), tf = factor(ifelse(transform == "log10", "Log10 predictors", "Raw predictors"), levels = c("Log10 predictors", "Raw predictors"))) |> add_panel()
  scol <- c("Spread only" = grey, "Spread + linear position" = petrol, "Spread + nonlinear position" = orange); dg <- position_dodge(width = .32)
  p1 <- ggplot(Pf, aes(order, heldout_median, colour = series, linetype = tf, group = interaction(series, tf))) + zero + geom_linerange(aes(ymin = heldout_min, ymax = heldout_max), alpha = .3, linewidth = .3, position = dg, show.legend = FALSE) +
    geom_line(linewidth = .6, position = dg) + geom_point(size = 1, position = dg, show.legend = FALSE) + facets() + order_scale() + scale_y_continuous(breaks = breaks_extended(5)) + scale_colour_manual(values = scol) + scale_linetype_manual(values = c("Log10 predictors" = "solid", "Raw predictors" = "42")) +
    guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = "solid")), linetype = guide_legend(order = 2, nrow = 1, override.aes = list(colour = ink))) + labs(tag = "A", x = "Zeta order, k", y = expression("Held-out " * D^2 * " (median; fold range)")) + base
  p2 <- ggplot(Pf |> filter(model != "spread_only"), aes(order, train_mean, colour = series, linetype = tf, group = interaction(series, tf))) + geom_line(linewidth = .6) + geom_point(size = 1, show.legend = FALSE) + facets() + order_scale() + scale_y_continuous(breaks = breaks_extended(4), limits = c(0, NA)) +
    scale_colour_manual(values = scol, drop = TRUE) + scale_linetype_manual(values = c("Log10 predictors" = "solid", "Raw predictors" = "42")) + labs(tag = "B", x = "Zeta order, k", y = expression("Training " * D^2 * " (mean)")) + base + no_strip + theme(legend.position = "none")
  Ft <- Fp |> filter(contrast == "delta_nonlinear_vs_linear") |> mutate(tf = factor(ifelse(transform == "log10", "Log10 predictors", "Raw predictors"), levels = c("Log10 predictors", "Raw predictors"))) |> add_panel(); St <- Sp |> filter(contrast == "delta_nonlinear_vs_linear") |> mutate(tf = factor(ifelse(transform == "log10", "Log10 predictors", "Raw predictors"), levels = c("Log10 predictors", "Raw predictors"))) |> add_panel()
  p3 <- ggplot() + zero + geom_point(data = Ft, aes(order, delta_D2, colour = tf), size = .5, alpha = .4, position = dg, show.legend = FALSE) + geom_line(data = St, aes(order, median, colour = tf, linetype = tf), linewidth = .7, position = dg) + facets() + order_scale() + scale_y_continuous(breaks = breaks_extended(4)) +
    scale_colour_manual(values = c("Log10 predictors" = ink, "Raw predictors" = purple)) + scale_linetype_manual(values = c("Log10 predictors" = "solid", "Raw predictors" = "42")) + labs(tag = "C", x = "Zeta order, k", y = expression("Nonlinear - linear position, " * Delta * D^2)) + base + no_strip
  export(stack(list(p1, p2, p3), c(1.2, 1, 1)), paste0(if (pool == "shared") "FigureS5r_1_performance_four_configurations" else "FigureS5r_full_pool_1_performance_four_configurations"), 215)
  dump(Pf, paste0(if (pool == "shared") "FigureS5r_1_performance_four_configurations" else "FigureS5r_full_pool_1_performance_four_configurations"), "values"); dump(Fp |> filter(contrast == "delta_nonlinear_vs_linear"), paste0(if (pool == "shared") "FigureS5r_1_performance_four_configurations" else "FigureS5r_full_pool_1_performance_four_configurations"), "paired_nonlinear_vs_linear_folds")
  # S: raw versus log importance (four configurations)
  I4 <- Ia |> mutate(cfg = factor(cfg_lab(transform, model), levels = c("Log10 + linear position", "Log10 + nonlinear position", "Raw + linear position", "Raw + nonlinear position")))
  nm2 <- paste0(if (pool == "shared") "FigureS5r_2_importance_raw_vs_log" else "FigureS5r_full_pool_2_importance_raw_vs_log")
  export(stack(list(make_importance(I4, "", rows = "cfg") + theme(plot.tag = element_blank())), 1), nm2, 228); dump(I4, nm2, "values")
  # S: conditional TP curves on response scale, and centred link change
  nm3 <- paste0(if (pool == "shared") "FigureS5r_3_TP_conditional_response_scale" else "FigureS5r_full_pool_3_TP_conditional_response_scale"); cr <- curve_block(Gd, "TP", "response", "A")
  export(stack(list(cr[[1]], cr[[2]]), c(1, 1.15)), nm3, 140); dump(Gd |> filter(predictor == "TP"), nm3, "values")
  nm3b <- sub("3_TP_conditional_response_scale", "3b_TP_conditional_centred_link", nm3); cc_ <- curve_block(Gd, "TP", "centred", "A"); export(stack(list(cc_[[1]], cc_[[2]]), c(1, 1.15)), nm3b, 140)
  # S: TN and Chl-a conditional curves (link scale)
  nm4 <- paste0(if (pool == "shared") "FigureS5r_4_TN_Chla_conditional" else "FigureS5r_full_pool_4_TN_Chla_conditional"); tn <- curve_block(Gd, "TN", "link", "A", legend_last = FALSE); ch <- curve_block(Gd, "Chla", "link", "B", legend_last = TRUE)
  export(stack(list(tn[[1]], tn[[2]], ch[[1]], ch[[2]]), c(1, 1, 1, 1.15)), nm4, 225); dump(Gd |> filter(predictor %in% c("TN", "Chla")), nm4, "values")
  # S: relative importance without clipping
  nm5 <- paste0(if (pool == "shared") "FigureS5r_5_relative_importance_unclipped" else "FigureS5r_full_pool_5_relative_importance_unclipped")
  export(stack(list(make_importance(Ibl, "", relative = TRUE) + theme(plot.tag = element_blank())), 1), nm5, 100); dump(Ibl, nm5, "values")
  # S: link versus response change, ceiling
  Sx <- Sm |> filter(transform == "log10", position_model == "nonlinear", predictor == "TP") |> mutate(axis_lab = factor(ifelse(axis == "spread", "TP spread", "TP position"), levels = c("TP spread", "TP position"))) |> add_panel()
  ax_cols <- c("TP spread" = petrol, "TP position" = orange); mk <- function(y, lab, tag) ggplot(Sx, aes(order, .data[[y]], colour = axis_lab, group = axis_lab)) + zero + geom_line(linewidth = .65) + geom_point(size = 1.05) + facets() + order_scale() + scale_colour_manual(values = ax_cols) + scale_y_continuous(breaks = breaks_extended(4)) + labs(tag = tag, x = "Zeta order, k", y = lab) + base + (if (tag != "A") no_strip else theme())
  nm6 <- paste0(if (pool == "shared") "FigureS5r_6_saturation_link_vs_response" else "FigureS5r_full_pool_6_saturation_link_vs_response")
  g6 <- mk("link_change", "Link-scale change across\ncommon TP range", "A"); g7 <- mk("response_change", "Dissimilarity change across\ncommon TP range", "B") + theme(legend.position = "none"); g8 <- mk("mean_dissimilarity", "Mean fitted dissimilarity\n(ceiling = 1)", "C") + theme(legend.position = "none")
  export(stack(list(g6 + theme(legend.position = "none"), g7, g8 + theme(legend.position = "bottom")), c(1, 1, 1)) + plot_layout(guides = "keep"), nm6, 190); dump(Sx, nm6, "values")
  # S: extreme-value sensitivity (raw vs log, with/without the 99th-percentile lakes)
  fT <- file.path(in_dir, sprintf("sens_extremes_trim_comparison_%s.csv", pool))
  if (file.exists(fT)) { E <- read.csv(fT) |> filter(model == "nonlinear_const") |> pivot_longer(c(heldout_median_all, heldout_median_trim99), names_to = "variant", values_to = "heldout") |> mutate(variant = factor(ifelse(grepl("trim99", variant), "Extreme lakes removed", "All lakes"), levels = c("All lakes", "Extreme lakes removed")), tf = factor(ifelse(transform == "log10", "Log10 predictors", "Raw predictors"), levels = c("Log10 predictors", "Raw predictors"))) |> add_panel()
    e1 <- ggplot(E, aes(order, heldout, colour = tf, linetype = variant, group = interaction(tf, variant))) + zero + geom_line(linewidth = .65) + geom_point(size = 1, show.legend = FALSE) + facets() + order_scale() + scale_y_continuous(breaks = breaks_extended(5)) +
      scale_colour_manual(values = c("Log10 predictors" = petrol, "Raw predictors" = orange)) + scale_linetype_manual(values = c("All lakes" = "solid", "Extreme lakes removed" = "42")) + guides(colour = guide_legend(order = 1, nrow = 1), linetype = guide_legend(order = 2, nrow = 1, override.aes = list(colour = ink))) +
      labs(tag = "A", x = "Zeta order, k", y = expression("Held-out " * D^2 * " (median; nonlinear position)")) + base
    fI <- file.path(in_dir, sprintf("sens_extremes_trim_importance_%s.csv", pool)); Ei <- read.csv(fI) |> filter(model == "nonlinear_const", transform == "raw") |> mutate(pred = pred_lab(predictor), term = factor(ifelse(component == "position", "Position", "Spread or distance"), levels = c("Spread or distance", "Position"))) |> add_panel()
    e2 <- ggplot(Ei, aes(order, change_abs_trim, colour = pred, linetype = term, group = interaction(pred, term))) + zero + geom_line(linewidth = .6, na.rm = TRUE) + facets() + order_scale() + scale_colour_manual(values = pred_cols, breaks = pred_levels) + scale_linetype_manual(values = c("Spread or distance" = "solid", "Position" = "42")) +
      guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = "solid")), linetype = guide_legend(order = 2, nrow = 1, override.aes = list(colour = ink, linetype = c("solid", "42")))) + labs(tag = "B", x = "Zeta order, k", y = expression("Change in raw-scale contribution (" * Delta * D^2 * ")")) + base + no_strip + theme(legend.box = "vertical")
    nm7 <- paste0(if (pool == "shared") "FigureS5r_7_extreme_nutrient_sensitivity" else "FigureS5r_full_pool_7_extreme_nutrient_sensitivity"); export(stack(list(e1, e2), c(1, 1)), nm7, 160); dump(E, nm7, "A_values"); dump(Ei, nm7, "B_values") }
  invisible(NULL) }
for (p in c("shared", "full")) { build_pool(p); message("figures written for pool ", p) }
