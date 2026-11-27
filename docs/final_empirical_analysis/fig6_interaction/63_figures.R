# 63_figures.R: diagnostic figures for the Figure 6 interaction analysis (plotting only; reads 60-62 tables). Run from final/: Rscript fig6_interaction/63_figures.R
# 183 mm wide; PDF (cairo), SVG, 600 dpi PNG. NOT final Figure 6 panels: diagnostic figures for inspection before the final figure is chosen.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(scales)}); stopifnot(requireNamespace("svglite", quietly = TRUE))
IN <- "fig6_interaction/tables"; OUTF <- "fig6_interaction/figures"; dir.create(OUTF, FALSE, TRUE)
ink <- "#252B2D"; petrol <- "#246878"; orange <- "#C45C19"; purple <- "#615A83"; grey <- "#848B8D"; sage <- "#79936C"; grid_col <- "#E5E9E9"
TXL <- c(BM = "EPA macroinvertebrates", ZP = "EPA zooplankton"); PEN_LAB <- c("0.01", "1e-04", "1e-06"); PEN_NUM <- c(.01, 1e-4, 1e-6); PEN_COL <- setNames(c(petrol, orange, purple), PEN_LAB)
VAR_LAB <- c(interaction_all = "Interaction: all predictors", interaction_TP = "Interaction: TP only"); VAR_TYPE <- c(interaction_all = "solid", interaction_TP = "42"); VAR_SHAPE <- c(interaction_all = 16, interaction_TP = 17)
COMP_LAB <- c(total = "Total", occupancy = "Occupancy-associated", heterogeneity = "Abundance heterogeneity"); STAGE_LAB <- c(total = "Total-dissimilarity stage", share = "Heterogeneity-share stage")
POS_COL <- c(Low = sage, Intermediate = petrol, High = orange)
base <- theme_minimal(base_size = 8.5, base_family = "sans") + theme(text = element_text(colour = ink), axis.text = element_text(size = 7, colour = ink), axis.title = element_text(size = 8, colour = ink),
  axis.line = element_line(colour = ink, linewidth = .3), axis.ticks = element_line(colour = ink, linewidth = .25), axis.ticks.length = unit(1, "mm"), panel.grid.major = element_line(colour = grid_col, linewidth = .25), panel.grid.minor = element_blank(),
  panel.spacing = unit(2.2, "mm"), strip.background = element_blank(), strip.text = element_text(size = 7.2, face = "bold", colour = ink), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7.2),
  legend.key.width = unit(7, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"), plot.margin = margin(2, 3, 1, 2), plot.title = element_text(size = 9, face = "bold"))
zero <- geom_hline(yintercept = 0, colour = "#9DA5A7", linewidth = .35)
exp3 <- function(p, name, h) { for (e in c("pdf", "svg", "png")) ggsave(file.path(OUTF, paste0(name, ".", e)), p, width = 183, height = h, units = "mm", dpi = 600, bg = "white", device = switch(e, pdf = grDevices::cairo_pdf, svg = svglite::svglite, png = "png")) }
pen_f <- function(x) factor(format(x, scientific = FALSE, trim = TRUE, digits = 6), levels = format(PEN_NUM, scientific = FALSE, trim = TRUE, digits = 6)); pen_lab_nice <- c("0.01", "0.0001", "0.000001")
pen_nice <- function(x) pen_lab_nice[match(signif(x, 3), signif(PEN_NUM, 3))]
pen_fac <- function(x) factor(pen_nice(x), levels = pen_lab_nice)
order_scale <- function() scale_x_continuous(breaks = function(l) if (max(l) <= 6) 2:5 else c(2, 4, 6, 8, 10, 12), expand = expansion(mult = .06))
for (pool in c("shared", "full")) {
  Ps <- read.csv(file.path(IN, sprintf("interaction_paired_summary_%s.csv", pool))); Pf <- read.csv(file.path(IN, sprintf("interaction_paired_folds_%s.csv", pool)))
  S <- read.csv(file.path(IN, sprintf("interaction_spline_values_%s.csv", pool))); L <- read.csv(file.path(IN, sprintf("interaction_link_decomposition_%s.csv", pool)))
  suf <- if (pool == "shared") "" else "_full_pool"
  # ---- 1. paired interaction gains (held-out and in-sample), all orders, all penalties, both interaction scopes
  rowlab <- function(d) d |> mutate(rowlab = factor(paste0(COMP_LAB[component], ", ", ifelse(split == "test", "held-out", "in-sample")), levels = c(paste0(COMP_LAB[c("total", "occupancy", "heterogeneity")], ", held-out"), paste0(COMP_LAB[c("total", "occupancy", "heterogeneity")], ", in-sample"))), pen = pen_fac(penalty))
  Ps1 <- rowlab(Ps); Pf1 <- rowlab(Pf) |> filter(split == "test")
  p1 <- ggplot() + zero + geom_point(data = Pf1, aes(order, delta_D2, colour = pen, shape = variant), size = .55, alpha = .30, position = position_dodge(width = .45), show.legend = FALSE) +
    geom_line(data = Ps1, aes(order, median, colour = pen, linetype = variant, group = interaction(pen, variant)), linewidth = .55) + geom_point(data = Ps1, aes(order, median, colour = pen, shape = variant), size = 1.1) +
    facet_grid(rows = vars(rowlab), cols = vars(group), scales = "free", labeller = labeller(group = TXL)) + order_scale() + scale_y_continuous(breaks = breaks_extended(4)) +
    scale_colour_manual(values = setNames(PEN_COL, pen_lab_nice), name = "Endpoint penalty") + scale_linetype_manual(values = VAR_TYPE, labels = VAR_LAB, name = NULL) + scale_shape_manual(values = VAR_SHAPE, labels = VAR_LAB, name = NULL) +
    guides(colour = guide_legend(order = 1, nrow = 1), linetype = guide_legend(order = 2, nrow = 1), shape = guide_legend(order = 2, nrow = 1)) +
    labs(x = "Zeta order, k", y = expression("Interaction gain, " * Delta * D^2 * " (interaction - additive)")) + base + theme(legend.title = element_text(size = 7.2), legend.box = "vertical", strip.text.y = element_text(angle = 0, hjust = 0, size = 6.6))
  exp3(p1, paste0("FigureS6i_1_paired_gains", suf), 215)
  # ---- 2. spline shapes: model scale, per taxon x variant x stage; rows = additive + three penalties; columns = orders; colour = TP position
  Sx <- S |> filter(position %in% c("Low", "Intermediate", "High")) |> mutate(position = factor(position, levels = names(POS_COL)))
  for (tx in c("BM", "ZP")) for (vr in c("interaction_all", "interaction_TP")) for (st in c("total", "share")) {
    d <- Sx |> filter(group == tx, stage == st, model == "additive" | variant == vr) |> mutate(mrow = factor(ifelse(model == "additive", "Additive", paste0("Explicit interaction, penalty ", pen_nice(penalty))), levels = c("Additive", paste0("Explicit interaction, penalty ", pen_lab_nice))), ok = TRUE, ord = factor(order))
    p <- ggplot(d, aes(TP_ugL, f_centred, colour = position)) + zero + geom_line(linewidth = .55) + facet_grid(rows = vars(mrow), cols = vars(order), labeller = labeller(order = function(x) paste0("k = ", x)), scales = "fixed") +
      scale_x_log10(breaks = c(10, 100), labels = c("10", "100")) + scale_y_continuous(breaks = breaks_extended(3)) + scale_colour_manual(values = POS_COL, name = "TP position") +
      labs(x = expression("TP (" * mu * g ~ L^-1 * ")"), y = paste0(STAGE_LAB[[st]], ": TP spread shape f(TP), logit scale, centred"), title = paste0(TXL[[tx]], ", ", sub("interaction_", "", sub("_all", " (all predictors)", sub("_TP", " (TP only)", vr))), if (pool != "shared") ", full pool" else "")) + base + theme(legend.title = element_text(size = 7.2), strip.text.y = element_text(angle = 0, hjust = 0, size = 6.4))
    exp3(p, sprintf("FigureS6i_2_spline_shapes_%s_%s_%s%s", tx, sub("interaction_", "", vr), st, suf), 150)
  }
  # ---- 3. explicit interaction contrast on the model scale (High minus Low position), additive = 0 by construction
  W <- S |> filter(position %in% c("Low", "High"), model != "additive") |> select(group, order, model, variant, penalty, stage, position, TP_ugL, f_centred) |> pivot_wider(names_from = position, values_from = f_centred) |> mutate(contrast = High - Low, pen = pen_fac(penalty), stagef = factor(STAGE_LAB[stage], levels = STAGE_LAB))
  for (tx in c("BM", "ZP")) {
    p <- ggplot(W |> filter(group == tx), aes(TP_ugL, contrast, colour = pen, linetype = variant)) + zero + geom_line(linewidth = .55) + facet_grid(rows = vars(stagef), cols = vars(order), labeller = labeller(order = function(x) paste0("k = ", x)), scales = "free_y") +
      scale_x_log10(breaks = c(10, 100), labels = c("10", "100")) + scale_y_continuous(breaks = breaks_extended(3)) + scale_colour_manual(values = setNames(PEN_COL, pen_lab_nice), name = "Endpoint penalty") + scale_linetype_manual(values = VAR_TYPE, labels = VAR_LAB, name = NULL) +
      guides(colour = guide_legend(order = 1, nrow = 1), linetype = guide_legend(order = 2, nrow = 1)) + labs(x = expression("TP (" * mu * g ~ L^-1 * ")"), y = "Explicit interaction: f(high TP position) - f(low TP position), logit scale", title = paste0(TXL[[tx]], if (pool != "shared") ", full pool" else "")) + base + theme(legend.title = element_text(size = 7.2), strip.text.y = element_text(angle = 0, hjust = 0, size = 6.4))
    exp3(p, sprintf("FigureS6i_3_interaction_contrast_%s%s", tx, suf), 110)
  }
  # ---- 4. link-scale decomposition: response-scale spread effect (pp) at the three TP positions; additive model = link-scale effect only
  Ld <- L |> pivot_longer(c(d_total_pp, d_het_pp, d_occ_pp), names_to = "component", values_to = "pp") |> mutate(component = factor(c(d_total_pp = "Total", d_occ_pp = "Occupancy-associated", d_het_pp = "Abundance heterogeneity")[component], levels = c("Total", "Occupancy-associated", "Abundance heterogeneity")),
    position = factor(position, levels = names(POS_COL)), vlab = factor(ifelse(variant == "additive", "Additive (link-scale effect only)", VAR_LAB[variant]), levels = c("Additive (link-scale effect only)", VAR_LAB)), pen = ifelse(is.na(penalty), "additive", pen_nice(penalty)))
  for (tx in c("BM", "ZP")) {
    d <- Ld |> filter(group == tx)
    p <- ggplot(d, aes(position, pp, colour = vlab, linetype = pen, group = model)) + zero + geom_line(linewidth = .55) + geom_point(size = .9) + facet_grid(rows = vars(component), cols = vars(order), labeller = labeller(order = function(x) paste0("k = ", x)), scales = "free_y") +
      scale_colour_manual(values = setNames(c(ink, petrol, orange), levels(d$vlab))) + scale_linetype_manual(values = c(additive = "solid", "0.01" = "solid", "0.0001" = "42", "0.000001" = "13"), breaks = c("0.01", "0.0001", "0.000001"), name = "Endpoint penalty") +
      scale_x_discrete(labels = c("L", "I", "H")) + labs(x = "TP position (L = ~10, I = ~35, H = ~100 µg/L)", y = "Response-scale change for the same TP spread increment (pp)", title = paste0(TXL[[tx]], if (pool != "shared") ", full pool" else "")) + base + theme(legend.title = element_text(size = 7.2), legend.box = "vertical") + guides(colour = guide_legend(order = 1, nrow = 1), linetype = guide_legend(order = 2, nrow = 1))
    exp3(p, sprintf("FigureS6i_4_link_decomposition_%s%s", tx, suf), 120)
  }
  message("figures for ", pool)
}
