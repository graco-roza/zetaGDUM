# Figure 2 | Rarity and species representation
# Run: Rscript figure2_rarity.R [output_directory] [summary_csv] [balanced_forms_csv]
# Dependencies: ggplot2, patchwork, tidyr, dplyr, scales (svglite optional)
# Default data: rounded medians transcribed from manuscript Table S4 (20 simulated communities per scenario; order 2). No simulated data are generated here.
# Submission version: run with the replicate-derived tables written by figure2_summary.R from results/rarity_forms (medians and 25th-75th percentiles over the 20 replicates):
#   Rscript figure2_rarity.R figures figure2_summary.csv figure2_balanced_forms.csv
# Optional summary CSV: one row per scenario with columns scenario, w_occ, w_ab, ratio; optional w_occ_lo/hi, w_ab_lo/hi, ratio_lo/hi.
# Optional balanced CSV: rows for forms Lr-Bn-La and Lr-Bn-Sa with w_occ, w_ab and optional *_lo/_hi.
# No panel titles or subtitles (the dashed reference lines and the bar colour rule belong in the caption: 70% of the species pool in A, 1 = equal effective numbers in B, 12.5% of species per form in C; bar orange when abundance > incidence).
# Dumbbell version: incidence and relative abundance share one row per form (no vertical offset); the bar between them is coloured by the direction of the change.
# Edits to the original script: (1) third argument for the balanced-pool intervals and interval bars in panel C; (2) the shaded band now covers the four large-range forms (ymin 5.5; it was 4.5 and also shaded 'Small / Broad / High'); (3) one legend instead of two; (4) wider right margin in A so '100%' does not touch panel B. Everything else is unchanged.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(dplyr)
  library(scales)
})
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Large/small range; broad/narrow niche; high/low local abundance.
forms <- c("Lr-Bn-La", "Lr-Bn-Sa", "Lr-Nn-La", "Lr-Nn-Sa",
           "Sr-Bn-La", "Sr-Bn-Sa", "Sr-Nn-La", "Sr-Nn-Sa", "Balanced")
labels <- c("Large  /  Broad  /  High", "Large  /  Broad  /  Low",
            "Large  /  Narrow  /  High", "Large  /  Narrow  /  Low",
            "Small  /  Broad  /  High", "Small  /  Broad  /  Low",
            "Small  /  Narrow  /  High", "Small  /  Narrow  /  Low",
            "Balanced pool")

s <- data.frame(
  scenario = forms,
  w_occ = c(.95, .87, .57, .23, .35, .24, .07, .02, NA),
  w_ab  = c(.99, .60, .39, .06, .43, .09, .04, .00, NA),
  ratio = c(.78, .35, .38, .37, .69, .39, .51, .50, .54)
)
if (length(args) >= 2) {
  s <- read.csv(args[2], stringsAsFactors = FALSE)
  stopifnot(all(c("scenario", "w_occ", "w_ab", "ratio") %in% names(s)),
            nrow(s) == 9L, !anyDuplicated(s$scenario), setequal(s$scenario, forms))
  s <- s[match(forms, s$scenario), ]
}
# Optional intervals: absent values are left missing, never manufactured.
for (v in c("w_occ_lo", "w_occ_hi", "w_ab_lo", "w_ab_hi", "ratio_lo", "ratio_hi"))
  if (!v %in% names(s)) s[[v]] <- NA_real_
stopifnot(all(is.finite(s$ratio)), all(s$ratio >= 0),
          all(s$w_occ[1:8] >= 0 & s$w_occ[1:8] <= 1),
          all(s$w_ab[1:8] >= 0 & s$w_ab[1:8] <= 1))
s$row <- 9:1
s$label <- labels
s$label_ratio <- sprintf("%.2f", s$ratio)
# Two selected forms in the balanced scenario (12.5% of species per form).
# Default values are rounded medians as reported in the S4 results text.
bal <- data.frame(form = c("Large range, broad niche, high abundance",
                           "Large range, broad niche, low abundance"),
                  row = c(2, 1), w_occ = c(.57, .28), w_ab = c(.82, .09))
if (length(args) >= 3) {
  b <- read.csv(args[3], stringsAsFactors = FALSE)
  stopifnot(all(c("form", "w_occ", "w_ab") %in% names(b)), all(c("Lr-Bn-La", "Lr-Bn-Sa") %in% b$form))
  b <- b[match(c("Lr-Bn-La", "Lr-Bn-Sa"), b$form), ]
  bal[c("w_occ", "w_ab")] <- b[c("w_occ", "w_ab")]
  for (v in c("w_occ_lo", "w_occ_hi", "w_ab_lo", "w_ab_hi")) bal[[v]] <- if (v %in% names(b)) b[[v]] else NA_real_
}
for (v in c("w_occ_lo", "w_occ_hi", "w_ab_lo", "w_ab_hi")) if (!v %in% names(bal)) bal[[v]] <- NA_real_

# Visual identity, matching the petrol/orange empirical figures.
col_inc <- "#246878"
col_ab  <- "#C45C19"
ink     <- "#252B2D"
muted   <- "#677378"
grid    <- "#E5E9E9"
link    <- "#B8C1C3"
ref     <- "#919DA0"
pal <- c("Incidence" = col_inc, "Relative abundance" = col_ab)
font <- "sans"

th <- theme_minimal(base_size = 9, base_family = font) +
  theme(text = element_text(colour = ink),
        plot.title = element_text(size = 10, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 8, colour = muted, lineheight = 1.08),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8, colour = ink),
        panel.grid.major.x = element_line(colour = grid, linewidth = .3),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = ink, linewidth = .25),
        axis.ticks.length = grid::unit(1.2, "mm"),
        axis.line.x = element_line(colour = ink, linewidth = .3),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin = margin(4, 7, 4, 5))

# Shared vertical scale for A and B: exact row alignment including balanced pool.
y_scale <- scale_y_continuous(breaks = 9:1, labels = labels,
                              limits = c(.55, 9.45), expand = c(0, 0))
pct <- label_percent(accuracy = 1)

# Convert paired representation summaries to long form with optional IQR.
paired <- function(d) {
  bind_rows(transmute(d, row, repr = "Incidence", weight = w_occ,
                      lo = w_occ_lo, hi = w_occ_hi),
            transmute(d, row, repr = "Relative abundance", weight = w_ab,
                      lo = w_ab_lo, hi = w_ab_hi)) |>
    mutate(repr = factor(repr, levels = names(pal)),
           # Slightly separate representations vertically to avoid overlap.
           yy = row)
}
wa <- paired(s[1:8, ])
# A. Share of expected shared composition attributed to dominant rarity form.
dumb <- function(d) mutate(d, dir = ifelse(w_ab >= w_occ, "Relative abundance", "Incidence"))
A <- ggplot() +
  annotate("rect", xmin = -.04, xmax = 1.03, ymin = 5.5, ymax = 9.45,
           fill = "#F5F7F7", colour = NA) +
  geom_vline(xintercept = .70, colour = ref, linewidth = .45, linetype = "22") +
  # dumbbell bar: coloured by the direction of the change (abundance minus incidence)
  geom_segment(data = dumb(s[1:8, ]), aes(x = w_occ, xend = w_ab, y = row, yend = row, colour = dir),
               linewidth = 2.6, alpha = .32, lineend = "round", show.legend = FALSE) +
  geom_segment(data = filter(wa, is.finite(lo), is.finite(hi)),
               aes(x = lo, xend = hi, y = row, yend = row, colour = repr),
               linewidth = .8, alpha = .9, show.legend = FALSE) +
  geom_point(data = wa, aes(weight, row, fill = repr, shape = repr), colour = "white",
             size = 3.3, stroke = .8) +
  annotate("text", x = .01, y = 1, label = "No dominant form",
           colour = muted, hjust = 0, size = 2.6, fontface = "italic") +
  y_scale + scale_x_continuous(breaks = seq(0, 1, .25), labels = pct,
                               limits = c(-.04, 1.03), expand = c(0, 0)) +
  scale_colour_manual(values = pal, drop = FALSE) + scale_fill_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = c("Incidence" = 21, "Relative abundance" = 24)) +
  labs(x = "Share of shared composition", y = "Range / niche / abundance") + th +
  guides(colour = "none", fill = "none", shape = "none") + theme(plot.margin = margin(4, 14, 4, 5))

# B. Effective number of contributing species (abundance/incidence).
B <- ggplot(s, aes(ratio, row)) +
  annotate("rect", xmin = 0, xmax = 1.12, ymin = 5.5, ymax = 9.45,
           fill = "#F5F7F7", colour = NA) +
  geom_vline(xintercept = 1, colour = ref, linewidth = .45, linetype = "22") +
  geom_segment(aes(x = 0, xend = ratio, yend = row), colour = link, linewidth = .65) +
  geom_segment(data = filter(s, is.finite(ratio_lo), is.finite(ratio_hi)),
               aes(x = ratio_lo, xend = ratio_hi, yend = row),
               colour = col_ab, linewidth = .5, alpha = .75) +
  geom_point(colour = col_ab, size = 2.4) +
  geom_text(aes(label = label_ratio), hjust = -.5, size = 2.6, colour = ink) +
  scale_y_continuous(limits = c(.55, 9.45), breaks = 9:1, expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1.12), breaks = c(0, .25, .5, .75, 1), expand = c(0, 0)) +
  labs(x = "Effective species number\n(abundance/incidence)", y = NULL) + th +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# C. Same range and niche, opposite local abundance in the balanced pool.
# This is explicitly a comparison of two selected forms, not a full partition.
bal_long <- bind_rows(transmute(bal, row, repr = "Incidence", weight = w_occ, lo = w_occ_lo, hi = w_occ_hi),
                      transmute(bal, row, repr = "Relative abundance", weight = w_ab, lo = w_ab_lo, hi = w_ab_hi)) |>
  mutate(repr = factor(repr, levels = names(pal)),
         yy = row)
C <- ggplot() +
  geom_vline(xintercept = .125, colour = ref, linewidth = .45, linetype = "22") +
  geom_segment(data = dumb(bal), aes(x = w_occ, xend = w_ab, y = row, yend = row, colour = dir),
               linewidth = 2.6, alpha = .32, lineend = "round", show.legend = FALSE) +
  geom_segment(data = filter(bal_long, is.finite(lo), is.finite(hi), hi > lo),
               aes(x = lo, xend = hi, y = row, yend = row, colour = repr),
               linewidth = .8, alpha = .9, show.legend = FALSE) +
  geom_point(data = bal_long, aes(weight, row, fill = repr, shape = repr), colour = "white", size = 3.3, stroke = .8) +
  geom_text(data = bal_long, aes(weight, row + .22, label = pct(weight)),
            hjust = .5, size = 2.6, colour = ink, show.legend = FALSE) +
  scale_y_continuous(breaks = 2:1, labels = bal$form, limits = c(.55, 2.45), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, .25), labels = pct, limits = c(0, 1.03), expand = c(0, 0)) +
  scale_colour_manual(values = pal, drop = FALSE) + scale_fill_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = c("Incidence" = 21, "Relative abundance" = 24)) +
  labs(x = "Share of shared composition", y = NULL) + th

fig <- (A | B) / C +
  plot_layout(widths = c(1.55, 1), heights = c(3.2, 1.25), guides = "collect") +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 11, face = "bold"))) &
  theme(legend.position = "bottom", plot.tag = element_text(size = 11, face = "bold"))

# Order 2 throughout. Bars show the 25th-75th percentiles across the 20 replicates when intervals are supplied.
ggsave(file.path(out, "Figure2_rarity.pdf"), fig, width = 183, height = 150,
       units = "mm", device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure2_rarity.png"), fig, width = 183, height = 150,
       units = "mm", dpi = 600, bg = "white")
if (requireNamespace("svglite", quietly = TRUE))
  ggsave(file.path(out, "Figure2_rarity.svg"), fig, width = 183, height = 150,
         units = "mm", device = svglite::svglite, bg = "white") else
  ggsave(file.path(out, "Figure2_rarity.svg"), fig, width = 183, height = 150, units = "mm", bg = "white")
write.csv(s, file.path(out, "Figure2_rarity_scenario_values.csv"), row.names = FALSE)
write.csv(bal, file.path(out, "Figure2_rarity_balanced_values.csv"), row.names = FALSE)
