# Figure 3 | Taxon representation across multisite orders
# Run from docs/final_empirical_analysis/:
#   Rscript figure3_species_representation.R tables figures shared
# Arguments: [input table directory] [output figure directory] [shared|full]
# Uses committed CSV files only. No zetaGDUM installation, raw data, or refits.
# Requires ggplot2, dplyr, tidyr, patchwork, scales, svglite.
# A: contribution-weighted occupancy, with site-bootstrap central 95% ranges.
# B: effective number of taxa contributing to similarity, observed vs shuffled.
# C: exact relative-abundance decomposition; dashed line = null median D_H.
# The three columns have different log-order ranges; scales match within columns.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L) args[1] else "tables"
out <- if (length(args) >= 2L) args[2] else "figures"
pool <- if (length(args) >= 3L) args[3] else "shared"
stopifnot(pool %in% c("shared", "full"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_input <- function(stem) {
  f <- file.path(input, paste0(stem, "_", pool, ".csv"))
  if (!file.exists(f)) stop("Missing committed input: ", f)
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
obs <- read_input("fig3_observed")
boot <- read_input("fig3_bootstrap_draws")
perm <- read_input("fig3_permutation_draws")

community_levels <- c("BCI trees (1985)", "EPA macroinvertebrates", "EPA zooplankton")
community_labels <- c("BCI trees", "EPA macroinvertebrates", "EPA zooplankton")
occ_cols <- c("occ_w_sim_inc", "occ_w_sim_rel", "occ_w_dis_inc", "occ_w_dis_rel")
required <- c("community", "order", occ_cols, "neff_sim_inc", "neff_sim_rel", "zeta_rel", "D_O", "D_H")
stopifnot(all(required %in% names(obs)), all(c("community", "order", "draw", occ_cols) %in% names(boot)),
          all(c("community", "order", "permutation", "neff_sim_rel", "D_H") %in% names(perm)),
          setequal(unique(obs$community), community_levels),
          !anyDuplicated(obs[c("community", "order")]),
          all(is.finite(obs$order) & obs$order >= 2),
          all(abs(obs$zeta_rel + obs$D_O + obs$D_H - 1) < 1e-8))

# Ensure the existing 100 site-resampling draws and 499 permutations are present.
keys <- obs |> distinct(community, order)
nb <- boot |> filter(draw > 0) |> group_by(community, order) |> summarise(n = n_distinct(draw), .groups = "drop")
np <- perm |> filter(permutation > 0) |> group_by(community, order) |> summarise(n = n_distinct(permutation), .groups = "drop")
stopifnot(nrow(anti_join(keys, nb, by = c("community", "order"))) == 0L,
          nrow(anti_join(keys, np, by = c("community", "order"))) == 0L,
          all(nb$n == 100L), all(np$n == 499L))
fac <- function(x) mutate(x, community = factor(community, levels = community_levels, labels = community_labels))
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

# A: mean occupancy, weighted by species-specific contributions.
# The incidence dissimilarity weights are the exact unnormalised richness deficits.
occ_long <- function(d) d |>
  select(community, order, any_of("draw"), all_of(occ_cols)) |>
  pivot_longer(all_of(occ_cols), names_to = "metric", values_to = "value")
A_data <- occ_long(obs) |>
  left_join(occ_long(filter(boot, draw > 0)) |>
              group_by(community, order, metric) |>
              summarise(lo = q(value, .025), hi = q(value, .975), .groups = "drop"),
            by = c("community", "order", "metric")) |>
  mutate(representation = ifelse(grepl("_inc$", metric), "Incidence", "Relative abundance"),
         contribution = ifelse(grepl("_sim_", metric), "Similarity", "Dissimilarity")) |>
  fac()
A_data$representation <- factor(A_data$representation, levels = c("Incidence", "Relative abundance"))
A_data$contribution <- factor(A_data$contribution, levels = c("Similarity", "Dissimilarity"))
stopifnot(all(is.finite(A_data$value)), all(A_data$value >= 0 & A_data$value <= 1),
          all(is.finite(A_data$lo) & is.finite(A_data$hi)))

# B: effective number 1/sum(w^2), where w is a normalised similarity weight.
# Within-site shuffling preserves taxa present and the local abundance distribution.
# Permutation zero is the observed dataset, not part of the null distribution.
null <- perm |> filter(permutation > 0) |>
  group_by(community, order) |>
  summarise(neff_med = median(neff_sim_rel), neff_lo = q(neff_sim_rel, .025),
            neff_hi = q(neff_sim_rel, .975), DH_med = median(D_H),
            DH_lo = q(D_H, .025), DH_hi = q(D_H, .975), .groups = "drop")
B_data <- obs |> select(community, order, neff_sim_inc, neff_sim_rel) |>
  pivot_longer(c(neff_sim_inc, neff_sim_rel), names_to = "metric", values_to = "value") |>
  mutate(representation = ifelse(metric == "neff_sim_inc", "Incidence", "Relative abundance")) |>
  fac()
B_data$representation <- factor(B_data$representation, levels = c("Incidence", "Relative abundance"))
N_data <- fac(null)
stopifnot(all(is.finite(B_data$value) & B_data$value > 0),
          all(is.finite(N_data$neff_lo) & N_data$neff_lo > 0),
          all(is.finite(N_data$neff_hi)))

# C: explicit polygon bounds ensure the dashed D_H null curve uses the
# same origin as the observed orange band. Exact closure is checked above.
C_data <- bind_rows(
  transmute(obs, community, order, component = "Abundance heterogeneity", lo = 0, hi = D_H),
  transmute(obs, community, order, component = "Occupancy-associated", lo = D_H, hi = D_H + D_O),
  transmute(obs, community, order, component = "Shared abundance", lo = D_H + D_O, hi = 1)
) |> fac() |> mutate(value = hi - lo,
                     component = factor(component, levels = c("Shared abundance", "Occupancy-associated", "Abundance heterogeneity")))

# Match Figure 2: white page, petrol-blue, burnt orange, restrained grids.
ink <- "#252B2D"; muted <- "#677378"; grid_col <- "#E5E9E9"
col_inc <- "#246878"; col_ab <- "#C45C19"; col_null <- "#777D81"; col_shared <- "#A8B89C"
repr_pal <- c("Incidence" = col_inc, "Relative abundance" = col_ab)
comp_pal <- c("Shared abundance" = col_shared, "Occupancy-associated" = col_inc,
              "Abundance heterogeneity" = col_ab)
base <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(text = element_text(colour = ink),
        axis.title = element_text(size = 8.5), axis.text = element_text(size = 8, colour = ink),
        panel.grid.major = element_line(colour = grid_col, linewidth = .25),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ink, linewidth = .3),
        axis.ticks = element_line(colour = ink, linewidth = .25),
        axis.ticks.length = grid::unit(1, "mm"),
        strip.background = element_blank(), strip.text = element_text(size = 8.5, face = "bold", colour = ink, margin = margin(0, 0, 5, 0)),
        panel.spacing.x = grid::unit(4, "mm"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7.5), legend.key.width = grid::unit(7, "mm"),
        legend.key.height = grid::unit(3.5, "mm"),
        legend.spacing.x = grid::unit(1.5, "mm"), legend.margin = margin(1, 0, 0, 0),
        plot.margin = margin(3, 4, 2, 3), plot.tag = element_text(size = 11, face = "bold"))
# Different informative orders require free x axes, not free y axes.
# All rows use the same data keys and the same log-order breaks.
x_breaks <- function(lim) {
  if (max(lim) <= 6) c(2, 3, 4, 5) else if (max(lim) <= 15) c(2, 3, 5, 8, 12) else c(2, 5, 20, 80, 480)
}
xscale <- function() scale_x_log10(breaks = x_breaks, labels = label_number(accuracy = 1),
                                   expand = expansion(mult = c(.025, .045)))
facets <- function() facet_wrap(vars(community), nrow = 1, scales = "free_x")
no_x <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
no_strip <- theme(strip.text = element_blank())

# A. Representation (colour), contribution (line type).
A <- ggplot(A_data, aes(order, value, colour = representation,
                        linetype = contribution, group = interaction(representation, contribution))) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = representation), colour = NA, alpha = .10, show.legend = FALSE) +
  geom_line(linewidth = .65) + geom_point(size = 1.3, show.legend = FALSE) +
  facets() + xscale() + scale_y_continuous(breaks = seq(0, 1, .25), labels = label_percent(accuracy = 1), limits = c(0, 1.02), expand = expansion(mult = c(.01, .02))) +
  scale_colour_manual(values = repr_pal) + scale_fill_manual(values = repr_pal, guide = "none") +
  scale_linetype_manual(values = c("Similarity" = "dashed", "Dissimilarity" = "solid")) +
  labs(x = NULL, y = "Contribution-weighted\noccupancy", colour = NULL, linetype = NULL) + base + no_x +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid")),
         linetype = guide_legend(order = 2, override.aes = list(colour = ink)))

# B. Show the permutation reference as a grey envelope and dashed median.
# Separate observed incidence and abundance with the same colors as A.
B <- ggplot() +
  geom_ribbon(data = N_data, aes(order, ymin = neff_lo, ymax = neff_hi), fill = col_null, alpha = .20) +
  geom_line(data = N_data, aes(order, neff_med, colour = "Within-site shuffle"), linewidth = .7, linetype = "dashed") +
  geom_line(data = B_data, aes(order, value, colour = representation), linewidth = .65) +
  geom_point(data = B_data, aes(order, value, colour = representation), size = 1.3, show.legend = FALSE) +
  facets() + xscale() + scale_y_log10(breaks = c(1, 3, 10, 30, 100), labels = label_number(accuracy = 1)) +
  scale_colour_manual(values = c(repr_pal, "Within-site shuffle" = col_null),
                      breaks = c("Incidence", "Relative abundance", "Within-site shuffle")) +
  labs(x = NULL, y = "Effective number of taxa\ncontributing to similarity", colour = NULL) + base + no_x + no_strip +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed"))))

# C. Exact decomposition of a unit of relative abundance.
C <- ggplot(C_data, aes(order, ymin = lo, ymax = hi, fill = component, group = component)) +
  geom_ribbon(colour = "white", linewidth = .20) +
  geom_line(data = N_data, aes(order, DH_med, linetype = "Shuffled heterogeneity"), inherit.aes = FALSE,
            colour = ink, linewidth = .6) +
  facets() + xscale() + scale_y_continuous(breaks = seq(0, 1, .25), labels = label_percent(accuracy = 1),
                                           limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = comp_pal) + scale_linetype_manual(values = c("Shuffled heterogeneity" = "dashed")) +
  labs(x = "Zeta order, k (number of sites)", y = "Partition of\nrelative abundance", fill = NULL, linetype = NULL) + base + no_strip +
  guides(fill = guide_legend(order = 1, nrow = 1), linetype = guide_legend(order = 2, nrow = 1)) +
  theme(legend.box = "vertical")

fig <- (A / B / C) + plot_layout(heights = c(1.13, 1, 1.12)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 11, face = "bold")))
name <- if (pool == "shared") "Figure3_main_redesign" else "FigureS3_full_pool_redesign"
ggsave(file.path(out, paste0(name, ".pdf")), fig, width = 183, height = 175, units = "mm", device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(out, paste0(name, ".png")), fig, width = 183, height = 175, units = "mm", dpi = 600, bg = "white")
if (!requireNamespace("svglite", quietly = TRUE)) stop("Install svglite to export editable SVG.")
ggsave(file.path(out, paste0(name, ".svg")), fig, width = 183, height = 175, units = "mm", device = svglite::svglite, bg = "white")

# Retain all plotted data, with exact observed values and the original intervals.
write.csv(A_data, file.path(out, paste0(name, "_A_occupancy.csv")), row.names = FALSE)
write.csv(B_data, file.path(out, paste0(name, "_B_effective_number.csv")), row.names = FALSE)
write.csv(C_data, file.path(out, paste0(name, "_C_components.csv")), row.names = FALSE)
write.csv(N_data, file.path(out, paste0(name, "_null.csv")), row.names = FALSE)
message("Saved ", name, " in ", normalizePath(out))
