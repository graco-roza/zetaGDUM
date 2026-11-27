# 05_figure4.R: main Figure 4 (2 x 2): Local compositional distinctiveness and regional representativeness; BCI 1985 census, tree records, 1250 quadrats.
source("fig_common.R"); d85 <- S[S$set == "1985_tree", ]; prof <- readRDS("out/retention_profiles.rds")[["1985_tree"]]
d85$LCBD3 <- d85$LCBD_hel * 1e3
pA <- map_layer(d85, "LCBD3", "Conventional LCBD (Hellinger)", sig = TRUE, option = "mako", dir = -1, lab = expression("LCBD" ~ (10^-3))) + labs(subtitle = "orange outline: swamp; dots: Holm-significant LCBD")
pB <- map_layer(d85, "T_occ", "Characteristic retention order T (occurrence)", option = "viridis", dir = -1, lab = "T (orders retained)") + labs(subtitle = "orange outline: swamp; dark = low retention")
# ---- C: retention profiles of two quadrats with near-identical LCBD in the swamp but contrasting retention, against the regional envelope
ex <- c("x014_y013", "x020_y014"); ord <- 1:640
mk <- function(P, type) { env <- data.frame(k = ord, med = apply(P[, ord], 2, median), lo = apply(P[, ord], 2, quantile, .05), hi = apply(P[, ord], 2, quantile, .95), type = type)
  q <- bind_rows(lapply(ex, function(e) data.frame(k = ord, v = P[e, ord], quad = e, type = type))); list(env = env, q = q) }
a <- mk(prof$relative, "Relative abundance"); b <- mk(prof$occurrence, "Occurrence"); env <- bind_rows(a$env, b$env); qq <- bind_rows(a$q, b$q); lab <- d85[match(ex, d85$quad), ]
qq$quad <- factor(qq$quad, ex, sprintf("%s (LCBD %.2f e-3, T_occ %.1f)", ex, lab$LCBD3, lab$T_occ)); env$type <- factor(env$type, c("Relative abundance", "Occurrence")); qq$type <- factor(qq$type, levels(env$type))
pC <- ggplot() + geom_ribbon(data = env, aes(k, ymin = lo, ymax = hi), fill = pal[["grey"]], alpha = .2) + geom_line(data = env, aes(k, med), colour = pal[["grey"]], linewidth = .5, linetype = 2) + geom_line(data = qq, aes(k, v, colour = quad), linewidth = .7) + facet_wrap(~type, nrow = 1) +
  scale_x_log10(breaks = c(1, 2, 5, 20, 100, 640)) + scale_colour_manual(values = c(pal[["orange"]], pal[["petrol"]]), name = NULL) + labs(x = "Zeta order k", y = "Shared composition (conditional zeta)", subtitle = "Grey: regional median and 5-95% range") + th() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5)) + guides(colour = guide_legend(nrow = 2))
# ---- D: LCBD versus characteristic retention order (abundance and occurrence)
dd <- bind_rows(d85 |> transmute(habitat, LCBD3, T = T_rel, type = "Relative abundance"), d85 |> transmute(habitat, LCBD3, T = T_occ, type = "Occurrence")); dd$type <- factor(dd$type, levels(env$type))
rr <- dd |> group_by(type) |> summarise(r = cor(LCBD3, T, method = "spearman"), .groups = "drop") |> mutate(lab = sprintf("Spearman rho = %.2f", r))
pD <- ggplot(dd, aes(LCBD3, T, colour = habitat)) + geom_point(size = .7, alpha = .75) + geom_point(data = dd[dd$habitat == "swamp", ], size = 1.2, colour = pal[["orange"]]) + geom_text(data = rr, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 2.2) + facet_wrap(~type, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = hab_col, labels = hab_lab, name = NULL) + labs(x = expression("Conventional LCBD" ~ (10^-3)), y = "Characteristic retention order T") + th() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5)) + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.5, alpha = 1)))
fig <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 9))
sv(fig, "Figure4_main", 7.4, 6.6); cat("done\n")
