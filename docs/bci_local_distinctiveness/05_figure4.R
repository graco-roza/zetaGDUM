# 05_figure4.R: main Figure 4 (2 x 2 blocks): Local compositional distinctiveness and regional representativeness; BCI 1985 census, tree records, 1250 quadrats.
# Every comparison pairs like with like: abundance LCBD (Hellinger, abundance) with relative-abundance retention; occurrence LCBD (Hellinger, presence-absence) with occurrence retention.
source("fig_common.R"); d85 <- S[S$set == "1985_tree", ]; prof <- readRDS("out/retention_profiles.rds")[["1985_tree"]]; d85$LCBDa <- d85$LCBD_hel * 1e3; d85$LCBDo <- d85$LCBD_hel_occ * 1e3
lab_l <- expression("LCBD" ~ (10^-3))
pA1 <- map_layer(d85, "LCBDa", "Abundance", sig = TRUE, option = "mako", dir = -1, lab = lab_l); pA2 <- map_layer(d85, "LCBDo", "Occurrence", sig = TRUE, option = "mako", dir = -1, lab = lab_l, occ = TRUE)
pB1 <- map_layer(d85, "T_rel", "Abundance", option = "viridis", dir = -1, lab = "T"); pB2 <- map_layer(d85, "T_occ", "Occurrence", option = "viridis", dir = -1, lab = "T")
pA1 <- pA1 + labs(tag = "A"); pA2 <- pA2 + labs(tag = ""); pB1 <- pB1 + labs(tag = "B"); pB2 <- pB2 + labs(tag = "")
# ---- C: retention profiles of two quadrats matched on BOTH conventional LCBDs (within 1.2%) and on abundance retention, with very different occurrence retention
ex <- c("x005_y024", "x014_y012"); ord <- 1:640
mk <- function(P, type) { env <- data.frame(k = ord, med = apply(P[, ord], 2, median), lo = apply(P[, ord], 2, quantile, .05), hi = apply(P[, ord], 2, quantile, .95), type = type); q <- bind_rows(lapply(ex, function(e) data.frame(k = ord, v = P[e, ord], quad = e, type = type))); list(env = env, q = q) }
a <- mk(prof$relative, "Relative abundance"); b <- mk(prof$occurrence, "Occurrence"); env <- bind_rows(a$env, b$env); qq <- bind_rows(a$q, b$q); L <- d85[match(ex, d85$quad), ]
qq$quad <- factor(qq$quad, ex, sprintf("%s (%s): LCBD %.2f abund., %.2f occ. (10^-3); T %.1f abund., %.1f occ.", ex, hab_lab[as.character(L$habitat)], L$LCBDa, L$LCBDo, L$T_rel, L$T_occ)); env$type <- factor(env$type, c("Relative abundance", "Occurrence")); qq$type <- factor(qq$type, levels(env$type))
pC <- ggplot() + geom_ribbon(data = env, aes(k, ymin = lo, ymax = hi), fill = pal[["grey"]], alpha = .2) + geom_line(data = env, aes(k, med), colour = pal[["grey"]], linewidth = .5, linetype = 2) + geom_line(data = qq, aes(k, v, colour = quad), linewidth = .7) + facet_wrap(~type, nrow = 1) +
  scale_x_log10(breaks = c(1, 2, 5, 20, 100, 640)) + scale_colour_manual(values = c(pal[["orange"]], pal[["petrol"]]), name = NULL) + labs(x = "Zeta order k", y = "Shared composition (conditional zeta)", subtitle = "Grey: regional median and 5-95% range") + th() + theme(legend.position = "bottom", legend.text = element_text(size = 5.3)) + guides(colour = guide_legend(nrow = 2))
# ---- D: LCBD versus characteristic retention order, matched representations
dd <- bind_rows(d85 |> transmute(habitat, LCBD = LCBDa, T = T_rel, type = "Abundance: LCBD vs T"), d85 |> transmute(habitat, LCBD = LCBDo, T = T_occ, type = "Occurrence: LCBD vs T")); dd$type <- factor(dd$type, unique(dd$type))
rr <- dd |> group_by(type) |> summarise(r = cor(LCBD, T, method = "spearman"), .groups = "drop") |> mutate(lab = sprintf("Spearman rho = %.2f", r))
pC <- pC + labs(tag = "C"); pD <- ggplot(dd, aes(LCBD, T, colour = habitat)) + geom_point(size = .7, alpha = .75) + geom_point(data = dd[dd$habitat == "swamp", ], size = 1.2, colour = pal[["orange"]]) + geom_text(data = rr, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 2.2) + facet_wrap(~type, nrow = 1, scales = "free") +
  scale_colour_manual(values = hab_col, labels = hab_lab, name = NULL) + labs(x = expression("Conventional LCBD" ~ (10^-3)), y = "Characteristic retention order T") + th() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5)) + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.5, alpha = 1)))
pD <- pD + labs(tag = "D")
fig <- wrap_plots(pA1, pB1, pA2, pB2, pC, pD, ncol = 2, heights = c(1, 1, 1.1)) & theme(plot.tag = element_text(face = "bold", size = 9))
sv(fig, "Figure4_main", 7.4, 8.8); cat("done\n")
