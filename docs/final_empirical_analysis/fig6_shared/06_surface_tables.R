# 06_surface_tables.R: panel-C quantities. Predicted dissimilarity vs TP spread at oligotrophic/mesotrophic/eutrophic positions, restricted to spread bins
# supported at ALL displayed positions (common spread range). Positions: bins nearest to TP = 10, 35, 100 ug/L (set-mean position, log10 scale), accepted only if
# a supported bin lies within 0.15 log10 units (about 40%) of the target. Models: additive (primary), explicit interaction (default penalty), explicit interaction (penalty 1e-4).
source("lib6.R"); TARGET <- c(Oligotrophic = 10, Mesotrophic = 35, Eutrophic = 100); TOL <- .15; BINS <- 12; KP <- c(3, 5)
S0 <- readRDS("out/surfaces.rds"); S1 <- readRDS("out/surfaces_pen1e-4.rds")
get <- function(SF, g, k, m) { z <- Filter(function(z) z$group == g & z$order == k, SF); if (!length(z)) return(NULL); z[[1]]$surfaces |> filter(model == m, bins == BINS) }
out <- list(); sel_tab <- list()
for (g in names(D)) for (k in sort(unique(S0 |> sapply(function(z) z$order)))) { add <- get(S0, g, k, "additive"); if (is.null(add)) next
  sup <- add |> filter(supported); pos <- sort(unique(sup$position)); if (!length(pos)) next
  chosen <- sapply(log10(TARGET), function(t) { p <- pos[which.min(abs(pos - t))]; if (abs(p - t) <= TOL) p else NA }); names(chosen) <- names(TARGET); chosen <- chosen[!is.na(chosen) & !duplicated(chosen)]
  if (length(chosen) < 2) { sel_tab[[length(sel_tab) + 1]] <- data.frame(group = g, order = k, positions = NA, note = "fewer than two supported trophic positions"); next }
  sp <- Reduce(intersect, lapply(chosen, function(p) round(sup$spread[abs(sup$position - p) < 1e-9], 6))); sel_tab[[length(sel_tab) + 1]] <- data.frame(group = g, order = k, positions = paste(names(chosen), round(10^chosen, 1), collapse = "; "), n_common_spread_bins = length(sp), common_spread_min = if (length(sp)) min(sp) else NA, common_spread_max = if (length(sp)) max(sp) else NA, note = if (length(sp) < 3) "common spread support <3 bins" else "")
  if (length(sp) < 2) next
  for (m in c("spread_only", "additive", "interaction", "interaction_pen1e-4")) { s <- if (m == "interaction_pen1e-4") get(S1, g, k, "interaction") else get(S0, g, k, m); if (is.null(s)) next
    for (nm in names(chosen)) { r <- s |> filter(abs(position - chosen[[nm]]) < 1e-9, round(spread, 6) %in% sp, supported); if (nrow(r)) out[[length(out) + 1]] <- data.frame(group = g, order = k, model = m, trophic = nm, TP_position = 10^chosen[[nm]], r[, c("spread", "total", "occupancy", "heterogeneity", "training_sets", "training_sites", "evaluated_sets")]) } } }
SUR <- bind_rows(out); write.csv(SUR, "tables/panelC_predictions.csv", row.names = FALSE); sel <- bind_rows(sel_tab); write.csv(sel, "tables/panelC_position_selection.csv", row.names = FALSE); print(sel)
# context dependence: change in prediction across the common spread range, per trophic position, and the ratio between highest and lowest displayed position
cd <- SUR |> pivot_longer(c(total, occupancy, heterogeneity), names_to = "component", values_to = "pred") |> group_by(group, order, model, trophic, component) |> arrange(spread, .by_group = TRUE) |>
  summarise(TP_position = first(TP_position), spread_range = diff(range(spread)), delta = last(pred) - first(pred), slope = (last(pred) - first(pred)) / (last(spread) - first(spread)), .groups = "drop")
write.csv(cd, "tables/panelC_context_dependence.csv", row.names = FALSE)
cmp <- cd |> group_by(group, order, model, component) |> arrange(TP_position, .by_group = TRUE) |> summarise(low = first(trophic), high = last(trophic), slope_low = first(slope), slope_high = last(slope), slope_ratio_high_over_low = last(slope) / first(slope), slope_diff = last(slope) - first(slope), .groups = "drop")
write.csv(cmp, "tables/panelC_slope_contrast.csv", row.names = FALSE); options(width = 220); print(cmp |> filter(model == "additive", order %in% KP) |> as.data.frame(), digits = 3)
