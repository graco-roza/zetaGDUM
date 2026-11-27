# 04_analysis.R: tables answering the six questions (tree-record data principal; stems as sensitivity). All quantities are exact/analytical zetaGDUM outputs.
source("bci_lib.R"); options(width = 220); S <- readRDS("out/site_metrics.rds"); B <- readRDS("out/bci_matrices.rds"); h <- B$habitat; K <- c(2, 3, 4, 5, 8, 12, 20, 40, 80, 160, 320, 640)
S$swamp <- S$habitat == "swamp"; pr <- function(v) rank(v, na.last = "keep") / sum(is.finite(v)); rho <- function(a, b) cor(a, b, method = "spearman", use = "complete.obs")
top <- function(v, p = .1) v >= quantile(v, 1 - p, na.rm = TRUE); jac <- function(a, b) sum(a & b) / sum(a | b)
cliff <- function(a, b) { m <- outer(a, b, "-"); (sum(m > 0) - sum(m < 0)) / length(m) }
Q1 <- list(); Q4 <- list(); Q3 <- list(); Q3b <- list(); Q2v <- list(); Q2p <- list(); Q5 <- list(); Q5oh <- list(); Q5sp <- list()
for (nm in unique(S$set)) { d <- S[S$set == nm, ]; x <- B$mat[[nm]]; x <- x[, colSums(x) > 0]; stopifnot(identical(rownames(x), d$quad))
  # ---- Q1/Q4: persistence of distinctiveness across orders; abundance vs incidence
  for (k in K) { lr <- d[[paste0("LCZD_rel_", k)]]; li <- d[[paste0("LCZD_inc_", k)]]
    Q1[[length(Q1) + 1]] <- data.frame(set = nm, order = k, rho_LCBD_LCZDrel = rho(d$LCBD_hel, lr), rho_LCBD_LCZDinc = rho(d$LCBD_hel, li), top10_LCBD_in_top10_LCZDrel = mean(top(lr)[top(d$LCBD_hel)]), top10_LCBD_in_top10_LCZDinc = mean(top(li)[top(d$LCBD_hel)]),
                                      swamp_rank_pct_LCBD = mean(pr(d$LCBD_hel)[d$swamp]), swamp_rank_pct_LCZDrel = mean(pr(lr)[d$swamp]), swamp_rank_pct_LCZDinc = mean(pr(li)[d$swamp]), swamp_in_top10_LCZDrel = mean(top(lr)[d$swamp]), swamp_in_top10_LCZDinc = mean(top(li)[d$swamp]))
    Q4[[length(Q4) + 1]] <- data.frame(set = nm, order = k, rho_rel_inc = rho(lr, li), jaccard_top10 = jac(top(lr), top(li)), jaccard_top10_swamp_excluded = jac(top(lr)[!d$swamp], top(li)[!d$swamp])) }
  # ---- Q3: habitats
  vars <- c("LCBD_hel", "LCZD_rel_2", "LCZD_rel_20", "LCZD_rel_80", "LCZD_inc_2", "LCZD_inc_20", "T_rel", "T_plus_rel", "T_occ", "T_plus_occ", "zeta2_rel", "zeta2_occ", "richness", "individuals")
  Q3[[nm]] <- d |> group_by(set, habitat) |> summarise(n = n(), across(all_of(vars), list(median = ~ median(., na.rm = TRUE), q25 = ~ quantile(., .25, na.rm = TRUE), q75 = ~ quantile(., .75, na.rm = TRUE))), .groups = "drop")
  Q3b[[nm]] <- bind_rows(lapply(vars, function(v) data.frame(set = nm, variable = v, swamp_median = median(d[[v]][d$swamp], na.rm = T), rest_median = median(d[[v]][!d$swamp], na.rm = T), cliffs_delta = cliff(d[[v]][d$swamp], d[[v]][!d$swamp]))))
  # ---- Q2: similar LCBD, different retention
  rk <- pr(d$LCBD_hel); for (v in c("T_rel", "T_plus_rel", "T_occ", "T_plus_occ", "zeta2_rel")) { f <- loess(d[[v]] ~ d$LCBD_hel, span = .5); r <- residuals(f); Q2v[[length(Q2v) + 1]] <- data.frame(set = nm, variable = v, rho_with_LCBD = rho(d$LCBD_hel, d[[v]]), R2_loess_on_LCBD = 1 - var(r) / var(d[[v]]), resid_sd = sd(r), range_q05_q95 = diff(quantile(d[[v]], c(.05, .95)))) }
  ord <- order(d$LCBD_hel); bins <- split(ord, ceiling(seq_along(ord) / 25))   # 50 bins of 25 quadrats with adjacent LCBD
  for (b in bins) for (v in c("T_rel", "T_plus_rel", "T_occ")) { vv <- d[[v]][b]; Q2p[[length(Q2p) + 1]] <- data.frame(set = nm, bin_LCBD_median = median(d$LCBD_hel[b]), bin_LCBD_range_pct = 100 * diff(range(d$LCBD_hel[b])) / median(d$LCBD_hel[b]), variable = v, bin_min = min(vv, na.rm = T), bin_max = max(vv, na.rm = T), bin_range = diff(range(vv, na.rm = T)), n_swamp = sum(d$swamp[b])) }
  # ---- Q5: species contributions (relative abundance and incidence), by habitat group and order
  for (k in c(2, 10, 40)) for (ty in c("relative", "incidence")) { sl <- species_local(x, k, ty)
    for (g in c("swamp", "rest")) { ii <- if (g == "swamp") d$swamp else !d$swamp; cm <- colMeans(sl$C[ii, , drop = FALSE]); om <- colMeans(sl$O[ii, , drop = FALSE]); hm <- colMeans(sl$H[ii, , drop = FALSE])
      Q5[[length(Q5) + 1]] <- data.frame(set = nm, order = k, type = ty, group = g, species = colnames(x), contribution = cm, occupancy_part = om, heterogeneity_part = hm, occ_all = colMeans(x > 0), occ_group = colMeans(x[ii, , drop = FALSE] > 0), mean_q_group = colMeans(sl$q[ii, , drop = FALSE]), mean_q_other = colMeans(sl$q[!ii, , drop = FALSE]), region_occupancy = colSums(x > 0)) }
    for (hb in unique(d$habitat)) { ii <- d$habitat == hb; Q5oh[[length(Q5oh) + 1]] <- data.frame(set = nm, order = k, type = ty, habitat = hb, expected_dissimilarity = mean(rowSums(sl$C[ii, , drop = FALSE])), occupancy_part = mean(rowSums(sl$O[ii, , drop = FALSE])), heterogeneity_part = mean(rowSums(sl$H[ii, , drop = FALSE]))) } }
}
w <- function(l, f) write.csv(bind_rows(l), file.path("tables", f), row.names = FALSE)
w(Q1, "q1_persistence_across_orders.csv"); w(Q4, "q4_abundance_vs_incidence.csv"); w(Q3, "q3_habitat_summary.csv"); w(Q3b, "q3_swamp_vs_rest.csv"); w(Q2v, "q2_variance_in_retention_explained_by_LCBD.csv"); w(Q2p, "q2_LCBD_bins_retention_range.csv"); w(Q5, "q5_species_contributions_swamp_vs_rest.csv"); w(Q5oh, "q5_occupancy_heterogeneity_by_habitat.csv")
# ---- Q2 concrete pairs (tree records): nearest-LCBD pairs with largest contrast in T_rel
pairs <- list(); for (nm in c("1985_tree", "2015_tree")) { d <- S[S$set == nm, ]; o <- order(d$LCBD_hel); for (i in seq_len(nrow(d) - 1)) for (j in (i + 1):min(i + 15, nrow(d))) { a <- o[i]; b <- o[j]; pairs[[length(pairs) + 1]] <- data.frame(set = nm, quad_a = d$quad[a], quad_b = d$quad[b], habitat_a = d$habitat[a], habitat_b = d$habitat[b], LCBD_a = d$LCBD_hel[a], LCBD_b = d$LCBD_hel[b], dLCBD_pct = 100 * abs(d$LCBD_hel[a] - d$LCBD_hel[b]) / mean(d$LCBD_hel[c(a, b)]),
      T_rel_a = d$T_rel[a], T_rel_b = d$T_rel[b], Tplus_rel_a = d$T_plus_rel[a], Tplus_rel_b = d$T_plus_rel[b], T_occ_a = d$T_occ[a], T_occ_b = d$T_occ[b], rich_a = d$richness[a], rich_b = d$richness[b], n_a = d$individuals[a], n_b = d$individuals[b]) } }
P <- bind_rows(pairs) |> mutate(dT_rel = abs(T_rel_a - T_rel_b), dTplus_rel = abs(Tplus_rel_a - Tplus_rel_b), dT_occ = abs(T_occ_a - T_occ_b)); write.csv(P |> filter(dLCBD_pct <= 2) |> arrange(set, desc(dT_rel)) |> group_by(set) |> slice_head(n = 200), "tables/q2_matched_pairs_top.csv", row.names = FALSE)
cat("\nQ1\n"); print(bind_rows(Q1) |> filter(set %in% c("1985_tree", "2015_tree")) |> mutate(across(where(is.numeric), ~ round(., 2))))
cat("\nQ4\n"); print(bind_rows(Q4) |> filter(set %in% c("1985_tree", "2015_tree")) |> mutate(across(where(is.numeric), ~ round(., 2))))
cat("\nQ2 variance\n"); print(bind_rows(Q2v) |> mutate(across(where(is.numeric), ~ round(., 2))))
cat("\nQ3 swamp vs rest\n"); print(bind_rows(Q3b) |> mutate(across(where(is.numeric), ~ round(., 3))) |> filter(set %in% c("1985_tree", "2015_tree")))
cat("\nQ2 pairs (dLCBD<=2%) summary: pairs", sum(P$dLCBD_pct <= 2 & P$set == "1985_tree"), " median |dT_rel|", round(median(P$dT_rel[P$dLCBD_pct <= 2 & P$set == "1985_tree"]), 2), " 95th", round(quantile(P$dT_rel[P$dLCBD_pct <= 2 & P$set == "1985_tree"], .95), 2), "\n")
cat("\nQ5 O/H share by habitat\n"); print(bind_rows(Q5oh) |> filter(set %in% c("1985_tree", "2015_tree"), type == "relative") |> mutate(H_share = heterogeneity_part / expected_dissimilarity) |> select(set, order, habitat, expected_dissimilarity, H_share) |> pivot_wider(names_from = order, values_from = c(expected_dissimilarity, H_share)) |> mutate(across(where(is.numeric), ~ round(., 3))))
