# 02_metrics.R: local metrics for each census (1985 principal, 2015 temporal) and table type (tree records = principal; all stems = sensitivity).
#  - conventional LCBD: Hellinger transformation, Legendre & De Caceres (2013); own computation cross-checked against adespatial::beta.div (with permutation tests, Holm-adjusted)
#  - LCZD (relative abundance and incidence) at orders K, exact analytical (zetaGDUM::zeta_lczd)
#  - characteristic retention order T and conditional retention order T_plus (zeta_local_retention, max_order = NULL, i.e. all orders), T vs T_full check
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr); library(adespatial)})
B <- readRDS("out/bci_matrices.rds"); h <- B$habitat; K <- c(2, 3, 4, 5, 8, 12, 20, 40, 80, 160, 320, 640); NPERM <- 999
lcbd_hel <- function(x) { y <- sqrt(x / rowSums(x)); ss <- scale(y, scale = FALSE)^2; list(LCBD = unname(rowSums(ss) / sum(ss)), SS_total = sum(ss), BDtotal = sum(ss) / (nrow(x) - 1)) }
res <- list(); prof <- list(); chk <- list()
for (nm in names(B$mat)) { x <- B$mat[[nm]]; x <- x[, colSums(x) > 0, drop = FALSE]; n <- nrow(x); cat("\n##", nm, n, "x", ncol(x), "\n"); stopifnot(all(rowSums(x) > 0))
  # --- LCBD (Hellinger) and checks
  hel <- lcbd_hel(x); set.seed(1); bd <- beta.div(x, method = "hellinger", sqrt.D = FALSE, nperm = NPERM, adj = TRUE, save.D = FALSE)
  chk[[nm]] <- data.frame(set = nm, own_vs_adespatial_max_abs_diff = max(abs(hel[["LCBD"]] - bd$LCBD)), BDtotal_own = hel$BDtotal, BDtotal_adespatial = bd$beta[["BDtotal"]], n_sig_p05 = sum(bd$p.LCBD <= .05), n_sig_holm = sum(bd$p.adj <= .05), swamp_sig_p05 = sum(bd$p.LCBD[h$habitat == "swamp"] <= .05), swamp_sig_holm = sum(bd$p.adj[h$habitat == "swamp"] <= .05))
  # incidence LCBD: Sorensen-based (Legendre & De Caceres 2013) on presence-absence
  bdi <- beta.div(1 * (x > 0), method = "sorensen", nperm = 0, save.D = FALSE)
  # --- exact LCZD
  lr <- zeta_lczd(x, K, "relative"); li <- zeta_lczd(x, K, "incidence")
  # --- retention over all orders
  t0 <- Sys.time(); rr <- zeta_local_retention(x, "relative", max_order = NULL); ro <- zeta_local_retention(x, "occurrence", max_order = NULL); cat("retention", format(Sys.time() - t0), "\n")
  chk[[nm]]$T_vs_Tfull_rel <- max(abs(rr$summary$T - rr$summary$T_full), na.rm = TRUE); chk[[nm]]$T_vs_Tfull_occ <- max(abs(ro$summary$T - ro$summary$T_full), na.rm = TRUE)
  site <- data.frame(set = nm, quad = rownames(x), x_m = h$x[match(rownames(x), h$quad)], y_m = h$y[match(rownames(x), h$quad)], habitat = h$habitat[match(rownames(x), h$quad)], individuals = rowSums(x), richness = rowSums(x > 0),
                     LCBD_hel = bd$LCBD, p_LCBD = bd$p.LCBD, p_LCBD_holm = bd$p.adj, LCBD_sorensen = bdi$LCBD, T_rel = rr$summary$T, T_plus_rel = rr$summary$T_plus, T_full_rel = rr$summary$T_full, T_occ = ro$summary$T, T_plus_occ = ro$summary$T_plus, T_full_occ = ro$summary$T_full,
                     zeta2_rel = rr$profile[, 2], zeta2_occ = ro$profile[, 2])
  for (k in K) { a <- lr[lr$order == k, ]; b <- li[li$order == k, ]; site[[paste0("LCZD_rel_", k)]] <- a$LCZD; site[[paste0("LCZD_inc_", k)]] <- b$LCZD; site[[paste0("ED_rel_", k)]] <- a$expected_dissimilarity; site[[paste0("ED_inc_", k)]] <- b$expected_dissimilarity }
  res[[nm]] <- site; prof[[nm]] <- list(relative = rr$profile, occurrence = ro$profile, loss_rel = rr$loss) }
CH <- bind_rows(chk); print(CH, width = 200); write.csv(CH, "out/metric_checks.csv", row.names = FALSE); S <- bind_rows(res); write.csv(S, "out/site_metrics.csv", row.names = FALSE); saveRDS(prof, "out/retention_profiles.rds"); saveRDS(S, "out/site_metrics.rds")
