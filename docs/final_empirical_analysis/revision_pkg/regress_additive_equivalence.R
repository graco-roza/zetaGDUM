args <- commandArgs(TRUE); lib <- args[1]; out <- args[2]; .libPaths(c(lib, .libPaths())); suppressMessages(library(zetaGDUM)); new <- exists("ugdm_position_design", asNamespace("zetaGDUM"))
setwd("/Users/graco-roza/Library/CloudStorage/OneDrive-UniversityofHelsinki/Admin/Academic-History/Post-graduation-UERJ/Doutorado/Capitulo IV/Brasil_das_aguas/Brasil_das_aguas/paper_mee/final")
P <- readRDS("out/pools.rds"); d <- P$shared$BM; set.seed(1); keep <- sort(sample(nrow(d$x), 260)); x <- d$x[keep, ]; x <- x[, colSums(x) > 0]; e <- data.frame(TP = log10(d$env$TP[keep]), TN = log10(d$env$TN[keep]), Chla = log10(d$env$Chla[keep])); xy <- as.matrix(d$env[keep, c("lon", "lat")]); k <- 3L
idx <- sample_sets(nrow(x), k, 700, seed = 5); DM <- coords_dist(xy, lonlat = TRUE); R <- list()
r <- ugdm_response(x, idx, "relative"); dd <- ugdm_design(e, idx, dist_matrix = DM); ok <- r$keep
pdesign <- function(...) if (new) ugdm_position_design(...) else ugdm_interaction_design(..., spread = "constant")
pfit <- function(...) if (new) ugdm_fit_position(...) else ugdm_fit_interactions(...)
R$ugdm_fit <- unclass(ugdm_fit(dd, r$y[ok], r$trials[ok])[c("beta", "lambda", "intercept", "D2")])
for (pos in c("linear", "nonlinear")) { dp <- pdesign(dd, e, position = pos); f <- pfit(dp, r$y, r$trials, spread_smooth = .002, position_smooth = .002, ridge = 1e-4, spread_null = .001, position_null = .001); R[[paste0("pos_", pos)]] <- list(coef = f$coef, D2 = f$D2, pred = ugdm_predict(f, dp), conv = f$converged)
  g <- pfit(dp, r$y, r$trials, ridge = 1e-3, penalty_scale = "range", distance = FALSE); R[[paste0("pos_range_", pos)]] <- list(coef = g$coef, D2 = g$D2)
  dg <- if (new) ugdm_position_diagnostics(f, dp, r$y, r$trials) else ugdm_interaction_diagnostics(f, dp, r$y, r$trials); R[[paste0("diag_", pos)]] <- list(rank = dg$rank, cond = dg$condition, alloc = dg$deviance_allocation[c("spread", "position")], support = dg$support, resid = dg$nonlinear_position_residual)
  im <- ugdm_within_order_importance(f, dp, r$y, r$trials, permutations = 60, seed = 3, exact_max = 3); R[[paste0("imp_", pos)]] <- im$hierarchy }
pen <- data.frame(spread_smooth = c(.001, .01), position_smooth = c(.001, .01), ridge = c(.001, .01)); folds <- rep(1:4, length.out = nrow(e)); resp <- function(i) { z <- ugdm_response(x, i, "relative"); list(y = z$y, trials = z$trials, keep = z$keep) }
if (new) sel <- ugdm_select_penalties(e, folds, 3L, resp, pen, dist_matrix = DM, M = 120L, seed = 9) else sel <- ugdm_select_interactions(e, folds, 3L, resp, transform(pen, interaction_shrink = 0), spread = "constant", dist_matrix = DM, M = 120L, seed = 9)
R$select <- list(mean_loss = sel$mean_loss, selected = sel$selected, fold_loss = sel$fold_loss)
cr <- ugdm_component_response(x, idx)
for (m in c("additive", "spread_only")) for (a in c("coherent", "independent")) { cd <- ugdm_component_design(e, idx, m, coordinates = xy, lonlat = TRUE); cf <- ugdm_fit_components(cd, cr, a); pr <- ugdm_predict_components(cf, cd); R[[paste(m, a)]] <- list(coef = lapply(cf$fits, function(z) z$coef), pred = as.matrix(pr[, c("total", "occupancy", "heterogeneity")]), conv = cf$converged, diag = ugdm_component_diagnostics(cf, cd, cr)$stages)
  if (m == "additive") { ci <- ugdm_component_importance(cf, cd, cr, permutations = 40, exact_max = 3, seed = 4); R[[paste("imp", m, a)]] <- ci$hierarchy[, c("response", "predictor", "term", "absolute", "D2")] } }
fold <- rep(1:3, length.out = nrow(x)); vc <- if (new) ugdm_validate_components(x, e, xy, 2L, fold = fold, combinations = 120, seed = 6, lonlat = TRUE, geographic = TRUE) else ugdm_validate_components(x, e, xy, 2L, fold = fold, models = c("spread_only", "additive"), combinations = 120, seed = 6, lonlat = TRUE, geographic = TRUE)
R$validate <- vc[, c("fold", "model", "architecture", "D2", "RMSE", "closure_RMSE")]; saveRDS(R, out); cat("saved", out, "new =", new, "\n")
