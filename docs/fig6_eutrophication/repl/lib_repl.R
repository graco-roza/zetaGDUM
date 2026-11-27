# lib_repl.R: replacement-oriented components with the original/extended MS-GDUM machinery of zetaGDUM f7e115f
# (ugdm_response -> ugdm_design -> ugdm_interaction_design -> ugdm_fit_interactions -> ugdm_within_order_importance). No package code is modified.
source("lib6.R"); NBLK <- 5; COMPS <- c(turnover = "Turnover (incidence)", nestedness = "Nestedness-resultant (incidence)", balanced = "Balanced variation (abundance)", gradient = "Abundance gradient (abundance)")
MODELS <- data.frame(model = c("linear_const", "nonlinear_const", "nonlinear_varying"), position = c("linear", "nonlinear", "nonlinear"), spread = c("constant", "constant", "varying"), stringsAsFactors = FALSE)
PEN <- list(spread_smooth = .002, position_smooth = .002, interaction_shrink = .005, ridge = 1e-4)
hav <- function(xy) { r <- pi / 180; lo <- xy[, 1] * r; la <- xy[, 2] * r; n <- nrow(xy); m <- matrix(0, n, n); for (i in seq_len(n)) { dl <- lo - lo[i]; a <- sin((la - la[i]) / 2)^2 + cos(la[i]) * cos(la) * sin(dl / 2)^2; m[i, ] <- 2 * 6371 * asin(pmin(1, sqrt(a))) }; m }
resp <- function(x, il, comp) { r <- ugdm_response(x, il, comp); r }   # x: community of the partition, il: local set indices
build <- function(e, il, dmat, model, ref = NULL, basis_sites = NULL) { m <- MODELS[MODELS$model == model, ]
  d0 <- if (is.null(ref)) ugdm_design(e, il, dist_matrix = dmat) else ugdm_design(e, il, dist_matrix = dmat, ref = ref$base)
  list(base = d0, int = ugdm_interaction_design(d0, e, position = m$position, spread = m$spread)) }
fit1 <- function(dsg, r, pen = PEN) { k <- r$keep; ugdm_fit_interactions(dsg$int, r$y[k], r$trials[k], spread_smooth = pen$spread_smooth, position_smooth = pen$position_smooth, interaction_shrink = if (is.null(pen$interaction_shrink)) 0 else pen$interaction_shrink, ridge = pen$ridge) }
