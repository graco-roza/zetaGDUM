# 00_ident_common.R: shared helpers for the decisive identifiability / interaction validation (all in simulation scripts; the zetaGDUM public API is not modified; package f7e115f, internal functions used read-only).
# Custom design builder: columns of the extended MS-GDUM layout are assembled from the package's own pieces (ugdm_design -> ugdm_interaction_design: base spread columns, nonlinear position columns, set-position coordinates u),
# so any interaction specification (all predictors = current package model, focal same-predictor, focal cross-predictor with the productivity position coordinate) can be built with the SAME endpoint-interpolation architecture,
# and penalties are built explicitly so that the equal-endpoint restriction of an interaction model has EXACTLY the additive penalty (the package's ridge doubles on spread coefficients in the varying parameterisation).
source("revision_sim/00_sim_common.R")
# spec: list(type = "spread_only"|"additive"|"interaction", focal = predictor name or NULL (NULL = all predictors, i.e. the package 'varying' layout), from = predictor whose set-position coordinate u interpolates the endpoints (default = focal))
build_X <- function(env, idx, spec, ref = NULL) { d0 <- if (is.null(ref)) ugdm_design(env, idx) else ugdm_design(env, idx, ref = ref$d0); di <- ugdm_interaction_design(d0, env, position = "nonlinear", spread = "constant"); vars <- d0$vars; base <- di$base; u <- di$u; colnames(u) <- vars
  cols <- list(intercept = rep(1, nrow(base))); nn <- character(); groups <- list()
  for (j in seq_along(vars)) { v <- vars[j]; B <- base[, 3 * (j - 1) + 1:3, drop = FALSE]; inter <- spec$type == "interaction" && (is.null(spec$focal) || v == spec$focal)
    if (inter) { src <- if (is.null(spec$from)) v else spec$from; uu <- u[, src]; for (i in 1:3) cols[[paste0(v, "_f", i, "_L")]] <- B[, i] * (1 - uu); for (i in 1:3) cols[[paste0(v, "_f", i, "_U")]] <- B[, i] * uu; groups[[v]] <- "inter" } else { for (i in 1:3) cols[[paste0(v, "_f", i)]] <- B[, i]; groups[[v]] <- "const" } }
  dirnames <- setdiff(names(cols), "intercept")
  if (spec$type != "spread_only") for (j in seq_along(vars)) { P <- di$pos[, 3 * (j - 1) + 1:3, drop = FALSE]; for (i in 1:3) cols[[colnames(P)[i]]] <- P[, i] }
  X <- do.call(cbind, cols); list(X = X, dir = which(colnames(X) %in% dirnames), vars = vars, groups = groups, di = di, d0 = d0, spec = spec, u = u) }
pen_matrix <- function(bx, s, matched = TRUE) { X <- bx$X; q <- ncol(X); P <- diag(c(0, rep(s$ridge, q - 1)), q); colnames(P) <- rownames(P) <- colnames(X); di <- bx$di
  for (j in seq_along(bx$vars)) { v <- bx$vars[j]; C <- di$curvature[[j]]
    if (bx$groups[[v]] == "const") { at <- match(paste0(v, "_f", 1:3), colnames(X)); P[at, at] <- P[at, at] + s$spread_smooth * C }
    else { at <- match(c(paste0(v, "_f", 1:3, "_L"), paste0(v, "_f", 1:3, "_U")), colnames(X)); I3 <- diag(3); A <- matrix(0, 6, 6); A[1:3, 1:3] <- A[4:6, 4:6] <- s$spread_smooth * C / 2; dl <- cbind(-I3, I3)
      if (matched) P[at, at] <- P[at, at] - diag(s$ridge / 2, 6)   # equal endpoints => total ridge = ridge * |b|^2, identical to the additive model
      P[at, at] <- P[at, at] + A + s$interaction_shrink * crossprod(dl, dl) } }
  if (any(grepl("_quadratic$", colnames(X)))) for (j in seq_along(bx$vars)) { v <- bx$vars[j]; at <- match(c(v, paste0(v, "_quadratic"), paste0(v, "_cubic")), colnames(X)); P[at, at] <- P[at, at] + s$position_smooth * di$position_curvature[j] * diag(c(0, 8, 24)) }
  P }
STAGES <- list(unpenalised = c(ridge = 0, spread_smooth = 0, position_smooth = 0, interaction_shrink = 0), ridge = c(ridge = 1e-4, spread_smooth = 0, position_smooth = 0, interaction_shrink = 0), curvature = c(ridge = 1e-4, spread_smooth = .002, position_smooth = .002, interaction_shrink = 0),
  shrink = c(ridge = 1e-4, spread_smooth = .002, position_smooth = .002, interaction_shrink = .005), weakshrink = c(ridge = 1e-4, spread_smooth = .002, position_smooth = .002, interaction_shrink = 1e-6))
# package optimiser (the package's own L-BFGS-B working-loss solver, called through its internal entry point) and an INDEPENDENT constrained optimiser (nlminb/PORT, unscaled parameterisation, analytic gradient) of the identical objective
objective_fun <- function(X, y, w, P) { w <- w / sum(w); function(b) { eta <- drop(X %*% b); sum(w * (log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta)) + drop(crossprod(b, P %*% b)) / 2 } }
fit_pkg <- function(bx, y, w, P) { f <- core_fn(bx$X, y, w, P, bx$dir); f$coef <- setNames(f$coef, colnames(bx$X)); f$obj <- objective_fun(bx$X, y, w, P)(f$coef); f }
fit_indep <- function(bx, y, w, P, start = NULL) { X <- bx$X; wn <- w / sum(w); q <- ncol(X); lower <- rep(-Inf, q); lower[bx$dir] <- 0; mu0 <- min(max(sum(wn * y), 1e-6), 1 - 1e-6); b0 <- c(qlogis(mu0), rep(0, q - 1)); if (!is.null(start)) b0 <- pmax(start, lower)
  fn <- objective_fun(X, y, w, P); gr <- function(b) { drop(crossprod(X, wn * (plogis(drop(X %*% b)) - y)) + P %*% b) }
  best <- NULL; for (s0 in list(b0, c(b0[1], rep(0.05, q - 1)))) { s0 <- pmax(s0, lower); o <- nlminb(s0, fn, gr, lower = lower, control = list(iter.max = 3000, eval.max = 6000, rel.tol = 1e-13, x.tol = 1e-12)); if (is.null(best) || o$objective < best$objective) best <- o }
  list(coef = setNames(best$par, colnames(X)), obj = best$objective, converged = best$convergence == 0) }
pred_eta <- function(bx, coef) drop(bx$X %*% coef[colnames(bx$X)])
d2_frac <- function(y, p, w, mu0) 1 - bdev(y, p, w) / bdev(y, rep(mu0, length(y)), w)
