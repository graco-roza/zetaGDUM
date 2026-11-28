# 00_sim_common.R: shared helpers for the spread-by-position interaction simulations (Simulation A: parameter recovery; Simulation B: guild mechanism).
# Package: zetaGDUM f7e115f, unmodified (ZG_LIB6). Framework: single-stage fractional-logit MS-GDUM (ugdm_design / ugdm_interaction_design / ugdm_fit_interactions), used for the Figure 5 models.
# Models (identical I-spline bases, logit link, nonnegativity constraints, penalties and sampled sets):
#   spread_only : position-independent spread curves only                  (spread = "constant", no position terms)
#   additive    : spread + nonlinear (cubic) position, position-independent spread    (spread = "constant", position = "nonlinear")
#   interaction : additive + explicit spread-by-position (endpoint curves interpolated by set position)  (spread = "varying", position = "nonlinear")
# Penalties (same in every model; the endpoint-difference shrinkage applies only to the interaction model): spread_smooth .002, position_smooth .002, interaction_shrink .005, ridge 1e-4 (as in the empirical Figure 5 fits).
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr); library(parallel)})
ns <- asNamespace("zetaGDUM"); lay_fn <- get(".interaction_layout", ns); pen_fn <- get(".interaction_penalties", ns); core_fn <- get(".interaction_fit_core", ns); wts_fn <- get(".model_weights", ns); bdev <- get(".binom_dev", ns)
PEN <- c(spread_smooth = .002, position_smooth = .002, interaction_shrink = .005, ridge = 1e-4); MODELS <- c("spread_only", "additive", "interaction")
mk_design <- function(env, idx, kind, ref = NULL) { d0 <- if (is.null(ref)) ugdm_design(env, idx) else ugdm_design(env, idx, ref = ref$base)
  list(base = d0, int = ugdm_interaction_design(d0, env, position = if (kind == "spread_only") "linear" else "nonlinear", spread = if (kind == "interaction") "varying" else "constant")) }
fit_model <- function(kind, dsg, y, trials, pen = PEN) { d <- dsg$int; a <- lay_fn(d, FALSE); w <- wts_fn(length(y), trials)
  if (kind == "spread_only") { keep <- c(1L, a$dir); P <- pen_fn(d, a, pen)[keep, keep]; f <- core_fn(a$X[, keep, drop = FALSE], y, w, P, a$dir); f$coef <- stats::setNames(f$coef, colnames(a$X)[keep]); f$kind <- kind; f$keep <- keep; f }
  else { f <- ugdm_fit_interactions(d, y, trials, distance = FALSE, spread_smooth = pen[["spread_smooth"]], position_smooth = pen[["position_smooth"]], interaction_shrink = pen[["interaction_shrink"]], ridge = pen[["ridge"]]); f$kind <- kind; f$keep <- NULL; f } }
eta_of <- function(f, dsg) { X <- lay_fn(dsg$int, FALSE)$X; if (!is.null(f$keep)) X <- X[, f$keep, drop = FALSE]; stopifnot(identical(colnames(X), names(f$coef))); drop(X %*% f$coef) }
d2_heldout <- function(y, p, trials, y_train, trials_train) { mu0 <- sum(y_train * trials_train) / sum(trials_train); 1 - bdev(y, p, trials) / bdev(y, rep(mu0, length(y)), trials) }
# synthetic matched profiles for conditional surfaces: sets of k sites whose focal predictor values are m + h * q (q fixed equally spaced with mean pairwise |difference| = 1), nuisance predictors at their means
prof_q <- function(k) { q <- seq(-1, 1, length.out = k); q / mean(dist(q)) }
surface_env <- function(env, focal, grid, k) { q <- prof_q(k); rows <- vector("list", nrow(grid)); for (i in seq_len(nrow(grid))) { e <- as.data.frame(lapply(env, function(v) rep(mean(v), k))); e[[focal]] <- grid$m[i] + grid$h[i] * q; rows[[i]] <- e }; E <- do.call(rbind, rows); idx <- matrix(seq_len(nrow(E)), nrow(grid), k, byrow = TRUE); list(env = E, idx = idx) }
set_stats <- function(v, idx) { k <- ncol(idx); pos <- rowMeans(matrix(v[idx], nrow(idx))); sp <- rep(0, nrow(idx)); for (a in 1:(k - 1)) for (b in (a + 1):k) sp <- sp + abs(v[idx[, a]] - v[idx[, b]]); data.frame(pos = pos, spread = sp / choose(k, 2)) }
support_cells <- function(v, idx, grid, k, wpos, wspr, minsets = 30) { s <- set_stats(v, idx); q <- prof_q(k); inside <- vapply(seq_len(nrow(grid)), function(i) { z <- grid$m[i] + grid$h[i] * q; all(z >= min(v) & z <= max(v)) }, TRUE)
  n <- vapply(seq_len(nrow(grid)), function(i) sum(abs(s$pos - grid$m[i]) <= wpos & abs(s$spread - grid$h[i]) <= wspr), 1L); data.frame(grid, n_sets = n, inside_range = inside, supported = inside & n >= minsets) }
# conditional contrasts from link-scale predictions on a (position x spread) grid; baseline = shortest spread at that position; response-scale via plogis
contrast_tab <- function(g, eta, m_lo, m_hi, h1, h2) { f <- function(m, h) eta[g$m == m & g$h == h]; d <- function(m) f(m, h2) - f(m, h1); r <- function(m) plogis(f(m, h2)) - plogis(f(m, h1)); data.frame(link_lo = d(m_lo), link_hi = d(m_hi), link_did = d(m_hi) - d(m_lo), resp_lo = r(m_lo), resp_hi = r(m_hi), resp_did = r(m_hi) - r(m_lo)) }
# explicit-interaction part of a conditional contrast = contrast of the model minus contrast of the same model with its two endpoint curves replaced by their average (i.e. forced position-independent spread; nonlinear I-spline sensitivity retained)
avg_endpoints <- function(cf) { for (v in unique(sub('_f[123]_[LU]$', '', grep('_f[123]_[LU]$', names(cf), value = TRUE)))) for (j in 1:3) { l <- paste0(v, '_f', j, '_L'); u <- paste0(v, '_f', j, '_U'); m <- (cf[[l]] + cf[[u]]) / 2; cf[[l]] <- m; cf[[u]] <- m }; cf }
# ---- generic evaluation of the competing models for one replicate / order (identical sets, bases, penalties) ----
CONFIGS <- data.frame(config = c("spread_only", "additive", "interaction_5e-03", "interaction_1e-04", "interaction_1e-06"), kind = c("spread_only", "additive", rep("interaction", 3)), shrink = c(NA, NA, 5e-3, 1e-4, 1e-6), stringsAsFactors = FALSE)
GRID_STD <- expand.grid(m = seq(-1.5, 1.5, by = .5), h = seq(.4, 2.4, by = .4)); M_LO <- -1; M_HI <- 1; WPOS <- .25; WSPR <- .25; MINSETS <- 30
# env: data.frame with focal predictor "P"; mu_P, sd_P: population mean/sd used to place the fixed grid. truth: list(eta = link-scale truth on GRID (raw units), resp = response-scale truth on GRID)
evaluate_order <- function(env, tr, te, k, il_tr, il_te, y_tr, y_te, w_tr, w_te, truth, mu_P, sd_P, keep_surface = FALSE) {
  grid <- data.frame(m = mu_P + sd_P * GRID_STD$m, h = sd_P * GRID_STD$h, m_std = GRID_STD$m, h_std = GRID_STD$h); sup <- support_cells(env$P[tr], il_tr, grid, k, WPOS * sd_P, WSPR * sd_P, MINSETS)
  # matched contrast: spreads supported at both positions (response-independent)
  common <- sort(intersect(sup$h_std[sup$supported & sup$m_std == M_LO], sup$h_std[sup$supported & sup$m_std == M_HI])); ok <- length(common) >= 2 && diff(range(common)) >= .8; h1 <- if (ok) min(common) else NA; h2 <- if (ok) max(common) else NA
  tl <- truth$eta; tr_resp <- truth$resp; tl_avg <- truth$eta_avg; out <- list(); surf <- list()
  tc <- if (ok) contrast_tab(transform(grid, hh = h_std), tl, M_LO, M_HI, h1, h2) else NULL
  gx <- grid; gx$m <- gx$m_std; gx$h <- gx$h_std; if (ok) tcx <- contrast_tab(transform(gx), tl, M_LO, M_HI, h1, h2) else tcx <- NULL; tcx_avg <- if (ok && !is.null(tl_avg)) contrast_tab(gx, tl_avg, M_LO, M_HI, h1, h2) else NULL
  syn <- surface_env(env[tr, , drop = FALSE], "P", grid, k)
  for (ci in seq_len(nrow(CONFIGS))) { cf <- CONFIGS[ci, ]; pen <- PEN; if (!is.na(cf$shrink)) pen[["interaction_shrink"]] <- cf$shrink
    r <- tryCatch({ dtr <- mk_design(env[tr, , drop = FALSE], il_tr, cf$kind); dte <- mk_design(env[te, , drop = FALSE], il_te, cf$kind, ref = dtr); f <- fit_model(cf$kind, dtr, y_tr, w_tr, pen); p <- plogis(eta_of(f, dte))
      dsyn <- mk_design(syn$env, syn$idx, cf$kind, ref = dtr); es <- eta_of(f, dsyn); s <- sup$supported
      cm <- if (ok) contrast_tab(gx, es, M_LO, M_HI, h1, h2) else NULL; ea <- if (cf$kind == 'interaction') { Xs <- lay_fn(dsyn$int, FALSE)$X; drop(Xs %*% avg_endpoints(f$coef)) } else es; cma <- if (ok) contrast_tab(gx, ea, M_LO, M_HI, h1, h2) else NULL
      surf_rows <- if (keep_surface) data.frame(config = cf$config, order = k, m_std = grid$m_std, h_std = grid$h_std, supported = s, eta_est = es, resp_est = plogis(es), eta_true = tl, resp_true = tr_resp) else NULL
      list(row = data.frame(config = cf$config, model = cf$kind, order = k, D2_test = d2_heldout(y_te, p, w_te, y_tr, w_tr), D2_train = f$D2, converged = isTRUE(f$converged), n_supported_cells = sum(s), matched_support = ok, h1 = h1, h2 = h2,
        rmse_link = if (any(s)) sqrt(mean((es[s] - tl[s])^2)) else NA, rmse_resp = if (any(s)) sqrt(mean((plogis(es[s]) - tr_resp[s])^2)) else NA, cor_link = if (sum(s) > 2) cor(es[s], tl[s]) else NA, cor_resp = if (sum(s) > 2) cor(plogis(es[s]), tr_resp[s]) else NA,
        est_link_did = if (ok) cm$link_did else NA, est_resp_did = if (ok) cm$resp_did else NA, true_link_did = if (ok) tcx$link_did else NA, true_resp_did = if (ok) tcx$resp_did else NA,
        est_link_did_explicit = if (ok) cm$link_did - cma$link_did else NA, est_resp_did_explicit = if (ok) cm$resp_did - cma$resp_did else NA, true_link_did_explicit = if (ok && !is.null(tcx_avg)) tcx$link_did - tcx_avg$link_did else NA, true_resp_did_explicit = if (ok && !is.null(tcx_avg)) tcx$resp_did - tcx_avg$resp_did else NA, est_link_lo = if (ok) cm$link_lo else NA, est_link_hi = if (ok) cm$link_hi else NA, true_link_lo = if (ok) tcx$link_lo else NA, true_link_hi = if (ok) tcx$link_hi else NA, est_resp_lo = if (ok) cm$resp_lo else NA, est_resp_hi = if (ok) cm$resp_hi else NA, true_resp_lo = if (ok) tcx$resp_lo else NA, true_resp_hi = if (ok) tcx$resp_hi else NA), surf = surf_rows, coef = f$coef) },
      error = function(e) { message("ERR ", cf$config, " k", k, ": ", conditionMessage(e)); NULL }); if (!is.null(r)) { out[[length(out) + 1]] <- r$row; if (!is.null(r$surf)) surf[[length(surf) + 1]] <- r$surf } }
  list(rows = do.call(rbind, out), surf = if (length(surf)) do.call(rbind, surf) else NULL) }
