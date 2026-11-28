# 00_sim_common.R: shared generators and evaluation for the sampling-requirements simulation (principal additive package; NO explicit interactions).
# Run from final/ with ZG_LIB6 = library built from branch principal/additive-msgdum.
# Environment is simulated directly on the log10 scale used by the empirical Figure 6, so the support rules (positions 10/35/100 ug/L, >= 10 common profiles, >= 3 grid points,
# >= 0.20 log10 spread width, >= 2 spread bins, per-predictor cell support and nearest-neighbour support) are applied UNCHANGED (revision_sampling/01_support.R, copied from the Figure 6 machinery).
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(parallel); library(RANN)})
source("revision_sampling/01_support.R"); source("revision_sampling/02_contrasts.R")
MAXPROF <- 40L                       # reference configurations used for predictions and Monte Carlo truth (support is judged on all reference configurations)
LO <- 0.4; HI <- 2.8                 # log10 TP limits (2.5 - 631 ug/L)
ENVS <- c("uniform", "normal", "clustered"); SAMPS <- c("random", "balanced"); KS <- c(2L, 3L, 5L, 10L, 20L)
# ---------------- environment ----------------
# TP (focal): uniform U(0.5, 2.5); normal N(1.5, 0.4) truncated to [0.4, 2.8]; clustered = three equally weighted clusters (0.9, 1.55, 2.2; sd 0.10).
# TN and Chl-a are correlated nuisance predictors (log10): TN = 2.5 + 0.45 (TP - 1.5) + N(0, .25); Chla = 0.3 + 0.60 (TP - 1.5) + N(0, .30). Coordinates uniform in a 20 x 15 degree box; the true effect of geographic distance is zero.
gen_env <- function(n, dist) {
  tp <- switch(dist, uniform = runif(n, .5, 2.5), normal = { x <- rnorm(n, 1.5, .4); while (any(b <- x < LO | x > HI)) x[b] <- rnorm(sum(b), 1.5, .4); x },
    clustered = c(.9, 1.55, 2.2)[sample.int(3, n, TRUE)] + rnorm(n, 0, .1))
  list(E = data.frame(TP = tp, TN = 2.5 + .45 * (tp - 1.5) + rnorm(n, 0, .25), Chla = .3 + .6 * (tp - 1.5) + rnorm(n, 0, .3)), xy = cbind(lon = runif(n, -100, -80), lat = runif(n, 30, 45))) }
# ---------------- combination sampling ----------------
rand_sets <- function(n, k, C) t(vapply(seq_len(C), function(i) sample.int(n, k), integer(k)))
pair_spread <- function(x) { k <- length(x); if (k == 2) return(abs(diff(x))); mean(dist(x)) }
shape_q <- function(k, type) { p <- (seq_len(k) - .5) / k
  v <- switch(type, even = seq(-1, 1, length.out = k), normal = qnorm(p), low = qbeta(p, 1, 2.5), up = qbeta(p, 2.5, 1), ends = c(rep(-1, ceiling(k / 2)), rep(1, floor(k / 2))) + seq(-.05, .05, length.out = k))
  v <- v - mean(v); v / pair_spread(v) }
SHAPES <- c("even", "normal", "low", "up", "ends")
# environmentally balanced combinations: target cells of (mean TP, TP spread) on a regular grid (mean 0.8-2.3 by 0.1; spread 0.05-1.30 by 0.05); a target profile m + s q (random arrangement q) is filled with a randomly chosen unused observed site within +-0.03 of each target value (nearest unused site if none); a set is accepted
# only when the achieved mean and spread are within 0.05 of the targets. Infeasible cells are skipped, not relaxed; the achieved number of sets is reported. Only sites already sampled are used (no additional field sampling).
bal_sets <- function(tp, k, C, tol = .05, w = .03, max_try = 15L * C) {
  n <- length(tp); ord <- order(tp); sv <- tp[ord]; out <- matrix(NA_integer_, C, k); got <- 0L; tries <- 0L; rg <- range(tp); mg <- seq(.8, 2.3, by = .1); sg <- seq(.05, 1.3, by = .05)
  while (got < C && tries < max_try) { tries <- tries + 1L; m <- sample(mg, 1); s <- sample(sg, 1); tx <- m + s * shape_q(k, sample(SHAPES, 1)); if (min(tx) < rg[1] || max(tx) > rg[2]) next
    used <- logical(n); sel <- integer(k); okk <- TRUE
    for (i in seq_len(k)) { a <- findInterval(tx[i] - w, sv) + 1L; b <- findInterval(tx[i] + w, sv); cand <- if (b >= a) a:b else integer(); cand <- cand[!used[ord[cand]]]
      if (length(cand)) j <- cand[sample.int(length(cand), 1L)] else { p <- findInterval(tx[i], sv); cand <- max(1L, p - k):min(n, p + k + 1L); cand <- cand[!used[ord[cand]]]; if (!length(cand)) { okk <- FALSE; break }; j <- cand[which.min(abs(sv[cand] - tx[i]))] }
      sel[i] <- ord[j]; used[ord[j]] <- TRUE }
    if (!okk) next; xs <- tp[sel]; if (abs(mean(xs) - m) > tol || abs(pair_spread(xs) - s) > tol) next
    got <- got + 1L; out[got, ] <- sel }
  list(idx = out[seq_len(got), , drop = FALSE], tries = tries) }
# ---------------- site-level metacommunity ----------------
# 200 taxa: 120 specialists (Gaussian niche on rescaled TP u = (TP - 0.4)/2.4, SD 0.08, optima evenly spaced -0.05..1.05) and 80 generalists (SD 0.35, optima 0..1); abundance factors lognormal(0, .5) rescaled within guild.
# Generalist dominance increases with TP (0.1 + 0.8 plogis(10 (u - .5))). Counts are generated per SITE (multinomial, N = 2000), never per sampled combination. TN, Chl-a and coordinates do not affect the community.
N_TOT <- 2000L
set.seed(4242); SP <- data.frame(guild = rep(c("specialist", "generalist"), c(120, 80)), opt = c(seq(-.05, 1.05, length.out = 120), seq(0, 1, length.out = 80)), sd = rep(c(.08, .35), c(120, 80)), a = exp(rnorm(200, 0, .5)))
for (g in c("specialist", "generalist")) { i <- SP$guild == g; SP$a[i] <- SP$a[i] / sum(SP$a[i] * exp(-.5 * ((.5 - SP$opt[i]) / SP$sd[i])^2)) }
meta_p <- function(tp) { u <- (tp - LO) / (HI - LO); N <- exp(-.5 * outer(u, SP$opt, "-")^2 / rep(SP$sd, each = length(u))^2) * rep(SP$a, each = length(u)); pg <- .1 + .8 * plogis(10 * (u - .5)); isg <- SP$guild == "generalist"; W <- N; W[, isg] <- W[, isg] * pg; W[, !isg] <- W[, !isg] * (1 - pg); W / rowSums(W) }
sim_counts <- function(tp) { p <- meta_p(tp); cnt <- t(vapply(seq_len(nrow(p)), function(i) rmultinom(1, N_TOT, p[i, ])[, 1], numeric(200))); colnames(cnt) <- sprintf("sp%03d", 1:200); cnt }
# independent Monte Carlo truth for a matched cell (position m, spread s) over reference configurations ids: fresh communities, R draws per reference configuration, package response functions
truth_cell <- function(p, ids, m, s, R = 25L) { k <- p$k; tp <- as.vector(t(m + s * p$shape[ids, , drop = FALSE])); P <- meta_p(tp); nr <- length(ids)
  cnt <- t(vapply(rep(seq_along(tp), each = R), function(i) rmultinom(1, N_TOT, P[i, ])[, 1], numeric(200))); colnames(cnt) <- sprintf("sp%03d", 1:200)   # rows: (site, draw)
  idx <- matrix(0L, nr * R, k); for (a in seq_len(nr)) for (r in seq_len(R)) idx[(a - 1) * R + r, ] <- ((a - 1) * k + seq_len(k) - 1L) * R + r
  y <- unclass(ugdm_component_response(cnt, idx)); data.frame(total = mean(y[, "total"]), occupancy = mean(y[, "occupancy"]), heterogeneity = mean(y[, "heterogeneity"])) }
# ---------------- helpers shared by both simulations ----------------
tab_support <- function(p) { iv <- names(p$valid); nv <- length(iv); best <- best_interval(iv); mx <- max(colSums(p$ok)); ij <- if (!is.na(best)) parse_iv(best) else c(NA, NA)
  # reason for unsupported: widest run of consecutive grid points with >= MINREF common profiles (before bin rule)
  good <- colSums(p$ok) >= MINREF; run <- 0L; mxrun <- 0L; for (g in good) { run <- if (g) run + 1L else 0L; mxrun <- max(mxrun, run) }
  data.frame(supported = !is.na(best), n_valid_intervals = nv, interval = best, s_from = if (!is.na(best)) GRID[ij[1]] else NA_real_, s_to = if (!is.na(best)) GRID[ij[2]] else NA_real_, interval_width = if (!is.na(best)) GRID[ij[2]] - GRID[ij[1]] else NA_real_, n_ref_supported = if (!is.na(best)) length(p$valid[[best]]) else 0L,
    max_common_profiles = mx, longest_run_grid_points = mxrun, n_ref = if (is.null(p$ref)) 0L else nrow(p$ref), nn_threshold = p$threshold) }
# contrasts from a set of predicted cells G (names "model position spread") for a model; returns A at three positions and B Low->High at mid spread
contrasts_for <- function(G, model, svals) { s_lo <- svals[1]; s_hi <- svals[length(svals)]; s_mid <- svals[ceiling(length(svals) / 2)]
  A <- do.call(rbind, lapply(names(PCOL), function(pn) { ca <- contrast_A(G, model, pn, s_lo, s_hi); data.frame(contrast = paste0("A_", pn), model = model, dT = ca$dT_pp, dO = ca$dO_pp, dH = ca$dH_pp, dEtaT = ca$dEtaT, dEtaS = ca$dEtaS) }))
  cb <- contrast_B(G, model, "Low", "High", s_mid); B <- data.frame(contrast = "B_LowHigh", model = model, dT = cb$dT_pp, dO = cb$dO_pp, dH = cb$dH_pp, dEtaT = cb$dEtaT, dEtaS = cb$dEtaS); rbind(A, B) }
truth_contrasts <- function(p, ids, svals, R = 25L) { s_lo <- svals[1]; s_hi <- svals[length(svals)]; s_mid <- svals[ceiling(length(svals) / 2)]; need <- unique(rbind(expand.grid(pn = names(PCOL), s = c(s_lo, s_hi), stringsAsFactors = FALSE), data.frame(pn = names(PCOL), s = s_mid)))
  V <- lapply(seq_len(nrow(need)), function(i) truth_cell(p, ids, log10(PCOL[[need$pn[i]]]), need$s[i], R)); names(V) <- paste(need$pn, need$s); g <- function(pn, s) V[[paste(pn, s)]]
  A <- do.call(rbind, lapply(names(PCOL), function(pn) { a <- g(pn, s_lo); b <- g(pn, s_hi); data.frame(contrast = paste0("A_", pn), dT = 100 * (b$total - a$total), dO = 100 * (b$occupancy - a$occupancy), dH = 100 * (b$heterogeneity - a$heterogeneity)) }))
  a <- g("Low", s_mid); b <- g("High", s_mid); rbind(A, data.frame(contrast = "B_LowHigh", dT = 100 * (b$total - a$total), dO = 100 * (b$occupancy - a$occupancy), dH = 100 * (b$heterogeneity - a$heterogeneity))) }
predict_cells <- function(p, dd, fit, model, ids, svals) { bT <- setNames(bvec(fit, "total"), colnames(layoutfn(dd, TRUE)$X)); bS <- setNames(bvec(fit, "share"), names(bT)); G <- list(); for (pn in names(PCOL)) for (sv in svals) G[[paste(model, pn, sv)]] <- cell_parts(p, dd, model, ids, log10(PCOL[[pn]]), sv, bT, bS); G }
dev_frac <- function(y, p, ybar) { f <- function(y, p) { p <- pmin(pmax(p, 1e-9), 1 - 1e-9); sum(ifelse(y > 0, y * log(y / p), 0) + ifelse(y < 1, (1 - y) * log((1 - y) / (1 - p)), 0)) }; 1 - f(y, p) / f(y, rep(ybar, length(y))) }
heldout <- function(fit, dtest, yte, ytr_means) { pr <- ugdm_predict_components(fit, dtest); y <- unclass(yte); c(D2_total = dev_frac(y[, "total"], pr$total, ytr_means[["total"]]), rmse_total_pp = 100 * sqrt(mean((pr$total - y[, "total"])^2)), mean_pred_total = mean(pr$total), mean_obs_total = mean(y[, "total"]),
  rmse_occupancy_pp = 100 * sqrt(mean((pr$occupancy - y[, "occupancy"])^2)), rmse_heterogeneity_pp = 100 * sqrt(mean((pr$heterogeneity - y[, "heterogeneity"])^2))) }
