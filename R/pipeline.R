#' Multi-site GDUM across zeta orders
#'
#' Convenience pipeline: computes the zeta retention curve and its breakpoint, chooses orders
#' (2, the midpoint and the breakpoint, or the largest order not above it at which at least
#' `min_share` of the sets share a species), builds sets, fits the multi-site GDUM at each order
#' and, with paired Bayesian-bootstrap draws (the same site weights at every order), returns
#' importance shares (and optionally Shapley contributions to D2 for groups of predictors).
#'
#' @param comm Sites x species matrix (counts or presence-absence).
#' @param env Sites x predictors data frame.
#' @param coords Optional two-column coordinates; adds distance as a directional predictor.
#' @param lonlat Coordinates are decimal degrees (haversine distance).
#' @param component Response: `"turnover"`, `"total"`, `"nestedness"`, or (abundance data)
#'   `"relative"` (relative-abundance dissimilarity), `"balanced"`,
#'   `"bray_curtis"`, `"gradient"` (raw-count responses).
#' @param orders Integer vector of orders; default `c(2, midpoint, break)`.
#' @param M Number of random sets per order (all pairs at order 2 when feasible).
#' @param ndraw Number of paired bootstrap draws.
#' @param groups Optional named list of predictor groups for the Shapley decomposition.
#' @param min_share Legacy minimum estimated sharing fraction used only in
#'   automatic order selection. To assess numerical precision separately from
#'   the ecological breakpoint, use [zeta_order_precision()] with explicit orders.
#' @param seed Random seed.
#' @param cores Number of cores for the draws.
#' @param paired If `TRUE` (default) the same site weights are used at every order within a draw; if `FALSE` each order uses independent draws (for comparison only).
#' @param free_sign Allow directional effects of either sign (see [ugdm_fit()]). The default (`NULL`) constrains them to be non-negative, as in the pairwise GDUM; setting `TRUE` is advisable for the nestedness and abundance-gradient components, for which there is no expectation about the direction of richness or abundance gradients along environmental gradients.
#' @param weighting `"direct"` (default) samples, for every draw, `M` sets of sites with probability proportional to the product of the site weights ([sample_sets_weighted()]), so that all sets have the same weight; `"product"` samples sets uniformly once and weights them by that product ([set_weights()]), which concentrates the weight on a few sets at high orders (the effective number of sets falls to a small fraction of `M`) and inflates the deviance explained. At order 2 with all pairs available, the pairs are always weighted.
#' @param sampling `"conventional"` (default) or `"species"` for the mixture
#'   in [sample_sets_species()]. Exhaustive pairs are always included when feasible.
#' @param alpha Targeted mixture probability for species-informed sampling.
#'   Supply one number, a numeric vector named by order (for example `k2`,
#'   `k5`), or `"adaptive"` to select the fraction by an independent pilot
#'   for each order. Adaptive selection requires `sampling = "species"`.
#' @param pilot_M,pilot_nrep Pilot budget and replicate count per fraction
#'   when `alpha = "adaptive"`. See [ugdm_sampling_pilot()].
#' @param working_weights Response weighting mode from [ugdm_response()].
#' @return An object of class `zeta_orders` with point fits, paired draws and
#'   sampling diagnostics. The diagnostics count retained draws and retained
#'   draws with shared species, and report Kish effective sample size based
#'   on the combined sampling and bootstrap weights before working weights.
#'   The selected fractions and any independent pilot results are included.
#' @export
#' @examples
#' sim <- simulate_metacommunity(60, 60, seed = 1)
#' res <- zeta_orders(sim$comm, sim$env, sim$coords, M = 2000, ndraw = 10,
#'                    groups = list(env1 = "env1", env2 = "env2", distance = "distance"))
#' res
#' compare_orders(res, "importance")
zeta_orders <- function(comm, env, coords = NULL, lonlat = FALSE, component = "turnover", orders = NULL, M = 1e4, ndraw = 100,
                        groups = NULL, min_share = 0.04, seed = 1, cores = 1, paired = TRUE, free_sign = NULL, weighting = c("direct", "product"),
                        sampling = c("conventional", "species"), alpha = 0.5, working_weights = c("default", "constant"),
                        pilot_M = 1000L, pilot_nrep = 8L) {
  weighting <- match.arg(weighting)
  sampling <- match.arg(sampling)
  working_weights <- match.arg(working_weights)
  adaptive <- identical(alpha, "adaptive")
  if (adaptive && sampling != "species") stop('Adaptive alpha requires `sampling = "species"`.')
  if (!adaptive) {
    if (!is.numeric(alpha) || !length(alpha)) stop("Supply numeric mixture fractions or `alpha = 'adaptive'`.")
    for (fraction in alpha) .check_alpha(fraction)
  }
  comm <- .as_abundance(comm)
  fs <- isTRUE(free_sign)
  n <- nrow(comm)
  if (n < 2L) stop("At least two sites are required.")
  if (nrow(env) != n) stop("`comm` and `env` must have the same number of sites.")
  cur <- zeta_curve(comm, kmax = max(2L, min(n - 1L, 80L)), type = "presence")
  kb <- zeta_break(cur)
  P <- (comm > 0) * 1L
  storage.mode(P) <- "integer"
  feasible <- function(k) {
    idx <- sample_sets(n, k, M = 4000, seed = seed)
    mean(shared_species_cpp(P, idx) > 0) >= min_share
  }
  if (is.null(orders)) {
    kbr <- if (is.na(kb)) 2L else kb
    while (kbr > 3 && !feasible(kbr)) kbr <- kbr - 1
    orders <- sort(unique(c(2L, as.integer(round((2 + kbr) / 2)), as.integer(kbr))))
  }
  dm <- if (!is.null(coords)) coords_dist(coords, lonlat) else NULL
  dr <- if (!is.null(dm)) .pair_distance_range(dm) else NULL
  keys <- paste0("k", orders)
  pilots <- stats::setNames(vector("list", length(orders)), keys)
  if (adaptive) {
    for (j in seq_along(orders)) {
      pilots[[j]] <- ugdm_sampling_pilot(comm, env, orders[j], M = pilot_M, nrep = pilot_nrep,
        component = component, working_weights = working_weights, free_sign = fs,
        seed = seed + 100000L * j, coords = coords, lonlat = lonlat)
    }
    alpha_by_order <- vapply(pilots, function(p) p$alpha, numeric(1L))
  } else if (length(alpha) == 1L) {
    alpha_by_order <- stats::setNames(rep(alpha, length(keys)), keys)
  } else {
    if (anyDuplicated(names(alpha)) || !all(keys %in% names(alpha))) {
      stop("Supply a named alpha value for every selected order.")
    }
    alpha_by_order <- alpha[keys]
  }
  if (sampling == "conventional") alpha_by_order[] <- 0
  W <- site_weights(n, ndraw, seed)
  terms <- c(colnames(env), if (!is.null(dm)) "distance")
  imp <- array(NA_real_, c(ndraw, length(orders), length(terms)), dimnames = list(NULL, paste0("k", orders), terms))
  shp <- if (!is.null(groups)) array(NA_real_, c(ndraw, length(orders), length(groups)), dimnames = list(NULL, paste0("k", orders), names(groups))) else NULL
  d2 <- matrix(NA_real_, ndraw, length(orders), dimnames = list(NULL, paste0("k", orders)))
  hgt <- array(NA_real_, c(ndraw, length(orders), length(terms)), dimnames = list(NULL, paste0("k", orders), terms))
  lam <- array(NA_real_, c(ndraw, length(orders), ncol(env)), dimnames = list(NULL, paste0("k", orders), colnames(env)))
  point <- list()
  metrics <- c("n_sets", "n_shared", "ess", "alpha")
  diagnostic <- array(NA_real_, c(ndraw, length(orders), length(metrics)),
                      dimnames = list(NULL, paste0("k", orders), metrics))
  point_diagnostic <- matrix(NA_real_, length(orders), length(metrics),
                             dimnames = list(paste0("k", orders), metrics))
  for (oi in seq_along(orders)) {
    k <- orders[oi]
    fraction <- alpha_by_order[[oi]]
    Wo <- if (paired) W else site_weights(n, ndraw, seed + 1000 * oi)
    direct <- weighting == "direct" && !(k == 2 && choose(n, 2) <= M)
    sets0 <- .model_sets(P, k, M, sampling, fraction, seed = seed + k)
    r0 <- ugdm_response(comm, sets0$idx, component, working_weights = working_weights)
    if (!any(r0$keep)) stop("No sets have a defined response at order ", k, ".")
    idx0 <- sets0$idx[r0$keep, , drop = FALSE]
    sw0 <- sets0$sampling_weights[r0$keep]
    y0 <- r0$y[r0$keep]
    tr0 <- r0$trials[r0$keep]
    des0 <- ugdm_design(env, idx0, dm, dist_range = dr)
    point[[paste0("k", k)]] <- ugdm_fit(des0, y0, tr0, free_sign = fs, sampling_weights = sw0)
    point_diagnostic[oi, ] <- .sampling_diagnostics(sets0, r0$keep, sw0)
    one <- function(b) {
      if (direct) {
        sets <- .model_sets(P, k, M, sampling, fraction, Wo[, b])
        rb <- ugdm_response(comm, sets$idx, component, working_weights = working_weights)
        if (!any(rb$keep)) return(NULL)
        sw <- sets$sampling_weights[rb$keep]
        y <- rb$y[rb$keep]
        tr <- rb$trials[rb$keep]
        des <- suppressWarnings(ugdm_design(env, sets$idx[rb$keep, , drop = FALSE], dm, ref = des0))
        diag <- .sampling_diagnostics(sets, rb$keep, sw)
      } else {
        y <- y0
        tr <- tr0
        des <- des0
        sw <- .product_sampling_weights(Wo[, b], idx0, sw0)
        diag <- .sampling_diagnostics(sets0, r0$keep, sw)
      }
      f <- tryCatch(ugdm_fit(des, y, tr, free_sign = fs, sampling_weights = sw), error = function(e) NULL)
      if (is.null(f)) return(NULL)
      list(imp = tryCatch(ugdm_importance(f, des, trials = tr, sampling_weights = sw), error = function(e) NULL),
           d2 = f$D2, hgt = c(f$height, if (f$has_distance) c(distance = f$dist_height)), lam = f$lambda,
           diagnostics = diag,
           shp = if (!is.null(groups)) tryCatch(
             ugdm_shapley(des, y, tr, groups, free_sign = fs, sampling_weights = sw),
             error = function(e) NULL))
    }
    res <- if (cores > 1) parallel::mclapply(seq_len(ndraw), one, mc.cores = cores) else lapply(seq_len(ndraw), one)
    for (b in seq_len(ndraw)) {
      if (is.null(res[[b]]) || inherits(res[[b]], "try-error")) next
      if (!is.null(res[[b]]$imp)) imp[b, oi, ] <- res[[b]]$imp[terms]
      diagnostic[b, oi, ] <- res[[b]]$diagnostics
      d2[b, oi] <- res[[b]]$d2
      hgt[b, oi, ] <- res[[b]]$hgt[terms]
      lam[b, oi, ] <- res[[b]]$lam[colnames(env)]
      if (!is.null(shp) && !is.null(res[[b]]$shp)) shp[b, oi, ] <- res[[b]]$shp
    }
  }
  structure(list(curve = cur, break_order = kb, orders = orders, component = component, point = point, importance = imp,
                 shapley = shp, D2 = d2, height = hgt, lambda = lam, n = n, sampling = sampling,
                 point_diagnostics = point_diagnostic, draw_diagnostics = diagnostic,
                 alpha_by_order = alpha_by_order, pilots = if (adaptive) pilots else NULL), class = "zeta_orders")
}

#' @export
print.zeta_orders <- function(x, ...) {
  cat("zeta_orders: ", x$n, " sites; retention break at order ", x$break_order, "; orders used: ", paste(x$orders, collapse = ", "),
      "; response: ", x$component, "\n", sep = "")
  cat("D2 (point fits): ", paste(sapply(x$point, function(f) round(f$D2, 3)), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Compare predictor importance between orders within paired draws
#'
#' For each term, the change between two orders in its importance share (or Shapley contribution),
#' computed draw by draw: mean change, 95% percentile interval, the proportion of draws with an
#' increase, and, for importance, the proportion of draws in which the term's rank improves,
#' worsens or is kept.
#'
#' @param x A `zeta_orders` object.
#' @param what `"importance"` or `"shapley"`.
#' @param from,to Orders to compare (default first and last used).
#' @return A data frame with one row per term.
#' @export
compare_orders <- function(x, what = c("importance", "shapley"), from = NULL, to = NULL) {
  what <- match.arg(what)
  A <- if (what == "importance") x$importance else x$shapley
  if (is.null(A)) stop("No `", what, "` results in this object (supply `groups` for Shapley).")
  if (is.null(from)) from <- x$orders[1]
  if (is.null(to)) to <- x$orders[length(x$orders)]
  terms <- dimnames(A)[[3]]
  order_draws <- function(k) {
    matrix(A[, paste0("k", k), , drop = FALSE], nrow = dim(A)[1],
           ncol = length(terms), dimnames = list(NULL, terms))
  }
  I1 <- order_draws(from)
  I2 <- order_draws(to)
  ok <- stats::complete.cases(I1) & stats::complete.cases(I2)
  I1 <- I1[ok, , drop = FALSE]
  I2 <- I2[ok, , drop = FALSE]
  D <- I2 - I1
  interval <- function(prob) {
    if (!nrow(D)) return(rep(NA_real_, length(terms)))
    apply(D, 2, stats::quantile, probs = prob)
  }
  out <- data.frame(term = terms, change = colMeans(D), lo = interval(0.025), hi = interval(0.975),
                    p_increase = colMeans(D > 0), draws = nrow(D), row.names = NULL)
  if (what == "importance") {
    ranks <- function(I) matrix(vapply(seq_len(nrow(I)), function(i) rank(-I[i, ]),
                                      numeric(ncol(I))), nrow(I), ncol(I), byrow = TRUE)
    R1 <- ranks(I1)
    R2 <- ranks(I2)
    out$p_rank_up <- colMeans(R2 < R1)
    out$p_rank_down <- colMeans(R2 > R1)
    out$p_rank_kept <- colMeans(R2 == R1)
  }
  out
}
