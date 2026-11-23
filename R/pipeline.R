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
#'   `"balanced"`, `"bray_curtis"`, `"gradient"`.
#' @param orders Integer vector of orders; default `c(2, midpoint, break)`.
#' @param M Number of random sets per order (all pairs at order 2 when feasible).
#' @param ndraw Number of paired bootstrap draws.
#' @param groups Optional named list of predictor groups for the Shapley decomposition.
#' @param min_share Minimum share of sets sharing at least one species for an order to be used.
#' @param seed Random seed.
#' @param cores Number of cores for the draws.
#' @param paired If `TRUE` (default) the same site weights are used at every order within a draw; if `FALSE` each order uses independent draws (for comparison only).
#' @param free_sign Allow directional effects of either sign (see [ugdm_fit()]). The default (`NULL`) constrains them to be non-negative, as in the pairwise GDUM; setting `TRUE` is advisable for the nestedness and abundance-gradient components, for which there is no expectation about the direction of richness or abundance gradients along environmental gradients.
#' @return An object of class `zeta_orders`.
#' @export
#' @examples
#' sim <- simulate_metacommunity(60, 60, seed = 1)
#' res <- zeta_orders(sim$comm, sim$env, sim$coords, M = 2000, ndraw = 10,
#'                    groups = list(env1 = "env1", env2 = "env2", distance = "distance"))
#' res
#' compare_orders(res, "importance")
zeta_orders <- function(comm, env, coords = NULL, lonlat = FALSE, component = "turnover", orders = NULL, M = 1e4, ndraw = 100,
                        groups = NULL, min_share = 0.04, seed = 1, cores = 1, paired = TRUE, free_sign = NULL) {
  comm <- .as_abundance(comm)
  fs <- isTRUE(free_sign)
  n <- nrow(comm)
  if (nrow(env) != n) stop("`comm` and `env` must have the same number of sites.")
  cur <- zeta_curve(comm, type = "presence")
  kb <- zeta_break(cur)
  P <- (comm > 0) * 1L
  storage.mode(P) <- "integer"
  feasible <- function(k) {
    idx <- sample_sets(n, k, M = 4000, seed = seed)
    mean(shared_species_cpp(P, idx) > 0) >= min_share
  }
  if (is.null(orders)) {
    kbr <- kb
    while (kbr > 3 && !feasible(kbr)) kbr <- kbr - 1
    orders <- sort(unique(c(2L, as.integer(round((2 + kbr) / 2)), as.integer(kbr))))
  }
  dm <- if (!is.null(coords)) coords_dist(coords, lonlat) else NULL
  W <- site_weights(n, ndraw, seed)
  terms <- c(colnames(env), if (!is.null(dm)) "distance")
  imp <- array(NA_real_, c(ndraw, length(orders), length(terms)), dimnames = list(NULL, paste0("k", orders), terms))
  shp <- if (!is.null(groups)) array(NA_real_, c(ndraw, length(orders), length(groups)), dimnames = list(NULL, paste0("k", orders), names(groups))) else NULL
  d2 <- matrix(NA_real_, ndraw, length(orders), dimnames = list(NULL, paste0("k", orders)))
  hgt <- array(NA_real_, c(ndraw, length(orders), length(terms)), dimnames = list(NULL, paste0("k", orders), terms))
  lam <- array(NA_real_, c(ndraw, length(orders), ncol(env)), dimnames = list(NULL, paste0("k", orders), colnames(env)))
  point <- list()
  for (oi in seq_along(orders)) {
    k <- orders[oi]
    idx <- sample_sets(n, k, M = M, seed = seed + k)
    r <- ugdm_response(comm, idx, component)
    idx <- idx[r$keep, , drop = FALSE]
    y <- r$y[r$keep]
    tr <- r$trials[r$keep]
    des <- ugdm_design(env, idx, dm)
    Wo <- if (paired) W else site_weights(n, ndraw, seed + 1000 * oi)
    point[[paste0("k", k)]] <- ugdm_fit(des, y, tr, free_sign = fs)
    one <- function(b) {
      w <- set_weights(Wo[, b], idx)
      f <- tryCatch(ugdm_fit(des, y, tr, w, free_sign = fs), error = function(e) NULL)
      if (is.null(f)) return(NULL)
      list(imp = tryCatch(ugdm_importance(f, des, w, tr), error = function(e) NULL), d2 = f$D2, hgt = c(f$height, if (f$has_distance) c(distance = f$dist_height)), lam = f$lambda,
           shp = if (!is.null(groups)) tryCatch(ugdm_shapley(des, y, tr, groups, w, free_sign = fs), error = function(e) NULL))
    }
    res <- if (cores > 1) parallel::mclapply(seq_len(ndraw), one, mc.cores = cores) else lapply(seq_len(ndraw), one)
    for (b in seq_len(ndraw)) {
      if (is.null(res[[b]])) next
      if (!is.null(res[[b]]$imp)) imp[b, oi, ] <- res[[b]]$imp[terms]
      d2[b, oi] <- res[[b]]$d2
      hgt[b, oi, ] <- res[[b]]$hgt[terms]
      lam[b, oi, ] <- res[[b]]$lam[colnames(env)]
      if (!is.null(shp) && !is.null(res[[b]]$shp)) shp[b, oi, ] <- res[[b]]$shp
    }
  }
  structure(list(curve = cur, break_order = kb, orders = orders, component = component, point = point, importance = imp,
                 shapley = shp, D2 = d2, height = hgt, lambda = lam, n = n), class = "zeta_orders")
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
  I1 <- A[, paste0("k", from), , drop = FALSE][, 1, ]
  I2 <- A[, paste0("k", to), , drop = FALSE][, 1, ]
  ok <- stats::complete.cases(I1) & stats::complete.cases(I2)
  I1 <- I1[ok, , drop = FALSE]
  I2 <- I2[ok, , drop = FALSE]
  D <- I2 - I1
  out <- data.frame(term = colnames(I1), change = colMeans(D), lo = apply(D, 2, stats::quantile, 0.025), hi = apply(D, 2, stats::quantile, 0.975),
                    p_increase = colMeans(D > 0), draws = nrow(D), row.names = NULL)
  if (what == "importance") {
    R1 <- t(apply(-I1, 1, rank))
    R2 <- t(apply(-I2, 1, rank))
    out$p_rank_up <- colMeans(R2 < R1)
    out$p_rank_down <- colMeans(R2 > R1)
    out$p_rank_kept <- colMeans(R2 == R1)
  }
  out
}
