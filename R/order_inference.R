.order_effects <- function(fit) {
  c(stats::setNames(fit$height, paste0("spread.", names(fit$height))),
    if (fit$has_distance) c(spread.distance = fit$dist_height),
    stats::setNames(fit$lambda, paste0("position.", names(fit$lambda))))
}

.order_terms <- function(design, fit) {
  d <- .design_matrix(design, fit$vars, fit$has_distance)
  cbind(intercept = 1, d$Xdir, d$Xpos)
}

.order_information <- function(design, fit, w) {
  X <- .order_terms(design, fit)
  mu <- ugdm_predict(fit, design)
  w <- w / sum(w)
  H <- crossprod(X * sqrt(w * mu * (1 - mu)))
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  tol <- max(ev) * ncol(X) * .Machine$double.eps
  rank <- sum(ev > tol)
  list(X = X, mu = mu, H = H, rank = rank, parameters = ncol(X),
       min_eigenvalue = min(ev), condition = if (rank == ncol(X)) max(ev) / min(ev) else Inf,
       active = sum(fit$beta <= 1e-7))
}

.order_effect_map <- function(design, fit) {
  p <- length(fit$beta)
  q <- length(fit$lambda)
  contrast <- matrix(0, length(.order_effects(fit)), 1L + p + q)
  width <- if (is.null(design$spread_width)) 3L else design$spread_width
  for (j in seq_len(q)) contrast[j, 1L + width * (j - 1L) + seq_len(width)] <- 1
  if (fit$has_distance) contrast[q + 1L, 1L + q * width + seq_len(3L)] <- 1
  contrast[q + as.integer(fit$has_distance) + seq_len(q), 1L + p + seq_len(q)] <- diag(q)
  rownames(contrast) <- names(.order_effects(fit))
  contrast
}

.inclusion_sums <- function(idx, a, n) {
  out <- matrix(0, n, ncol(a))
  for (j in seq_len(ncol(idx))) {
    totals <- rowsum(a, idx[, j], reorder = FALSE)
    sites <- as.integer(rownames(totals))
    out[sites, ] <- out[sites, , drop = FALSE] + totals
  }
  out
}

.subset_order_design <- function(d, keep) {
  for (term in c("idx", "Xdir", "Xpos", "Xdist", "Xlinear")) {
    if (!is.null(d[[term]])) d[[term]] <- d[[term]][keep, , drop = FALSE]
  }
  d
}

.kish <- function(w) {
  if (!any(w > 0)) return(0)
  w <- w / max(w)
  sum(w)^2 / sum(w^2)
}

.jackknife_sd <- function(delete) {
  n <- nrow(delete)
  centered <- sweep(delete, 2L, colMeans(delete), "-")
  sqrt((n - 1) / n * colSums(centered^2))
}

#' Site-level uncertainty of spread and position effects
#'
#' Refit after each site deletion, holding the reference predictor basis fixed.
#' Use all sets for exact conditional deletion fits, or an existing probability
#' sample with its proposal weights. In the latter case the deletion fits
#' retain only sampled sets excluding that site and can contain Monte Carlo
#' error. For regular interior fits, a site-score multiplier estimates the
#' first-order projection of the estimating equation. It uses the
#' linearised delete-one estimating equation and also returns the raw
#' inclusion-score multiplier for comparison. These approximations assume independent sites and a nondegenerate projection; it is not
#' calibrated for large order-to-site ratios, spatial dependence, or active
#' monotonicity constraints.
#'
#' @param design,y,trials,sampling_weights,free_sign As in [ugdm_fit()].
#' @param ndraw Number of Gaussian site multiplier draws. Use zero to omit draws.
#' @param seed Optional seed for multiplier draws.
#' @param sites Eligible site indices. Defaults to every site in `design`.
#'   Supply the original training site indices when the design retains rows
#'   for a larger set of sites.
#' @return A list of point effects, delete-one effects and Kish sizes,
#'   jackknife standard errors, site influence values and multiplier effects
#'   (when regular), and information diagnostics. `multiplier_raw` uses
#'   the inclusion-score projection without the deletion normalisation.
#'   Position effects are per
#'   training-standardised unit; spread heights use the common spline bounds.
#' @export
ugdm_site_uncertainty <- function(design, y, trials, sampling_weights = NULL,
                                  free_sign = FALSE, ndraw = 399L, seed = NULL,
                                  sites = seq_len(design$n)) {
  if (length(ndraw) != 1L || !is.finite(ndraw) || ndraw < 0 || ndraw != floor(ndraw)) {
    stop("`ndraw` must be a non-negative integer.")
  }
  if (!is.numeric(sites) || !length(sites) || anyNA(sites) || anyDuplicated(sites) ||
      any(sites != floor(sites) | sites < 1L | sites > design$n) ||
      !all(design$idx %in% sites)) {
    stop("`sites` must contain distinct valid indices and include every site used in the design.")
  }
  n <- length(sites)
  k <- design$k
  weights <- if (is.null(sampling_weights)) rep(1, length(y)) else sampling_weights
  w <- .model_weights(length(y), trials, sampling_weights = weights)
  fit <- ugdm_fit(design, y, trials, free_sign = free_sign, sampling_weights = weights)
  effect <- .order_effects(fit)
  inf <- .order_information(design, fit, w)
  deletes <- matrix(NA_real_, n, length(effect), dimnames = list(NULL, names(effect)))
  deleted_ess <- numeric(n)
  regular <- inf$rank == inf$parameters && fit$converged && is.finite(fit$D2)
  if (k < n && regular) {
    excludes <- split(rep(seq_len(nrow(design$idx)), k), as.vector(design$idx))
    for (i in seq_len(n)) {
      keep <- rep(TRUE, length(y))
      keep[excludes[[as.character(sites[i])]]] <- FALSE
      deleted_ess[i] <- .kish(w[keep])
      if (sum(keep & w > 0) < inf$parameters) next
      d <- .subset_order_design(design, keep)
      f <- ugdm_fit(d, y[keep], trials[keep], free_sign = free_sign, sampling_weights = weights[keep])
      info <- .order_information(d, f, w[keep])
      if (f$converged && is.finite(f$D2) && info$rank == info$parameters) {
        deletes[i, ] <- .order_effects(f)
      }
    }
  }
  jk <- .jackknife_sd(deletes)
  influence <- multiplier <- multiplier_raw <- NULL
  reason <- if (k >= n) "no sets after deletion" else if (!regular) "unidentified or nonconverged fit" else
    if (!free_sign && inf$active > 0) "active directional constraint" else NULL
  if (is.null(reason)) {
    weights_normalized <- w / sum(w)
    score <- (weights_normalized * (y - inf$mu)) * inf$X
    included <- .inclusion_sums(design$idx, cbind(weights_normalized, score), design$n)[sites, , drop = FALSE]
    mass <- included[, 1L]
    score_in <- included[, -1L, drop = FALSE]
    score_centered <- score_in - outer(mass, colSums(score))
    if (any(1 - mass <= .Machine$double.eps)) {
      reason <- "no effective sets after at least one deletion"
    } else {
      map <- solve(inf$H, t(.order_effect_map(design, fit)))
      raw <- n * sweep(score_centered, 2L, colMeans(score_centered), "-") %*% map
      pseudo <- (n - 1L) * score_centered / (1 - mass)
      influence <- sweep(pseudo, 2L, colMeans(pseudo), "-") %*% map
      if (ndraw > 0) {
        if (!is.null(seed)) set.seed(seed)
        multipliers <- matrix(stats::rnorm(ndraw * n), ndraw, n) / sqrt(n * (n - 1))
        multiplier <- sweep(multipliers %*% influence, 2, effect, "+")
        multiplier_raw <- sweep(multipliers %*% raw, 2, effect, "+")
      }
    }
  }
  list(effect = effect, sites = sites, delete = deletes, delete_ess = deleted_ess, jackknife_se = jk,
       influence = influence, multiplier = multiplier, multiplier_raw = multiplier_raw,
       multiplier_unavailable = reason,
       diagnostics = inf[c("rank", "parameters", "min_eigenvalue", "condition", "active")],
       converged = fit$converged)
}

.inference_sets <- function(P, k, M, alpha, exact_limit = 0, omega = NULL, seed = NULL) {
  if (choose(nrow(P), k) <= exact_limit) {
    idx <- t(utils::combn(nrow(P), k))
    weights <- if (is.null(omega)) rep(1, nrow(idx)) else .product_sampling_weights(omega, idx, rep(1, nrow(idx)))
    return(list(idx = idx, sampling_weights = weights, shared = shared_species_cpp(P, idx),
                alpha = 0, exact = TRUE))
  }
  s <- .model_sets(P, k, M, if (alpha > 0) "species" else "conventional", alpha, omega, seed)
  s$exact <- k == 2L && choose(nrow(P), k) <= M
  s
}

.inference_fit <- function(comm, env, sets, component, working_weights, free_sign, groups = NULL, ref = NULL,
                           dist_matrix = NULL, dist_range = NULL) {
  r <- ugdm_response(comm, sets$idx, component, working_weights = working_weights)
  keep <- r$keep
  if (!any(keep)) stop("No sampled sets have a defined response.")
  d <- ugdm_design(env, sets$idx[keep, , drop = FALSE], dist_matrix, ref = ref, dist_range = dist_range)
  w <- sets$sampling_weights[keep]
  y <- r$y[keep]
  tr <- r$trials[keep]
  f <- ugdm_fit(d, y, tr, free_sign = free_sign, sampling_weights = w)
  imp <- ugdm_importance(f, d, sampling_weights = w, trials = tr)
  shp <- if (is.null(groups)) numeric() else
    ugdm_shapley(d, y, tr, groups, free_sign = free_sign, sampling_weights = w)
  names(imp) <- paste0("importance.", names(imp))
  if (length(shp)) names(shp) <- paste0("shapley.", names(shp))
  list(fit = f, design = d, y = y, trials = tr, weights = w,
       effects = c(.order_effects(f), imp, stats::setNames(as.numeric(shp), names(shp))),
       information = .order_information(d, f, .model_weights(length(y), tr, sampling_weights = w)),
       keep = keep)
}
