.score_input <- function(d, y, trials, sampling_weights, distance, position) {
  if (nrow(d$idx) != length(y) || any(!is.finite(y) | y < 0 | y > 1)) {
    stop("Supply a finite fractional response for every design row.")
  }
  idx <- .species_indices(d$idx, d$n)
  if (d$k >= d$n) stop("Site score deletions require k < n.")
  if (is.null(trials)) trials <- rep(1, length(y))
  if (is.null(sampling_weights)) sampling_weights <- rep(1, length(y))
  eff <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  x <- .design_matrix(d, distance = distance)
  pos <- if (position == "mean") x$Xpos / d$k else x$Xpos
  colnames(pos) <- paste0(colnames(pos), ".position")
  X <- cbind(intercept = 1, x$Xdir, pos)
  W <- sampling_weights / sum(sampling_weights)
  mass <- score_inclusion_cpp(idx, X[, 1L, drop = FALSE], W * 0, W, d$n)[seq_len(d$n), 1L]
  if (any(mass >= 1 - 1e-12)) stop("At least one deletion has no proposal support.")
  list(idx = idx, X = X, y = y, w = W, effective = W * trials, raw = eff / length(y),
       mass = mass, n = d$n, k = d$k)
}

.site_score <- function(p, theta, method) {
  mu <- stats::plogis(drop(p$X %*% theta))
  inc <- score_inclusion_cpp(p$idx, p$X, p$effective * (p$y - mu), p$w, p$n)
  full <- inc[p$n + 1L, -1L]
  incl <- inc[seq_len(p$n), -1L, drop = FALSE]
  deleted <- sweep(-incl, 2L, full, "+") / (1 - p$mass)
  adjusted <- switch(method, deletion = sweep(-(p$n - 1) * deleted, 2L, p$n * full, "+"),
    centered = {
      value <- -(p$n - 1) * deleted
      sweep(value, 2L, colMeans(value) - full)
    },
    inclusion = {
      # Horvitz--Thompson inclusion averages use the raw correction ratios.
      raw <- score_inclusion_cpp(p$idx, p$X, p$raw * (p$y - mu), p$w, p$n)
      full <- raw[p$n + 1L, -1L]
      cond <- raw[seq_len(p$n), -1L, drop = FALSE] * (p$n / p$k)
      sweep(sweep(cond, 2L, full) * (p$k * (p$n - 1) / (p$n - p$k)), 2L, full, "+")
    })
  colnames(adjusted) <- colnames(deleted) <- colnames(p$X)
  names(full) <- colnames(p$X)
  list(pseudovalues = adjusted, full = full, deleted = deleted,
       conditional_ratio = incl / ifelse(p$mass > 0, p$mass, NA_real_),
       inclusion_mass = p$mass, mean_discrepancy = colMeans(adjusted) - full)
}

#' Site score pseudovalues for a fixed MS-GDUM parameter
#'
#' Computes genuine matched set-retention deletion ratios, centred
#' versions, or inverse-inclusion estimates using the raw proposal ratios.
#' Each agrees for complete uniform sampling with fixed model bases.
#' @param design,y,trials,sampling_weights,distance As in [ugdm_fit_variant()].
#' @param theta Intercept, spread and geographic slopes, then position slopes.
#' @param position Scale of the supplied position slopes.
#' @param method Deletion ratios, centred deletion ratios or raw inclusion correction.
#' @return Full and deleted scores, site pseudovalues and numerical diagnostics.
#'   Deleted scores always use ratio weights.
#' @export
ugdm_score_pseudovalues <- function(design, y, theta, trials = NULL,
                                    sampling_weights = NULL, distance = !is.null(design$Xdist),
                                    position = c("mean", "sum"),
                                    method = c("deletion", "centered", "inclusion")) {
  p <- .score_input(design, y, trials, sampling_weights, distance, match.arg(position))
  if (length(theta) != ncol(p$X) || any(!is.finite(theta))) stop("Invalid theta for this design.")
  .site_score(p, theta, match.arg(method))
}

#' Dirichlet bootstrap of matched site pseudovalues
#'
#' Allows any scalar or vector of fitted summaries, including paired
#' effects and Shapley summaries. Preserve site row alignment across
#' columns. No truncation is applied at coefficient boundaries.
#' @param estimate Full-sample estimates, one per deletion column.
#' @param deletions Sites by estimates matrix of matched deletion fits.
#' @param draws Number of Dirichlet draws without supplied weights.
#' @param weights Optional sites by draws matrix, e.g. [site_weights()].
#'   Columns are normalised to sum to one.
#' @param centered Centre the draws at the full-sample estimate.
#' @param seed Optional random seed.
#' @return Bootstrap draws, pseudovalues and exact bootstrap and jackknife
#'   covariance matrices. The draws approximate site-based uncertainty.
#' @export
ugdm_pseudovalue_bootstrap <- function(estimate, deletions, draws = 1000L, weights = NULL,
                                       centered = TRUE, seed = NULL) {
  a <- as.matrix(deletions)
  n <- nrow(a)
  if (n < 2L || length(estimate) != ncol(a) || any(!is.finite(c(estimate, a)))) {
    stop("Supply finite matched deletions for every estimate.")
  }
  if (is.null(weights)) weights <- site_weights(n, .set_design_integer(draws, "draws"), seed)
  weights <- as.matrix(weights)
  if (nrow(weights) != n || !ncol(weights) || any(!is.finite(weights) | weights < 0) ||
      any(colSums(weights) <= 0)) stop("Invalid site weights.")
  weights <- sweep(weights, 2L, colSums(weights), "/")
  V <- sweep(-(n - 1) * a, 2L, n * estimate, "+")
  delta <- sweep(V, 2L, colMeans(V))
  dr <- if (centered) sweep(crossprod(weights, delta), 2L, estimate, "+") else crossprod(weights, V)
  cp <- crossprod(delta)
  list(draws = dr, pseudovalues = V, bootstrap_covariance = cp / (n * (n + 1)),
       jackknife_covariance = cp / (n * (n - 1)))
}

#' Empirical likelihood under zero site moments
#'
#' Solves the ordinary dual or the exponentially tilted alternative.
#' A failure to find a full-rank interior solution is reported explicitly.
#' @param moments Sites by moments matrix.
#' @param method Ordinary empirical likelihood or exponential tilting.
#' @return Log likelihood relative to uniform weights, estimated weights,
#'   and an interior-solution diagnostic.
#' @export
ugdm_empirical_likelihood <- function(moments, method = c("EL", "ETEL")) {
  method <- match.arg(method)
  x <- as.matrix(moments)
  n <- nrow(x)
  p <- ncol(x)
  if (!is.numeric(x) || any(!is.finite(x)) || p < 1L || n <= p) {
    stop("Need finite moments with more sites than moment columns.")
  }
  fail <- function(reason) list(loglik = -Inf, weights = rep(NA_real_, n), ok = FALSE, reason = reason)
  U <- tryCatch(chol(crossprod(x) / n), error = function(e) NULL)
  if (is.null(U) || rcond(U) < 1e-10) return(fail("rank-deficient moments"))
  z <- x %*% backsolve(U, diag(p))
  t <- numeric(p)
  evaluate <- function(t) {
    h <- drop(z %*% t)
    if (method == "EL") {
      den <- 1 + h
      if (any(!is.finite(den) | den <= 0)) return(list(value = -Inf))
      w <- 1 / den
      list(value = mean(log(den)), grad = colMeans(z * w), hess = crossprod(z * w) / n,
           w = w / n, loglik = -sum(log(den)))
    } else {
      expo <- exp(h - max(h))
      w <- expo / sum(expo)
      grad <- colSums(z * w)
      list(value = -max(h) - log(mean(expo)), grad = -grad,
           hess = crossprod(sweep(z, 2L, grad) * sqrt(w)), w = w,
           loglik = sum(h - max(h) - log(mean(expo))))
    }
  }
  a <- evaluate(t)
  for (iter in seq_len(80L)) {
    if (max(abs(a$grad)) < 1e-9 && abs(sum(a$w) - 1) < 1e-6 &&
        max(abs(colSums(z * a$w))) < 1e-7) {
      return(list(loglik = min(a$loglik, 0), weights = a$w / sum(a$w), ok = TRUE, reason = NA_character_))
    }
    step <- tryCatch(solve(a$hess, a$grad), error = function(e) NULL)
    if (is.null(step)) break
    size <- 1
    repeat {
      b <- evaluate(t + size * step)
      if (is.finite(b$value) && b$value >= a$value - 1e-14) break
      if (size < 1e-10) return(fail("no certified interior solution"))
      size <- size / 2
    }
    t <- t + size * step
    if (max(abs(t)) > 1e5 || min(b$w) < 1e-14) break
    a <- b
  }
  fail("no certified interior solution")
}

#' Joint low-dimensional jackknife empirical-likelihood posterior
#'
#' Stacks scores at the site level before solving the empirical-likelihood
#' constraints. Uses a Gaussian prior in mean-position units and a
#' multivariate Student importance proposal. This experimental method
#' requires unconstrained regular fits and a common set of site rows.
#' @param data List of order-specific lists with design, y, and optional trials
#'   and sampling_weights. All designs must share a site index.
#' @param prior_sd Positive prior SD, scalar or one per coefficient.
#' @param draws Number of independent importance proposals.
#' @param method Ordinary or exponentially tilted empirical likelihood.
#' @param score Pseudovalue construction used in the likelihood.
#' @param seed Optional random seed.
#' @return Joint draws and normalised weights, point estimate and numerical
#'   diagnostics (importance ESS, rejected likelihood fraction and timing).
#' @export
ugdm_jel <- function(data, prior_sd = 2, draws = 1500L, method = c("EL", "ETEL"),
                      score = c("centered", "deletion", "inclusion"), seed = NULL) {
  method <- match.arg(method)
  score <- match.arg(score)
  start <- proc.time()[["elapsed"]]
  draws <- .set_design_integer(draws, "draws", 2L)
  if (!length(data)) stop("Supply at least one regular design.")
  if (any(vapply(data, function(d) !is.null(d$fit) && !isTRUE(d$fit$free_sign), TRUE))) {
    stop("JEL currently requires regular unconstrained models. Active constraints use KKT conditions.")
  }
  prep <- lapply(data, function(a) .score_input(a$design, a$y, a$trials, a$sampling_weights,
     !is.null(a$design$Xdist), "mean"))
  n <- prep[[1L]]$n
  blocks <- lengths(lapply(prep, function(p) colnames(p$X)))
  d <- sum(blocks)
  if (any(vapply(prep, function(p) p$n != n, TRUE)) || d > 12L || n <= d + 1L) {
    stop("Use the same sites and at most 12 joint coefficients for this experimental method.")
  }
  if (length(prior_sd) == 1L) prior_sd <- rep(prior_sd, d)
  if (length(prior_sd) != d || any(!is.finite(prior_sd) | prior_sd <= 0)) stop("Invalid prior SD.")
  fits <- lapply(prep, function(p) .variant_solution(p$X, p$y, p$effective,
                                                   matrix(0, ncol(p$X), ncol(p$X)), integer()))
  b <- unlist(lapply(fits, function(f) f$coef), use.names = FALSE)
  if (!all(vapply(fits, function(x) x$converged && !x$constant && x$null_deviance > 1e-15, TRUE))) {
    stop("Posterior requires converged nonconstant unconstrained fits.")
  }
  index <- split(seq_len(d), rep(seq_along(blocks), blocks))
  influence <- lapply(seq_along(prep), function(j) {
    p <- prep[[j]]
    v <- b[index[[j]]]
    mu <- stats::plogis(drop(p$X %*% v))
    H <- crossprod(p$X * sqrt(p$effective * mu * (1 - mu)))
    score <- .site_score(p, v, score)$pseudovalues
    sweep(score, 2L, colMeans(score)) %*% solve(H)
  })
  V <- crossprod(do.call(cbind, influence)) / (n * (n - 1))
  U <- tryCatch(chol(1.5 * V), error = function(e) NULL)
  if (is.null(U) || rcond(U) < 1e-10) stop("Joint site influence is numerically rank deficient.")
  if (!is.null(seed)) set.seed(seed)
  df <- 5
  z <- matrix(stats::rnorm(draws * d), draws, d) / sqrt(stats::rchisq(draws, df) / df)
  sample <- sweep(z %*% U, 2L, b, "+")
  lp <- numeric(draws)
  for (i in seq_len(draws)) {
    candidate <- sample[i, ]
    J <- do.call(cbind, lapply(seq_along(prep), function(j) {
      .site_score(prep[[j]], candidate[index[[j]]], score)$pseudovalues
    }))
    ll <- ugdm_empirical_likelihood(J, method)$loglik
    lp[i] <- ll + sum(stats::dnorm(candidate, sd = prior_sd, log = TRUE)) +
      (d + df) / 2 * log1p(sum(z[i, ]^2) / df)
  }
  valid <- is.finite(lp)
  if (!any(valid)) stop("No posterior proposal passed the moment constraints.")
  w <- exp(lp - max(lp))
  w <- w / sum(w)
  nm <- unlist(lapply(seq_along(prep), function(j) paste0("order", prep[[j]]$k, ".", colnames(prep[[j]]$X))))
  colnames(sample) <- make.unique(nm)
  names(b) <- colnames(sample)
  list(draws = sample, weights = w, estimate = b, ess = 1 / sum(w^2),
       rejected = mean(!valid), max_weight = max(w), method = method, score = score,
       seconds = proc.time()[["elapsed"]] - start)
}
