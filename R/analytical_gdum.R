.environment_for_moments <- function(env, design) {
  env <- as.matrix(as.data.frame(env))
  b <- design$environment_basis
  if (is.null(b) || nrow(env) != design$n || !identical(colnames(env), design$vars) ||
      !is.numeric(env) || any(!is.finite(env))) {
    stop("Supply the original numeric `env` and a design with a stored environmental basis.")
  }
  z <- scale(env, center = b$center, scale = b$scale)
  spline <- lapply(seq_len(ncol(z)), function(j) {
    limits <- attr(b$splines[[j]], "Boundary.knots")
    stats::predict(b$splines[[j]], pmin(pmax(z[, j], limits[1L]), limits[2L]))
  })
  list(z = z, spline = do.call(cbind, spline))
}

.environment_moments <- function(z, spline, k, overall, site, pair) {
  n <- nrow(z)
  spread <- vapply(seq_len(ncol(spline)), function(j) {
    d <- abs(outer(spline[, j], spline[, j], "-"))
    sum(d[lower.tri(d)] * pair) / choose(n, 2)
  }, numeric(1))
  c(overall, spread, k / n * drop(crossprod(z, site)))
}

#' Exact response moments for multisite GDUM
#'
#' Calculates the intercept, environmental spread and environmental position
#' entries of the uniform-set expectation of the working weight times the
#' response times the model design. The companion predictor moments exclude
#' the response. The existing geographic distance splines require sampled
#' moments, indicated by `NA`.
#'
#' @param x Sites by species matrix in the same order as `env`.
#' @param env Original sites by environmental predictors used in the design.
#' @param design A design from [ugdm_design()] with a stored environmental basis.
#' @param type `"relative"` for unit-weight relative-abundance dissimilarity
#'   or `"incidence"` for mean-richness-weighted total incidence dissimilarity.
#' @return List of named `response` and `predictor` moment vectors in
#'   intercept, spread (including any distance columns), and position order.
#' @export
ugdm_crossproducts <- function(x, env, design, type = c("relative", "incidence")) {
  type <- match.arg(type)
  e <- .environment_for_moments(env, design)
  k <- design$k
  x <- .conditional_community(x)
  n <- nrow(x)
  if (n != design$n) stop("The community matrix and design must have the same sites.")
  z <- conditional_zeta_dissimilarity(x, k, type, pairs = TRUE)
  response <- .environment_moments(e$z, e$spline, k, z$mean, z$site, z$pair[lower.tri(z$pair)])
  if (type == "incidence") {
    predictor <- .environment_moments(e$z, e$spline, k,
      z$mean + z$shared$mean, z$site + z$shared$site,
      (z$pair + z$shared$pair)[lower.tri(z$pair)])
  } else {
    predictor <- .environment_moments(e$z, e$spline, k, 1, rep(1, n), rep(1, choose(n, 2)))
  }
  nm <- c("intercept", colnames(design$Xdir), colnames(design$Xdist), colnames(design$Xpos))
  insert_distance <- function(v) {
    if (!is.null(design$Xdist)) v <- append(v, rep(NA_real_, ncol(design$Xdist)),
                                           after = 1L + ncol(design$Xdir))
    stats::setNames(v, nm)
  }
  list(response = insert_distance(response), predictor = insert_distance(predictor),
       type = type, order = k)
}

.softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

.moment_certificate <- function(X, w, b, anchor) {
  keep <- w > 0
  x <- X[keep, , drop = FALSE]
  w <- w[keep]
  if (nrow(x) < ncol(x) || qr(x * sqrt(w))$rank < ncol(x)) {
    return(list(ok = FALSE, reason = "rank-deficient positive-weight design"))
  }
  u <- pmin(pmax(stats::plogis(anchor[keep]), 1e-6), 1 - 1e-6)
  variance <- u * (1 - u)
  H <- crossprod(x * sqrt(w * variance))
  g <- b - drop(crossprod(x, w * u))
  shift <- tryCatch(solve(H, g), error = function(e) NULL)
  if (is.null(shift)) return(list(ok = FALSE, reason = "could not solve the moment certificate"))
  pseudo <- u + variance * drop(x %*% shift)
  err <- max(abs(drop(crossprod(x, w * pseudo)) - b))
  margin <- min(pseudo, 1 - pseudo)
  ok <- is.finite(err) && err <= 1e-9 * max(1, max(abs(b))) && margin > 1e-9
  list(ok = ok, reason = if (ok) NULL else "strict interior moment certificate failed",
       residual = err, margin = margin)
}

#' Experimental analytical–Monte Carlo GDUM fit
#'
#' Uses exact response moments for relative abundance or richness-weighted
#' incidence. The control-variate variant centres the nonlinear integral
#' at an intercept-only reference and uses the exact predictor moments.
#' The current geographic spline columns use sampled moments in either
#' variant. Returns the conventional fit if a strict compatibility
#' certificate or optimisation fails, unless `fallback = FALSE`.
#'
#' @inheritParams ugdm_crossproducts
#' @param sampling_weights Raw base-to-proposal probability ratios as
#'   returned by [sample_sets_species()]. Use unit weights for uniform designs.
#' @param free_sign As in [ugdm_fit()].
#' @param method Exact linear term alone (`"hybrid"`) or an additional
#'   intercept-centred control variate (`"control_variate"`).
#' @param fallback Return the conventional fit when the hybrid fails.
#' @param moments Optional cached moments from [ugdm_crossproducts()] for
#'   exactly the same sites, order, predictor basis, and distance columns.
#' @return A `ugdm_fit` object with an additional `analytical` diagnostic.
#'   Deviance summaries use the sampled sets and applied working weights.
#' @export
ugdm_fit_analytical <- function(x, env, design, type = c("relative", "incidence"),
                                sampling_weights = NULL, free_sign = FALSE,
                                method = c("hybrid", "control_variate"), fallback = TRUE,
                                moments = NULL) {
  type <- match.arg(type)
  method <- match.arg(method)
  x <- .conditional_community(x)
  if (nrow(x) != design$n) stop("Community and design site counts differ.")
  if (type == "relative" && any(!is.finite(rowSums(x)) | rowSums(x) <= 0)) {
    stop("Relative abundance requires positive finite site totals.")
  }
  check <- suppressWarnings(ugdm_design(env, design$idx, ref = design))
  if (!isTRUE(all.equal(check$Xdir, design$Xdir)) || !isTRUE(all.equal(check$Xpos, design$Xpos))) {
    stop("`env` and the environmental design columns differ.")
  }
  response <- ugdm_response(x, design$idx, if (type == "relative") "relative" else "total")
  y <- response$y
  trials <- response$trials
  # Zero-richness incidence sets contribute zero to the weighted objective.
  trials[!response$keep] <- 0
  y[!response$keep] <- 0
  weights <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  base <- ugdm_fit(design, y, trials, free_sign = free_sign, sampling_weights = sampling_weights)
  w <- weights / length(y)
  X <- .order_terms(design, base)
  if (is.null(moments)) moments <- ugdm_crossproducts(x, env, design, type)
  if (!identical(moments$type, type) || moments$order != design$k ||
      !identical(names(moments$response), colnames(X)) ||
      !identical(names(moments$predictor), colnames(X))) {
    stop("Cached moments are incompatible with the design, order, or response.")
  }
  exact <- is.finite(moments$response)
  b <- drop(crossprod(X, w * y))
  b[exact] <- moments$response[exact]
  if (method == "control_variate") {
    mean_x <- drop(crossprod(X, w))
    mu0 <- moments$response[[1L]] / moments$predictor[[1L]]
    b[exact] <- b[exact] + mu0 * (mean_x[exact] - moments$predictor[exact])
  } else mu0 <- 0
  old <- c(base$intercept, base$beta, base$lambda)
  certificate <- .moment_certificate(X, w, b, drop(X %*% old))
  fail <- function(reason, convergence = NA_integer_) {
    if (!isTRUE(fallback)) stop(reason)
    base$analytical <- list(requested = method, used = "conventional", reason = reason,
                            certificate = certificate, optimisation_code = convergence)
    base
  }
  if (!certificate$ok) return(fail(certificate$reason))
  lower <- rep(-Inf, ncol(X))
  nd <- length(base$beta)
  if (!free_sign) lower[1L + seq_len(nd)] <- 0
  anchor <- if (mu0 > 0 && mu0 < 1) stats::qlogis(mu0) else 0
  # Subtract centred terms before summation to reduce cancellation. Add back
  # the matching linear term so this remains the certified objective.
  linear <- b - mu0 * drop(crossprod(X, w))
  fn <- function(p) {
    eta <- drop(X %*% p)
    sum(w * (.softplus(eta) - .softplus(anchor) - mu0 * (eta - anchor))) - sum(p * linear)
  }
  gr <- function(p) drop(crossprod(X, w * (stats::plogis(drop(X %*% p)) - mu0))) - linear
  scale <- sqrt(colSums(w * X^2) / sum(w))
  opt <- tryCatch(stats::optim(pmax(old, lower), fn, gr, method = "L-BFGS-B", lower = lower,
    control = list(maxit = 1500L, factr = 1e3, pgtol = 1e-9, parscale = 1 / scale)),
    error = function(e) NULL)
  if (is.null(opt) || opt$convergence != 0L || any(!is.finite(opt$par))) {
    return(fail("analytical optimisation failed", if (is.null(opt)) NA_integer_ else opt$convergence))
  }
  co <- opt$par
  projected <- gr(co)
  at_lower <- is.finite(lower) & co <= lower + 1e-9
  projected[at_lower] <- pmin(projected[at_lower], 0)
  if (max(abs(projected / scale)) > 1e-6 * max(1, sum(w))) {
    return(fail("analytical optimisation did not reach stationarity", opt$convergence))
  }
  base$intercept <- co[1L]
  base$beta <- unname(co[1L + seq_len(nd)])
  base$lambda[] <- co[1L + nd + seq_along(base$lambda)]
  base$height[] <- rowSums(matrix(base$beta[seq_len(3L * length(base$vars))], ncol = 3L, byrow = TRUE))
  if (base$has_distance) base$dist_height <- sum(utils::tail(base$beta, 3L))
  wt <- weights / mean(weights)
  base$deviance <- .binom_dev(y, stats::plogis(drop(X %*% co)), wt)
  base$D2 <- if (base$null_deviance > 0 && is.finite(base$D2)) 1 - base$deviance / base$null_deviance else NA_real_
  base$converged <- TRUE
  base$analytical <- list(requested = method, used = method, reason = NULL,
                          certificate = certificate, projected_gradient = max(abs(projected / scale)),
                          optimisation_code = opt$convergence)
  base
}

#' Observed, fitted and residual local distinctiveness
#'
#' Fitted values are conditional sample estimates. Spread and position
#' share the change from the intercept-only response via two-block Shapley
#' values, with distance counted as spread.
#'
#' @param x Sites by species matrix.
#' @param fit Model from [ugdm_fit()] or [ugdm_fit_analytical()].
#' @param design Set-level design used to obtain predictions.
#' @param type Response used by the fit. For a richness-weighted incidence
#'   fit use `"normalised"` to describe its response on the original scale.
#' @param sampling_weights Optional base-to-proposal ratios. Site summaries
#'   target the base distribution independently of fitting working weights.
#' @return Site table with observed, predicted and residual dissimilarity,
#'   conditional spread and position contributions and intercept baseline.
#'   Unrepresented sites have missing conditional predictions.
#' @export
ugdm_local_distinctiveness <- function(x, fit, design, type = c("relative", "normalised"),
                                       sampling_weights = NULL) {
  type <- match.arg(type)
  n <- nrow(x)
  if (n != design$n) stop("Community and design site counts differ.")
  idx <- design$idx
  w <- .model_weights(nrow(idx), sampling_weights = sampling_weights)
  obs <- if (type == "relative") zeta_lczd(x, design$k, type) else
    zeta_lczd(x, design$k, type, idx, w)
  if (type == "normalised") w[!ugdm_response(x, idx, "total")$keep] <- 0
  d <- .design_matrix(design, fit$vars, fit$has_distance)
  a <- fit$intercept
  s <- drop(d$Xdir %*% fit$beta)
  p <- drop(d$Xpos %*% fit$lambda)
  empty <- stats::plogis(a)
  both <- stats::plogis(a + s + p)
  spread <- .5 * (stats::plogis(a + s) - empty + both - stats::plogis(a + p))
  vals <- .conditional_sample_mean(idx, cbind(both, spread, both - empty - spread), w, n)
  data.frame(site = obs$site, order = design$k, observed = obs$expected_dissimilarity,
             predicted = vals[, 1L], residual = obs$expected_dissimilarity - vals[, 1L],
             baseline = empty, spread = vals[, 2L], position = vals[, 3L])
}
