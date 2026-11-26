#' Environmental spread design for MS-GDUM variants
#'
#' Model B uses mean pairwise absolute differences of standardised site
#' predictors. Models A and C retain the existing environmental splines.
#' The geographic spline is shared by all three designs.
#' @param design Design from [ugdm_path_data()] or [ugdm_design()].
#' @param model Original spline A, linear spread B, or penalised spline C.
#' @param env Original environment, needed for B if the design has no cached
#'   linear spread columns.
#' @return A design for [ugdm_fit_variant()] and [ugdm_predict()].
#' @export
ugdm_model_design <- function(design, model = c("A", "B", "C"), env = NULL) {
  model <- match.arg(model)
  if (!is.null(design$spread_width)) {
    if (model == "B") return(design)
    stop("Use an original spline design as the starting point for A or C.")
  }
  if (model != "B") return(design)
  if (is.null(design$Xlinear)) {
    if (is.null(env)) stop("Supply env to calculate linear spread.")
    e <- .environment_for_moments(env, design)$z
    mat <- vapply(seq_len(ncol(e)), function(j) {
      meanpair_dist_cpp(abs(outer(e[, j], e[, j], "-")), design$idx)
    }, numeric(nrow(design$idx)))
    design$Xlinear <- matrix(mat, nrow(design$idx), ncol(e))
  }
  design$Xdir <- design$Xlinear
  colnames(design$Xdir) <- paste0(design$vars, "_linear")
  design$spread_width <- 1L
  design
}

#' Integrated curvature of the environmental GDUM basis
#'
#' Returns the quadratic form for the integrated squared second derivative
#' on each stored standardised environmental range. Two-point Gaussian
#' quadrature is exact for the degree-two I-splines used by the design.
#' @param design Original environmental spline design from [ugdm_design()].
#' @return Named list of positive semidefinite matrices, one per predictor.
#' @export
ugdm_curvature_penalty <- function(design) {
  basis <- design$environment_basis$splines
  if (is.null(basis) || !is.null(design$spread_width)) stop("Supply an environmental I-spline design.")
  out <- lapply(basis, function(b) {
    limits <- attr(b, "Boundary.knots")
    knots <- attr(b, "knots")
    degree <- attr(b, "degree")
    if (!identical(as.integer(degree), 2L)) stop("Curvature currently requires degree-two I-splines.")
    breaks <- sort(unique(c(limits, knots)))
    half <- diff(breaks) / 2
    center <- utils::head(breaks, -1L) + half
    nodes <- c(center - half / sqrt(3), center + half / sqrt(3))
    d <- splines2::iSpline(nodes, knots = knots, degree = degree,
      Boundary.knots = limits, intercept = attr(b, "intercept"), derivs = 2L)
    crossprod(d * sqrt(rep(half, 2L)))
  })
  stats::setNames(out, design$vars)
}

.variant_penalty <- function(d, D, smooth, ridge, regularization = "adjacent") {
  p <- ncol(D) + ncol(d$Xpos)
  P <- diag(c(0, rep(ridge, p)), p + 1L)
  if (smooth > 0) {
    blocks <- if (regularization == "curvature") ugdm_curvature_penalty(d) else
      rep(list(crossprod(diff(diag(3L)))), length(d$vars))
    for (j in seq_along(d$vars)) {
      idx <- 1L + 3L * (j - 1L) + seq_len(3L)
      P[idx, idx] <- P[idx, idx] + smooth * blocks[[j]]
    }
  }
  P
}

.variant_solution <- function(X, y, w, P, nonnegative, start = NULL, score_shift = NULL) {
  w <- w / sum(w)
  mu0 <- pmin(pmax(sum(w * y), 1e-10), 1 - 1e-10)
  null <- .binom_dev(y, rep(mu0, length(y)), w)
  a <- stats::qlogis(mu0)
  if (is.null(score_shift) && all(y[w > 0] == y[which(w > 0)[1L]])) {
    return(list(coef = c(a, rep(0, ncol(X) - 1L)), converged = TRUE,
                deviance = null, null_deviance = null, projected_gradient = 0, constant = TRUE))
  }
  m <- colSums(X * w)
  m[1L] <- 0
  S <- sweep(X, 2L, m)
  scale <- sqrt(colSums(S^2 * w))
  scale[!is.finite(scale) | scale < 1e-12] <- 1
  Z <- sweep(S, 2L, scale, "/")
  Pscaled <- P / outer(scale, scale)
  shift <- if (is.null(score_shift)) rep(0, ncol(X)) else
    (score_shift - score_shift[1L] * m) / scale
  lower <- rep(-Inf, ncol(X))
  lower[nonnegative] <- 0
  fn <- function(b) {
    eta <- drop(Z %*% b)
    sum(w * (.softplus(eta) - y * eta)) + drop(crossprod(b, Pscaled %*% b)) / 2 - sum(shift * b)
  }
  gr <- function(b) {
    drop(crossprod(Z, w * (stats::plogis(drop(Z %*% b)) - y)) + Pscaled %*% b - shift)
  }
  initial <- if (is.null(start)) c(a, rep(0, ncol(X) - 1L)) else {
    c(sum(start * c(1, m[-1L])), start[-1L] * scale[-1L])
  }
  fit <- stats::optim(pmax(initial, lower), fn, gr, method = "L-BFGS-B",
    lower = lower, control = list(maxit = 1200L, factr = 1e3, pgtol = 1e-9))
  g <- gr(fit$par)
  active <- which(is.finite(lower) & fit$par <= lower + 1e-8)
  g[active] <- pmin(g[active], 0)
  b <- fit$par / scale
  b[1L] <- b[1L] - sum(m * b)
  list(coef = b, converged = fit$convergence == 0L && max(abs(g)) <= 1e-5,
       projected_gradient = max(abs(g)), constant = FALSE,
       deviance = .binom_dev(y, stats::plogis(drop(X %*% b)), w), null_deviance = null)
}

#' Fit original, linear spread or penalised multisite GDUM
#'
#' Fitting is done in canonical mean-position units, with scaling for the
#' numerical optimiser. The standard prediction interface stores position
#' coefficients in sum units. The requested representation is also returned.
#' Model C penalises adjacent coefficients or integrated curvature and
#' optionally shrinks all slopes in canonical units.
#' @param design Matching design from [ugdm_model_design()].
#' @param y,trials,sampling_weights,distance,free_sign As in [ugdm_fit()].
#' @param model Original A, reduced B or regularised C.
#' @param position Representation of the reported position coefficients.
#' @param smooth,ridge Nonnegative strengths in the average working loss for C.
#' @param start Optional fit with matching model and coefficient layout, for warm starts.
#' @param regularization Environmental smoothness penalty for C. The
#'   integrated curvature leaves linear functions unpenalised.
#' @return A fit usable with [ugdm_predict()] on the matching variant design.
#' @export
ugdm_fit_variant <- function(design, y, trials = rep(1, length(y)),
                              model = c("A", "B", "C"), position = c("sum", "mean"),
                              smooth = .001, ridge = .0001, sampling_weights = NULL,
                              distance = TRUE, free_sign = FALSE, start = NULL,
                              regularization = c("adjacent", "curvature")) {
  regularization <- match.arg(regularization)
  model <- match.arg(model)
  position <- match.arg(position)
  if ((model == "B") != !is.null(design$spread_width)) stop("Use the matching variant design.")
  if (length(y) != nrow(design$idx) || any(!is.finite(y) | y < 0 | y > 1)) {
    stop("Supply one finite fractional response per design row.")
  }
  if (length(smooth) != 1L || length(ridge) != 1L ||
      any(!is.finite(c(smooth, ridge)) | c(smooth, ridge) < 0)) stop("Invalid penalty strengths.")
  w <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  d <- .design_matrix(design, distance = distance)
  X <- cbind(intercept = 1, d$Xdir, d$Xpos / design$k)
  nd <- ncol(d$Xdir)
  nc <- if (free_sign) integer() else 1L + seq_len(nd)
  P <- .variant_penalty(d = design, D = d$Xdir,
    smooth = if (model == "C") smooth else 0, ridge = if (model == "C") ridge else 0, regularization = regularization)
  initial <- if (!is.null(start)) c(start$intercept, start$beta, design$k * start$lambda)
  if (!is.null(initial) && (length(initial) != ncol(X) || any(!is.finite(initial)))) {
    stop("Start must match this coefficient layout.")
  }
  f <- .variant_solution(X, y, w, P, nc, initial)
  beta <- f$coef[1L + seq_len(nd)]
  mean_coef <- stats::setNames(utils::tail(f$coef, length(design$vars)), design$vars)
  width <- if (model == "B") 1L else 3L
  nl <- width * length(mean_coef)
  h <- stats::setNames(rowSums(matrix(beta[seq_len(nl)], ncol = width, byrow = TRUE)), design$vars)
  has_geo <- distance && !is.null(design$Xdist)
  structure(list(beta = unname(beta), lambda = mean_coef / design$k, height = h,
    dist_height = if (has_geo) sum(utils::tail(beta, 3L)) else NA_real_,
    intercept = unname(f$coef[1L]), D2 = if (!f$constant && f$null_deviance > 1e-15) 1 - f$deviance / f$null_deviance else NA_real_,
    deviance = f$deviance, null_deviance = f$null_deviance, converged = f$converged,
    projected_gradient = f$projected_gradient, vars = design$vars, has_distance = has_geo,
    model = model, free_sign = free_sign, position = position, penalty = P, regularization = regularization,
    position_coefficients = if (position == "mean") mean_coef else mean_coef / design$k),
    class = "ugdm_fit")
}

#' Scaled information and spread-position allocations
#' @param fit Fitted variant.
#' @param design Matching variant design.
#' @param y,trials,sampling_weights As in [ugdm_fit_variant()].
#' @param contrasts Optional rows of coefficient contrasts in canonical
#'   mean-position units. Projection into near-null information directions
#'   and local inverse curvature are reported per contrast.
#' @return Working information rank and scaled condition, unpenalised fitted
#'   deviance and two-block Shapley allocations. For C each coalition refits
#'   the same restriction of the original penalty.
#' @export
ugdm_variant_diagnostics <- function(fit, design, y, trials = rep(1, length(y)),
                                     sampling_weights = NULL, contrasts = NULL) {
  w <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  d <- .design_matrix(design, distance = fit$has_distance)
  X <- cbind(intercept = 1, d$Xdir, d$Xpos / design$k)
  mu <- ugdm_predict(fit, design)
  info <- crossprod(X * sqrt(w * mu * (1 - mu) / sum(w)))
  scale <- sqrt(diag(info))
  scale[scale == 0] <- 1
  eig <- eigen(info / outer(scale, scale), symmetric = TRUE, only.values = is.null(contrasts))
  ev <- eig$values
  rank <- sum(ev > max(ev) * length(ev) * .Machine$double.eps)
  directions <- NULL
  if (!is.null(contrasts)) {
    L <- as.matrix(contrasts)
    if (ncol(L) != ncol(X) || any(!is.finite(L))) stop("Supply finite contrasts matching the coefficient layout.")
    projections <- sweep(L, 2L, scale, "/") %*% eig$vectors
    good <- ev > max(ev) * 1e-8
    magnitude <- rowSums(projections^2)
    weak <- rowSums(projections[, !good, drop = FALSE]^2)
    inverse <- rowSums(sweep(projections[, good, drop = FALSE]^2, 2L, ev[good], "/"))
    directions <- data.frame(weak_direction_fraction = ifelse(magnitude > 0, weak / magnitude, 0),
                             inverse_curvature = inverse)
  }
  spread <- 1L + seq_len(ncol(d$Xdir))
  pos <- 1L + ncol(d$Xdir) + seq_len(ncol(d$Xpos))
  coalition <- function(idx, constrained) {
    keep <- c(1L, idx)
    f <- .variant_solution(X[, keep, drop = FALSE], y, w, fit$penalty[keep, keep, drop = FALSE],
                           if (constrained) 1L + seq_along(idx) else integer())
    if (f$converged && !f$constant && f$null_deviance > 1e-15) 1 - f$deviance / f$null_deviance else NA_real_
  }
  e <- coalition(spread, !fit$free_sign)
  p <- coalition(pos, FALSE)
  shapley <- c(spread = (fit$D2 + e - p) / 2, position = (fit$D2 + p - e) / 2)
  list(rank = rank, parameters = ncol(X), condition = if (rank == ncol(X) && length(ev))
         max(ev) / min(ev) else Inf, active = if (fit$free_sign) 0L else sum(fit$beta <= 1e-8),
       contrast_information = directions, shapley = shapley, relative_shapley = if (is.finite(fit$D2) && fit$D2 > 1e-12)
         shapley / fit$D2 else shapley * NA_real_)
}
