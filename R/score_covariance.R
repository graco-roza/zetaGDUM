.score_covariance_parts <- function(total, products, R, n) {
  average <- products / (R * (n - 1))
  pooled <- crossprod(total / R) / (n - 1)
  cross <- (R * pooled - average) / (R - 1)
  list(average = average, pooled = pooled, cross = cross, numerical = average - cross)
}

#' Independent-pool covariance of joint site-score pseudovalues
#'
#' Estimate average within-pool, pooled and distinct cross-pool covariance at
#' a fixed parameter and basis. Each pool must use independent combination
#' sampling conditional on the same sites. Align site rows and coefficient
#' columns. The cross estimate removes expected numerical error only under
#' conditional unbiasedness and independence, and can be indefinite.
#'
#' Regularisations operate in coordinates scaled by the average within-pool
#' score SDs, even if the cross covariance has negative diagonal entries.
#' Diagonal shrinkage targets the identity. Block shrinkage targets blocks of
#' the standardised within-pool average, so matching effect blocks can retain
#' cross-order covariance. The estimated shrinkage uses leave-one-pool-out
#' numerical variability, conditional on the observed sites. This is a
#' sensitivity estimator inspired by linear covariance shrinkage and
#' does not use combinations as independent ecological observations.
#' @param pools List of two or more sites-by-joint-coefficients matrices.
#' @param method Covariance to regularise or report. Raw uses the first matrix.
#' @param regularization None, diagonal shrinkage, block shrinkage, an explicit
#'   eigenvalue floor, or positive-semidefinite projection in score coordinates.
#' @param blocks Block labels, one per coefficient, for block shrinkage.
#' @param alpha Shrinkage strength in [0,1]. If missing, estimate strength
#'   from at least three independent pools.
#' @param eigen_floor Strictly positive standardised eigenvalue floor for
#'   eigenvalue regularisation. The PSD projection instead uses zero.
#' @param contrasts Optional matrix with one column per joint coefficient.
#' @return Covariance, average, cross, numerical, pooled and raw covariances,
#'   before/after standardised spectra, adjustment diagnostics and numerical
#'   standard errors across pools. Contrast covariances are returned if requested.
#' @export
ugdm_score_covariance <- function(pools, method = c("cross", "pooled", "average", "raw"),
    regularization = c("none", "diagonal", "block", "eigen", "psd"),
    blocks = NULL, alpha = NULL, eigen_floor = NULL, contrasts = NULL) {
  method <- match.arg(method)
  regularization <- match.arg(regularization)
  R <- length(pools)
  if (R < 2L) stop("Supply at least two independent pools.")
  pools <- lapply(pools, as.matrix)
  first <- pools[[1L]]
  n <- nrow(first)
  d <- ncol(first)
  if (n < 2L || d < 1L || any(!vapply(pools, function(a)
    is.numeric(a) && identical(dim(a), dim(first)) && identical(dimnames(a), dimnames(first)) &&
      all(is.finite(a)), logical(1L)))) {
    stop("Supply finite, equally sized matrices with aligned rows and columns.")
  }
  centres <- lapply(pools, function(a) sweep(a, 2L, colMeans(a)))
  products <- lapply(centres, crossprod)
  total <- Reduce("+", centres)
  sum_products <- Reduce("+", products)
  parts <- .score_covariance_parts(total, sum_products, R, n)
  parts$raw <- products[[1L]] / (n - 1)
  per_pool <- lapply(products, function(a) a / (n - 1))
  se <- matrix(NA_real_, d, d, dimnames = dimnames(parts$raw))
  leave_one_out <- NULL
  L <- .conditioned_contrast(contrasts, d)
  if (method == "raw") {
    deviations <- lapply(per_pool, function(a) a - parts$average)
    se <- sqrt(Reduce("+", lapply(deviations, function(a) a^2)) / (R - 1))
  } else if (R >= 3L) {
    leave_one_out <- lapply(seq_len(R), function(i)
      .score_covariance_parts(total - centres[[i]], sum_products - products[[i]], R - 1, n)[[method]])
    mean_leave_one_out <- Reduce("+", leave_one_out) / R
    deviations <- lapply(leave_one_out, function(a) a - mean_leave_one_out)
    se <- sqrt((R - 1) / R * Reduce("+", lapply(deviations, function(a) a^2)))
  }
  omega <- parts[[method]]
  scale <- sqrt(diag(parts$average))
  zero_scale <- which(scale == 0)
  scale[zero_scale] <- 1
  units <- tcrossprod(scale)
  S <- omega / units
  raw_spectrum <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  A <- S
  if (regularization %in% c("diagonal", "block")) {
    target <- diag(as.numeric(seq_len(d) %in% setdiff(seq_len(d), zero_scale)), d)
    if (regularization == "block") {
      if (length(blocks) != d || anyNA(blocks)) stop("Supply one block label per coefficient.")
      target <- parts$average / units * outer(blocks, blocks, "==")
    }
    if (is.null(alpha)) {
      if (R < 3L) stop("Supply alpha explicitly for fewer than three pools.")
      denominator <- sum((S - target)^2)
      alpha <- if (denominator > 0) min(1, sum((se / units)^2) / denominator) else 0
    }
    if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha < 0 || alpha > 1) {
      stop("alpha must lie in [0,1].")
    }
    A <- (1 - alpha) * S + alpha * target
  }
  floor <- NULL
  if (regularization %in% c("psd", "eigen")) {
    floor <- if (regularization == "psd") 0 else eigen_floor
    if (length(floor) != 1L || !is.numeric(floor) || !is.finite(floor) ||
        (regularization == "eigen" && floor <= 0)) {
      stop("Specify a strictly positive eigen_floor for eigenvalue regularisation.")
    }
    ed <- eigen(S, symmetric = TRUE)
    A <- tcrossprod(sweep(ed$vectors, 2L, sqrt(pmax(ed$values, floor)), "*"))
  }
  covariance <- A * units
  dimnames(covariance) <- dimnames(omega)
  contrast_covariance <- contrast_se <- NULL
  if (!is.null(L)) {
    transform <- function(a) L %*% a %*% t(L)
    contrast_covariance <- c(lapply(parts, transform), list(regularized = transform(covariance)))
    if (!is.null(leave_one_out) || method == "raw") {
      values <- lapply(if (method == "raw") per_pool else leave_one_out, transform)
      other <- Reduce("+", values) / R
      factor <- if (method == "raw") 1 / (R - 1) else (R - 1) / R
      contrast_se <- sqrt(factor * Reduce("+", lapply(values, function(a) (a - other)^2)))
    }
  }
  list(covariance = covariance, cross = parts$cross, average = parts$average,
       pooled = parts$pooled, raw = parts$raw, numerical = parts$numerical,
       spectrum = raw_spectrum, adjusted_spectrum = eigen(A, symmetric = TRUE, only.values = TRUE)$values,
       adjustment = sqrt(sum((A - S)^2)), scale = scale, zero_scale = zero_scale,
       alpha = alpha, floor = floor, numerical_se = se, contrasts = contrast_covariance,
       contrast_numerical_se = contrast_se, method = method, regularization = regularization,
       pools = R, sites = n)
}
