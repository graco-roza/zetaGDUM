.prepare_score_pools <- function(pools) {
  if (!length(pools)) stop("Provide at least one independent set pool.")
  prepared <- lapply(pools, function(a) .score_input(a$design, a$y, a$trials,
    a$sampling_weights, !is.null(a$design$Xdist), "mean"))
  first <- prepared[[1L]]
  if (!all(vapply(prepared, function(p) p$n == first$n && p$k == first$k &&
                  identical(colnames(p$X), colnames(first$X)), TRUE))) {
    stop("Set pools must share sites, order, and the coefficient layout.")
  }
  first_basis <- pools[[1L]]$design
  if (!all(vapply(pools, function(p)
    identical(p$design$environment_basis, first_basis$environment_basis) &&
      identical(p$design$dist_range, first_basis$dist_range), TRUE))) {
    stop("Set pools must use fixed common environmental and distance bases.")
  }
  prepared
}

.pool_scores <- function(prepared, theta, method) {
  count <- length(prepared)
  if (method != "conditional") {
    if (count == 1L) return(.site_score(prepared[[1L]], theta, method)$pseudovalues)
    return(Reduce("+", lapply(prepared, function(p) .site_score(p, theta, method)$pseudovalues)) / count)
  }
  p <- prepared[[1L]]
  scores <- lapply(prepared, .site_score, theta = theta, method = "deletion")
  full <- Reduce("+", lapply(scores, function(s) s$full)) / count
  cond <- Reduce("+", lapply(scores, function(s) s$conditional_ratio)) / count
  if (any(!is.finite(cond))) stop("A conditional-score pool has a site with no inclusion support.")
  z <- sweep(cond, 2L, colMeans(cond)) * p$k * (p$n - 1) / (p$n - p$k)
  sweep(z, 2L, full, "+")
}

#' Site score pseudovalues from independent set pools
#'
#' Pools share sites, order and predictor bases. Deletion, centred deletion
#' and analytical inclusion average their existing pseudovalues. Conditional
#' averages the site inclusion ratios before centring and applying the exact
#' complete-set identity. It can differ from deletion ratios at finite budgets.
#' @param pools List of set pools with design, y, and optional trials and
#'   sampling_weights. Site indices must identify the same sites in each pool.
#' @param theta Coefficients in canonical mean-position units.
#' @param method Site score construction.
#' @return A sites by coefficients matrix.
#' @export
ugdm_pool_scores <- function(pools, theta,
    method = c("conditional", "deletion", "centered", "inclusion")) {
  p <- .prepare_score_pools(pools)
  if (length(theta) != ncol(p[[1L]]$X) || any(!is.finite(theta))) stop("Incorrect coefficient layout.")
  .pool_scores(p, theta, match.arg(method))
}

#' Numerical covariance diagnostic for paired site pseudovalues
#'
#' Centres each matrix over the same sites. The difference covariance
#' estimates numerical variation from equally precise independent set
#' pools on the same observed community. No covariance correction is applied.
#' @param first,second Finite sites by contrasts matrices with aligned rows.
#' @return Numerical and average jackknife covariance, their diagonal ratios,
#'   plus jackknife covariance of the mean pseudovalues.
#' @export
ugdm_pseudovalue_noise <- function(first, second) {
  a <- as.matrix(first)
  b <- as.matrix(second)
  n <- nrow(a)
  if (!identical(dim(a), dim(b)) || n < 2L || ncol(a) < 1L ||
      !is.numeric(a) || !is.numeric(b) || any(!is.finite(c(a, b)))) {
    stop("Supply finite, equally sized site pseudovalue matrices.")
  }
  a <- sweep(a, 2L, colMeans(a))
  b <- sweep(b, 2L, colMeans(b))
  den <- n * (n - 1)
  mc <- crossprod(a - b) / (2 * den)
  jk <- (crossprod(a) + crossprod(b)) / (2 * den)
  list(numerical_covariance = mc, average_jackknife = jk,
       numerical_fraction = ifelse(diag(jk) > 0, diag(mc) / diag(jk), NA_real_),
       pooled_jackknife = crossprod((a + b) / 2) / den)
}
