#' Fixed order-two kernel covariance and its Hoeffding components
#'
#' Estimates covariance of a complete, symmetric order-two U-statistic
#' from products of edges having two, one or zero shared sites. Each
#' unordered pair is included once. Independent identically distributed
#' sites and a fixed square-integrable kernel justify unbiasedness.
#' The resulting matrices can be indefinite in a finite sample.
#' @param scores Numeric matrix of fixed vector kernel evaluations, one
#'   row per pair. A numeric vector represents a scalar kernel.
#' @param pairs Matrix with two integer site indices per row.
#' @param n Number of sites. Requires at least four.
#' @return Covariance estimates for the first and second Hoeffding
#'   components, their sum, and the delete-one-site jackknife covariance.
#' @export
ugdm_pair_covariance <- function(scores, pairs, n) {
  n <- .set_design_integer(n, "n", 4L)
  pairs <- .species_indices(pairs, n)
  scores <- as.matrix(scores)
  m <- choose(n, 2L)
  if (ncol(pairs) != 2L || nrow(pairs) != m || nrow(scores) != m ||
      !ncol(scores) || !is.numeric(scores) || any(!is.finite(scores))) {
    stop("Supply one kernel evaluation per unordered site pair.")
  }
  low <- pmin(pairs[, 1L], pairs[, 2L])
  high <- pmax(pairs[, 1L], pairs[, 2L])
  if (any(low == high) || anyDuplicated(low + as.double(n) * (high - 1L))) {
    stop("Pairs must be complete and unique.")
  }

  # A constant shift cancels in every covariance component. Centring
  # protects the disjoint-edge product from numerical cancellation.
  x <- sweep(scores, 2L, colMeans(scores))
  endpoints <- c(pairs[, 1L], pairs[, 2L])
  rows <- rowsum(rbind(x, x), endpoints)
  same <- crossprod(x)
  shared <- crossprod(rows) - 2 * same
  disjoint <- tcrossprod(colSums(x)) - same - shared
  A <- same / m
  B <- shared / (6 * choose(n, 3L))
  C <- disjoint / (6 * choose(n, 4L))
  first <- 4 * (B - C) / n
  second <- (A - 2 * B + C) / m
  list(first = first, second = second, covariance = first + second,
    jackknife = crossprod(rows) * 4 / (n * (n - 1) * (n - 2)^2))
}
