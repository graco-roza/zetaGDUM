.set_design_integer <- function(x, name, minimum = 1L, maximum = .Machine$integer.max) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) ||
      x < minimum || x > maximum) {
    stop("`", name, "` must be an integer in [", minimum, ", ", maximum, "].")
  }
  as.integer(x)
}

.disjoint_block <- function(n, k, m) {
  matrix(sample.int(n, k * m), m, k, byrow = TRUE)
}

.is_prime_size <- function(n) {
  limit <- floor(sqrt(n))
  n >= 2L && (limit < 2L || !any(n %% seq.int(2L, limit) == 0))
}

#' Randomised balanced incomplete set designs
#'
#' Each row has the uniform law of a distinct-site set, with dependence
#' between rows. Site balancing uses disjoint blocks. Pair balancing
#' uses complete pair cycles at order two and random relabelling of
#' affine cycles for prime site counts at higher orders. At higher orders
#' a complete prime cycle has `n * (n - 1)` rows. A truncated cycle
#' retains uniform individual rows and can have imperfect balance.
#' Environmental stratification draws stratum counts from the uniform
#' set law and uses cyclic random permutations within each stratum.
#'
#' @inheritParams sample_sets
#' @param design Balance sites, pairs, or sites within environmental strata.
#' @param strata Fixed environmental stratum labels, one per site.
#' @return An integer matrix with one set per row. Unit weights target
#'   the uniform distribution. Repeated sets across rows are possible.
#' @export
sample_sets_balanced <- function(n, k, M = 1e4,
                                 design = c("site", "pair", "environment"),
                                 strata = NULL, seed = NULL) {
  design <- match.arg(design)
  n <- .set_design_integer(n, "n", 2L)
  k <- .zeta_order_integer(k, n, 2L)
  M <- .set_design_integer(M, "M")
  if (design == "environment" && (length(strata) != n || anyNA(strata))) {
    stop("Environmental balance requires one nonmissing stratum per site.")
  }
  if (!is.null(seed)) set.seed(seed)
  if (k == n) return(matrix(rep(seq_len(n), each = M), M, k))
  out <- matrix(0L, M, k)
  if (design == "site") {
    size <- n %/% k
    for (first in seq.int(1L, M, by = size)) {
      m <- min(size, M - first + 1L)
      out[seq.int(first, length.out = m), ] <- .disjoint_block(n, k, m)
    }
    return(out)
  }
  if (design == "pair") {
    if (k > 2L && !.is_prime_size(n)) stop("Higher-order affine pair balance needs a prime number of sites.")
    offsets <- if (k == 2L) c(0, cumsum(as.double(seq_len(n - 1L))))
    size <- as.double(n) * (n - 1) / if (k == 2L) 2 else 1
    for (first in seq(1, M, by = size)) {
      m <- min(size, M - first + 1)
      ids <- sample.int(size, m) - 1
      if (k == 2L) {
        last <- findInterval(ids, offsets)
        block <- cbind(as.integer(ids - offsets[last] + 1), last + 1L)
      } else {
        labels <- sample.int(n)
        a <- ids %/% n + 1
        b <- ids %% n
        block <- matrix(labels[((outer(a, seq_len(k) - 1, "*") + b) %% n) + 1], m, k)
      }
      out[seq.int(first, length.out = m), ] <- block
    }
    return(out)
  }
  codes <- match(strata, unique(strata))
  sites <- split(seq_len(n), codes)
  size <- max(1L, n %/% k)
  sizes <- lengths(sites)
  for (first in seq.int(1L, M, by = size)) {
    permutations <- lapply(sites, function(i) i[sample.int(length(i))])
    cursor <- integer(length(sites))
    for (j in seq.int(first, min(first + size - 1L, M))) {
      counts <- tabulate(codes[sample.int(n, k)], length(sites))
      block <- lapply(which(counts > 0), function(h) {
        permutations[[h]][(cursor[h] + seq_len(counts[h]) - 1L) %% sizes[h] + 1L]
      })
      out[j, ] <- unlist(block, use.names = FALSE)
      cursor <- (cursor + counts) %% sizes
    }
  }
  out
}

#' Uniform matched sets for overlap-aware covariance
#'
#' Generates independent groups of randomly selected mutually disjoint
#' sets. Each set marginally has the uniform law. Groups reuse the site
#' population. This is the matched sampling scheme of Xu, Zhu and Shao
#' (2024) for orders at most half the number of sites.
#'
#' @inheritParams sample_sets
#' @param groups Number of matched groups (at least two).
#' @param sets_per_group Number of disjoint sets in each group.
#' @return A list with site indices `idx` and matched group labels `group`.
#' @export
sample_sets_matched <- function(n, k, groups = 100L, sets_per_group = n %/% k, seed = NULL) {
  n <- .set_design_integer(n, "n", 4L)
  k <- .zeta_order_integer(k, n %/% 2L, 2L)
  m <- .set_design_integer(sets_per_group, "sets_per_group", 2L, n %/% k)
  B <- .set_design_integer(groups, "groups", 2L, .Machine$integer.max %/% m)
  if (!is.null(seed)) set.seed(seed)
  idx <- do.call(rbind, lapply(seq_len(B), function(j) .disjoint_block(n, k, m)))
  list(idx = idx, group = rep(seq_len(B), each = m))
}

.kernel_matrix <- function(values, n) {
  values <- as.matrix(values)
  if (!is.numeric(values) || nrow(values) != n || !ncol(values) || any(!is.finite(values))) {
    stop("Kernel values must be finite numeric vectors or matrices with one row per set.")
  }
  values
}

.sandwich_covariance <- function(covariance, hessian) {
  if (is.null(hessian)) return(NULL)
  if (!is.matrix(hessian) || !identical(dim(hessian), dim(covariance)) ||
      any(!is.finite(hessian)) || !isSymmetric(hessian)) stop("Supply a compatible symmetric Hessian.")
  Hinv <- chol2inv(chol(hessian))
  res <- Hinv %*% covariance %*% Hinv
  (res + t(res)) / 2
}

.minimum_eigenvalue <- function(x) min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)

#' Matched-sample covariance of a fixed set kernel
#'
#' Matrix version of the matched-sample variance estimator in Xu, Zhu and
#' Shao (2024). Reports the covariance of the incomplete matched mean
#' and its bias-corrected complete-set counterpart. These finite-sample
#' results require independent sites and a fixed kernel. A score evaluated
#' at a fitted parameter yields an experimental plug-in covariance.
#' Negative eigenvalues are reported and preserved.
#'
#' @param scores Numeric kernel values, one row per matched set.
#' @param matched Matched sets from [sample_sets_matched()].
#' @param hessian Optional positive-definite symmetric Hessian for
#'   sandwich propagation to a regular interior parameter estimate.
#' @return Kernel covariance matrices, within and total covariance,
#'   minimum eigenvalues, and optional sandwich matrices.
#' @export
u_matched_covariance <- function(scores, matched, hessian = NULL) {
  idx <- as.matrix(matched$idx)
  if (!is.numeric(idx) || !length(idx) || ncol(idx) < 2L || anyNA(idx) || any(!is.finite(idx) | idx < 1 | idx != floor(idx)) ||
      length(matched$group) != nrow(idx) || anyNA(matched$group)) {
    stop("Supply valid matched sets and group labels.")
  }
  scores <- .kernel_matrix(scores, nrow(idx))
  parts <- split(seq_len(nrow(idx)), matched$group)
  sizes <- lengths(parts)
  m <- sizes[1L]
  B <- length(parts)
  if (B < 2L || m < 2L || any(sizes != m) ||
      any(vapply(parts, function(i) anyDuplicated(as.vector(idx[i, , drop = FALSE])) > 0L, logical(1)))) {
    stop("Use at least two equal-size groups of at least two mutually disjoint sets.")
  }
  N <- m * B
  within <- Reduce("+", lapply(parts, function(i) stats::cov(scores[i, , drop = FALSE]))) / B
  total <- stats::cov(scores)
  delta <- (m - 1) / (N - 1)
  complete <- (within - total) / (1 - delta)
  incomplete <- within - (N - 1) / N * total
  list(complete = complete, incomplete = incomplete, within = within, total = total,
       min_eigenvalue = c(complete = .minimum_eigenvalue(complete),
                          incomplete = .minimum_eigenvalue(incomplete)),
       sandwich_complete = .sandwich_covariance(complete, hessian),
       sandwich_incomplete = .sandwich_covariance(incomplete, hessian))
}

#' Overlap-resolved covariance of a fixed set kernel
#'
#' Independently draws pairs of sets at each overlap. Half the outer
#' product of the two kernel differences estimates a covariance
#' contrast without plugging in a mean score. The complete covariance
#' uses exact hypergeometric overlap probabilities and the disjoint
#' contrast. This experimental estimator needs independent sites,
#' a fixed kernel and `2 * k <= n`.
#'
#' @inheritParams sample_sets
#' @param kernel Function mapping an integer matrix of sets to one numeric
#'   kernel row per set.
#' @param pairs_per_overlap Number of independently drawn pairs per overlap.
#' @param hessian Optional Hessian as in [u_matched_covariance()].
#' @return Overlap probabilities, covariance contrasts relative to zero
#'   overlap, complete covariance, and optional sandwich covariance.
#'   Estimates can have negative eigenvalues.
#' @export
u_overlap_covariance <- function(n, k, kernel, pairs_per_overlap = 1000L, hessian = NULL, seed = NULL) {
  n <- .set_design_integer(n, "n", 4L)
  k <- .zeta_order_integer(k, n %/% 2L, 2L)
  M <- .set_design_integer(pairs_per_overlap, "pairs_per_overlap")
  if (!is.null(seed)) set.seed(seed)
  pi <- stats::dhyper(0:k, k, n - k, k)
  differences <- lapply(seq_len(k) - 1L, function(j) {
    size <- 2L * k - j
    u <- matrix(vapply(seq_len(M), function(i) sample.int(n, size), integer(size)),
                M, size, byrow = TRUE)
    first <- u[, seq_len(k), drop = FALSE]
    second <- u[, c(seq_len(j), seq.int(k + 1L, size)), drop = FALSE]
    a <- .kernel_matrix(kernel(first), M)
    b <- .kernel_matrix(kernel(second), M)
    if (ncol(a) != ncol(b)) stop("Kernel dimensions vary between set evaluations.")
    crossprod(a - b) / (2 * M)
  })
  zero <- differences[[1L]]
  differences[[k + 1L]] <- zero * 0
  contrasts <- lapply(differences, function(x) zero - x)
  complete <- Reduce("+", Map("*", pi, contrasts))
  list(overlap = 0:k, probability = pi, contrast = contrasts, complete = complete,
       min_eigenvalue = .minimum_eigenvalue(complete),
       sandwich = .sandwich_covariance(complete, hessian))
}
