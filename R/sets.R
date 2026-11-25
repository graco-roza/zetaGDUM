#' Sample sets of sites
#'
#' All pairs when `k = 2` and the number of pairs does not exceed `M`; otherwise `M`
#' random sets of `k` distinct sites.
#'
#' @param n Number of sites.
#' @param k Set size (zeta order).
#' @param M Maximum number of sets.
#' @param seed Optional random seed.
#' @return Integer matrix with one set per row (1-based site indices).
#' @export
#' @examples
#' dim(sample_sets(20, 2))
#' dim(sample_sets(100, 5, M = 1000, seed = 1))
sample_sets <- function(n, k, M = 2e4, seed = NULL) {
  if (k < 2 || k > n) stop("`k` must be between 2 and `n`.")
  if (!is.null(seed)) set.seed(seed)
  if (k == 2 && choose(n, 2) <= M) {
    idx <- t(utils::combn(n, 2))
  } else {
    idx <- t(vapply(seq_len(M), function(j) sample.int(n, k), integer(k)))
  }
  storage.mode(idx) <- "integer"
  idx
}

#' Set-wise similarity components
#'
#' For each set of sites, computes shared species, minimum and mean richness, the
#' Simpson-type similarity (`shared / min richness`; replacement, or turnover), the
#' Sorensen-type similarity (`shared / mean richness`; total) and the nestedness-resultant
#' component (their difference), so that `1 - total = (1 - turnover) + nestedness`. With
#' abundance data it returns the abundance-based counterparts: the summed minimum abundance
#' `a`, balanced variation (`a / min total abundance`), the Bray-Curtis-type similarity
#' (`a / mean total abundance`) and the abundance-gradient component (their difference).
#' For pairs these are the Baselga (2010, 2013) components.
#'
#' @param x Numeric sites x species matrix. If it contains values other than 0/1 the
#'   abundance components are returned as well.
#' @param idx Integer matrix of sets (rows) from [sample_sets()].
#' @param abundance Logical; compute abundance components (default: if `x` is not binary).
#' @return A data frame with one row per set.
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rpois(300, 1), 30, 10)
#' idx <- sample_sets(30, 3, M = 500, seed = 2)
#' head(zeta_sets(x, idx))
zeta_sets <- function(x, idx, abundance = NULL) {
  x <- .as_abundance(x)
  if (is.null(abundance)) abundance <- any(x > 1 | (x > 0 & x < 1))
  P <- (x > 0) * 1L
  storage.mode(P) <- "integer"
  storage.mode(idx) <- "integer"
  rich <- rowSums(P)
  rs <- matrix(rich[idx], nrow(idx))
  m <- apply(rs, 1, min)
  rb <- rowMeans(rs)
  sh <- shared_species_cpp(P, idx)
  out <- data.frame(shared = sh, richness_min = m, richness_mean = rb,
                    turnover = ifelse(m > 0, sh / m, NA_real_), total = ifelse(rb > 0, sh / rb, NA_real_))
  out$nestedness <- out$turnover - out$total
  out$richness_heterogeneity <- ifelse(rb > 0, 1 - m / rb, NA_real_)
  if (abundance) {
    a <- min_abundance_cpp(x, idx)
    tot <- rowSums(x)
    ns <- matrix(tot[idx], nrow(idx))
    nmin <- apply(ns, 1, min)
    nbar <- rowMeans(ns)
    out$abundance_shared <- a
    out$abundance_min_total <- nmin
    out$abundance_mean_total <- nbar
    out$balanced <- ifelse(nmin > 0, a / nmin, NA_real_)
    out$bray_curtis <- ifelse(nbar > 0, a / nbar, NA_real_)
    out$gradient <- out$balanced - out$bray_curtis
    out$abundance_heterogeneity <- ifelse(nbar > 0, 1 - nmin / nbar, NA_real_)
  }
  out
}

#' Pairwise partitions of dissimilarity
#'
#' For every pair of sites returns (a) the Baselga (2010) partition of Sorensen dissimilarity
#' into turnover (Simpson, `min(b, c) / (a + min(b, c))`) and nestedness-resultant, and (b) the
#' Podani-Schmera / Legendre (2014) Jaccard-based partition into replacement
#' (`2 min(b, c) / (a + b + c)`) and richness difference (`|b - c| / (a + b + c)`), where
#' `a` is the number of shared species and `b`, `c` the species unique to each site.
#' These are defined for pairs only; for sets use [zeta_sets()].
#'
#' @param x Numeric sites x species matrix (only presence is used).
#' @return A data frame with one row per pair and columns `i`, `j`, `a`, `b`, `c`,
#'   `sorensen`, `simpson`, `nestedness`, `jaccard`, `replacement`, `richness_difference`.
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rbinom(300, 1, 0.3), 30, 10)
#' head(pair_partition(x))
pair_partition <- function(x) {
  x <- .as_abundance(x)
  P <- (x > 0) * 1
  idx <- t(utils::combn(nrow(P), 2))
  a <- rowSums(P[idx[, 1], , drop = FALSE] * P[idx[, 2], , drop = FALSE])
  r <- rowSums(P)
  b <- r[idx[, 1]] - a
  cc <- r[idx[, 2]] - a
  sor <- (b + cc) / (2 * a + b + cc)
  sim <- ifelse(a + pmin(b, cc) > 0, pmin(b, cc) / (a + pmin(b, cc)), NA_real_)
  jac <- (b + cc) / (a + b + cc)
  data.frame(i = idx[, 1], j = idx[, 2], a = a, b = b, c = cc, sorensen = sor, simpson = sim, nestedness = sor - sim,
             jaccard = jac, replacement = 2 * pmin(b, cc) / (a + b + cc), richness_difference = abs(b - cc) / (a + b + cc))
}

#' Sample sets of sites from a weighted distribution
#'
#' Draws `M` sets of `k` distinct sites with probability proportional to the product of the site weights
#' (conditional Poisson sampling). This is the distribution that the set weights of [set_weights()] define over
#' all sets, so every sampled set has the same weight and the effective number of sets equals `M`; weighting uniformly
#' sampled sets instead concentrates the weight on a few sets at high orders.
#'
#' @param omega Positive site weights, for example one column of [site_weights()].
#' @param k Number of sites per set.
#' @param M Number of sets.
#' @param seed Optional random seed.
#' @return An integer matrix with `M` rows and `k` columns of site indices.
#' @export
#' @examples
#' w <- site_weights(50, 1, seed = 1)[, 1]
#' idx <- sample_sets_weighted(w, k = 5, M = 100, seed = 1)
#' dim(idx)
sample_sets_weighted <- function(omega, k, M = 1e4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (k < 2) stop("`k` must be at least 2.")
  cps_sets_cpp(as.numeric(omega), as.integer(k), as.integer(M))
}
