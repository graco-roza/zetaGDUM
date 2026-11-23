#' Dominance matrix
#'
#' Converts an abundance matrix into a presence-absence matrix of local
#' dominance: in each site, the smallest set of species whose cumulative relative
#' abundance reaches `p` is marked as dominant. Zeta diversity of this matrix is
#' the zeta diversity of dominance (species that are locally dominant in all
#' sites of a set).
#'
#' @param x Numeric sites x species abundance matrix.
#' @param p Cumulative share of the abundance of a site that defines dominance (0 < p <= 1).
#' @return Integer sites x species matrix of 0/1.
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rpois(60, 3), 10, 6)
#' dominance_matrix(x, 0.8)
dominance_matrix <- function(x, p = 0.8) {
  x <- .as_abundance(x)
  if (p <= 0 || p > 1) stop("`p` must be in (0, 1].")
  D <- matrix(0L, nrow(x), ncol(x), dimnames = dimnames(x))
  for (i in seq_len(nrow(x))) {
    v <- x[i, ]
    tot <- sum(v)
    if (tot <= 0) next
    o <- order(-v)
    cs <- cumsum(v[o]) / tot
    m <- which(cs >= p - 1e-12)[1]
    D[i, o[seq_len(m)]] <- 1L
  }
  D
}

#' Exact species weights at a zeta order
#'
#' The share of each species in the zeta diversity of order `k`. For presence-absence
#' data a species with occupancy `o` contributes `choose(o, k)` shared sets, so species
#' with occupancy below `k` have zero weight and raising `k` removes species in increasing
#' order of occupancy. For abundance data the weight is the expected minimum relative
#' abundance over all sets of `k` sites, computed exactly from the sorted abundances.
#'
#' @param x Numeric sites x species matrix (presence-absence or abundance).
#' @param k Zeta order (integer >= 1).
#' @param abundance If `TRUE`, weights follow the abundance (multi-site Czekanowski) definition
#'   based on relative abundances within sites; if `FALSE`, presence-absence.
#' @return Numeric vector of species weights summing to one (zero vector if no species is shared).
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rpois(200, 1), 20, 10)
#' zeta_weights(x, 3)
#' zeta_weights(x, 3, abundance = TRUE)
zeta_weights <- function(x, k, abundance = FALSE) {
  x <- .as_abundance(x)
  n <- nrow(x)
  if (k < 1 || k > n) stop("`k` must be between 1 and the number of sites.")
  v <- if (abundance) .species_zeta_abundance(x, k) else .species_zeta_presence(x, k)
  s <- sum(v)
  if (s <= 0) return(rep(0, length(v)))
  v / s
}

# per-species contribution to zeta of order k, presence-absence
.species_zeta_presence <- function(x, k) {
  n <- nrow(x)
  o <- colSums(x > 0)
  l <- lchoose(o, k) - lchoose(n, k)
  l[!is.finite(l)] <- -Inf
  exp(l)
}

# per-species contribution to zeta of order k, abundance: expected minimum relative abundance
.species_zeta_abundance <- function(x, k) {
  n <- nrow(x)
  rs <- rowSums(x)
  rs[rs == 0] <- 1
  Q <- x / rs
  Qs <- apply(Q, 2, sort, decreasing = FALSE)
  if (is.null(dim(Qs))) Qs <- matrix(Qs, nrow = n)
  j <- seq_len(n - k + 1)
  w <- exp(lchoose(n - j, k - 1) - lchoose(n, k))
  colSums(Qs[j, , drop = FALSE] * w)
}

#' Zeta diversity curve and retention rate
#'
#' Mean number of species (or mean summed minimum relative abundance, or mean number
#' of locally dominant species) shared by `k` sites, for `k = 1, ..., kmax`, computed exactly,
#' and the retention rate `zeta_k / zeta_(k-1)`.
#'
#' @param x Numeric sites x species matrix.
#' @param kmax Largest order (default `min(n - 1, 80)`).
#' @param type `"presence"` (shared species), `"abundance"` (summed minimum relative
#'   abundance) or `"dominance"` (species dominant in all sites; uses `p`).
#' @param p Cumulative abundance share defining dominance when `type = "dominance"`.
#' @return A data frame with columns `k`, `zeta` and `retention`.
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rpois(400, 1), 40, 10)
#' head(zeta_curve(x, kmax = 10))
zeta_curve <- function(x, kmax = NULL, type = c("presence", "abundance", "dominance"), p = 0.8) {
  type <- match.arg(type)
  x <- .as_abundance(x)
  n <- nrow(x)
  if (is.null(kmax)) kmax <- min(n - 1, 80)
  kmax <- min(kmax, n)
  if (type == "dominance") x <- dominance_matrix(x, p)
  ks <- seq_len(kmax)
  z <- vapply(ks, function(k) {
    if (type == "abundance") sum(.species_zeta_abundance(x, k)) else sum(.species_zeta_presence(x, k))
  }, numeric(1))
  data.frame(k = ks, zeta = z, retention = c(NA_real_, z[-1] / z[-length(z)]))
}

#' Breakpoint of the zeta retention rate
#'
#' The first order at which the retention rate reaches a given fraction (default 97%) of
#' its maximum. Beyond it the species still shared are, essentially, those that occupy
#' most sites.
#'
#' @param curve A data frame from [zeta_curve()].
#' @param level Fraction of the maximum retention rate (default 0.97).
#' @return An integer order.
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rbinom(2000, 1, runif(2000, 0.05, 0.9)), 100, 20)
#' zeta_break(zeta_curve(x))
zeta_break <- function(curve, level = 0.97) {
  r <- curve$retention
  ok <- is.finite(r) & curve$zeta > 1e-12
  if (!any(ok)) return(NA_integer_)
  ks <- curve$k[ok]
  rr <- r[ok]
  as.integer(ks[which(rr >= level * max(rr))[1]])
}

.as_abundance <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be numeric.")
  x[is.na(x)] <- 0
  if (any(x < 0)) stop("`x` cannot contain negative values.")
  x
}
