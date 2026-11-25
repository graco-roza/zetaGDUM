#' Exact species contributions to set-level compositional dissimilarity
#'
#' For each set, partitions total dissimilarity into non-negative species terms.
#' With `type = "presence"`, the response is one minus shared richness divided
#' by mean site richness. With `type = "relative"`, counts are first divided by
#' each site's total, and the response is one minus summed minimum relative
#' abundance. With `type = "absolute"`, the response is one minus summed minimum
#' count divided by mean site total. This last response includes variation in
#' total abundance and requires comparable sampling effort. The relative response
#' does not. In every case the species terms sum exactly to the response.
#'
#' @param x Non-negative numeric sites by species matrix.
#' @param idx Integer matrix of site sets, one set per row.
#' @param type One of `"presence"`, `"relative"`, or `"absolute"`.
#' @return A list with matrices `similarity` and `dissimilarity` (sets by
#'   species), and vectors `response` and `defined`. Undefined sets have `NA`
#'   contributions and response (all sites in a set have zero total, or any
#'   site has zero total for the relative-abundance response).
#' @export
species_contributions <- function(x, idx, type = c("presence", "relative", "absolute")) {
  type <- match.arg(type)
  x <- .as_abundance(x)
  idx <- .species_indices(idx, nrow(x))
  if (type == "presence") x <- (x > 0) * 1
  if (type == "relative") {
    totals <- rowSums(x)
    x <- x / ifelse(totals > 0, totals, 1)
  }
  nr <- nrow(idx)
  ns <- ncol(x)
  sim <- matrix(NA_real_, nr, ns, dimnames = list(NULL, colnames(x)))
  dis <- sim
  for (j in seq_len(nr)) {
    if (type == "relative" && any(totals[idx[j, ]] == 0)) next
    a <- x[idx[j, ], , drop = FALSE]
    mean_sp <- colMeans(a)
    shared_sp <- apply(a, 2L, min)
    denom <- sum(mean_sp)
    if (denom <= 0) next
    sim[j, ] <- shared_sp / denom
    dis[j, ] <- (mean_sp - shared_sp) / denom
  }
  defined <- rowSums(is.na(dis)) == 0L
  list(similarity = sim, dissimilarity = dis,
       response = rowSums(dis), defined = defined)
}

.species_indices <- function(idx, n) {
  idx <- as.matrix(idx)
  if (!is.numeric(idx) || ncol(idx) < 2L || anyNA(idx) ||
      any(idx != floor(idx)) || any(idx < 1L | idx > n) ||
      any(apply(idx, 1L, anyDuplicated) > 0L)) {
    stop("`idx` must contain distinct, valid integer site indices in each set.")
  }
  idx
}

#' Expected species contributions over uniform site sets
#'
#' Computes the exact expectation over all sets of `k` distinct sites. For
#' incidence, the contribution is marginal occurrence minus the probability
#' of sharing the species. It is in richness units. The expected contribution
#' to richness-normalised incidence dissimilarity is a separate quantity
#' obtained by averaging the set-level ratios from [species_contributions()].
#' For relative abundance the contribution is mean site relative abundance
#' minus the expected minimum relative abundance. Its sum is expected
#' relative-abundance dissimilarity.
#'
#' @param x Non-negative numeric sites by species matrix.
#' @param k Set size, an integer from 1 to the number of sites.
#' @param type `"presence"` for unnormalised incidence or `"relative"` for
#'   relative abundance. Relative abundance requires positive totals at all sites.
#' @return A data frame with species names, marginal means, sharing
#'   probabilities, conditional mean minima (`NA` when no set can share the
#'   species), expected minima and expected dissimilarity contributions.
#' @export
expected_species_contributions <- function(x, k, type = c("presence", "relative")) {
  type <- match.arg(type)
  x <- .as_abundance(x)
  n <- nrow(x)
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k != floor(k) || k < 1 || k > n) {
    stop("`k` must be an integer between 1 and the number of sites.")
  }
  if (!ncol(x)) stop("`x` must contain at least one species.")
  occupied <- colSums(x > 0)
  P <- .species_zeta_presence(x, k)
  if (type == "presence") {
    marginal <- occupied / n
    conditional <- ifelse(occupied >= k, 1, NA_real_)
    minimum <- P
  } else {
    totals <- rowSums(x)
    if (any(totals == 0)) stop("Relative abundance requires positive totals at all sites.")
    q <- x / totals
    marginal <- colMeans(q)
    conditional <- vapply(seq_len(ncol(q)), function(s) {
      o <- occupied[s]
      if (o < k) return(NA_real_)
      v <- sort(q[x[, s] > 0, s])
      j <- seq_len(o - k + 1L)
      sum(v[j] * exp(lchoose(o - j, k - 1L) - lchoose(o, k)))
    }, numeric(1))
    minimum <- P * replace(conditional, is.na(conditional), 0)
  }
  data.frame(species = if (is.null(colnames(x))) seq_len(ncol(x)) else colnames(x),
             marginal_mean = marginal, sharing_probability = P,
             conditional_minimum = conditional, expected_minimum = minimum,
             dissimilarity = if (k == 1) 0 else pmax(marginal - minimum, 0),
             row.names = NULL)
}
