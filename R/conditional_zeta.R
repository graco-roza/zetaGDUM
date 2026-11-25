.zeta_order_integer <- function(k, n, minimum = 1L) {
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) ||
      k != floor(k) || k < minimum || k > n) {
    stop("`k` must be a single integer between ", minimum, " and the number of sites.")
  }
  as.integer(k)
}

.conditional_community <- function(x, relative = FALSE) {
  x <- .as_abundance(x)
  if (!nrow(x) || !ncol(x) || any(!is.finite(x))) {
    stop("Supply a finite matrix with at least one site and one species.")
  }
  if (relative) {
    totals <- rowSums(x)
    if (any(!is.finite(totals) | totals <= 0)) {
      stop("Relative abundance requires positive finite site totals.")
    }
    x <- x / totals
  }
  x
}

#' Exact conditional shared composition
#'
#' Expectations over uniform sets of distinct sites, conditional on a
#' specified set of sites being included. The relative-abundance calculation
#' adds threshold increments of within-site species proportions. Site and
#' pair expectations use cumulative threshold lookups.
#'
#' @param x Numeric sites by species matrix.
#' @param k Number of sites in a set.
#' @param type Incidence shared richness or shared relative abundance.
#' @param include Optional integer vector of forced site indices, or a list
#'   of such vectors. A zero-length vector refers to the unconditional mean.
#' @param pairs Return the symmetric matrix for conditioning on pairs of
#'   distinct sites. Requires `k >= 2` and uses quadratic output memory.
#' @return A list with `mean`, `site`, optional `pair` (diagonal `NA`),
#'   and `forced` (one entry per requested subset).
#' @export
conditional_zeta <- function(x, k, type = c("incidence", "relative"), include = NULL, pairs = FALSE) {
  type <- match.arg(type)
  x <- .conditional_community(x, type == "relative")
  n <- nrow(x)
  k <- .zeta_order_integer(k, n, if (isTRUE(pairs)) 2L else 1L)
  if (is.null(include)) include <- list()
  if (!is.list(include)) include <- list(include)
  a <- lapply(include, function(i) {
    if (length(i) && (!is.numeric(i) || anyNA(i) || anyDuplicated(i) ||
                      any(i < 1 | i > n | i != floor(i)) || length(i) > k)) {
      stop("Each forced subset must contain at most `k` distinct valid site indices.")
    }
    as.integer(i)
  })
  out <- conditional_zeta_cpp(x, k, type == "relative", isTRUE(pairs), a)
  names(out$site) <- rownames(x)
  names(out$forced) <- names(include)
  if (isTRUE(pairs)) dimnames(out$pair) <- list(rownames(x), rownames(x)) else out$pair <- NULL
  out
}

.mean_richness_given <- function(richness, k, sites) {
  a <- length(sites)
  n <- length(richness)
  if (a == n) return(mean(richness))
  forced <- sum(richness[sites])
  (forced + (k - a) / (n - a) * (sum(richness) - forced)) / k
}

#' Conditional expected dissimilarity
#'
#' Uses exact relative-abundance or unnormalised incidence dissimilarity.
#' The latter is mean site richness minus shared richness.
#'
#' @inheritParams conditional_zeta
#' @return The same structure as [conditional_zeta()], with expected
#'   dissimilarities in place of expected shared composition and a
#'   `shared` element containing the original zeta expectations.
#' @export
conditional_zeta_dissimilarity <- function(x, k, type = c("incidence", "relative"),
                                           include = NULL, pairs = FALSE) {
  type <- match.arg(type)
  x <- .conditional_community(x)
  shared <- conditional_zeta(x, k, type, include, pairs)
  out <- shared
  if (type == "relative") {
    for (term in names(shared)) out[[term]] <- pmax(1 - shared[[term]], 0)
  } else {
    r <- rowSums(x > 0)
    n <- length(r)
    out$mean <- mean(r) - shared$mean
    out$site <- (r + (k - 1) * (sum(r) - r) / max(1, n - 1)) / k - shared$site
    if (isTRUE(pairs)) {
      sums <- outer(r, r, "+")
      out$pair <- (sums + (k - 2) * (sum(r) - sums) / max(1, n - 2)) / k - shared$pair
    }
    if (length(shared$forced)) {
      if (!is.list(include)) include <- list(include)
      out$forced <- vapply(include, function(i) .mean_richness_given(r, k, i), numeric(1L)) - shared$forced
    }
    for (j in names(out)) out[[j]] <- pmax(out[[j]], 0)
  }
  out$shared <- shared
  out
}

.conditional_sample_mean <- function(idx, values, w, n) {
  values <- as.matrix(values)
  if (!any(w > 0)) return(matrix(NA_real_, n, ncol(values)))
  w <- w / max(w)
  sums <- .inclusion_sums(idx, cbind(w, w * values), n)
  out <- sums[, -1L, drop = FALSE] / pmax(sums[, 1L], .Machine$double.xmin)
  out[sums[, 1L] == 0, ] <- NA_real_
  out
}

.lczd_summary <- function(site, k, shared, dissimilarity, overall) {
  total <- sum(dissimilarity)
  tol <- 64 * .Machine$double.eps * length(site)
  if (is.finite(total) && total <= tol) {
    dissimilarity[abs(dissimilarity) <= tol] <- 0
    total <- sum(dissimilarity)
  }
  lc <- if (is.finite(total) && total > 0) dissimilarity / total else rep(NA_real_, length(site))
  data.frame(site = site, order = k, expected_shared = shared,
             expected_dissimilarity = dissimilarity, LCZD = lc,
             relative_distinctiveness = length(site) * lc,
             centred_distinctiveness = dissimilarity - overall,
             rank = rank(-dissimilarity, na.last = "keep"), row.names = NULL)
}

#' Local contributions to zeta dissimilarity
#'
#' The exact versions are relative-abundance dissimilarity and
#' unnormalised incidence dissimilarity. Richness-normalised incidence
#' averages set-level ratios within each site using a supplied sample.
#' Its centred values use the sample mean from the same proposal.
#'
#' @param x Numeric sites by species matrix.
#' @param k One or more distinct orders. The sampled version uses one order.
#' @param type Response type.
#' @param idx Optional matrix of sampled sets for richness-normalised incidence.
#' @param sampling_weights Optional target-to-proposal ratios for `idx`.
#' @return One row per site and order, including the descending rank
#'   (averaged ties). LCZD is undefined for a zero total or missing conditional
#'   estimates. In the exact version each site's expected dissimilarity is the
#'   mean over all sets containing it.
#' @export
zeta_lczd <- function(x, k, type = c("relative", "incidence", "normalised"),
                      idx = NULL, sampling_weights = NULL) {
  type <- match.arg(type)
  x <- .conditional_community(x)
  n <- nrow(x)
  if (!length(k) || anyDuplicated(k)) stop("Provide distinct orders.")
  site <- if (is.null(rownames(x))) seq_len(n) else rownames(x)
  if (type != "normalised") {
    if (!is.null(idx) || !is.null(sampling_weights)) stop("Sample arguments are used only with normalised incidence.")
    return(do.call(rbind, lapply(k, function(order) {
      z <- conditional_zeta_dissimilarity(x, order, type)
      .lczd_summary(site, order, z$shared$site, z$site, z$mean)
    })))
  }
  k <- .zeta_order_integer(k, n, 2L)
  idx <- .species_indices(idx, n)
  if (ncol(idx) != k) stop("The number of columns in `idx` must equal `k`.")
  w <- .model_weights(nrow(idx), sampling_weights = sampling_weights)
  z <- zeta_sets(x, idx, abundance = FALSE)
  ok <- is.finite(z$total) & w > 0
  if (!any(ok)) stop("No sampled sets have defined normalised incidence dissimilarity.")
  vals <- .conditional_sample_mean(idx[ok, , drop = FALSE],
                                  cbind(z$shared[ok], 1 - z$total[ok]), w[ok], n)
  .lczd_summary(site, k, vals[, 1L], vals[, 2L],
                stats::weighted.mean(1 - z$total[ok], w[ok]))
}
