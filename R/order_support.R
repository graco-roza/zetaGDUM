#' Exact ecological support by zeta order
#'
#' Uses species occupancy and pairwise co-occupancy to compute the first
#' two moments of shared richness under uniform site-set sampling. Relative
#' abundance and the effective species counts use exact species minima.
#'
#' @param x Numeric sites by species matrix.
#' @param orders Distinct site-set sizes between one and the number of sites.
#' @param abundance Compute relative-abundance diagnostics. Positive finite
#'   totals at every site are required when enabled.
#' @return An order table with moments, analytical sharing bounds,
#'   incidence and abundance retention, and inverse concentration of
#'   expected shared species contributions. Effective species counts are
#'   missing when there is no expected shared contribution.
#' @export
zeta_order_support <- function(x, orders = seq_len(nrow(x)), abundance = TRUE) {
  x <- .conditional_community(x, relative = FALSE)
  if (isTRUE(abundance) && any(!is.finite(rowSums(x)) | rowSums(x) <= 0)) {
    stop("Relative abundance requires positive finite site totals.")
  }
  n <- nrow(x)
  if (!length(orders) || anyDuplicated(orders)) stop("Provide distinct orders.")
  orders <- vapply(orders, .zeta_order_integer, integer(1), n = n)
  P <- x > 0
  o <- colSums(P)
  common <- crossprod(P)
  pairs <- tabulate(common[lower.tri(common)] + 1L, n + 1L)
  moment <- function(k) {
    center <- lchoose(n, k)
    prob <- exp(lchoose(o, k) - center)
    mu <- sum(prob)
    twice_joint <- 2 * sum(pairs * exp(lchoose(0:n, k) - center))
    second <- mu + twice_joint
    m <- sum(o >= k)
    positive <- mu > 0
    lower <- if (positive) max(max(prob), mu * (mu / second), mu - twice_joint / 2, mu / m) else 0
    upper <- if (positive) min(1, mu, mu - (second - mu) / m) else 0
    if (any(o == n)) lower <- upper <- 1
    c(expected_richness = mu, richness_second_moment = second,
      richness_variance = max(0, second - mu^2),
      sharing_lower = min(max(lower, 0), 1),
      sharing_upper = min(max(upper, lower), 1),
      effective_incidence_species = if (positive) 1 / sum((prob / mu)^2) else NA_real_)
  }
  ab <- function(k) {
    if (!isTRUE(abundance)) return(c(NA_real_, NA_real_))
    w <- expected_species_contributions(x, k, "relative")$expected_minimum
    total <- sum(w)
    c(total, if (total > 0) 1 / sum((w / total)^2) else NA_real_)
  }
  # Include predecessor orders even when the requested grid has gaps.
  grid <- sort(unique(c(orders, pmax(orders - 1L, 1L))))
  rich <- t(vapply(grid, moment, numeric(6)))
  rel <- t(vapply(grid, ab, numeric(2)))
  pick <- match(orders, grid)
  prev <- match(pmax(orders - 1L, 1L), grid)
  retention <- function(y) {
    value <- ifelse(y[prev] > 0, y[pick] / y[prev], NA_real_)
    value[orders == 1L] <- NA_real_
    value
  }
  data.frame(order = orders, ratio = orders / n, rich[pick, , drop = FALSE],
    abundance_zeta = rel[pick, 1L], effective_abundance_species = rel[pick, 2L],
    incidence_retention = retention(rich[, 1L]), abundance_retention = retention(rel[, 1L]),
    row.names = NULL)
}

.iid_uniform_sets <- function(n, k, M) {
  t(vapply(seq_len(M), function(i) sort.int(sample.int(n, k)), integer(k)))
}

#' Stratify uniform site-set integration by species sharing
#'
#' Samples uniformly within each sharing stratum. Exact stratum mass is
#' used for small populations and analytically determined cases. Otherwise
#' a separate IID uniform pilot estimates the mass. In that case weighted
#' set averages are unbiased before conditioning on the mass pilot, and
#' final variance assessments must repeat the pilot. Rejection can be
#' expensive for imbalanced strata.
#'
#' @param x,k,M,seed As in [sample_sets_species()].
#' @param positive_fraction Allocation fraction for sets sharing a species.
#' @param mass_draws Independent uniform draws for estimating stratum mass.
#'   The default grows with the integration budget.
#' @param exact_limit Maximum number of sets to enumerate for exact mass.
#' @param max_proposals Maximum candidate draws per stratum. If exhausted,
#'   return a fresh IID uniform draw with unit sampling weights.
#' @return A list with sampled indices, sampling weights, sharing counts,
#'   estimated positive mass and its pilot standard error, mass source,
#'   proposal counts and fallback reason when applicable.
#' @export
sample_sets_stratified <- function(x, k, M = 1e4, positive_fraction = .5, mass_draws = max(10000L, M),
                                   exact_limit = 50000L, max_proposals = 20 * M, seed = NULL) {
  x <- .conditional_community(x)
  n <- nrow(x)
  P <- (x > 0) * 1L
  k <- .zeta_order_integer(k, n, 2L)
  M <- .set_design_integer(M, "M", 2L)
  mass_draws <- .set_design_integer(mass_draws, "mass_draws", 2L)
  exact_limit <- .set_design_integer(exact_limit, "exact_limit", 0L)
  max_proposals <- .set_design_integer(max_proposals, "max_proposals")
  if (!is.numeric(positive_fraction) || length(positive_fraction) != 1L ||
      !is.finite(positive_fraction) || positive_fraction <= 0 || positive_fraction >= 1) {
    stop("Supply a positive allocation fraction strictly below one.")
  }
  if (!is.null(seed)) set.seed(seed)
  uniform <- function(m) .iid_uniform_sets(n, k, m)
  sharing <- function(idx) shared_species_cpp(P, idx)
  plan <- .sampling_plan(P, k, NULL)
  proposals <- c(zero = 0L, positive = 0L)
  mass <- se <- NA_real_
  source <- "independent_pilot"
  plain <- function(reason = NULL) {
    idx <- uniform(M)
    list(idx = idx, sampling_weights = rep(1, M), shared = sharing(idx),
      positive_mass = mass, mass_se = se, mass_source = source, proposals = proposals, fallback = reason)
  }
  all_idx <- if (lchoose(n, k) <= log(exact_limit) + 1e-12) t(utils::combn(n, k))
  if (!is.finite(plan$log_expected_shared) || any(colSums(P) == n) || !is.null(all_idx)) {
    source <- "exact"
    all_z <- if (!is.null(all_idx)) sharing(all_idx)
    mass <- if (!is.finite(plan$log_expected_shared)) 0 else
      if (any(colSums(P) == n)) 1 else mean(all_z > 0)
    se <- 0
    if (mass %in% c(0, 1)) return(plain())
  }
  target <- function(m) {
    s <- sample.int(ncol(P), m, TRUE, plan$species_probability)
    idx <- matrix(0L, m, k)
    by_s <- split(seq_len(m), s)
    for (key in names(by_s)) {
      sites <- plan$occupied[[as.integer(key)]]
      idx[by_s[[key]], ] <- matrix(sites[.iid_uniform_sets(length(sites), k, length(by_s[[key]]))],
                                  length(by_s[[key]]), k)
    }
    idx
  }
  rejection <- function(m, positive) {
    count <- tried <- 0L
    idx <- matrix(0L, m, k)
    z <- integer(m)
    while (count < m && tried < max_proposals) {
      batch <- min(max_proposals - tried, max(64L, m - count))
      fresh <- if (positive) target(batch) else uniform(batch)
      shares <- sharing(fresh)
      accept <- if (positive) stats::runif(batch) * shares < 1 else shares == 0
      rows <- utils::head(which(accept), m - count)
      if (length(rows)) {
        take <- count + seq_along(rows)
        idx[take, ] <- fresh[rows, , drop = FALSE]
        z[take] <- shares[rows]
      }
      count <- count + length(rows)
      tried <- tried + batch
    }
    list(idx = idx, z = z, tried = tried, full = count == m)
  }
  np <- min(M - 1L, max(1L, round(M * positive_fraction)))
  nn <- M - np
  if (!is.null(all_idx)) {
    zero <- which(all_z == 0)
    pos <- which(all_z > 0)
    rows <- c(zero[sample.int(length(zero), nn, TRUE)], pos[sample.int(length(pos), np, TRUE)])
    idx <- all_idx[rows, , drop = FALSE]
    z <- all_z[rows]
  } else {
    # The pilot is independent of both conditional draws and fallback selection.
    if (source == "independent_pilot") {
      mass <- mean(sharing(uniform(mass_draws)) > 0)
      se <- sqrt(mass * (1 - mass) / (mass_draws - 1))
    }
    zero <- rejection(nn, FALSE)
    pos <- if (zero$full) rejection(np, TRUE)
    proposals <- c(zero = zero$tried, positive = if (is.null(pos)) 0L else pos$tried)
    if (!zero$full || !pos$full) return(plain("proposal limit before completing quotas"))
    idx <- rbind(zero$idx, pos$idx)
    z <- c(zero$z, pos$z)
  }
  list(idx = idx, sampling_weights = rep(c(M * (1 - mass) / nn, M * mass / np), c(nn, np)),
       shared = z, positive_mass = mass, mass_se = se, mass_source = source,
       proposals = proposals, fallback = NULL)
}
