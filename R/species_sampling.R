#' Sample site sets with a species-informed mixture
#'
#' Mixes conventional site-set draws with draws guaranteed to share a species.
#' The base distribution is uniform, or proportional to the product of site
#' weights if `omega` is supplied. The targeted branch selects a species in
#' proportion to its probability of being shared under the base distribution,
#' then draws from that distribution conditional on sharing the selected species.
#' Draws are independent, and the same combination can occur more than once.
#'
#' @param x Numeric sites by species matrix. Positive entries indicate presence.
#' @param k Set size, an integer from 2 to the number of eligible sites.
#' @param M Positive integer number of draws.
#' @param alpha Probability of using the targeted branch, in `[0, 1)`.
#' @param omega Optional finite non-negative site weights. Zero-weight sites
#'   are ineligible. By default all sites have equal weight.
#' @param seed Optional random seed.
#' @return A list with `idx` (site indices), `sampling_weights` (base-to-mixture
#'   probability ratios), `shared` (species shared per set), `targeted` (branch
#'   indicators), `expected_shared` (base expectation), `species_probability`
#'   (targeted species probabilities), and `alpha` (effective mixture probability,
#'   zero when no eligible combination can share a species).
#' @export
sample_sets_species <- function(x, k, M = 1e4, alpha = 0.5, omega = NULL, seed = NULL) {
  x <- .as_abundance(x)
  .sample_sets_mixture((x > 0) * 1L, k, M, alpha, omega, seed)
}

.check_alpha <- function(alpha) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha < 0 || alpha >= 1) {
    stop("`alpha` must be a number in [0, 1).")
  }
}

.sampling_plan <- function(P, k, omega) {
  eligible <- if (is.null(omega)) rep(TRUE, nrow(P)) else omega > 0
  occupied <- lapply(seq_len(ncol(P)), function(s) which(eligible & P[, s] > 0L))
  sizes <- lengths(occupied)
  if (is.null(omega)) {
    log_mass <- lchoose(sizes, k)
    log_base <- lchoose(sum(eligible), k)
  } else {
    log_omega <- log(omega)
    log_omega <- log_omega - max(log_omega)
    log_mass <- vapply(occupied, function(sites) {
      if (length(sites) < k) -Inf else log_esym_cpp(log_omega[sites], k)
    }, numeric(1))
    log_base <- log_esym_cpp(log_omega, k)
  }
  top <- if (length(log_mass)) max(log_mass) else -Inf
  mass <- if (is.finite(top)) exp(log_mass - top) else numeric(length(occupied))
  total <- sum(mass)
  list(occupied = occupied, species_probability = if (total > 0) mass / total else mass,
       log_expected_shared = if (total > 0) (top - log_base) + log(total) else -Inf)
}

.mixture_weights <- function(shared, alpha, log_expected_shared) {
  if (alpha == 0 || log_expected_shared == -Inf) return(rep(1, length(shared)))
  log_target <- log(alpha) + log(shared) - log_expected_shared
  log_base <- log1p(-alpha)
  scale <- pmax(log_base, log_target)
  exp(-scale) / (exp(log_base - scale) + exp(log_target - scale))
}

.sample_sets_mixture <- function(P, k, M, alpha, omega = NULL, seed = NULL) {
  .check_alpha(alpha)
  n <- nrow(P)
  if (!is.null(omega) && (!is.numeric(omega) || length(omega) != n ||
      any(!is.finite(omega) | omega < 0))) {
    stop("`omega` must contain one finite, non-negative weight per site.")
  }
  n_eligible <- if (is.null(omega)) n else sum(omega > 0)
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k != floor(k) || k < 2 || k > n_eligible) {
    stop("`k` must be an integer between 2 and the number of eligible sites.")
  }
  if (length(M) != 1L || !is.numeric(M) || !is.finite(M) || M != floor(M) ||
      M < 1 || M > .Machine$integer.max) {
    stop("`M` must be a positive integer representable in R.")
  }
  if (!is.null(seed)) set.seed(seed)
  plan <- .sampling_plan(P, k, omega)
  if (!is.finite(plan$log_expected_shared)) alpha <- 0
  draw <- function(sites, m) {
    if (length(sites) == k) return(matrix(rep(sites, each = m), m, k))
    local <- if (is.null(omega)) {
      t(vapply(seq_len(m), function(i) sort.int(sample.int(length(sites), k)), integer(k)))
    } else {
      sample_sets_weighted(omega[sites], k, m)
    }
    matrix(sites[local], m, k)
  }
  targeted <- if (alpha > 0) stats::runif(M) < alpha else rep(FALSE, M)
  idx <- matrix(0L, M, k)
  if (any(!targeted)) idx[!targeted, ] <- draw(seq_len(n), sum(!targeted))
  if (any(targeted)) {
    species <- sample.int(ncol(P), sum(targeted), replace = TRUE, prob = plan$species_probability)
    by_species <- split(which(targeted), species)
    for (s in names(by_species)) {
      rows <- by_species[[s]]
      idx[rows, ] <- draw(plan$occupied[[as.integer(s)]], length(rows))
    }
  }
  shared <- shared_species_cpp(P, idx)
  list(idx = idx, sampling_weights = .mixture_weights(shared, alpha, plan$log_expected_shared),
       shared = shared, targeted = targeted, expected_shared = exp(plan$log_expected_shared),
       species_probability = plan$species_probability, alpha = alpha)
}

.model_sets <- function(P, k, M, sampling, alpha, omega = NULL, seed = NULL) {
  n <- nrow(P)
  exhaustive <- k == 2 && choose(n, 2) <= M
  if (sampling == "species" && alpha > 0 && !exhaustive) {
    return(.sample_sets_mixture(P, k, M, alpha, omega, seed))
  }
  idx <- if (is.null(omega) || exhaustive) sample_sets(n, k, M, seed) else
    sample_sets_weighted(omega, k, M, seed)
  list(idx = idx, sampling_weights = rep(1, nrow(idx)),
       shared = shared_species_cpp(P, idx), alpha = 0)
}

.sampling_diagnostics <- function(sets, keep, weights) {
  n <- sum(keep)
  c(n_sets = n, n_shared = sum(sets$shared[keep] > 0),
    ess = if (any(weights > 0)) sum(weights / max(weights))^2 / sum((weights / max(weights))^2) else 0,
    alpha = sets$alpha)
}

.product_sampling_weights <- function(omega, idx, sampling_weights) {
  log_w <- rowSums(matrix(log(omega)[idx], nrow(idx))) + log(sampling_weights)
  exp(log_w - max(log_w))
}
