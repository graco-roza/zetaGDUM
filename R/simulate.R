#' Simulate a metacommunity along environmental gradients
#'
#' Sites are placed at random in a unit square; the first environmental gradient is spatially
#' smooth, the second is independent noise. Species have Gaussian niches in the gradients
#' (optima uniform, breadths log-normal) and a maximum local abundance drawn from a log-normal
#' distribution, so that locally dominant, regionally rare and widespread species all occur.
#' Counts are Poisson. The broader the niches, the more widespread the species. Optionally,
#' a site-level covariate `u` thins the individuals of a site at random (a non-directional
#' effect on richness and nestedness), and the second gradient can be made irrelevant.
#'
#' @param n_sites,n_species Numbers of sites and species.
#' @param breadth Median niche breadth (standard deviation in standardised units) of the species.
#' @param breadth_sd Log-scale standard deviation of niche breadths.
#' @param abundance_mean,abundance_sd Log-scale mean and sd of the maximum local abundance.
#' @param spatial Strength (0-1) of the spatial smoothness of the first gradient.
#' @param use_env Logical vector of length two: which gradients shape the niches.
#' @param thin Strength of the site-level thinning by `u` (0 = none); the retained fraction of
#'   individuals is `plogis(2 - thin * u)`.
#' @param species Optional data frame from a previous call's `truth`, to simulate new sites for
#'   the same species pool.
#' @param seed Optional random seed.
#' @return A list with `comm` (counts), `env` (data frame with `env1`, `env2` and, if `thin > 0`, `u`),
#'   `coords`, and `truth` (species optima, breadths and maximum abundances).
#' @export
#' @examples
#' sim <- simulate_metacommunity(60, 50, seed = 1)
#' dim(sim$comm)
simulate_metacommunity <- function(n_sites = 100, n_species = 100, breadth = 0.8, breadth_sd = 0.5,
                                   abundance_mean = 1.5, abundance_sd = 1, spatial = 0.7, use_env = c(TRUE, TRUE),
                                   thin = 0, species = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(species)) {
    opt <- cbind(stats::runif(n_species, -1.5, 1.5), stats::runif(n_species, -1.5, 1.5))
    br <- exp(stats::rnorm(n_species, log(breadth), breadth_sd))
    amax <- exp(stats::rnorm(n_species, abundance_mean, abundance_sd))
  } else {
    opt <- as.matrix(species[, c("opt1", "opt2")])
    br <- species$breadth
    amax <- species$amax
    n_species <- nrow(species)
  }
  coords <- cbind(x = stats::runif(n_sites), y = stats::runif(n_sites))
  smooth <- sin(2.5 * coords[, 1]) + cos(2 * coords[, 2])
  e1 <- spatial * as.numeric(scale(smooth)) + sqrt(1 - spatial^2) * stats::rnorm(n_sites)
  e2 <- stats::rnorm(n_sites)
  env <- data.frame(env1 = as.numeric(scale(e1)), env2 = as.numeric(scale(e2)))
  mu <- sapply(seq_len(n_species), function(s) {
    d <- (if (use_env[1]) (env$env1 - opt[s, 1])^2 else 0) + (if (use_env[2]) (env$env2 - opt[s, 2])^2 else 0)
    amax[s] * exp(-d / (2 * br[s]^2))
  })
  if (thin != 0) {
    u <- stats::rnorm(n_sites)
    keep <- stats::plogis(2 - thin * u)
    mu <- mu * keep
    env$u <- u
  }
  comm <- matrix(stats::rpois(length(mu), mu), n_sites, n_species)
  colnames(comm) <- paste0("sp", seq_len(n_species))
  rownames(comm) <- paste0("site", seq_len(n_sites))
  rownames(env) <- rownames(comm)
  list(comm = comm, env = env, coords = coords, truth = data.frame(species = colnames(comm), opt1 = opt[, 1], opt2 = opt[, 2], breadth = br, amax = amax))
}
