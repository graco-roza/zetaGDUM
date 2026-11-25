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
#'   the same species pool. If it has columns `cx`, `cy` and `rad`, each species is also restricted
#'   to the sites within distance `rad` of the point (`cx`, `cy`) (a geographic range window).
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
  has_range <- !is.null(species) && all(c("cx", "cy", "rad") %in% names(species))
  mu <- sapply(seq_len(n_species), function(s) {
    d <- rep(0, n_sites)
    if (use_env[1]) d <- d + (env$env1 - opt[s, 1])^2
    if (use_env[2]) d <- d + (env$env2 - opt[s, 2])^2
    m <- amax[s] * exp(-d / (2 * br[s]^2))
    if (has_range) m <- m * (sqrt((coords[, 1] - species$cx[s])^2 + (coords[, 2] - species$cy[s])^2) <= species$rad[s])
    m
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
  truth <- data.frame(species = colnames(comm), opt1 = opt[, 1], opt2 = opt[, 2], breadth = br, amax = amax)
  if (has_range) truth <- cbind(truth, species[, setdiff(names(species), names(truth)), drop = FALSE])
  list(comm = comm, env = env, coords = coords, truth = truth)
}

#' Species pool built from the forms of rarity of Rabinowitz
#'
#' Species differ in range size (a geographic window: large or small), habitat specificity (niche
#' breadth: broad or narrow) and local abundance (maximum abundance: large or small), giving eight
#' forms, one of them common (large range, broad niche, large abundance) and seven rare
#' (Rabinowitz 1981). A community can be dominated by one form (a share `p_dominant` of its
#' species) with the others sharing the rest equally, or be balanced across forms.
#'
#' @param n_species Number of species in the pool.
#' @param dominant Form that dominates: an integer 1-8 or a label from the returned `forms` table; `NULL` for a balanced pool.
#' @param p_dominant Share of the species belonging to the dominant form.
#' @param range_radius,breadths,abundances Length-two vectors (large, small) with the radius of the range window
#'   (unit square), the median niche breadth and the maximum abundance of the two levels of each axis.
#' @param seed Optional random seed.
#' @return A data frame of species (`opt1`, `opt2`, `breadth`, `amax`, `cx`, `cy`, `rad`, `form`) to be passed to
#'   [simulate_metacommunity()] as `species`; the table of forms is returned as the attribute `forms`.
#' @export
#' @examples
#' sp <- rarity_species_pool(100, dominant = 4, seed = 1)
#' table(sp$form)
rarity_species_pool <- function(n_species = 200, dominant = NULL, p_dominant = 0.7, range_radius = c(0.9, 0.3),
                                breadths = c(1.6, 0.45), abundances = c(20, 3), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  forms <- expand.grid(abundance = c("large", "small"), habitat = c("broad", "narrow"), range = c("large", "small"), stringsAsFactors = FALSE)[, c("range", "habitat", "abundance")]
  forms$label <- paste0(c(large = "Lr", small = "Sr")[forms$range], "-", c(broad = "Bn", narrow = "Nn")[forms$habitat], "-", c(large = "La", small = "Sa")[forms$abundance])
  if (is.character(dominant)) dominant <- match(dominant, forms$label)
  pr <- if (is.null(dominant)) rep(1 / 8, 8) else { p <- rep((1 - p_dominant) / 7, 8); p[dominant] <- p_dominant; p }
  f <- sample.int(8, n_species, replace = TRUE, prob = pr)
  lvl <- function(x, v) ifelse(x == "large", v[1], v[2])
  sp <- data.frame(opt1 = stats::runif(n_species, -1.5, 1.5), opt2 = stats::runif(n_species, -1.5, 1.5),
                   breadth = exp(stats::rnorm(n_species, log(ifelse(forms$habitat[f] == "broad", breadths[1], breadths[2])), 0.2)),
                   amax = exp(stats::rnorm(n_species, log(lvl(forms$abundance[f], abundances)), 0.3)),
                   cx = stats::runif(n_species), cy = stats::runif(n_species),
                   rad = lvl(forms$range[f], range_radius), form = forms$label[f], stringsAsFactors = FALSE)
  attr(sp, "forms") <- forms
  sp
}
