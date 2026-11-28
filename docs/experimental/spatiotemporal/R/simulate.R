st_scenarios <- c("spatial_filtering", "gains_only", "losses_only", "expansion", "contraction",
                  "abundance_only", "redistribution", "homogenisation", "differentiation",
                  "uncommon_effect", "widespread_effect", "persistent_effect", "interaction", "null")

simulate_st_gdum <- function(scenario = "null", n_sites = 30L, n_species = 48L,
                             times = 0:4, abundance = c("constant", "count", "gamma"),
                             detection = 1, variable_effort = FALSE, missing = 0, seed = NULL) {
  abundance <- match.arg(abundance)
  if (!length(scenario) || anyDuplicated(scenario) || !all(scenario %in% st_scenarios)) {
    stop("Supply distinct scenario names from st_scenarios.")
  }
  n <- .st_one(n_sites, 1e5, "n_sites")
  s <- .st_one(n_species, 1e5, "n_species")
  nt <- length(times)
  if (n < 6L || s < 8L || nt < 2L || !is.numeric(times) ||
      any(!is.finite(times)) || any(diff(times) <= 0)) stop("Require >= 6 sites, >= 8 species and increasing times.")
  for (v in list(detection, missing)) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0 || v > 1) {
      stop("detection and missing are probabilities.")
    }
  }
  if (!is.logical(variable_effort) || length(variable_effort) != 1L || is.na(variable_effort)) {
    stop("variable_effort must be TRUE or FALSE.")
  }
  if (!is.null(seed)) set.seed(seed)
  items <- rep(seq_len(2L), length.out = s)
  species <- data.frame(species = paste0("sp", seq_len(s)), widespread = items == 2L,
                        persistent = rep(c(FALSE, FALSE, TRUE, TRUE), length.out = s))
  species$occupancy_probability <- ifelse(species$widespread, 0.72, 0.22)
  ids <- list(site = paste0("site", seq_len(n)), time = as.character(times), species = species$species)
  x <- array(0, c(n, nt, s), ids)
  env <- array(0, c(n, nt, 1), list(ids$site, ids$time, "gradient"))
  env[, 1, 1] <- stats::rnorm(n)
  for (t in seq.int(2L, nt)) env[, t, 1] <- env[, t - 1L, 1] + stats::rnorm(n, 0, 0.75)
  baseline <- matrix(rep(species$occupancy_probability, each = n), n, s)
  if ("spatial_filtering" %in% scenario) {
    preferences <- rep(c(-1, -1, 1, 1), length.out = s)
    baseline <- stats::plogis(stats::qlogis(baseline) + 1.8 * outer(env[, 1, 1], preferences))
  }
  if ("differentiation" %in% scenario) baseline[] <- 0.8
  x[, 1, ] <- matrix(stats::runif(n * s), n, s) < baseline
  truth <- c(gain_change = 0, loss_change = 0, gain_interaction = 0, abundance_change = 0)
  if ("gains_only" %in% scenario) truth["gain_change"] <- 1
  if ("losses_only" %in% scenario) truth["loss_change"] <- -1
  if ("interaction" %in% scenario) truth["gain_interaction"] <- 1.25
  if ("abundance_only" %in% scenario) truth["abundance_change"] <- 0.65
  targeted <- list(uncommon_effect = !species$widespread, widespread_effect = species$widespread,
                   persistent_effect = species$persistent)
  effects <- intersect(names(targeted), scenario)
  targets <- rep(FALSE, s)
  if (length(effects)) targets <- Reduce(`|`, targeted[effects])
  beta_gain <- truth["gain_change"] + 1.2 * targets
  transition_lag <- diff(times)
  for (t in seq.int(2L, nt)) {
    initial <- env[, t - 1L, 1]
    delta <- env[, t, 1] - initial
    prev <- x[, t - 1L, ] > 0
    now <- prev
    if (any(c("gains_only", "expansion", "interaction", names(targeted), "null") %in% scenario)) {
      gp <- stats::plogis(-2.1 + outer(delta, beta_gain) + truth["gain_interaction"] * initial * delta)
      eligible <- if ("expansion" %in% scenario && !any(c("gains_only", "interaction", names(targeted), "null") %in% scenario)) {
        matrix(rep(colSums(prev) > 0, each = n), n, s)
      } else {
        matrix(TRUE, n, s)
      }
      now <- now | (eligible & matrix(stats::runif(n * s), n, s) < gp)
    }
    if (any(c("losses_only", "contraction", names(targeted), "interaction", "null") %in% scenario)) {
      lp <- stats::plogis(-2 + truth["loss_change"] * delta -
                          outer(rep(1, n), 2.2 * species$persistent))
      now <- now & !(prev & matrix(stats::runif(n * s), n, s) < lp)
      if ("contraction" %in% scenario && !any(c("losses_only", names(targeted), "interaction", "null") %in% scenario)) {
        for (j in which(colSums(prev) > 0 & colSums(now) == 0)) {
          occupied <- which(prev[, j])
          now[occupied[sample.int(length(occupied), 1L)], j] <- TRUE
        }
      }
    }
    if ("redistribution" %in% scenario) {
      for (j in seq_len(s)) {
        old <- which(now[, j])
        new <- which(!now[, j])
        swaps <- min(length(old), length(new), max(1L, floor(length(old) / 4)))
        if (swaps > 0L) {
          now[old[sample.int(length(old), swaps)], j] <- FALSE
          now[new[sample.int(length(new), swaps)], j] <- TRUE
        }
      }
    }
    if ("homogenisation" %in% scenario) {
      core <- which.max(colSums(now))
      for (j in seq_len(s)) {
        goal <- if (j == core) TRUE else colSums(prev)[j] > n / 2
        update <- stats::runif(n) < 0.75
        now[update, j] <- goal
      }
      now[, core] <- TRUE
    }
    if ("differentiation" %in% scenario) {
      target <- outer(env[, t, 1] > stats::median(env[, t, 1]), species$widespread, `==`)
      if (t == nt) now <- target else now <- now & target | now & matrix(stats::runif(n * s) > 0.4, n, s)
    }
    x[, t, ] <- now
  }
  positive <- switch(abundance, constant = array(8, dim(x)), count = array(1 + stats::rpois(length(x), 7), dim(x)),
                     gamma = array(stats::rgamma(length(x), 8, rate = 1), dim(x)))
  if ("abundance_only" %in% scenario) {
    for (t in seq.int(2L, nt)) {
      multiplier <- exp(truth["abundance_change"] * (env[, t, 1] - env[, t - 1L, 1]))
      positive[, t, ] <- positive[, t, ] * multiplier
    }
  }
  if (abundance == "count") positive <- pmax(1, round(positive))
  x <- x * positive
  effort <- matrix(1, n, nt, dimnames = ids[1:2])
  if (variable_effort) effort[] <- exp(-0.5 + 0.65 * env[, , 1])
  p_detect <- 1 - (1 - detection)^effort
  detected <- array(stats::runif(length(x)), dim(x)) < as.vector(p_detect)
  observed <- x * detected
  observed_sites <- matrix(stats::runif(n * nt) >= missing, n, nt, dimnames = ids[1:2])
  observed[!rep(as.vector(observed_sites), s)] <- NA_real_
  list(latent = st_panel(x, times), observed = observed, env = env, effort = effort,
       p_detect = p_detect, observed_sites = observed_sites, species = species,
       truth = truth, gain_slope_by_species = setNames(beta_gain, species$species),
       scenario = scenario, lags = transition_lag)
}
