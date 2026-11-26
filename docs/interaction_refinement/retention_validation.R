library(zetaGDUM)
set.seed(30157)
root <- "docs/interaction_refinement/results"
dir.create(root, recursive = TRUE, showWarnings = FALSE)

exhaustive <- function(x, type) {
  q <- if (type == "relative") x / rowSums(x) else 1 * (x > 0)
  n <- nrow(q)
  enumerated <- lapply(seq_len(n), function(k) {
    sets <- utils::combn(n, k, simplify = FALSE)
    common <- vapply(sets, function(members) {
      sum(apply(q[members, , drop = FALSE], 2L, min))
    }, numeric(1))
    vapply(seq_len(n), function(i) {
      include <- vapply(sets, function(members) i %in% members, logical(1))
      mean(common[include]) / sum(q[i, ])
    }, numeric(1))
  })
  do.call(cbind, enumerated)
}

validation <- lapply(seq_len(40L), function(j) {
  n <- sample(3:8, 1)
  x <- matrix(sample(0:3, n * 7L, replace = TRUE), n)
  x[cbind(seq_len(n), sample.int(ncol(x), n, replace = TRUE))] <- 1
  do.call(rbind, lapply(c("relative", "occurrence"), function(type) {
    exact <- exhaustive(x, type)
    local <- zeta_local_retention(x, type)
    trunc <- zeta_local_retention(x, type, 3L)
    data.frame(replicate = j, type, n,
               profile_max_error = max(abs(exact - local$profile)),
               full_T_max_error = max(abs(rowSums(exact) - local$summary$T_full)),
               loss_mass_max_error = max(abs(rowSums(local$loss) - 1)),
               truncated_T_max_error = max(abs(rowSums(exact[, 1:3]) - trunc$summary$T)))
  }))
})
utils::write.csv(do.call(rbind, validation), file.path(root, "retention_exact.csv"), row.names = FALSE)

site_correlation <- function(x, y) {
  if (!(stats::sd(x) > 1e-12 && stats::sd(y) > 1e-12)) return(NA_real_)
  stats::cor(x, y, method = "spearman")
}

n <- 70L
p <- 30L
k <- 8L
simulate <- function(n, p, occupancy_saturated) {
  g <- as.numeric(stats::arima.sim(list(ar = .7), n))
  g <- as.numeric(scale(g + seq(-.7, .7, length.out = n)))
  env <- data.frame(gradient = g)
  niche_centres <- sort(stats::runif(p, -2, 1.4))
  breadth <- stats::runif(p, .3, .75)
  spine <- outer(g, niche_centres, "-")
  lambda <- 6 * exp(-.5 * (sweep(spine, 2L, breadth, "/"))^2)
  x <- matrix(stats::rpois(n * p, lambda), n, p)
  if (occupancy_saturated) x <- 1 + x
  x[, 1L] <- x[, 1L] + 1
  rownames(x) <- paste0("site", seq_len(n))
  list(x = x, env = env)
}

ecological <- lapply(c("niche_gradient", "occupancy_saturated"), function(scenario) {
  do.call(rbind, lapply(seq_len(12L), function(j) {
    d <- simulate(n, p, scenario == "occupancy_saturated")
    idx <- sample_sets(n, k, M = 12000L)
    design <- ugdm_design(d$env, idx)
    y <- ugdm_response(d$x, idx, "relative", working_weights = "constant")
    fit <- ugdm_fit(design, y$y, y$trials, distance = FALSE)
    # Separate prediction paths for site summaries, retaining the training scale.
    reference <- ugdm_design(d$env, sample_sets(n, k, M = 50000L), ref = design)
    components <- ugdm_local_distinctiveness(d$x, fit, reference, "relative")
    relative <- zeta_local_retention(d$x, "relative")
    occurrence <- zeta_local_retention(d$x, "occurrence")
    R <- relative$summary$T
    O <- occurrence$summary$T
    loss <- unname(relative$loss[, n])
    data.frame(scenario, replicate = j, order = k, converged = fit$converged,
      occurrence_median = stats::median(O), occurrence_sd = stats::sd(O),
      relative_median = stats::median(R), relative_sd = stats::sd(R),
      censored_mass = mean(loss), observed = site_correlation(R, components$observed),
      fitted = site_correlation(R, components$predicted),
      spread = site_correlation(R, components$spread),
      position = site_correlation(R, components$position),
      residual = site_correlation(R, components$residual))
  }))
})
utils::write.csv(do.call(rbind, ecological), file.path(root, "retention_environment.csv"),
                 row.names = FALSE)

first <- simulate(n, p, FALSE)
second <- simulate(n, p, FALSE)
third <- simulate(n, p, FALSE)
occasion <- list(early1 = first$x[-sample.int(n, 6L), ],
                 early2 = second$x[-sample.int(n, 8L), ],
                 late1 = third$x[-sample.int(n, 14L), ])
seasons <- c("early", "early", "late")
temporal <- lapply(c("common_sites", "truncated"), function(comparison) {
  fits <- zeta_temporal_retention(occasion, seasons, comparison, "relative", max_order = k)
  values <- fits$seasonal
  values$comparison <- comparison
  values
})
utils::write.csv(do.call(rbind, temporal), file.path(root, "retention_temporal.csv"),
                 row.names = FALSE)
