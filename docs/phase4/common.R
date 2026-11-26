source("docs/manuscript_validation/common.R")

phase4_write <- function(x, study, key) {
  sharded <- study %in% c("profile", "precision", "pilot", "bootstrap", "comparison",
                          "sensitivity", "sensitivity_precision", "composition")
  directory <- if (sharded) "docs/phase4/runs" else "docs/phase4/results"
  write_study(x, directory, paste(study, key, sep = "_"))
}

phase4_input <- function(x, env, dist = NULL) {
  p <- list(comm = x, P = 1L * (x > 0), env = env, dm = dist,
            dist_range = if (!is.null(dist)) internal(".pair_distance_range")(dist))
  n <- nrow(x)
  p$basis <- ugdm_design(env, cbind(seq_len(n - 1L), seq.int(2L, n)),
                         dist_matrix = dist, dist_range = p$dist_range)
  p
}

support_case <- function(n) {
  group <- ceiling(seq_len(n) / 4)
  x <- cbind(core = 1, 1e7 * outer(group, seq_len(max(group)), "=="))
  list(comm = x, env = data.frame(env1 = seq_len(n) / n))
}

safe_community <- function(x) {
  empty <- which(rowSums(x) == 0)
  if (!length(empty)) return(x)
  species <- matrix(0, nrow(x), length(empty))
  species[cbind(empty, seq_along(empty))] <- 1
  cbind(x, species)
}

phase4_profile <- function(x, env, orders, component, reference_M = 12000L,
                           min_response_sd = .005, seed = 191) {
  ugdm_order_limits(x, env, orders, names(env)[1L], c(spread = .10, position = .35),
    min_response_sd = min_response_sd, sampler = "conventional", pilot_budgets = c(1000, 4000),
    reference_M = reference_M, max_budget = 50000, nrep = 10, exact_limit = 50000,
    min_delete_ess = 100, epsilon = .25, ecological_noise = .35, require_change = TRUE,
    shapley = TRUE, component = component, working_weights = "constant", seed = seed)
}

paired_bootstrap <- function(x, env, orders, component, M, B, seed) {
  set.seed(seed)
  n <- nrow(x)
  p <- phase4_input(x, env)
  effect <- internal(".limit_effect_spec")(env, p$basis, names(env)[1L])
  omega <- site_weights(n, B)
  fits <- lapply(orders, function(k) t(vapply(seq_len(B), function(j) {
    idx <- sample_sets_weighted(omega[, j], k, M)
    sets <- list(idx = idx, sampling_weights = rep(1, nrow(idx)),
                 shared = integer(nrow(idx)), exact = FALSE)
    f <- internal(".limit_fit")(p, sets, component, "constant", FALSE, effect, TRUE)
    if (f$regular) unname(f$effects) else rep(NA_real_, 4)
  }, numeric(4))))
  stats <- function(draws, k, previous) {
    s <- vapply(seq_len(ncol(draws)), function(j) {
      d <- draws[, j]
      if (mean(is.finite(d)) < .95) return(c(sum(is.finite(d)), NA_real_, NA_real_, NA_real_))
      c(sum(is.finite(d)), sd(d, na.rm = TRUE), quantile(d, c(.025, .975), na.rm = TRUE, names = FALSE))
    }, numeric(4))
    data.frame(order = k, preceding = previous, effect = c("spread", "position", "shapley_spread", "shapley_position"),
               valid = s[1, ], bootstrap_se = s[2, ], lo = s[3, ], hi = s[4, ], B, M)
  }
  do.call(rbind, lapply(seq_along(orders), function(j) {
    level <- stats(fits[[j]], orders[j], NA_integer_)
    if (j == 1L) return(level)
    rbind(level, stats(fits[[j]] - fits[[j - 1L]], orders[j], orders[j - 1L]))
  }))
}
