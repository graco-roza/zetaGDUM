source("docs/phase4/common.R")
# Numerical budget verification when exhaustive evaluation is intentionally disabled.
n <- 24L
s <- simulate_metacommunity(n, 50, abundance_mean = 2, seed = 8834)
s$comm <- safe_community(s$comm)
env <- s$env["env1"]
ans <- ugdm_order_limits(s$comm, env, c(3L, 5L), "env1", c(spread = .10, position = .35),
  min_response_sd = .005, sampler = "site", pilot_budgets = c(700L, 2800L, 8400L),
  reference_M = 22000L, max_budget = 60000L, nrep = 16L, exact_limit = 0,
  min_delete_ess = 100, epsilon = .30, ecological_noise = .35, shapley = TRUE,
  component = "relative", seed = 10496)
phase4_write(cbind(ans$profile, highest_supported_order = ans$highest_supported_order),
             "verification", "profile")
phase4_write(ans$precision, "verification", "precision")
phase4_write(ans$pilots, "verification", "pilot")
