test_that("sharing moments and bounds agree with exhaustive populations", {
  x <- rbind(c(0,2,0,0,0), c(1,2,3,0,0), c(1,0,3,0,0),
             c(0,0,4,1,0), c(0,2,1,0,1), c(0,0,3,0,0))
  for (comm in list(x, cbind(x, 1), diag(6), matrix(1, 6, 2))) {
    n <- nrow(comm)
    ans <- zeta_order_support(comm)
    for (k in seq_len(n)) {
      idx <- t(combn(n, k))
      contributions <- t(apply(idx, 1, function(s) colSums(comm[s, , drop = FALSE] > 0) == k))
      z <- rowSums(contributions)
      ps <- colMeans(contributions)
      expect_equal(ans$expected_richness[k], mean(z), tolerance = 1e-10)
      expect_equal(ans$richness_second_moment[k], mean(z^2), tolerance = 1e-10)
      expect_equal(ans$richness_variance[k], mean((z - mean(z))^2), tolerance = 1e-10)
      expect_lte(ans$sharing_lower[k] - 1e-10, mean(z > 0))
      expect_gte(ans$sharing_upper[k] + 1e-10, mean(z > 0))
      expect_equal(ans$abundance_zeta[k], conditional_zeta(comm, k, "relative")$mean, tolerance = 1e-10)
      if (sum(ps) > 0) {
        expect_equal(ans$effective_incidence_species[k], sum(ps)^2 / sum(ps^2), tolerance = 1e-10)
      } else expect_true(is.na(ans$effective_incidence_species[k]))
    }
    expect_equal(zeta_order_support(comm, c(5, 2)), ans[c(5, 2), ], ignore_attr = "row.names")
  }
  expect_error(zeta_order_support(rbind(0, x)), "positive finite")
})

test_that("stratified conditional draws respect zero and multi-species strata", {
  n <- 7L
  k <- 3L
  x <- 1L * cbind(seq_len(n) <= 4, seq_len(n) %in% c(1, 2, 3, 5, 7),
             seq_len(n) %in% c(2, 3, 5, 6))
  all <- t(combn(n, k))
  z <- zeta_sets(x, all)$shared
  key <- function(idx) apply(idx, 1L, paste, collapse = "|")
  for (exact in c(TRUE, FALSE)) {
    s <- sample_sets_stratified(x, k, 14000L, positive_fraction = .35, mass_draws = 30000L,
                                exact_limit = if (exact) 1000L else 0L, seed = 9183)
    expect_equal(mean(s$shared > 0), .35)
    expect_equal(mean(s$sampling_weights), 1)
    expect_equal(s$positive_mass, mean(z > 0), tolerance = .012)
    expect_equal(s$mass_source, if (exact) "exact" else "independent_pilot")
    id <- match(key(s$idx), key(all))
    expect_false(anyNA(id))
    expect_equal(s$shared, z[id], ignore_attr = TRUE)
    w <- rowsum(s$sampling_weights, id, reorder = TRUE) / nrow(s$idx)
    expect_equal(nrow(w), nrow(all))
    expect_lt(max(abs(w[,1] - 1 / nrow(all))), .009)
  }
  expect_equal(sample_sets_stratified(diag(n), k, 60)$positive_mass, 0)
  expect_equal(sample_sets_stratified(cbind(x, 1), k, 60)$positive_mass, 1)
  fail <- sample_sets_stratified(x, k, 60, exact_limit = 0, max_proposals = 1, seed = 3)
  expect_match(fail$fallback, "proposal limit")
  expect_true(all(fail$sampling_weights == 1))
})

test_that("order limits require ecological stability and handle exhaustive and boundary orders", {
  set.seed(18)
  n <- 18L
  sim <- simulate_metacommunity(n, 30, abundance_mean = 2, seed = 8834)
  ans <- ugdm_order_limits(sim$comm, sim$env["env1"], c(2, 3, n), "env1",
    c(spread = .2, position = .5), .005, pilot_budgets = c(250, 1000), reference_M = 1000,
    max_budget = 3000, nrep = 4, min_delete_ess = 10, shapley = TRUE, seed = 333)
  expect_identical(ans$profile$pass, c(TRUE, TRUE, FALSE))
  expect_equal(ans$profile$confirmed_budget, c(choose(n, 2), choose(n, 3), NA))
  expect_identical(ans$highest_supported_order, 3L)
  expect_equal(ans$profile$minimum_delete_ess[3], 0)
  expect_true(all(is.na(ans$precision$confirmation_mc_sd[ans$precision$order == n])))
  expect_identical(ans$precision$preceding[ans$precision$order == 3 & ans$precision$effect == "spread"],
                   c(NA_integer_, 2L))
  pairdiff <- ans$deletes[[2]]$first - ans$deletes[[1]]$first
  expect_equal(unname(zetaGDUM:::.jackknife_sd(pairdiff)),
    ans$precision$jackknife_1[ans$precision$order == 3 & !is.na(ans$precision$preceding)])
  unstable <- zetaGDUM:::.limit_eco(matrix(1:8, 4), -matrix(1:8, 4), .35)
  expect_false(any(unstable$reliable))
  expect_true(all(unstable$noise_ratio > 1))
})

test_that("undefined incidence becomes a diagnostic failure", {
  x <- matrix(0, 8, 3)
  env <- data.frame(g = seq_len(nrow(x)))
  f <- ugdm_order_limits(x, env, 2L, "g", c(spread = .1, position = .3), .005,
    pilot_budgets = c(16, 32), reference_M = 32, max_budget = 64, nrep = 4, seed = 21)
  expect_identical(f$profile$pass, FALSE)
  expect_equal(f$profile$response_defined, 0)
  expect_identical(f$highest_supported_order, NA_integer_)
  expect_identical(f$retention_break, NA_integer_)
})
