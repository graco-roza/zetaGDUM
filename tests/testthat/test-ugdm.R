test_that("design at order 2 equals the pairwise absolute difference of I-spline bases", {
  sim <- simulate_metacommunity(30, 20, seed = 3)
  idx <- sample_sets(30, 2)
  d <- ugdm_design(sim$env[, "env1", drop = FALSE], idx)
  z <- as.numeric(scale(sim$env$env1))
  B <- unclass(splines2::isp(z, df = 3, degree = 2))[, 1:3]
  expect_equal(unname(d$Xdir[, 1]), abs(B[idx[, 1], 1] - B[idx[, 2], 1]), tolerance = 1e-10)
  expect_equal(unname(d$Xpos[, 1]), z[idx[, 1]] + z[idx[, 2]], tolerance = 1e-10)
})

test_that("the design for sets is the mean pairwise transformed distance", {
  sim <- simulate_metacommunity(30, 20, seed = 4)
  idx <- sample_sets(30, 4, M = 20, seed = 1)
  d <- ugdm_design(sim$env[, "env2", drop = FALSE], idx)
  z <- as.numeric(scale(sim$env$env2))
  B <- unclass(splines2::isp(z, df = 3, degree = 2))[, 2]
  brute <- apply(idx, 1, function(i) mean(dist(B[i])))
  expect_equal(unname(d$Xdir[, 2]), brute, tolerance = 1e-10)
})

test_that("fitted directional coefficients are non-negative and D2 is in [0, 1]", {
  sim <- simulate_metacommunity(50, 60, seed = 5)
  idx <- sample_sets(50, 2)
  d <- ugdm_design(sim$env, idx, coords_dist(sim$coords))
  r <- ugdm_response(sim$comm, idx, "turnover")
  f <- ugdm_fit(d, r$y[r$keep], r$trials[r$keep])
  expect_true(all(f$beta >= -1e-8))
  expect_true(f$D2 >= 0 && f$D2 <= 1)
  imp <- ugdm_importance(f, d, trials = r$trials[r$keep])
  expect_equal(sum(imp), 1, tolerance = 1e-8)
})

test_that("Shapley contributions sum to the D2 of the full model", {
  sim <- simulate_metacommunity(40, 50, seed = 6)
  idx <- sample_sets(40, 2)
  d <- ugdm_design(sim$env, idx, coords_dist(sim$coords))
  r <- ugdm_response(sim$comm, idx, "turnover")
  s <- ugdm_shapley(d, r$y, r$trials, list(env1 = "env1", env2 = "env2", distance = "distance"))
  full <- ugdm_fit(d, r$y, r$trials)$D2
  expect_equal(sum(s), full, tolerance = 1e-8)
  expect_equal(attr(s, "D2"), full, tolerance = 1e-8)
})

test_that("paired draws give the same site weights at every order", {
  w <- site_weights(20, 4, seed = 1)
  expect_equal(colMeans(w), rep(1, 4), tolerance = 1e-12)
  expect_equal(site_weights(20, 4, seed = 1), w)
})

test_that("zeta_orders runs end to end and compares orders", {
  sim <- simulate_metacommunity(50, 60, seed = 7)
  res <- zeta_orders(sim$comm, sim$env, sim$coords, M = 1500, ndraw = 6, groups = list(env = c("env1", "env2"), distance = "distance"), seed = 2)
  expect_s3_class(res, "zeta_orders")
  cmp <- compare_orders(res, "importance")
  expect_true(all(c("change", "lo", "hi", "p_increase", "p_rank_up") %in% names(cmp)))
  expect_equal(sum(cmp$change), 0, tolerance = 1e-8)
  expect_equal(nrow(compare_orders(res, "shapley")), 2)
})

test_that("at order 2 the model reproduces the published pairwise GDUM", {
  skip_on_cran()
  skip_if_not_installed("gdmmTMB")
  sim <- simulate_metacommunity(40, 50, seed = 8)
  P <- (sim$comm > 0) * 1
  idx <- sample_sets(40, 2)
  d <- ugdm_design(sim$env, idx)
  r <- ugdm_response(sim$comm, idx, "turnover")
  mine <- ugdm_fit(d, r$y, r$trials)
  a <- rowSums(P[idx[, 1], ] * P[idx[, 2], ])
  rich <- rowSums(P)
  mn <- pmin(rich[idx[, 1]] - a, rich[idx[, 2]] - a)
  den <- a + mn
  keep <- den > 0
  theirs <- suppressWarnings(suppressMessages(gdmmTMB::gdmm(Y_diss = mn[keep], Y_den = den[keep], D = as.data.frame(idx[keep, ]), X = sim$env,
    diss_formula = ~ gdmmTMB::isp(env1) + gdmmTMB::isp(env2), uniq_formula = ~ scale(env1) + scale(env2), mono = TRUE, family = "binomial", bboot = FALSE)))
  rp <- theirs$obj$report(theirs$obj$env$last.par.best)
  expect_gt(cor(unname(mine$lambda), rp$lambda), 0.9)
})

test_that("linearly dependent predictors trigger a warning", {
  sim <- simulate_metacommunity(30, 20, seed = 9)
  env <- sim$env[, c("env1", "env2")]
  env$env3 <- 1 - env$env1 - env$env2
  idx <- sample_sets(30, 2)
  expect_warning(ugdm_design(env, idx), "linearly dependent")
})

test_that("free-sign directional effects can recover a decline of dissimilarity with spread", {
  set.seed(3)
  n <- 60
  env <- data.frame(a = rnorm(n), b = rnorm(n))
  idx <- sample_sets(n, 2)
  d <- ugdm_design(env, idx)
  # response in which dissimilarity DEcreases with the spread of predictor a
  spread <- abs(env$a[idx[, 1]] - env$a[idx[, 2]])
  y <- plogis(0.5 - 1.2 * spread)
  tr <- rep(50, length(y))
  fc <- ugdm_fit(d, y, tr)
  ff <- ugdm_fit(d, y, tr, free_sign = TRUE)
  expect_gte(min(fc$height), -1e-8)
  expect_lt(ff$height["a"], -0.3)
  expect_gt(ff$D2, fc$D2)
})
