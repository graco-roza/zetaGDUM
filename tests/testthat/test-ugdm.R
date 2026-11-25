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


test_that("ugdm_design warns about strongly skewed predictors", {
  set.seed(1)
  env <- data.frame(a = rnorm(60), b = c(rep(0, 55), runif(5, 0.3, 0.5)))
  idx <- sample_sets(60, 2, M = 200, seed = 1)
  expect_warning(ugdm_design(env, idx), "skewed")
  env$b <- log1p(100 * env$b) + rnorm(60, sd = 1)
  expect_silent(ugdm_design(env, idx))
})

test_that("the C++ design matches the pairwise definition of the directional terms", {
  set.seed(2)
  env <- data.frame(a = rnorm(30), b = runif(30))
  idx <- sample_sets(30, 4, M = 50, seed = 3)
  d <- ugdm_design(env, idx)
  Z <- scale(as.matrix(env)); B <- unclass(splines2::isp(Z[, 1], df = 3, degree = 2))[, 1:3]
  ref <- t(apply(idx, 1, function(s) sapply(1:3, function(b) mean(abs(outer(B[s, b], B[s, b], "-"))[upper.tri(diag(4))]))))
  expect_equal(unname(d$Xdir[, 1:3]), ref, tolerance = 1e-10)
})

test_that("sample_sets_weighted draws sets with probability proportional to the product of site weights", {
  w <- c(0.3, 1, 2, 0.5, 1.5, 3)
  idx <- sample_sets_weighted(w, k = 3, M = 2e5, seed = 1)
  expect_true(all(apply(idx, 1, function(s) length(unique(s)) == 3)))
  key <- apply(idx, 1, paste, collapse = "-")
  all_sets <- t(utils::combn(6, 3)); p <- apply(all_sets, 1, function(s) prod(w[s])); p <- p / sum(p)
  obs <- as.numeric(table(factor(key, levels = apply(all_sets, 1, paste, collapse = "-")))) / length(key)
  expect_lt(max(abs(obs - p)), 0.004)
})

test_that("ugdm_predict reproduces the fitted deviance and ugdm_validate returns held-out metrics", {
  sim <- simulate_metacommunity(n_sites = 80, n_species = 60, seed = 1)
  idx <- sample_sets(80, 3, M = 500, seed = 1)
  r <- ugdm_response(sim$comm, idx, "turnover"); idk <- idx[r$keep, , drop = FALSE]
  d <- ugdm_design(sim$env, idk, coords_dist(sim$coords)); f <- ugdm_fit(d, r$y[r$keep], r$trials[r$keep])
  p <- ugdm_predict(f, d)
  d2 <- 1 - zetaGDUM:::.binom_dev(r$y[r$keep], p, r$trials[r$keep]) / zetaGDUM:::.binom_dev(r$y[r$keep], rep(stats::weighted.mean(r$y[r$keep], r$trials[r$keep]), length(p)), r$trials[r$keep])
  expect_equal(d2, f$D2, tolerance = 1e-6)
  v <- ugdm_validate(sim$comm, sim$env, sim$coords, order = 3, M = 500, nrep = 2, keep_predictions = TRUE)
  expect_true(all(c("D2_test", "r2_test", "slope_test") %in% names(v$summary)))
  expect_equal(nrow(v$summary), 2); expect_true(all(is.finite(v$predictions$predicted)))
})

test_that("zeta_decompose factorises the abundance zeta into occurrence times minimum abundance when shared", {
  set.seed(2); x <- matrix(rpois(10 * 9, 0.8), 10, 9); Q <- x / pmax(rowSums(x), 1); S <- t(utils::combn(10, 3))
  d <- zeta_decompose(x, 3)
  brute_occ <- colMeans(sapply(1:9, function(s) apply(S, 1, function(i) all(x[i, s] > 0))))
  brute_ab <- colMeans(sapply(1:9, function(s) apply(S, 1, function(i) min(Q[i, s]))))
  expect_equal(d$P, brute_occ, tolerance = 1e-12)
  expect_equal(d$term_abundance, brute_ab, tolerance = 1e-12)
  expect_equal(d$term_abundance, d$P * ifelse(is.na(d$m), 0, d$m), tolerance = 1e-12)
  s <- which(d$o >= 3)[1]; occ <- which(x[, s] > 0)
  expect_equal(d$m[s], mean(apply(utils::combn(occ, 3), 2, function(i) min(Q[i, s]))), tolerance = 1e-12)
  expect_equal(attr(d, "ratio"), sum(d$w_occ * ifelse(is.na(d$m), 0, d$m)), tolerance = 1e-12)
  expect_equal(d$w_abu, zeta_weights(x, 3, abundance = TRUE), tolerance = 1e-12)
})

test_that("rarity_species_pool and the range window build the eight forms", {
  sp <- rarity_species_pool(400, dominant = "Sr-Nn-La", p_dominant = 0.7, seed = 1)
  expect_true(all(c("opt1", "opt2", "breadth", "amax", "cx", "cy", "rad", "form") %in% names(sp)))
  expect_equal(length(unique(attr(sp, "forms")$label)), 8)
  expect_gt(mean(sp$form == "Sr-Nn-La"), 0.6)
  expect_true(all(sp$rad[sp$form == "Sr-Nn-La"] < 0.5))
  sim <- simulate_metacommunity(80, species = sp, thin = 1, seed = 2)
  expect_equal(dim(sim$comm), c(80, 400))
  occ <- colMeans(sim$comm > 0); expect_lt(median(occ[sp$form == "Sr-Nn-La"]), median(colMeans(sim$comm[, sp$form != "Sr-Nn-La", drop = FALSE] > 0)) + 1e-9)
  expect_true("form" %in% names(sim$truth))
})

test_that("simulator keeps a site-by-species mean when both gradients are disabled", {
  set.seed(10)
  x <- simulate_metacommunity(n_sites = 17, n_species = 200,
                              use_env = c(FALSE, FALSE), thin = 2, seed = 10)
  expect_equal(dim(x$comm), c(17, 200))
  expect_true(all(is.finite(x$comm)))
  # With no niche or range effects, a site's thinning must act on every species.
  u <- x$env$u
  row_total <- rowSums(x$comm)
  expect_lt(stats::cor(row_total, u), -0.5)
})

test_that("deviance is finite and silent at exact boundary responses", {
  expect_silent(v <- zetaGDUM:::.binom_dev(c(0, 1, 0.4), c(0.2, 0.8, 0.5), c(1, 1, 1)))
  expect_true(is.finite(v))
  expect_gt(v, 0)
})
