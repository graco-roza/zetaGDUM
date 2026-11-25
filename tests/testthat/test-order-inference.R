test_that("predictor transformations use only the basis sites and carry to new sites", {
  env <- data.frame(a = c(1:8, 400, 900), b = c(3, 6, 2, 5, 9, 10, 8, 7, -100, 90))
  train <- 1:8
  it <- t(utils::combn(train, 3))
  test <- rbind(c(1L, 9L, 10L), c(5L, 9L, 10L))
  fit_design <- ugdm_design(env, it, basis_sites = train)
  d <- ugdm_design(env, test, ref = fit_design)
  train_only <- ugdm_design(env[train, ], it)
  expect_equal(fit_design$Xdir, train_only$Xdir)
  expect_equal(fit_design$Xpos, train_only$Xpos)
  expect_equal(fit_design$environment_basis, train_only$environment_basis)
  expect_equal(d$Xpos[, "a"], rowSums(matrix((env$a[test] - mean(env$a[train])) / sd(env$a[train]), nrow(test))))
  altered <- env
  altered$a[9:10] <- c(5000, 9000)
  expect_equal(ugdm_design(altered, test, ref = fit_design)$Xdir, d$Xdir)
  expect_false(isTRUE(all.equal(ugdm_design(altered, test, ref = fit_design)$Xpos, d$Xpos)))
})

test_that("held-out predictor extremes cannot change training fit metrics", {
  sim <- simulate_metacommunity(32, 48, seed = 22)
  env <- sim$env
  set.seed(16)
  test <- sample.int(nrow(env), round(0.3 * nrow(env)))
  altered <- env
  altered$env1[test] <- altered$env1[test] + 50
  original <- ugdm_validate(sim$comm, env, order = 3, M = 300, nrep = 1, seed = 15)
  extreme <- ugdm_validate(sim$comm, altered, order = 3, M = 300, nrep = 1, seed = 15)
  expect_identical(original$summary$D2_train, extreme$summary$D2_train)
  expect_false(isTRUE(all.equal(original$summary$rmse_test, extreme$summary$rmse_test)))
})

test_that("delete-one and multiplier calculations use sites and the correct omission mass", {
  set.seed(24)
  n <- 12L
  k <- 4L
  idx <- t(utils::combn(n, k))
  d <- ugdm_design(data.frame(g = rnorm(n)), idx)
  eta <- .5 + .12 * d$Xpos[, 1] + drop(d$Xdir %*% c(.4, .3, .6))
  mu <- plogis(eta)
  y <- mu + .4 * mu * (1 - mu) * rowMeans(matrix(stats::runif(n, -1, 1)[idx], nrow(idx)))
  t <- rep(1, length(y))
  s <- ugdm_site_uncertainty(d, y, t, free_sign = TRUE, seed = 12, ndraw = 99)
  expect_equal(dim(s$delete), c(n, 2))
  expect_true(all(is.finite(s$jackknife_se)))
  expect_true(all(s$delete_ess == choose(n - 1L, k)))
  correction <- (n - 1) / (n - k)
  expect_equal(sweep(s$multiplier, 2, s$effect, "-"),
               correction * sweep(s$multiplier_raw, 2, s$effect, "-"), tolerance = 1e-10)
  expect_equal(colSums(s$influence), setNames(c(0, 0), names(s$effect)), tolerance = 1e-10)
  keep <- rowSums(idx == 1L) == 0
  deleted_fit <- ugdm_fit(.subset_order_design(d, keep), y[keep], t[keep], free_sign = TRUE)
  expect_equal(s$delete[1, ], .order_effects(deleted_fit))
})

test_that("site inference diagnoses an active constraint and a single possible set", {
  env <- data.frame(g = seq(-1, 1, length.out = 12))
  idx <- t(utils::combn(12, 3))
  d <- ugdm_design(env, idx)
  y <- plogis(.1 + .3 * d$Xpos[, 1] - 2 * rowSums(d$Xdir))
  s <- ugdm_site_uncertainty(d, y, rep(1, length(y)))
  expect_null(s$multiplier)
  expect_match(s$multiplier_unavailable, "active")
  d <- ugdm_design(env, matrix(1:12, 1))
  s <- ugdm_site_uncertainty(d, .6, 1)
  expect_null(s$multiplier)
  expect_true(all(is.na(s$jackknife_se)))
  expect_match(s$multiplier_unavailable, "no sets")
})

test_that("precision curves retain paired changes, sharing and identification diagnostics", {
  sim <- simulate_metacommunity(24, 60, seed = 19)
  z <- zeta_order_precision(sim$comm, sim$env, c(2, 3), budgets = c(300, 800),
                            reference_M = 2000, nrep = 3,
                            groups = list(a = "env1", b = "env2"), mc_fraction = 10, seed = 9)
  expect_true("k3 - k2" %in% z$precision$contrast)
  expect_true(all(c("spread.env1", "position.env2", "importance.env1", "shapley.a") %in% z$precision$effect))
  expect_true(all(z$diagnostics$exact[z$diagnostics$order == 2]))
  expect_true(all(z$diagnostics$sharing_lower <= z$diagnostics$sharing_upper + 1e-12))
  expect_true(all(z$diagnostics$min_delete_ess > 0))
  if (!is.na(z$confirmed_budget)) expect_true(any(z$precision$confirmation))
  p <- ugdm_sampling_pilot(sim$comm, sim$env, 3, fractions = c(0, .5), M = 400, nrep = 4)
  expect_equal(p$candidates$alpha, c(0, .5))
  expect_true(p$alpha %in% p$candidates$alpha)
})

test_that("site uncertainty can identify a training population inside a larger design", {
  env <- data.frame(g = c(seq(-1, 1, length.out = 12), 20, 50))
  sites <- 1:12
  idx <- t(utils::combn(sites, 3))
  full_design <- ugdm_design(env, idx, basis_sites = sites)
  train_design <- ugdm_design(env[sites, , drop = FALSE], idx)
  y <- plogis(.3 + drop(full_design$Xdir %*% c(.25, .35, .45)) - .1 * full_design$Xpos)
  a <- ugdm_site_uncertainty(full_design, drop(y), rep(1, length(y)), sites = sites, seed = 2)
  b <- ugdm_site_uncertainty(train_design, drop(y), rep(1, length(y)), seed = 2)
  expect_equal(a$delete, b$delete)
  expect_false(is.null(a$multiplier))
  expect_equal(a$multiplier, b$multiplier)
  expect_error(ugdm_site_uncertainty(full_design, drop(y), rep(1, length(y)), sites = 1:6), "include every site")
})

test_that("precision diagnostics can include a geographic distance predictor", {
  sim <- simulate_metacommunity(22, 50, seed = 85)
  a <- zeta_order_precision(sim$comm, sim$env, 3, budgets = 300, nrep = 3, reference_M = 600,
                            coords = sim$coords, groups = list(environment = names(sim$env), distance = "distance"))
  expect_true(all(c("spread.distance", "importance.distance", "shapley.distance") %in% a$precision$effect))
})

test_that("paired order fits can use pilot-selected or prespecified fractions", {
  sim <- simulate_metacommunity(26, 42, seed = 9)
  a <- zeta_orders(sim$comm, sim$env, orders = c(2, 4), M = 400, ndraw = 2,
                  sampling = "species", alpha = "adaptive", pilot_M = 400, pilot_nrep = 3, seed = 93)
  expect_equal(a$alpha_by_order[[1]], 0)
  expect_equal(a$alpha_by_order[[2]], a$pilots$k4$alpha)
  b <- zeta_orders(sim$comm, sim$env, orders = c(2, 4), M = 400, ndraw = 2,
                  sampling = "species", alpha = a$alpha_by_order, seed = 93)
  expect_identical(a$point, b$point)
  expect_identical(a$D2, b$D2)
  expect_error(zeta_orders(sim$comm, sim$env, alpha = "adaptive"), "requires")
})
