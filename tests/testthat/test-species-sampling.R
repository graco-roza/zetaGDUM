sampling_community <- function() {
  cbind(a = c(1, 1, 1, 1, 0, 0), b = c(2, 1, 0, 0, 1, 0),
        c = c(0, 3, 2, 0, 1, 2), d = c(0, 0, 0, 1, 0, 0),
        e = c(1, 2, 3, 0, 0, 0))
}

test_that("species selection and mixture probabilities agree with enumeration", {
  x <- sampling_community()
  k <- 3
  alpha <- 0.6
  all <- t(utils::combn(nrow(x), k))
  keys <- apply(all, 1, paste, collapse = "-")
  shared <- t(apply(all, 1, function(i) colSums(x[i, ] > 0) == k))
  z <- rowSums(shared)
  expect_true(any(z > 1))
  for (omega in list(NULL, c(0.2, 1.7, 2.3, 0.8, 1, 0),
                      c(0.3, 1, 2, 0.5, 1.5, 3) * 1e270)) {
    log_p <- if (is.null(omega)) rep(0, nrow(all)) else apply(all, 1, function(i) sum(log(omega[i])))
    p <- exp(log_p - max(log_p))
    p <- p / sum(p)
    share_prob <- colSums(p * shared)
    ez <- sum(p * z)
    q <- p * (1 - alpha + alpha * z / ez)
    draw <- sample_sets_species(x, k, 30000, alpha, omega, seed = 47)
    ids <- match(apply(draw$idx, 1, paste, collapse = "-"), keys)
    freq <- tabulate(ids, nrow(all)) / length(ids)
    expect_false(anyNA(ids))
    expect_equal(draw$expected_shared, ez, tolerance = 1e-12)
    expect_equal(draw$species_probability, unname(share_prob / ez), tolerance = 1e-12)
    expect_equal(draw$shared, z[ids])
    expect_equal(draw$sampling_weights, (p / q)[ids])
    expect_true(all(draw$shared[draw$targeted] > 0))
    expect_true(any(draw$shared == 0))
    expect_lt(max(abs(freq - q)), 0.012)
    expect_equal(sum(q * (1 / (1 - alpha + alpha * z / ez)) * z), ez)
  }
})

test_that("mixture retains base support above maximum occupancy", {
  x <- diag(6)
  d <- sample_sets_species(x, 4, 50, omega = c(0, 1:5), seed = 2)
  expect_equal(d$alpha, 0)
  expect_equal(d$expected_shared, 0)
  expect_equal(d$sampling_weights, rep(1, 50))
  expect_false(any(d$targeted))
  expect_true(all(d$idx > 1))
  expect_equal(sample_sets_species(x, 3, 1, omega = c(1, 0, 1, 0, 1, 0))$idx,
               matrix(c(1L, 3L, 5L), 1))
  expect_error(sample_sets_species(x, 2, 5, alpha = 1), "alpha")
  expect_error(sample_sets_species(x, 2.5, 5), "integer")
  expect_error(sample_sets_species(x, 2, 3, omega = c(Inf, 1:5)), "finite")
})

test_that("log symmetric sums remain finite when products overflow or underflow", {
  w <- c(1e-250, 1e-250, 1e-260, 1e250, 1e250, 1e250)
  for (k in c(2, 3, 6)) {
    logs <- utils::combn(log(w), k, sum)
    expected <- max(logs) + log(sum(exp(logs - max(logs))))
    expect_equal(log_esym_cpp(log(w), k), expected, tolerance = 1e-12)
  }
  d <- sample_sets_weighted(w, 3, 10, seed = 20)
  expect_equal(d, matrix(rep(4:6, each = 10), 10))
})

test_that("proposal correction recovers fit, importance and Shapley for both base measures", {
  x <- sampling_community()
  idx <- t(utils::combn(nrow(x), 3))
  des <- ugdm_design(data.frame(a = 1:6, b = c(1, 4, 2, 6, 5, 3)), idx)
  r <- ugdm_response(x, idx, "total")
  z <- zeta_sets(x, idx)$shared
  groups <- list(a = "a", b = "b")
  for (base in list(rep(1, nrow(idx)), set_weights(c(1, 2, 3, 1, 0.5, 2), idx))) {
    correction <- 1 / (0.5 + 0.5 * z / weighted.mean(z, base))
    proposal <- base / correction
    ref <- ugdm_fit(des, r$y, r$trials, base)
    fit <- ugdm_fit(des, r$y, r$trials, proposal, sampling_weights = correction)
    expect_equal(fit, ref, tolerance = 1e-8)
    expect_equal(ugdm_importance(fit, des, proposal, r$trials, correction),
                 ugdm_importance(ref, des, base, r$trials), tolerance = 1e-8)
    expect_equal(ugdm_shapley(des, r$y, r$trials, groups, proposal, sampling_weights = correction),
                 ugdm_shapley(des, r$y, r$trials, groups, base), tolerance = 1e-8)
  }
})

test_that("working weights preserve response definitions and legacy defaults", {
  x <- sampling_community()
  x[1, ] <- 0
  idx <- t(utils::combn(nrow(x), 2))
  for (component in c("turnover", "total", "nestedness", "relative", "balanced")) {
    a <- ugdm_response(x, idx, component, abundance_trials = 3)
    b <- ugdm_response(x, idx, component, abundance_trials = 3, working_weights = "constant")
    expect_identical(a$y, b$y)
    expect_identical(a$keep, b$keep)
    expect_true(all(b$trials == if (component == "balanced") 3 else 1))
  }
  expect_true(any(ugdm_response(x, idx, "total")$trials > 1))
})

test_that("order pipeline supports both bootstrap schemes and explicit no-sharing orders", {
  sim <- simulate_metacommunity(45, 55, seed = 4)
  groups <- list(a = "env1", b = "env2")
  for (weighting in c("direct", "product")) {
    res <- zeta_orders(sim$comm, sim$env, orders = c(2, 4), M = 1000, ndraw = 3, groups = groups,
                       sampling = "species", weighting = weighting, working_weights = "constant")
    expect_equal(unname(res$point_diagnostics[, "alpha"]), c(0, 0.5))
    expect_equal(unname(res$draw_diagnostics[, "k4", "alpha"]), rep(0.5, 3))
    expect_equal(apply(res$shapley, c(1, 2), sum), res$D2, tolerance = 1e-8)
    expect_equal(unname(apply(res$importance, c(1, 2), sum)), matrix(1, 3, 2), tolerance = 1e-8)
  }
  base <- zeta_orders(sim$comm, sim$env, orders = 4, M = 200, ndraw = 2)
  zero <- zeta_orders(sim$comm, sim$env, orders = 4, M = 200, ndraw = 2,
                      sampling = "species", alpha = 0)
  expect_identical(base$point, zero$point)
  expect_identical(base$D2, zero$D2)
  x <- diag(12)
  res <- zeta_orders(x, data.frame(a = 1:12), orders = 4, M = 100, ndraw = 2, sampling = "species")
  expect_equal(unname(res$point_diagnostics[1, c("n_shared", "alpha")]), c(0, 0))
  expect_true(all(is.na(res$D2)))
  expect_true(is.na(res$point$k4$D2))
  expect_equal(compare_orders(res)$draws, 0)
})

test_that("held-out scoring corrects proposal bias and can reuse conventional test draws", {
  sim <- simulate_metacommunity(50, 35, seed = 1)
  v <- ugdm_validate(sim$comm, sim$env, order = 4, M = 600, nrep = 1,
                     sampling = "species", keep_predictions = TRUE)
  p <- v$predictions
  w <- p$trials * p$sampling_weights
  expect_equal(v$summary$rmse_test, sqrt(weighted.mean((p$observed - p$predicted)^2, w)))
  expect_equal(v$summary$slope_test, sum(p$sampling_weights * p$observed * p$predicted) /
                 sum(p$sampling_weights * p$predicted^2))
  expect_true(any(p$sampling_weights != 1))
  a <- ugdm_validate(sim$comm, sim$env, order = 4, M = 300, nrep = 1, keep_predictions = TRUE)
  b <- ugdm_validate(sim$comm, sim$env, order = 4, M = 300, nrep = 1, keep_predictions = TRUE,
                     sampling = "species", test_sampling = "conventional")
  expect_identical(a$predictions$observed, b$predictions$observed)
})
