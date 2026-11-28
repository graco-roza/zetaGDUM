source("../phase2/factorial.R", local = TRUE)

test_that("factorial worlds share innovations and the stationary parameters obey the Markov identities", {
  s <- st2_factorial_world(25, "gain", 1, n = 10, T = 6, group_sizes = rep(2,8))
  expect_identical(s$enabled$x[, 1, ], s$disabled$x[, 1, ])
  expect_true(any(s$enabled$x != s$disabled$x))
  expect_identical(s$enabled$x[, , s$group != 1], s$disabled$x[, , s$group != 1])
  g <- with(s$species, pi*(1-rho))
  l <- with(s$species, (1-pi)*(1-rho))
  expect_equal(g/(g+l), s$species$pi)
  expect_equal(1-g-l, s$species$rho)
  s <- st2_factorial_world(25, "null", n = 10, T = 6, group_sizes = rep(2,8))
  expect_identical(s$enabled$x, s$disabled$x)
})

test_that("uncertainty methods preserve complete affected histories", {
  s <- st2_factorial_world(27, "null", n = 14, T = 6, group_sizes = rep(2,8))
  d <- st2_transition_design(s$enabled, s$env, 2, sets = st2_draw_sets(14, 2, 21))
  fit <- st2_fit_transition(s$enabled, d)
  u <- st2_transition_uncertainty(fit, "gain", s$blocks, draws = 5, terms = "position_change.climate")
  ref <- st2_transition_uncertainty(fit, "gain", s$blocks, draws = 5, exact_delete = TRUE, terms = "position_change.climate")
  expect_true(all(is.finite(u$se)))
  expect_lt(abs(u$se[1, "delete_site"] - ref$se[1, "delete_site"]), .07)
  expect_equal(sum(.st2_moving_counts(10, 3)), 10)
})

source("../phase2/prediction_validation.R", local = TRUE)
source("../phase2/effect_validation.R", local = TRUE)

test_that("M2 window probabilities and response-specific order diagnostics agree with their domains", {
  skip_if_not_installed("lme4")
  s <- st2_factorial_world(19, "additive", n = 14, T = 6, group_sizes = rep(3,8))
  p <- s$enabled
  fit <- st2_fit_dynamic(p, st2_transition_design(p, s$env, 1))
  d <- st2_transition_design(p, s$env, 2, 2, sets = matrix(c(1L,2L),1L), basis = fit$design$basis)
  pr <- st2_m2_window(fit, p, s$env, d)
  expect_equal(pr$gain + pr$loss, matrix(1, nrow(pr$gain), ncol(pr$gain)))
  two <- st2_order_diagnostic(p, s$env, 2, 2, budget = 15, mc_replicates = 2)
  expect_equal(nrow(two$orders), 24L)
  expect_equal(sum(two$weights$weight[two$weights$spatial_order==1L & two$weights$temporal_order==1L]), 3)
  expect_false(identical(two$orders$events[two$orders$response == "presence"],
                         two$orders$events[two$orders$response == "dominance"]))
})

test_that("common-random-number factorial contrasts sum across true classes", {
  skip_if_not_installed("lme4")
  s <- st2_factorial_world(29, "gain", 6, n = 14, T = 6, group_sizes = rep(3,8))
  r <- st2_factorial_effect(s, 1, 1L, 1L, budget = 14L, draws = 2)
  a <- subset(r, estimand == "trajectory_surface" & response == "presence")
  expect_equal(a$truth[a$group == 0], sum(a$truth[a$group > 0]))
  expect_equal(a$estimate[a$group == 0], sum(a$estimate[a$group > 0]))
})
