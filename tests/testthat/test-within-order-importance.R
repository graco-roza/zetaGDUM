test_that("within-order grouped attribution has the declared full total", {
  set.seed(12)
  n <- 25L
  env <- data.frame(a = rnorm(n), b = rnorm(n))
  idx <- sample_sets(n, 3L, M = 250L)
  base <- ugdm_design(env, idx)
  d <- ugdm_interaction_design(base, env, "nonlinear", "constant")
  X <- cbind(1, d$dir, d$pos)
  y <- plogis(-.7 + .25 * d$dir[, 1L] + .2 * d$pos[, 1L] +
                rnorm(nrow(idx), sd = .13))
  fit <- ugdm_fit_interactions(d, y, distance = FALSE, ridge = .001)
  a <- ugdm_within_order_importance(fit, d, y)
  expect_equal(sum(a$predictors$absolute), a$predictors$total, tolerance = 1e-7)
  expect_equal(sum(a$components$absolute), a$components$total, tolerance = 1e-7)
  expect_setequal(names(a$predictors$absolute), c("a", "b"))
  expect_setequal(names(a$components$absolute), c("spread", "position"))
  b <- ugdm_within_order_importance(fit, d, y, evaluation = d, evaluation_y = y)
  expect_equal(a$predictors$absolute, b$predictors$absolute, tolerance = 1e-7)
})
