test_that("within-order support requires calibrated intervals and separate numerical precision", {
  set.seed(14)
  env <- data.frame(a = runif(20), b = runif(20))
  original <- ugdm_design(env, t(combn(20L, 2L)))
  design <- ugdm_interaction_design(original, env, "nonlinear", "constant")
  y <- plogis(original$Xpos[, 1L] / 2)
  fit <- ugdm_fit_interactions(design, y, ridge = 0.01)
  numeric <- cbind(a = c(.025, .03, .02), b = c(.01, .014, .008))
  ecological <- cbind(a = seq(.025, .035, len = 30), b = seq(.008, .012, len = 30))
  check <- function(calibrated, mc, ...) ugdm_within_order_diagnostic(
    fit, design, y, numeric, ecological, heldout_D2 = .1,
    min_response_sd = .05, min_heldout_D2 = .01, max_mc_se = mc,
    max_width = .05, max_condition = 1e8, calibrated = calibrated, ...)
  expect_identical(check(FALSE, .02)$status, "prediction_only")
  expect_identical(check(TRUE, .001)$status, "prediction_only")
  expect_identical(check(TRUE, .02)$status, "predictor_inference_supported")
  expect_equal(unname(check(TRUE, .02)$mc_se["a - b"]), sd(numeric[, 1] - numeric[, 2]))
  weights <- seq_along(y)
  shares <- y > median(y)
  diagnostic <- check(FALSE, .02, sharing = shares, sampling_weights = weights)
  expect_equal(diagnostic$sharing, weighted.mean(shares, weights))
  expect_equal(diagnostic$sampling_ESS, sum(weights)^2 / sum(weights^2))
})
