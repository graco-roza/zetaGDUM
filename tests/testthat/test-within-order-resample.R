test_that("within-order draws use a common null and coherent predictor rankings", {
  set.seed(721)
  n <- 27L
  env <- data.frame(temperature = rnorm(n), moisture = rnorm(n))
  make_y <- function(d) plogis(-.4 + 1.2 * rowSums(d$dir[, 1:3]) +
    .5 * d$pos[, 1] + .25 * rowMeans(matrix(site_noise[d$idx], nrow(d$idx))))
  site_noise <- rnorm(n)
  idx <- sample_sets(n, 3L, M = 340)
  base <- ugdm_design(env, idx)
  d <- ugdm_position_design(base, env, "nonlinear")
  y <- make_y(d)
  fit <- ugdm_fit_position(d, y, distance = FALSE, ridge = .002)
  build <- function(idx, sites) {
    dd <- ugdm_position_design(ugdm_design(env[sites, ], idx, ref = base),
                                  env[sites, ], "nonlinear")
    list(design = dd, y = make_y(dd))
  }
  methods <- c("bayesian_bootstrap", "gaussian_gmm", "multiplier", "site_bootstrap")
  for (m in methods) {
    z <- ugdm_within_order_resample(fit, d, y, evaluation = d, evaluation_y = y,
      method = m, draws = 9L, M = 340L, build = build, threshold = .03, seed = 918L)
    expect_true(all(z$valid), info = paste(m, paste(z$errors, collapse = " ")))
    expect_equal(sum(z$predictors$p_leading), 1)
    expect_equal(drop(z$difference_draws), z$predictor_draws[, 1] - z$predictor_draws[, 2])
    expect_equal(rowSums(z$component_draws), rowSums(z$predictor_draws), tolerance = 1e-8)
    expect_lte(z$differences$p_practical, z$differences$p_greater)
  }
  broken <- function(idx, sites) stop("response unavailable")
  bad <- ugdm_within_order_resample(fit, d, y, draws = 2L, build = broken)
  expect_false(any(bad$valid))
  expect_equal(bad$errors, rep("response unavailable", 2L))
  expect_true(all(is.na(bad$predictors$median)))
})

test_that("exact perturbation solves the original scaled loss", {
  X <- cbind(1, c(-1, 0, 2))
  y <- c(.2, .4, .55)
  w <- c(.2, .5, .3)
  P <- diag(c(0, .1))
  y2 <- y + c(.04, -.01, .015)
  shift <- drop(crossprod(X, w * (y2 - y)))
  a <- .variant_solution(X, y2, w, P, 2L)
  b <- .variant_solution(X, y, w, P, 2L, score_shift = shift)
  expect_true(a$converged && b$converged)
  expect_equal(a$coef, b$coef, tolerance = 1e-6)
})
