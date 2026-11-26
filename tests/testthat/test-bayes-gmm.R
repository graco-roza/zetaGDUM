test_that("frozen moments give the sandwich and a prior-adjusted Gaussian approximation", {
  set.seed(1081)
  n <- 22L
  z <- runif(n, -2, 2)
  make <- function(k) {
    idx <- t(combn(n, k))
    Z <- matrix(z[idx], nrow(idx))
    sp <- rowMeans(abs(Z[, 1L] - Z[, -1L, drop = FALSE]))
    y <- plogis(-.7 - .3 * sp + .2 * rowMeans(Z) + rnorm(nrow(idx), 0, .35))
    design <- list(idx = idx, n = n, k = k, Xdir = cbind(g = sp),
      Xpos = cbind(g = rowSums(Z)), vars = "g", spread_width = 1L)
    list(design = design, y = y)
  }
  a <- lapply(c(2L, 3L), make)
  out <- ugdm_bayes_gmm(a, draws = 100L, seed = 13, project = TRUE)
  V <- solve(out$derivative, out$omega) %*% solve(out$derivative) / n
  G <- solve(solve(V) + diag(1 / 4, 6L))
  expect_equal(unname(out$sandwich_covariance), unname(V), tolerance = 1e-8)
  expect_equal(unname(out$gaussian_covariance), unname(G), tolerance = 1e-8)
  expect_equal(unname(out$gaussian_mean), drop(G %*% solve(V, unname(out$estimate))),
    tolerance = 1e-8)
  expect_true(all(out$projected[, c(2L, 5L)] >= 0))
  expect_true(any(out$draws[, c(2L, 5L)] < 0))
  expect_true(any(abs(out$projected[, c(1L, 3L, 4L, 6L)] - out$draws[, c(1L, 3L, 4L, 6L)]) > 0))
  expect_error(ugdm_bayes_gmm(a, site_scores = function(theta) matrix(1, n, length(theta))),
    "zero or invalid marginal")
  J <- do.call(cbind, lapply(seq_along(a), function(j) {
    ugdm_pool_scores(list(a[[j]]), out$pilot_estimate[(j - 1L) * 3L + seq_len(3L)], "centered")
  }))
  same <- ugdm_bayes_gmm(a, site_scores = function(theta) J, draws = 100L, seed = 13)
  expect_equal(same$draws, out$draws)
  design_matrix <- function(a) cbind(intercept = 1, g = a$design$Xdir[, 1L],
                                    g.position = a$design$Xpos[, 1L] / a$design$k)
  custom <- ugdm_bayes_gmm(a, design_matrix = design_matrix, draws = 100L, seed = 13)
  expect_equal(custom$draws, out$draws)
  expect_error(ugdm_bayes_gmm(a, design_matrix = design_matrix, project = TRUE),
    "Supply bounded")
  nonlinear <- ugdm_bayes_gmm(a, integration = "mcmc", draws = 100L, warmup = 100L,
    chains = 2L, seed = 20)
  expect_equal(dim(nonlinear$draws), c(200L, 6L))
  expect_true(all(is.finite(nonlinear$draws)))
})
