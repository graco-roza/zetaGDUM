test_that("overlap strata recover each fixed-kernel covariance component", {
  n <- 6L
  pair <- t(combn(n, 2))
  z <- seq_len(n) / 5
  x <- cbind(rowMeans(matrix(z[pair], ncol = 2)), z[pair[, 1]] * z[pair[, 2]])
  fit <- ugdm_pair_covariance(x, pair, n)
  overlap <- tcrossprod(vapply(seq_len(n), function(i) rowSums(pair == i), numeric(nrow(pair))))
  brute <- lapply(c(2, 1, 0), function(j) {
    take <- which(overlap == j, arr.ind = TRUE)
    crossprod(x[take[, 1], ], x[take[, 2], ]) / nrow(take)
  })
  expect_equal(fit$first, 4 * (brute[[2]] - brute[[3]]) / n, tolerance = 1e-12)
  expect_equal(fit$second, (brute[[1]] - 2 * brute[[2]] + brute[[3]]) / nrow(pair),
    tolerance = 1e-12)
  del <- t(vapply(seq_len(n), function(i) colMeans(x[rowSums(pair == i) == 0, ]), numeric(2)))
  expect_equal(fit$jackknife, stats::cov(del) * (n - 1)^2 / n)
  expect_error(ugdm_pair_covariance(x, pair[c(2, 2:nrow(pair)), ], n), "unique")
})

test_that("penalised GMM aligns the score equation, derivative and stochastic covariance", {
  set.seed(571)
  n <- 24L
  pair <- t(combn(n, 2))
  z <- rnorm(n)
  Z <- matrix(z[pair], ncol = 2L)
  spread <- abs(Z[, 1L] - Z[, 2L])
  y <- plogis(-.5 + .2 * spread + rowMeans(Z) * .3 + rnorm(nrow(pair), sd = .4))
  design <- list(idx = pair, n = n, k = 2L, Xdir = cbind(z = spread),
                 Xpos = cbind(z = rowSums(Z)), vars = "z", spread_width = 1L)
  a <- list(list(design = design, y = y))
  X <- cbind(1, spread, rowMeans(Z))
  K <- diag(c(0, .12, .2))
  omega <- diag(c(1, 2, 3))
  out <- ugdm_bayes_gmm(a, penalty = list(K), prior_sd = Inf, draws = 25L,
    moment_covariance = function(theta) list(covariance = omega, raw = omega))
  mu <- plogis(drop(X %*% out$estimate))
  expect_lt(max(abs(crossprod(X, y - mu) / length(y) - K %*% out$estimate)), 1e-5)
  expect_equal(unname(out$derivative),
    unname(crossprod(X * sqrt(mu * (1 - mu) / length(y))) + K))
  expect_equal(unname(out$gaussian_mean), unname(out$estimate), tolerance = 1e-10)
  expect_equal(out$gaussian_covariance, out$sandwich_covariance, ignore_attr = TRUE)
  expect_equal(out$covariance_diagnostics$raw, omega)
  tangent <- ugdm_gmm_tangent(out, 2L, active = 2L, draws = 50L, seed = 12)
  expect_true(all(tangent$increments[, 2L] >= 0))
  expect_true(any(tangent$increments[, 2L] == 0))
  interior <- ugdm_gmm_tangent(out, 2L, active = integer(), draws = 50L, seed = 12)
  expect_true(any(interior$increments[, 2L] < 0))
})
