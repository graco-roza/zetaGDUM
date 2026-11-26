test_that("interaction design nests the original and averages transformed values", {
  z <- seq(-2, 2, length.out = 12)
  env <- data.frame(climate = z)
  d <- ugdm_design(env, t(combn(length(z), 3L)))
  a <- ugdm_interaction_design(d, env)
  b <- ugdm_interaction_design(d, env, "nonlinear", "varying")
  z <- drop(scale(z))
  t <- z / max(z)
  expect_equal(unname(b$pos[, 2L]), rowMeans(matrix(t[d$idx]^2, nrow(d$idx))))
  expect_equal(unname(b$pos[, 2L] - rowMeans(matrix(t[d$idx], nrow(d$idx)))^2),
               unname(b$position_variance[, 1L] / max(z)^2))
  expect_equal(unname(b$dir[, 1:3] + b$dir[, 4:6]), unname(a$dir))
  expect_true(all(b$dir >= 0))
  y <- plogis(-.3 + .5 * a$pos[, 1L] + drop(a$dir %*% c(.2, .5, .8)))
  f <- ugdm_fit_variant(d, y, distance = FALSE, position = "mean")
  g <- ugdm_fit_interactions(a, y)
  expect_equal(ugdm_predict(f, d), ugdm_predict(g, a), tolerance = 1e-8)
  h <- ugdm_fit_interactions(b, y, interaction_shrink = .1, position_smooth = .01)
  diag <- ugdm_interaction_diagnostics(h, b, y)
  expect_true(h$converged)
  expect_gte(min(h$coef[2:7]), 0)
  expect_equal(sum(diag$deviance_allocation), h$D2)
  expect_equal(nrow(ugdm_interaction_contrast(h, b, rep(1 / length(y), length(y)))), 1L)
})

test_that("linear varying spread at order two is aliased by curved turnover", {
  env <- data.frame(z = seq(-2, 2, length.out = 12))
  d <- ugdm_design(env, t(combn(12L, 2L)))
  e <- drop(scale(env$z))
  a <- ugdm_interaction_design(d, env, spread = "varying")
  width <- abs(e[d$idx[, 1L]] - e[d$idx[, 2L]])
  midpoint <- rowMeans(matrix(e[d$idx], nrow(d$idx)))
  f <- e + .15 * e^2
  expect_equal(abs(f[d$idx[, 1L]] - f[d$idx[, 2L]]), width * (1 + .3 * midpoint))
  expect_lt(max(abs(lm.fit(cbind(1, a$base), width * midpoint)$residuals)), 1e-10)
  y <- plogis(-1 + rowSums(a$base))
  fit <- ugdm_fit_interactions(a, y, ridge = .001)
  diag <- ugdm_interaction_diagnostics(fit, a, y)
  expect_lt(diag$rank, diag$parameters)
  expect_identical(diag$condition, Inf)
})

test_that("smoothing selection holds out sites and compares the same pools", {
  env <- data.frame(z = seq(-2, 2, length.out = 28))
  folds <- rep(1:4, each = 7L)
  eligible <- function(idx) {
    expect_true(all(apply(idx, 1L, function(i) length(unique(folds[i])) == 1L)) ||
                length(unique(folds[as.vector(idx)])) == 3L)
    plogis(-1 + rowMeans(matrix(env$z[idx]^2, nrow(idx))))
  }
  p <- data.frame(spread_smooth = c(.001, .01), position_smooth = c(.001, .01),
                  interaction_shrink = c(.001, .01), ridge = c(.001, .01))
  selection <- ugdm_select_interactions(env, folds, 3L, eligible, p,
    dist_matrix = coords_dist(cbind(env$z, 0)), M = 80L, seed = 442L)
  expect_equal(dim(selection$fold_loss), c(4L, 2L))
  expect_true(selection$selected %in% 1:2)
  expect_length(selection$selected_without_fold, 4L)
})
