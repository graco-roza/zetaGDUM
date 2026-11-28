test_that("position design nests the original and averages transformed values", {
  z <- seq(-2, 2, length.out = 12)
  env <- data.frame(climate = z)
  d <- ugdm_design(env, t(combn(length(z), 3L)))
  a <- ugdm_position_design(d, env)
  b <- ugdm_position_design(d, env, "nonlinear")
  z <- drop(scale(z))
  t <- z / max(z)
  expect_equal(unname(b$pos[, 2L]), rowMeans(matrix(t[d$idx]^2, nrow(d$idx))))
  expect_equal(unname(b$pos[, 2L] - rowMeans(matrix(t[d$idx], nrow(d$idx)))^2),
               unname(b$position_variance[, 1L] / max(z)^2))
  expect_equal(unname(b$dir), unname(a$dir))
  expect_true(all(b$dir >= 0))
  y <- plogis(-.3 + .5 * a$pos[, 1L] + drop(a$dir %*% c(.2, .5, .8)))
  f <- ugdm_fit_variant(d, y, distance = FALSE, position = "mean")
  g <- ugdm_fit_position(a, y)
  expect_equal(ugdm_predict(f, d), ugdm_predict(g, a), tolerance = 1e-8)
  h <- ugdm_fit_position(b, y, position_smooth = .01)
  diag <- ugdm_position_diagnostics(h, b, y)
  expect_true(h$converged)
  expect_gte(min(h$coef[2:4]), 0)
  expect_equal(sum(diag$deviance_allocation), h$D2, tolerance = 0.05)
  expect_equal(nrow(ugdm_position_contrast(h, b, rep(1 / length(y), length(y)))), 1L)
})

test_that("the spread-only comparator is the same design without position columns", {
  z <- seq(-2, 2, length.out = 14)
  env <- data.frame(climate = z)
  d <- ugdm_design(env, t(combn(length(z), 3L)))
  a <- ugdm_position_design(d, env, "nonlinear")
  y <- plogis(-.3 + .4 * a$pos[, 1L] + drop(a$dir %*% c(.2, .5, .8)))
  fit <- ugdm_fit_position(a, y, ridge = .001, spread_smooth = .002, position_smooth = .002)
  diag <- ugdm_position_diagnostics(fit, a, y)
  expect_named(diag$deviance_allocation, c("spread", "position"))
  expect_true(diag$spread_only_D2 <= diag$full_D2 + 1e-8)
  expect_true(is.finite(diag$position_only_D2))
  expect_false(any(grepl("interaction", names(diag))))
})

test_that("spread_only fitting matches the position model with position columns removed", {
  z <- seq(-2, 2, length.out = 14)
  env <- data.frame(climate = z)
  d <- ugdm_design(env, t(combn(length(z), 3L)))
  a <- ugdm_position_design(d, env, "nonlinear")
  y <- plogis(-.3 + .4 * a$pos[, 1L] + drop(a$dir %*% c(.2, .5, .8)))
  s <- ugdm_fit_position(a, y, ridge = .001, spread_smooth = .002, position_smooth = .002, spread_only = TRUE)
  expect_true(s$converged)
  expect_equal(length(s$coef), 1L + ncol(a$dir))
  expect_equal(length(ugdm_predict(s, a)), length(y))
  f <- ugdm_fit_position(a, y, ridge = .001, spread_smooth = .002, position_smooth = .002)
  expect_gte(f$D2, s$D2 - 1e-8)
})

test_that("penalty selection holds out sites and compares the same pools", {
  env <- data.frame(z = seq(-2, 2, length.out = 28))
  folds <- rep(1:4, each = 7L)
  eligible <- function(idx) {
    expect_true(all(apply(idx, 1L, function(i) length(unique(folds[i])) == 1L)) ||
                length(unique(folds[as.vector(idx)])) == 3L)
    plogis(-1 + rowMeans(matrix(env$z[idx]^2, nrow(idx))))
  }
  p <- data.frame(spread_smooth = c(.001, .01), position_smooth = c(.001, .01),
                  ridge = c(.001, .01))
  selection <- ugdm_select_penalties(env, folds, 3L, eligible, p,
    dist_matrix = coords_dist(cbind(env$z, 0)), M = 80L, seed = 442L)
  expect_equal(dim(selection$fold_loss), c(4L, 2L))
  expect_true(selection$selected %in% 1:2)
  expect_length(selection$selected_without_fold, 4L)
})

test_that("explicit interactions are not part of the public interface", {
  exports <- getNamespaceExports("zetaGDUM")
  expect_false(any(grepl("interaction|contextual", exports)))
  expect_false("interaction" %in% eval(formals(ugdm_component_design)$model))
})
