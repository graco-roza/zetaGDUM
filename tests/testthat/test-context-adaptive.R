test_that("contextual spread follows set versus pair midpoint algebra", {
  env <- data.frame(z = c(-2, -1, -.6, 0, 1, 2.2, 2.7, 4))
  idx <- t(combn(nrow(env), 3L))
  base <- ugdm_design(env, idx)
  d <- ugdm_interaction_design(base, env, "nonlinear", "contextual")
  z <- as.numeric(scale(env$z))
  basis <- base$environment_basis$splines[[1L]]
  B <- predict(basis, z)
  direct <- t(apply(idx, 1L, function(s) {
    edges <- combn(s, 2L)
    widths <- abs(B[edges[1L, ], ] - B[edges[2L, ], ])
    colMeans((mean(z[s]) - colMeans(matrix(z[edges], nrow(edges)))) * widths)
  }))
  limits <- attr(basis, "Boundary.knots")
  expect_equal(unname(d$contextual), unname(direct), tolerance = 1e-12)
  expect_equal(unname(d$bounded_context) * diff(limits), unname(direct), tolerance = 1e-12)
  expect_equal(unname(d$dir[, 1:3] + d$dir[, 4:6]), unname(base$Xdir))
  expect_equal(unname(d$dir[, 7:9] + d$dir[, 10:12]), unname(base$Xdir))
  expect_true(all(d$dir >= 0))
  policy <- ugdm_contextual_diagnostics(d, seq_len(nrow(idx)))
  expect_gt(unname(policy$rank_increment["context"]), 0)
  expect_lt(unname(policy$rank_increment["midpoint"]), ncol(d$pair_midpoint))
  expect_true(all(policy$context_residual_fraction <= 1 + 1e-12))
  y <- plogis(-.3 + rowSums(d$dir) + d$pos[, 1L] + 2 * d$contextual[, 1L])
  fit <- ugdm_fit_interactions(d, y, interaction_penalty = "derivative",
    penalty_scale = "range", interaction_shrink = 1e-5, ridge = 1e-5)
  expect_true(fit$converged)
  expect_true(all(fit$coef[2:13] >= 0))
  expect_equal(length(ugdm_predict(fit, d)), nrow(idx))
  expect_equal(ugdm_interaction_diagnostics(fit, d, y)$contextual$rank, policy$rank)
  pair <- ugdm_interaction_design(ugdm_design(env, t(combn(nrow(env), 2L)), ref = base),
                                  env, "nonlinear", "contextual")
  expect_true(all(pair$contextual == 0))
  expect_equal(unname(ugdm_contextual_diagnostics(pair)$rank_increment["context"]), 0L)
})

test_that("function-space penalty measures actual slope differences", {
  env <- data.frame(z = seq(-3, 3, length.out = 24L))
  d <- ugdm_interaction_design(ugdm_design(env, sample_sets(24L, 4L, M = 300L)),
                               env, "nonlinear", "varying")
  b <- d$basis$splines[[1L]]
  limits <- attr(b, "Boundary.knots")
  beta <- c(.2, -.24, .13)
  integral <- integrate(function(x) drop(predict(b, x, derivs = 1L) %*% beta)^2,
                        lower = limits[1L], upper = limits[2L])$value
  expected <- drop(crossprod(beta, zetaGDUM:::.context_penalty_basis(b)$derivative %*% beta))
  expect_equal(expected, integral, tolerance = 1e-12)
  a <- zetaGDUM:::.interaction_layout(d, FALSE)
  P <- function(..., penalty_scale = "native") {
    zetaGDUM:::.interaction_penalties(d, a, c(...), "derivative", penalty_scale)
  }
  quadratic <- function(b, P) drop(crossprod(b, P %*% b))
  basis_difference <- c(0, rep(.4, 3L), .4 + beta, rep(0, 3L))
  expect_equal(quadratic(basis_difference, P(interaction_shrink = 1)), integral,
               tolerance = 1e-12)
  expect_equal(quadratic(basis_difference, P(interaction_shrink = 1, penalty_scale = "range")),
               integral * diff(limits), tolerance = 1e-12)
  common_linear <- c(0, rep(1, 6), .5, 0, 0)
  expect_equal(quadratic(common_linear, P(interaction_shrink = 1, spread_smooth = 1,
    position_smooth = 1)), 0, tolerance = 1e-12)
  expect_gt(quadratic(common_linear, P(spread_null = 1, position_null = 1)), 0)
})

test_that("repeated site blocks allow a predeclared contrast criterion", {
  n <- 24L
  env <- data.frame(z = seq(-1, 1, length.out = n))
  coords <- cbind(env$z, sin(env$z * 8))
  folds <- ugdm_site_partitions(n, 3L, 2L, coords, seed = 95)
  expect_equal(dim(folds), c(n, 2L))
  expect_true(all(apply(folds, 2L, function(x) all(table(x) == n / 3L))))
  p <- data.frame(spread_smooth = c(0, .001), position_smooth = c(0, .002),
                  interaction_shrink = c(.001, .02), ridge = c(.001, .01))
  response <- function(idx) plogis(-1 + rowMeans(matrix(env$z[idx]^2, nrow(idx))))
  contrast <- function(fit, design, y, weights) as.numeric(fit$strengths[["ridge"]] < .01)
  out <- ugdm_select_interactions(env, folds, 3L, response, p, M = 80L, seed = 45,
    contrast_loss = contrast, selection = "contrast", interaction_penalty = "derivative")
  expect_equal(out$selected, 2L)
  expect_equal(dim(out$fold_loss), c(6L, 2L))
  expect_equal(out$mean_contrast_loss, c(1, 0))
  expect_true(all(out$selected_without_fold == 2L))
})

test_that("symmetric sets do not gain spurious rank from tiny columns", {
  env <- data.frame(z = seq(-2, 2, length.out = 20L))
  first <- t(combn(10L, 3L))
  d <- ugdm_interaction_design(ugdm_design(env, cbind(first, 21L - first)), env,
                               "nonlinear", "contextual")
  audit <- ugdm_contextual_diagnostics(d)
  expect_equal(unname(audit$rank_increment["context"]), 0L)
  expect_lt(max(abs(rowSums(d$contextual))), 1e-12)
  expect_gt(max(abs(d$contextual)), .005)
  layout <- zetaGDUM:::.interaction_layout(d, FALSE)
  strength <- c(spread_smooth = 1, position_smooth = 1, interaction_shrink = 1,
                spread_null = 1, position_null = 1, ridge = 1)
  P <- zetaGDUM:::.interaction_penalties(d, layout, strength, "derivative", "range")
  same_model <- c(0, rep(.5, 6L), rep(-.5, 6L), rep(0, 3L))
  expect_lt(max(abs(P %*% same_model)), 1e-10)
  expect_lt(max(abs(layout$X %*% same_model)), 1e-12)
})
