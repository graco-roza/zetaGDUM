test_that("conditional paths share prefixes and count every evaluation", {
  a <- sample_sets_conditioned(7, c(4, 2, 5), seq_len(7), seed = 19)
  expect_identical(a[["5"]]$idx[, 1L], rep(seq_len(7), seq_len(7)))
  expect_identical(a[["4"]]$idx, a[["5"]]$idx[, 1:4])
  expect_identical(a[["2"]]$idx, a[["5"]]$idx[, 1:2])
  expect_false(any(apply(a[["5"]]$idx, 1L, anyDuplicated) > 0))
  expect_equal(attr(a, "evaluations"), 3 * sum(1:7))
  low <- sample_sets_conditioned(7, 2, 2000, seed = 19)
  first <- low[[1L]]$idx
  expect_lt(max(abs(as.numeric(table(first[first[, 1L] == 1L, 2L])) / 2000 - 1/6)), .035)
})

test_that("conditional means and score identities agree with complete enumeration", {
  n <- 6L
  env <- data.frame(z = c(-2, -.7, .1, .2, .9, 2))
  x <- rbind(c(4, 0, 2), c(0, 2, 1), c(1, 4, 2), c(4, 0, 1), c(1, 2, 0), c(2, 3, 1))
  theta <- list(c(-.4, .12, .03, .21, -.31), c(.3, .04, .05, .11, .27))
  L <- c(0, -1, -1, -1, 0, 0, 1, 1, 1, 0)
  a <- ugdm_conditioned_score_plan(x, env, c(2, 4), 20, theta,
    controls = c("response", "predictors"), pilot_draws = 12, adaptive = TRUE,
    allocation_contrast = L, minimum_draws = 4, seed = 28)
  expect_equal(sum(a$counts), 120)
  expect_gte(min(a$counts), 4)
  expect_gt(length(unique(a$counts)), 1)
  r <- ugdm_conditioned_scores(a, theta, L)
  expect_equal(unname(r$contrast), unname(r$pseudovalues %*% L))
  expect_equal(colMeans(r$pseudovalues), r$score)

  # Replace production draws with every focal permutation prefix exactly once.
  orders <- c(2L, 4L)
  for (j in seq_along(orders)) {
    k <- orders[j]
    idx <- t(combn(n, k))
    d <- ugdm_design(env, idx)
    y <- ugdm_response(x, idx, "relative")$y
    exact <- ugdm_score_pseudovalues(d, y, theta[[j]])
    focal <- rep(1:n, each = choose(n - 1, k - 1))
    sampled <- do.call(rbind, lapply(1:n, function(i) {
      values <- idx[apply(idx == i, 1L, any), , drop = FALSE]
      cbind(i, matrix(t(values)[t(values) != i], ncol = k - 1L, byrow = TRUE))
    }))
    mode(sampled) <- "integer"
    sets <- structure(list(list(idx = sampled, sampling_weights = rep(1, nrow(sampled)))),
      focal = focal, draws = as.integer(tabulate(focal)), evaluations = nrow(sampled))
    control <- zetaGDUM:::.conditioned_data(x, env, sets, NULL, d, "relative", "constant",
      c("response", "predictors"), NULL)
    control$coefficients <- a$coefficients[j]
    control$strata <- split(seq_along(focal), focal)
    control$pilot_evaluations <- 0L
    p <- control$orders[[1L]]
    expect_equal(unname(rowsum(p$control, focal) / control$counts), unname(p$expectation))
    checked <- ugdm_conditioned_scores(control, theta[[j]])
    expect_equal(unname(checked$pseudovalues), unname(exact$pseudovalues), tolerance = 1e-14)
    expect_equal(unname(checked$score), unname(exact$full), tolerance = 1e-14)
  }
})

test_that("incidence controls use distinct, exact focal and weighted targets", {
  n <- 6L
  x <- rbind(c(0, 0), c(0, 0), c(1, 1), c(0, 1), c(1, 0), c(1, 1))
  env <- data.frame(z = c(-2, -.7, .1, .2, .9, 2))
  expect_error(ugdm_conditioned_score_plan(x, env, 2, type = "incidence"), "All-empty")
  a <- ugdm_conditioned_score_plan(x, env, 2, 10, rep(0, 5),
    type = "incidence", weighting = "richness", controls = c("response", "focal_incidence"),
    pilot_draws = 12, seed = 28)
  expect_true(all(is.finite(ugdm_conditioned_scores(a, rep(0, 5))$pseudovalues)))
  targets <- a$orders[[1]]$expectation
  expect_equal(unname(targets[, "response"]), conditional_zeta_dissimilarity(x, 2, "incidence")$site)
  expect_equal(unname(targets[1:2, "focal_incidence"]), c(0, 0))
  expect_false(identical(targets[, "response"], targets[, "focal_incidence"]))
})
