test_that("order information reproduces exhaustive compositional moments", {
  x <- rbind(c(3, 0, 0, 2), c(0, 4, 1, 2), c(5, 1, 0, 0),
             c(1, 1, 2, 0), c(0, 1, 0, 4), c(3, 0, 2, 2))
  f <- zeta_order_information(x, orders = c(4L, 2L, 6L), response = "all", method = "exact")
  variance <- function(y) mean((y - mean(y))^2)
  exact <- lapply(c(4L, 2L, 6L), function(k) {
    idx <- t(combn(nrow(x), k))
    z <- species_contributions(x, idx, "relative")
    occ <- rowSums(z$dissimilarity * (z$similarity == 0))
    het <- rowSums(z$dissimilarity * (z$similarity > 0))
    expect_equal(occ + het, ugdm_response(x, idx, "relative")$y)
    y <- cbind(ugdm_response(x, idx, "total")$y, z$response, occ, het)
    rows <- f[f$order == k, ]
    expect_equal(rows$response_mean, unname(colMeans(y)), tolerance = 1e-13)
    expect_equal(rows$sharing_probability, rep(mean(zeta_sets(x, idx)$shared > 0), 4L))
    apply(y, 2L, variance)
  })
  expect_identical(unique(f$order), c(4L, 2L, 6L))
  expect_equal(f$response_variance, unname(unlist(exact)), tolerance = 1e-13)
  expect_equal(f$relative_information, unname(unlist(lapply(exact, function(v) v / exact[[2L]]))))
  expect_true(all(f$variance_se == 0))
  expect_true(all(f$recommended_status == "ecology_only"))
})

test_that("an exact flat pairwise response leaves the index undefined", {
  x <- cbind(c(1,1,1,0), c(1,0,0,1), c(0,1,0,1), c(0,0,1,1),
             c(1,0,0,0), c(0,1,0,0), c(0,0,1,0))
  f <- zeta_order_information(x, orders = 2:4, response = "incidence", method = "exact")
  expect_equal(f$response_variance[1L], 0)
  expect_gt(f$response_variance[2L], 0)
  expect_true(all(is.na(f$relative_information)))
  g <- zeta_order_information(cbind(x, c(1,0,0,0)), orders = 2:4, response = "incidence")
  expect_gt(g$relative_information[2L], 1)

  pairs <- t(combn(4L, 2L))
  q <- matrix(0, 4L, nrow(pairs))
  for (j in seq_len(nrow(pairs))) q[pairs[j, ], j] <- 1
  a <- zeta_order_information(q, orders = 2:3, response = "all")
  expect_true(all(a$response_variance == 0))
  expect_true(all(is.na(a$relative_information)))
  expect_identical(a$variance_method[a$order == 2 & a$response == "relative"], "enumerated")
})

test_that("position and mean pairwise spread variance have exact finite-site forms", {
  x <- cbind(1, diag(6))
  env <- data.frame(temperature = c(-2,-.3,.1,.5,.7,3), water = c(0,3,1,4,2,7))
  coord <- cbind(env$temperature, env$water)
  f <- zeta_order_information(x, env, coord, 2:6, response = "relative")
  detail <- attr(f, "environment_detail")
  pop <- function(y) mean((y - mean(y))^2)
  dist <- coords_dist(coord)
  for (k in 2:6) {
    idx <- t(combn(nrow(x), k))
    base <- ugdm_design(env, idx, dist)
    actual <- ugdm_position_design(base, env, "nonlinear")
    d <- detail[detail$order == k, ]
    expect_equal(d$position_variance, unname(apply(actual$mean_position, 2, pop)))
    expect_equal(d$quadratic_variance, unname(apply(actual$pos[, c(2,5), drop=FALSE], 2, pop)))
    expect_equal(d$cubic_variance, unname(apply(actual$pos[, c(3,6), drop=FALSE], 2, pop)))
    expect_equal(as.matrix(d[grep("spread_basis",names(d))]),
      t(matrix(apply(actual$dir, 2, pop), nrow=3)), ignore_attr=TRUE)
    expect_equal(f$geographic_spread_variance[f$order == k], pop(zetaGDUM:::meanpair_dist_cpp(dist,idx)))
  }
  expect_identical(tail(f$recommended_status,1), "unsupported")
  bad <- zeta_order_information(x, data.frame(flat = 1, water = env$water), orders = 2)
  expect_identical(bad$unsupported_predictors, 1L)
  expect_identical(bad$recommended_status, "unsupported")
})

test_that("independent numerical pools expose error for a single response", {
  set.seed(3)
  x <- cbind(1, matrix(rpois(25 * 8, 1.2), 25))
  g <- function(seed) zeta_order_information(x, orders = c(3,5), method = "sample",
    combinations = 100, pools = 4, seed = seed, threshold_limit = 0)
  a <- g(91)
  expect_identical(a, g(91))
  expect_true(all(a$variance_se > 0))
  expect_true(all(a$relative_information_se > 0))
  expect_true(all(a$evaluated_sets == 400))
  expect_false(isTRUE(all.equal(a$response_variance, g(42)$response_variance)))
  expect_error(zeta_order_information(x, orders = 3, method = "exact", exact_limit=6, threshold_limit=0), "enumeration")
})
