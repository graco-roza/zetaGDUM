test_that("species terms exactly reconstruct incidence and abundance responses", {
  x <- rbind(c(2, 0, 1), c(0, 3, 1), c(2, 3, 1))
  for (idx in list(rbind(c(1L, 2L), c(1L, 3L)), matrix(c(1L, 2L, 3L), nrow = 1)))
  for (type in c("presence", "relative", "absolute")) {
    z <- species_contributions(x, idx, type)
    expect_true(all(z$defined))
    expect_true(all(z$dissimilarity >= -1e-12))
    expect_equal(rowSums(z$similarity) + z$response, rep(1, nrow(idx)))
    if (type == "presence") {
      s <- zeta_sets(x, idx, abundance = FALSE)
      expect_equal(z$response, 1 - s$total)
    } else if (type == "absolute") {
      s <- zeta_sets(x, idx, abundance = TRUE)
      expect_equal(z$response, 1 - s$bray_curtis)
    } else {
      q <- x / rowSums(x)
      brute <- apply(idx, 1, function(i) 1 - sum(apply(q[i, , drop = FALSE], 2, min)))
      expect_equal(z$response, brute)
    }
  }
  # A species absent from some sites contributes no similarity but can cause dissimilarity.
  pair <- matrix(c(1L, 2L), nrow = 1)
  expect_equal(species_contributions(x, pair, "presence")$similarity[1, 1], 0)
  expect_gt(species_contributions(x, pair, "presence")$dissimilarity[1, 1], 0)
})

test_that("relative and absolute responses answer different abundance questions", {
  x <- rbind(c(2, 2), c(20, 20))
  idx <- matrix(c(1L, 2L), nrow = 1)
  expect_equal(species_contributions(x, idx, "relative")$response, 0)
  expect_gt(species_contributions(x, idx, "absolute")$response, 0)
  r <- ugdm_response(x, idx, "relative")
  expect_equal(r$y, 0)
  expect_equal(r$trials, 1)
  expect_true(r$keep)
})

test_that("relative-abundance zeta rejects empty sites", {
  x <- rbind(c(1, 2), c(0, 0), c(2, 1))
  expect_error(zeta_curve(x, kmax = 2, type = "abundance"), "undefined")
  expect_error(zeta_weights(x, 2, abundance = TRUE), "undefined")
  expect_error(zeta_decompose(x, 2), "undefined")
  z <- species_contributions(x, matrix(c(1L, 2L), nrow = 1), "relative")
  expect_false(z$defined)
  expect_true(is.na(z$response))
})

test_that("relative species contributions are invariant to each site's positive total", {
  x <- rbind(c(1, 2, 0), c(2, 0, 1), c(1, 1, 3))
  scaled <- x * c(0.04, 2, 0.001)
  for (idx in list(rbind(c(1L, 2L), c(2L, 3L)), matrix(1:3, nrow = 1L))) {
    expected <- species_contributions(x, idx, "relative")
    expect_equal(species_contributions(scaled, idx, "relative"), expected)
    expect_equal(ugdm_response(scaled, idx, "relative")$y, expected$response)
  }
})

test_that("exact expected contributions match exhaustive enumeration at every order", {
  x <- cbind(a = c(4, 1, 3, 1, 0), b = c(1, 1, 0, 0, 2), c = c(0, 3, 0, 0, 0), d = 0)
  for (k in 1:nrow(x)) {
    idx <- t(utils::combn(nrow(x), k))
    for (type in c("presence", "relative")) {
      q <- if (type == "presence") (x > 0) * 1 else x / rowSums(x)
      means <- t(apply(idx, 1, function(i) colMeans(q[i, , drop = FALSE])))
      mins <- t(apply(idx, 1, function(i) apply(q[i, , drop = FALSE], 2, min)))
      analytical <- expected_species_contributions(x, k, type)
      expect_equal(analytical$dissimilarity, unname(colMeans(means - mins)), tolerance = 1e-12)
      expect_equal(analytical$expected_minimum, unname(colMeans(mins)), tolerance = 1e-12)
      expect_equal(analytical$sharing_probability, unname(colMeans(mins > 0)), tolerance = 1e-12)
      if (k > 1 && type == "relative") {
        expect_equal(analytical$dissimilarity,
                     unname(colMeans(species_contributions(x, idx, type)$dissimilarity)), tolerance = 1e-12)
      }
    }
  }
  idx <- t(utils::combn(nrow(x), 2))
  normalised <- colMeans(species_contributions(x, idx, "presence")$dissimilarity)
  expect_gt(max(abs(normalised - expected_species_contributions(x, 2, "presence")$dissimilarity)), 0.01)
})
