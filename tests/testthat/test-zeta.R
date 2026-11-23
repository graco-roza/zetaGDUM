set.seed(11)
x <- matrix(rpois(8 * 7, 1.2), 8, 7)
x[1, ] <- x[1, ] + 1L
P <- (x > 0) * 1

test_that("presence zeta equals brute-force enumeration", {
  cur <- zeta_curve(x, kmax = 6)
  for (k in 1:6) {
    cm <- utils::combn(nrow(P), k)
    brute <- mean(apply(cm, 2, function(i) sum(colSums(P[i, , drop = FALSE]) == k)))
    expect_equal(cur$zeta[k], brute, tolerance = 1e-10)
  }
})

test_that("abundance zeta equals brute-force enumeration of summed minimum relative abundance", {
  Q <- x / rowSums(x)
  cur <- zeta_curve(x, kmax = 6, type = "abundance")
  for (k in 1:6) {
    cm <- utils::combn(nrow(Q), k)
    brute <- mean(apply(cm, 2, function(i) sum(apply(Q[i, , drop = FALSE], 2, min))))
    expect_equal(cur$zeta[k], brute, tolerance = 1e-10)
  }
})

test_that("species weights sum to one, are zero above occupancy, and match enumeration", {
  o <- colSums(P > 0)
  for (k in 2:5) {
    w <- zeta_weights(x, k)
    expect_equal(sum(w), 1, tolerance = 1e-10)
    expect_true(all(w[o < k] == 0))
    cm <- utils::combn(nrow(P), k)
    brute <- rowMeans(apply(cm, 2, function(i) as.numeric(colSums(P[i, , drop = FALSE]) == k)))
    expect_equal(w, brute / sum(brute), tolerance = 1e-10)
    wa <- zeta_weights(x, k, abundance = TRUE)
    Q <- x / rowSums(x)
    bruteA <- rowMeans(apply(cm, 2, function(i) apply(Q[i, , drop = FALSE], 2, min)))
    expect_equal(wa, bruteA / sum(bruteA), tolerance = 1e-10)
  }
})

test_that("dominance matrix marks the smallest set reaching the cumulative share", {
  D <- dominance_matrix(x, 0.8)
  for (i in 1:nrow(x)) {
    v <- x[i, ]
    sel <- D[i, ] == 1
    expect_gte(sum(v[sel]) / sum(v), 0.8 - 1e-12)
    if (sum(sel) > 1) expect_lt((sum(v[sel]) - min(v[sel])) / sum(v), 0.8 + 1e-12)
  }
  expect_equal(dominance_matrix(x, 1)[x > 0], rep(1L, sum(x > 0)))
})

test_that("zeta_break returns an order within the curve", {
  set.seed(2)
  big <- matrix(rbinom(3000, 1, rep(c(0.05, 0.3, 0.9), each = 1000)), 100, 30)
  kb <- zeta_break(zeta_curve(big))
  expect_true(is.integer(kb) && kb >= 2 && kb <= 80)
})

test_that("set components satisfy the partition identities", {
  idx <- sample_sets(nrow(x), 3, M = 40, seed = 1)
  z <- zeta_sets(x, idx)
  expect_equal(1 - z$total, (1 - z$turnover) + z$nestedness, tolerance = 1e-12)
  expect_equal(z$balanced - z$bray_curtis, z$gradient, tolerance = 1e-12)
  expect_true(all(z$nestedness >= -1e-12, na.rm = TRUE))
})

test_that("pairs reproduce the Baselga and Podani-Legendre components", {
  pp <- pair_partition(x)
  idx <- as.matrix(pp[, c("i", "j")])
  storage.mode(idx) <- "integer"
  z <- zeta_sets(x, idx)
  expect_equal(1 - z$turnover, pp$simpson, tolerance = 1e-12)
  expect_equal(1 - z$total, pp$sorensen, tolerance = 1e-12)
  expect_equal(pp$sorensen - pp$simpson, pp$nestedness, tolerance = 1e-12)
  expect_equal(pp$replacement + pp$richness_difference, pp$jaccard, tolerance = 1e-12)
})

test_that("abundance components match the pairwise Bray-Curtis partition", {
  idx <- matrix(c(1L, 2L), 1)
  z <- zeta_sets(x, idx, abundance = TRUE)
  A <- sum(pmin(x[1, ], x[2, ])); B <- sum(pmax(x[1, ] - x[2, ], 0)); C <- sum(pmax(x[2, ] - x[1, ], 0))
  expect_equal(z$balanced, A / (A + min(B, C)), tolerance = 1e-12)
  expect_equal(z$bray_curtis, 2 * A / (2 * A + B + C), tolerance = 1e-12)
})
