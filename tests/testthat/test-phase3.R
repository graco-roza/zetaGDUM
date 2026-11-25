test_that("conditional zeta and dissimilarity agree with enumeration", {
  x <- rbind(c(3,0,3,0,0,0), c(3,1,2,0,0,0), c(3,0,0,2,0,0),
             c(3,1,0,0,0,0), c(3,0,0,0,4,0), c(3,0,3,0,0,0))
  examples <- list(x, 1 * (x > 0), diag(6), matrix(c(2, 2, 0), 6, 3, byrow = TRUE))
  for (comm in examples) for (kind in c("incidence", "relative")) {
    response <- if (kind == "incidence") 1 * (comm > 0) else comm / rowSums(comm)
    for (k in seq_len(nrow(comm))) {
      idx <- t(combn(nrow(comm), k))
      z <- apply(idx, 1, function(s) sum(apply(response[s, , drop = FALSE], 2, min)))
      d <- if (kind == "incidence") rowMeans(matrix(rowSums(response)[idx], nrow(idx))) - z else 1 - z
      take <- function(a) apply(idx, 1, function(s) all(a %in% s))
      forced <- list(none = integer(), first = 1L, last = tail(seq_len(nrow(comm)), k))
      result <- conditional_zeta(comm, k, kind, forced, pairs = k > 1)
      dis <- conditional_zeta_dissimilarity(comm, k, kind, forced, pairs = k > 1)
      expect_equal(result$mean, mean(z), tolerance = 1e-10)
      expect_equal(mean(result$site), mean(z), tolerance = 1e-10)
      expect_equal(dis$mean, mean(d), tolerance = 1e-10)
      expect_equal(mean(dis$site), mean(d), tolerance = 1e-10)
      expect_equal(result$site, vapply(seq_len(nrow(comm)), function(i) mean(z[take(i)]), 0.0),
                   tolerance = 1e-10)
      expect_equal(dis$site, vapply(seq_len(nrow(comm)), function(i) mean(d[take(i)]), 0.0),
                   tolerance = 1e-10)
      expect_equal(result$forced, vapply(forced, function(a) mean(z[take(a)]), 0.0), tolerance = 1e-10)
      expect_equal(dis$forced, vapply(forced, function(a) mean(d[take(a)]), 0.0), tolerance = 1e-10)
      if (k > 1) {
        pair <- combn(nrow(comm), 2)
        ij <- t(pair)
        expect_true(all(is.na(diag(dis$pair))))
        expect_equal(result$pair, t(result$pair))
        expect_equal(dis$pair[ij], apply(pair, 2, function(a) mean(d[take(a)])), tolerance = 1e-10)
        expect_equal(result$pair[ij], apply(pair, 2, function(a) mean(z[take(a)])), tolerance = 1e-10)
        expect_equal(mean(result$pair[ij]), mean(z), tolerance = 1e-10)
      }
      lc <- zeta_lczd(comm, k, kind)
      expect_equal(lc$expected_dissimilarity, dis$site, tolerance = 1e-10)
      if (mean(d) > 1e-12) expect_equal(sum(lc$LCZD), 1, tolerance = 1e-10)
      else expect_true(all(is.na(lc$LCZD)))
    }
  }
  expect_error(conditional_zeta(rbind(x, 0), 3, "relative"), "positive")
  expect_error(conditional_zeta(x, 3, include = c(1, 1)), "distinct")
})

test_that("sampled normalised incidence averages set ratios", {
  x <- rbind(c(1,0,0), c(1,1,1), c(0,1,0), c(0,0,0))
  idx <- t(combn(4, 2))
  w <- seq_len(nrow(idx))
  r <- ugdm_response(x, idx, "total")
  lc <- zeta_lczd(x, 2, "normalised", idx, w)
  expect_equal(lc$expected_dissimilarity, vapply(seq_len(4), function(i) {
    has <- idx[,1] == i | idx[,2] == i
    weighted.mean(r$y[has], w[has])
  }, 0.0))
  expect_false(isTRUE(all.equal(lc$expected_dissimilarity,
    conditional_zeta_dissimilarity(x, 2)$site / mean(rowSums(x)))))
})

test_that("analytical GDUM moments use the stored environmental basis and fit the same full objective", {
  set.seed(103)
  n <- 12L
  env <- data.frame(warm = rnorm(n), rain = rnorm(n))
  x <- cbind(runif(n, 4, 8), 7, matrix(sample(c(0, 0, 1:5), n * 10, TRUE), n))
  idx <- t(combn(n, 4))
  geo <- as.matrix(dist(matrix(rnorm(2 * n), n)))
  for (distance in list(NULL, geo)) {
    design <- ugdm_design(env, idx, dist_matrix = distance)
    for (kind in c("incidence", "relative")) {
      r <- ugdm_response(x, idx, if (kind == "relative") kind else "total")
      ref <- ugdm_fit(design, r$y, r$trials)
      X <- zetaGDUM:::.order_terms(design, ref)
      m <- ugdm_crossproducts(x, env, design, kind)
      exact <- !is.na(m$response)
      expect_equal(m$response[exact], colMeans(X * r$y * r$trials)[exact], tolerance = 1e-10)
      expect_equal(m$predictor[exact], colMeans(X * r$trials)[exact], tolerance = 1e-10)
      expect_identical(unname(!exact), grepl("geo|dist", names(m$response)))
      for (method in c("hybrid", "control_variate")) {
        f <- ugdm_fit_analytical(x, env, design, kind, method = method)
        expect_identical(f$analytical$used, method, info = f$analytical$reason)
        expect_equal(f$deviance, ref$deviance, tolerance = 1e-7)
        expect_lt(max(abs(ugdm_predict(f, design) - ugdm_predict(ref, design))), 1e-4)
      }
      local <- ugdm_local_distinctiveness(x, ref, design, if (kind == "relative") kind else "normalised")
      expect_equal(local$spread + local$position + local$baseline, local$predicted)
      expect_equal(local$observed - local$predicted, local$residual)
      expect_equal(mean(local$predicted), mean(ugdm_predict(ref, design)))
      m$response[1L] <- 3 * m$predictor[1L]
      f <- ugdm_fit_analytical(x, env, design, kind, moments = m)
      expect_identical(f$analytical$used, "conventional")
      expect_error(ugdm_fit_analytical(x, env, design, kind, moments = m, fallback = FALSE),
                   "certificate")
    }
  }
})

test_that("balanced designs and matched groups retain distinct sites", {
  p <- sample_sets_balanced(7, 3, 42, "pair", seed = 10)
  count <- crossprod(1 * outer(seq_len(nrow(p)), seq_len(7), Vectorize(function(i, j) j %in% p[i, ])))
  expect_true(all(diag(count) == 18L))
  expect_true(all(count[lower.tri(count)] == 6L))
  pp <- sample_sets_balanced(8, 2, 28, "pair", seed = 1)
  expect_equal(anyDuplicated(pp, MARGIN = 1), 0L)
  expect_true(all(pp[, 1] < pp[, 2]))
  expect_error(sample_sets_balanced(8, 3, 24, "pair"), "prime")
  matched <- sample_sets_matched(12, 3, 40, seed = 2)
  expect_true(all(vapply(split(seq_len(nrow(matched$idx)), matched$group), function(i) {
    anyDuplicated(as.vector(matched$idx[i, ])) == 0L
  }, TRUE)))
  for (method in c("site", "environment")) {
    p <- sample_sets_balanced(9, 4, 100, method, strata = rep(1:3, 3), seed = 98)
    expect_true(all(apply(p, 1, anyDuplicated) == 0L))
  }
})

test_that("matched and overlap covariance preserve signed estimates and sandwich algebra", {
  matched <- sample_sets_matched(20, 4, 25, seed = 84)
  scores <- cbind(sin(seq_len(nrow(matched$idx))), cos(seq_len(nrow(matched$idx))))
  fit <- u_matched_covariance(scores, matched, diag(c(2, 4)))
  n <- nrow(scores)
  m <- 5L
  expect_equal(fit$incomplete, fit$within - (n - 1) / n * fit$total)
  expect_equal(fit$complete, (fit$within - fit$total) / (1 - (m - 1) / (n - 1)))
  expect_equal(fit$sandwich_complete, diag(c(.5, .25)) %*% fit$complete %*% diag(c(.5, .25)))
  invalid <- matched
  invalid$idx[1, 1] <- invalid$idx[2, 1]
  expect_error(u_matched_covariance(scores, invalid), "disjoint")
  evaluate <- function(idx) cbind(rowMeans(idx), rowMeans(idx^2))
  fit <- u_overlap_covariance(10, 4, evaluate, 500, hessian = diag(2), seed = 3)
  expect_equal(sum(fit$probability), 1)
  expect_equal(fit$probability, dhyper(0:4, 4, 6, 4))
  expect_equal(fit$contrast[[1]], matrix(0, 2, 2))
  expect_equal(fit$complete, Reduce("+", Map("*", fit$probability, fit$contrast)))
  expect_equal(fit$sandwich, fit$complete)
  expect_error(u_overlap_covariance(9, 5, evaluate), "integer")
})

test_that("raw proposal correction recovers the exhaustive weighted loss on a replicated design", {
  set.seed(51)
  n <- 11
  env <- data.frame(g = rnorm(n))
  comm <- cbind(4, matrix(rpois(n * 8, rlnorm(n * 8)), n))
  all <- t(combn(n, 3))
  r <- ugdm_response(comm, all, "total")
  d <- ugdm_design(env, all)
  target <- ugdm_fit(d, r$y, r$trials, free_sign = TRUE)
  reps <- 1L + (apply(all, 1, sum) %% 5L)
  subset <- rep(seq_len(nrow(all)), reps)
  w <- length(subset) / (nrow(all) * reps[subset])
  d_subset <- ugdm_design(env, all[subset, ], ref = d)
  for (method in c("hybrid", "control_variate")) {
    h <- ugdm_fit_analytical(comm, env, d_subset, "incidence", w, free_sign = TRUE,
                             method = method, fallback = FALSE)
    expect_lt(max(abs(ugdm_predict(h, d) - ugdm_predict(target, d))), 1e-4)
  }
})
