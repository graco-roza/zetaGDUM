testthat::test_that("rectangular sharing is exact and agrees with spatial zeta", {
  s <- simulate_st_gdum("interaction", n_sites = 10, n_species = 12, seed = 3)
  p <- s$latent
  d <- st_zeta_retention(st_zeta_decline(p, 4, 5))
  for (k in 1:4) for (m in 1:5) {
    z <- st_zeta(p, k, m)
    i <- d$spatial_order == k & d$temporal_order == m
    testthat::expect_equal(d$zeta[i], mean(z$zeta))
    testthat::expect_equal(z$zeta, lengths(z$species))
  }
  for (t in seq_along(p$times)) {
    native <- zetaGDUM::zeta_curve(.st_slice(p, t), 4)
    actual <- vapply(1:4, function(k) mean(st_zeta(p, k, 1, time_sets = matrix(t))$zeta), numeric(1))
    testthat::expect_equal(actual, native$zeta)
  }
  testthat::expect_true(all(d$spatial_retention <= 1 + 1e-12, na.rm = TRUE))
  testthat::expect_true(all(d$temporal_retention <= 1 + 1e-12, na.rm = TRUE))
  one <- array(c(1, 0, 1), c(1, 3, 1))
  testthat::expect_equal(st_zeta_decline(one)$zeta, c(2/3, 1/3, 0))
  testthat::expect_equal(st_zeta(one, 1, 3)$zeta, 0L)
  testthat::expect_error(st_zeta(p, 4, 4, limit = 40), "budget")
})

testthat::test_that("endpoint events and sampled recolonisation retain species identity", {
  a <- array(0, c(2, 3, 3), list(c("A", "B"), c("a", "b", "c"), c("x", "y", "z")))
  a[, c(1, 3), "x"] <- 1
  a[, 1:2, "y"] <- 1
  a[, 2:3, "z"] <- 1
  d <- st_decompose(st_panel(a, c(0, 1, 5)), 2)
  testthat::expect_equal(d$persistent_species[[1]], "x")
  testthat::expect_equal(d$recolonised_species[[1]], "x")
  testthat::expect_equal(d$interrupted_species[[1]], "x")
  testthat::expect_equal(d$disappeared_during_species[[1]], c("x", "y"))
  testthat::expect_equal(d$gained_during_species[[1]], c("x", "z"))
  testthat::expect_equal(d$continuous[[1]], 0)
  testthat::expect_equal(d$gains_species[[1]], "z")
  testthat::expect_equal(d$losses_species[[1]], "y")
  testthat::expect_identical(d$before, d$persistent + d$losses)
  testthat::expect_identical(d$after, d$persistent + d$gains)
  testthat::expect_identical(d$after - d$before, d$change)
  testthat::expect_equal(st_decompose(a, 1, c(1, 3))$continuous, c(1, 1))
  testthat::expect_error(st_decompose(a, 1, 3:1), "increasing")
})

testthat::test_that("analytical occupancy contributions reproduce enumerated events", {
  sim <- simulate_st_gdum(c("expansion", "redistribution", "abundance_only"),
                         n_sites = 8, n_species = 12, seed = 2)
  x <- sim$latent
  o <- st_occupancy(x, 1, 5)
  for (k in 1:8) {
    d <- st_decompose(x, k, c(1, 5))
    for (event in c("before", "after", "persistent", "gains", "losses", "change")) {
      testthat::expect_equal(mean(d[[event]]), sum(o$contributions[[event]][o$contributions$spatial_order == k]))
    }
  }
  sp <- o$species
  testthat::expect_true(any(sp$expansion & sp$redistribution))
  for (suffix in c("before", "after")) {
    occ <- sp[[paste0("occupancy_", suffix)]]
    conditional <- sp[[paste0("conditional_", suffix)]]
    testthat::expect_equal(sp[[paste0("total_", suffix)]], ifelse(occ == 0, 0, occ * conditional))
  }
  testthat::expect_equal(sp$total_change, sp$abundance_gained_sites - sp$abundance_lost_sites +
                          sp$abundance_change_persistent_sites)
  z <- st_occupancy(array(c(0, 3, 4, 0, 0, 0), c(1, 2, 3)))$species
  testthat::expect_equal(z$regional_appearance, c(TRUE, FALSE, FALSE))
  testthat::expect_equal(z$regional_disappearance, c(FALSE, TRUE, FALSE))
  testthat::expect_true(is.na(z$conditional_before[1]))
})

testthat::test_that("pair and established multisite beta obey partition properties", {
  s <- simulate_st_gdum("homogenisation", n_sites = 8, n_species = 12, seed = 2)
  p <- s$latent
  pair <- st_beta(p)
  ok <- is.finite(pair$nestedness)
  testthat::expect_equal(pair$sorensen[ok], pair$simpson[ok] + pair$nestedness[ok])
  for (time in seq_along(p$times)) {
    beta <- st_beta_multisite(p, c(1, 4), time)
    actual <- pair[pair$time == p$times[time] & pair$i == 1 & pair$j == 4, names(beta)]
    testthat::expect_equal(as.numeric(beta), unname(as.numeric(actual)))
  }
  testthat::expect_equal(st_beta_multisite(p), st_beta_multisite(p, 8:1))
  x <- array(c(1, 1, 1, 0, 1, 1, 0, 0, 1), c(3, 1, 3))
  testthat::expect_equal(unname(st_beta_multisite(x)["simpson"]), 0)
  x <- array(c(1, 1, 0, 1, 0, 1, 0, 1, 1), c(3, 1, 3))
  testthat::expect_equal(unname(st_beta_multisite(x)["nestedness"]), 0)
  testthat::expect_equal(unname(st_beta_multisite(x)["sorensen"]), 0.5)
  testthat::expect_true(all(is.na(st_beta_multisite(array(0, c(3, 1, 2))))))
  change <- st_beta_change(p)
  ok <- is.finite(change$gain_contribution + change$loss_contribution)
  testthat::expect_equal(change$change[ok], (change$gain_contribution + change$loss_contribution)[ok])
  if (requireNamespace("betapart", quietly = TRUE)) {
    ref <- betapart::beta.multi((.st_slice(p, 1) > 0) * 1L)
    testthat::expect_equal(unname(st_beta_multisite(p)), unname(unlist(ref[c("beta.SOR", "beta.SIM", "beta.SNE")])))
  }
})

testthat::test_that("spatial-only fitting exactly preserves native response and basis", {
  sim <- simulate_st_gdum("spatial_filtering", n_sites = 16, seed = 4)
  p <- sim$latent$x[, 1, , drop = FALSE]
  env <- sim$env[, 1, , drop = FALSE]
  idx <- t(combn(16, 2))
  for (component in c("turnover", "total", "nestedness", "relative", "balanced", "gradient")) {
    native <- zetaGDUM::ugdm_response(.st_slice(st_panel(p), 1), idx, component)
    design <- zetaGDUM::ugdm_design(.st_flat_env(env), idx[native$keep, , drop = FALSE])
    fit <- zetaGDUM::ugdm_fit(design, native$y[native$keep], native$trials[native$keep])
    experimental <- st_spatial_gdum(p, env, component = component)[[1]]
    testthat::expect_equal(experimental$response, native)
    testthat::expect_equal(experimental$fit, fit)
  }
})

testthat::test_that("conditional likelihood recovers effects and reports structural failures", {
  sim <- simulate_st_gdum("interaction", n_sites = 54, n_species = 100, seed = 14)
  x <- sim$latent
  d <- st_design(x, sim$env)
  terms <- c("position_initial.gradient", "position_change.gradient", "position_by_change.gradient")
  f <- st_gdum(x, d, terms, loss_terms = terms[1], blocks = seq_len(54))
  testthat::expect_equal(f$fit$gain$status, "estimated")
  expected <- 1.25 * d$ref$scale^2
  testthat::expect_lt(abs(f$fit$gain$coefficients[4] - expected), 0.35)
  testthat::expect_true(all(is.finite(f$fit$gain$jackknife_se)))
  y <- st_predict(f)
  testthat::expect_equal(y$before, y$persistent + y$losses)
  testthat::expect_equal(y$after, y$persistent + y$gains)
  testthat::expect_equal(f$counts$gain_opportunities + f$counts$loss_opportunities, rep(100, nrow(d$X)))
  bad <- st_gdum(x, d, c("position_initial.gradient", "position_final.gradient", "position_change.gradient"))
  testthat::expect_equal(bad$fit$gain$status, "aliased")
  testthat::expect_equal(st_gdum(x, d, "lag")$fit$gain$status, "aliased")
  sim <- simulate_st_gdum("gains_only", seed = 14)
  f <- st_gdum(sim$latent, st_design(sim$latent, sim$env))
  testthat::expect_equal(f$fit$loss$status, "boundary")
  testthat::expect_true(all(st_predict(f)$losses == 0))
  sim$latent$x[] <- 0
  empty <- st_gdum(sim$latent, st_design(sim$latent, sim$env))
  testthat::expect_equal(empty$fit$loss$status, "no_opportunities")
  testthat::expect_true(all(st_predict(empty)$persistent == 0))
})

testthat::test_that("typed panels handle incomplete sites and effort explicitly", {
  s <- simulate_st_gdum("null", missing = 0.15, detection = 0.8, variable_effort = TRUE, seed = 91)
  testthat::expect_error(st_panel(s$observed), "Balanced")
  p <- st_panel(s$observed, s$latent$times, missing = "complete_sites")
  testthat::expect_equal(length(p$dropped_sites), sum(rowSums(!s$observed_sites) > 0))
  testthat::expect_true(all(is.finite(st_design(p, s$env)$X)))
  testthat::expect_error(st_panel(s$observed, 5:1), "strictly increasing")
  testthat::expect_equal(dim(s$latent$x), dim(s$observed))
  s <- simulate_st_gdum("abundance_only", abundance = "gamma", seed = 4)
  d <- st_design(s$latent, s$env)
  a <- st_abundance_gdum(s$latent, d, "position_change.gradient", blocks = seq_len(30))
  testthat::expect_lt(abs(a$fit$coefficients[2] / d$ref$scale - 0.65), 0.07)
  paired <- st_abundance_gdum(s$latent, d, "position_change.gradient", response = "log_ratio")
  testthat::expect_true(all(paired$data$initial > 0 & paired$data$final > 0))
})

testthat::test_that("order diagnostics expose vanishing support on both axes", {
  x <- array(0, c(3, 2, 1))
  x[1, 1, 1] <- 1
  d <- st_order_diagnostic(x, 3, 2, matrix(1:6, 3, 2), budget = 10)
  testthat::expect_equal(d$zero_fraction, c(5/6, 1, 1, 1, 1, 1))
  testthat::expect_equal(d$relative_position_information, c(1, 0, 0, 0, 0, 0))
  testthat::expect_true(all(is.na(d$temporal_retention[d$temporal_order == 1])))
})

testthat::test_that("held-out order-two predictions reuse the training basis", {
  sim <- simulate_st_gdum("interaction", n_sites = 18, seed = 77)
  train <- sim$latent
  test <- sim$latent
  train$x <- train$x[1:12, , , drop = FALSE]
  test$x <- test$x[13:18, , , drop = FALSE]
  a <- st_design(train, sim$env, 2)
  b <- st_design(test, sim$env, 2, ref = a$ref)
  f <- st_gdum(train, a, "position_change.gradient")
  testthat::expect_identical(a$ref, b$ref)
  testthat::expect_error(st_predict(f, b), "baseline counts")
  actual <- st_gdum(test, b)$counts
  pred <- st_predict(f, b, actual$before)
  testthat::expect_true(all(is.finite(pred$after)))
  testthat::expect_equal(pred$before, pred$persistent + pred$losses)
})
