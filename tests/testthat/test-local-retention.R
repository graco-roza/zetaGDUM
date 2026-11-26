test_that("local retention and loss agree with exhaustive conditional zeta", {
  x <- rbind(a = c(1, 1, 0, 0), b = c(4, 0, 2, 0), c = c(1, 1, 0, 0),
             d = c(1, 0, 0, 2), e = c(0, 0, 0, 1))
  n <- nrow(x)
  for (type in c("occurrence", "relative")) {
    q <- if (type == "occurrence") 1L * (x > 0) else x / rowSums(x)
    den <- rowSums(q)
    exhaustive <- vapply(seq_len(n), function(k) {
      vapply(seq_len(n), function(i) {
        sets <- if (k == 1L) list(integer()) else
          utils::combn(setdiff(seq_len(n), i), k - 1L, simplify = FALSE)
        mean(vapply(sets, function(s) sum(apply(q[c(i, s), , drop = FALSE], 2L, min)),
                    numeric(1))) / den[i]
      }, numeric(1))
    }, numeric(n))
    a <- zeta_local_retention(x, type)
    expect_equal(unname(a$profile), exhaustive, tolerance = 1e-13)
    expect_equal(a$summary$T, rowSums(exhaustive), tolerance = 1e-13)
    expect_equal(a$summary$T_full, a$summary$T, tolerance = 1e-13)
    expect_equal(a$summary$T_plus, 1 + (rowSums(exhaustive) - 1) / exhaustive[, 2L],
                 tolerance = 1e-13)
    expect_true(all(a$loss >= 0))
    expect_equal(unname(rowSums(a$loss)), rep(1, n))
    expect_equal(unname(drop(a$loss %*% seq_len(n))), a$summary$T, tolerance = 1e-13)
    b <- zeta_local_retention(x, type, 3)
    expect_equal(b$summary$T, rowSums(exhaustive[, 1:3]), tolerance = 1e-13)
    expect_equal(b$summary$T_full, a$summary$T, tolerance = 1e-13)
    expect_equal(b$loss[, 3L], a$profile[, 3L], tolerance = 1e-13)
  }
})

test_that("retention distinguishes no sharing, empty and right-censored profiles", {
  x <- rbind(a = c(1, 0), b = c(1, 0), c = c(0, 1), d = c(0, 0))
  a <- zeta_local_retention(x)
  expect_identical(a$summary$T[3], 1)
  expect_true(is.na(a$summary$T_plus[3]))
  expect_true(all(is.na(a$profile[4, ])))
  expect_true(all(is.na(a$loss[4, ])))
  expect_error(zeta_local_retention(x, "relative"), "positive finite site totals")
  b <- zeta_local_retention(x[1:2, ], "relative")
  expect_equal(b$summary$T, c(2, 2))
  expect_equal(unname(b$loss), cbind(c(0, 0), c(1, 1)))
  c <- zeta_local_retention(x[1, , drop = FALSE], "relative")
  expect_equal(c$summary$T, 1)
  expect_true(is.na(c$summary$T_plus))
  tiny <- zeta_local_retention(diag(3) + cbind(rep(1e-22, 3), rep(0, 3), rep(0, 3)), "relative")
  expect_equal(tiny$summary$T_plus[2:3], rep(3, 2))
})

test_that("temporal retention uses observed sites and a common comparison domain", {
  a <- rbind(A = c(1, 0), B = c(1, 2), C = c(2, 0), D = c(0, 2))
  b <- rbind(D = c(1, 3), A = c(1, 1), C = c(1, 0))
  occasions <- list(may_a = a, may_b = b, june = a)
  season <- c("May", "May", "June")
  common <- zeta_temporal_retention(occasions, season)
  expect_identical(rownames(common$occasions[[2]]$profile), c("A", "C", "D"))
  expect_true(all(common$seasonal$occasions[common$seasonal$season == "May"] == 2L))
  expect_equal(common$max_order, 3)
  expect_error(zeta_local_retention(rbind(b, E = c(NA, NA))), "Omit unobserved sites")
  expect_error(zeta_temporal_retention(occasions, comparison = "truncated"), "Declare a common")
  truncated <- zeta_temporal_retention(occasions, season, "truncated", "relative", 2L)
  expect_equal(truncated$max_order, 2L)
  expect_equal(truncated$occasions[[2]]$summary$n, rep(3, 3))
  may <- subset(truncated$seasonal, season == "May")
  expect_equal(may$occasions[may$site == "B"], 1L)
  expect_true(is.na(may$T_sd[may$site == "B"]))
  ra <- zeta_local_retention(a, "relative", 2L)
  rb <- zeta_local_retention(b, "relative", 2L)
  av <- c(ra$summary$T[ra$summary$site == "A"], rb$summary$T[rb$summary$site == "A"])
  expect_equal(may$T_mean[may$site == "A"], mean(av))
  expect_equal(may$T_sd[may$site == "A"], sd(av))
  expect_equal(may$T_mean[may$site == "A"],
               sum(subset(truncated$seasonal_profile, site == "A" & season == "May")$zeta_mean))
  loss <- subset(truncated$seasonal_loss, site == "A" & season == "May")
  expect_equal(sum(loss$mass_mean), 1)
  expect_identical(loss$order, c("2", ">2"))
})
