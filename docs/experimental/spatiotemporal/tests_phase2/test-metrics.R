test_that("exact species contributions reproduce native results at every spatial order", {
  set.seed(74)
  a <- array(rpois(8 * 4 * 6, 2), c(8, 4, 6))
  for (t in 1:4) {
    x <- a[, t, ]
    v <- st2_zeta(a, times = t)
    for (type in c("presence", "abundance", "dominance")) {
      surface <- subset(v$surface, response == type)
      expect_equal(surface$zeta, zetaGDUM::zeta_curve(x, 8, type)$zeta)
    }
    for (k in 1:8) {
      z <- zetaGDUM::zeta_decompose(x, k)
      expect_equal(unname(v$contributions$presence[k, ]), z$P)
      expect_equal(unname(v$contributions$abundance[k, ]), z$term_abundance)
      expect_equal(unname(v$conditional_minimum[k, ]), z$m)
      expect_equal(unname(v$weights$presence[k, ]), z$w_occ)
      expect_equal(unname(v$weights$abundance[k, ]), z$w_abu)
    }
  }
})

test_that("parallel surfaces match enumeration and obey invariants", {
  set.seed(1)
  a <- array(rpois(5 * 4 * 4, 2), c(5, 4, 4))
  a[1, 1, 1] <- 1
  for (mode in c("all", "consecutive")) {
    v <- st2_zeta(a, mode = mode)
    for (m in 1:4) for (k in 1:5) {
      ts <- .st2_time_sets(1:4, m, mode, 100)
      xs <- t(combn(5, k))
      all <- lapply(seq_len(nrow(xs)), function(i) {
        lapply(seq_len(nrow(ts)), function(j) st2_zeta_rectangle(a, xs[i, ], ts[j, ])[2:4])
      })
      all <- do.call(rbind, unlist(all, recursive = FALSE))
      for (type in names(v$contributions)) {
        expect_equal(sum(all[[type]]) / (nrow(xs) * nrow(ts)),
                     subset(v$surface, response == type & spatial_order == k & temporal_order == m)$zeta)
      }
    }
    expect_true(all(v$surface$temporal_retention <= 1 + 1e-12, na.rm = TRUE))
    expect_true(all(v$surface$spatial_retention <= 1 + 1e-12, na.rm = TRUE))
    expect_true(all(v$contributions$dominance <= v$contributions$presence))
  }
  b <- sweep(a, c(1, 2), matrix(rexp(20), 5, 4), "*")
  expect_equal(st2_zeta(a)$surface$zeta, st2_zeta(b)$surface$zeta)
  b <- a * array(exp(rnorm(length(a))), dim(a))
  expect_equal(st2_zeta(a)$contributions$presence, st2_zeta(b)$contributions$presence)
  a[1, 1, ] <- 0
  expect_error(st2_zeta(a), "undefined")
  expect_true(all(is.na(st2_zeta(a, empty = "undefined")$contributions$abundance)))
  z <- array(1, c(1, 3, 1))
  expect_equal(st2_zeta(z)$surface$zeta, rep(1, 9))
})

test_that("the directional species ledger satisfies the exact accounting identities", {
  a <- array(c(1, 3, 2, 4, 4, 0, 0, 8, 0, 2, 3, 1, 4, 4, 3, 8), c(2, 2, 4))
  for (sites in list(1L, 2L, 1:2)) {
    z <- st2_directional(a, sites, 1, 2)
    for (r in c("presence", "dominance")) {
      q <- z[[r]]
      expect_equal(q[["before"]], q[["P"]] + q[["L"]])
      expect_equal(q[["after"]], q[["P"]] + q[["G"]])
      expect_equal(q[["D_temporal"]], q[["D_gain"]] + q[["D_loss"]])
    }
    q <- z$abundance
    expect_equal(q[["change"]], q[["abundance_gain"]] - q[["abundance_loss"]])
    expect_equal(q[["change"]], q[["abundance_gain_occ"]] - q[["abundance_loss_occ"]] +
                   q[["abundance_gain_num"]] - q[["abundance_loss_num"]])
    expect_equal(z$dominance[["G"]], z$dominance[["dominance_gain_occ"]] + z$dominance[["dominance_gain_num"]])
    expect_equal(z$dominance[["L"]], z$dominance[["dominance_loss_occ"]] + z$dominance[["dominance_loss_num"]])
  }
  a[1, 1, ] <- 0
  expect_true(is.na(st2_directional(a, 1:2, 1, 2)$abundance[["change"]]))
})

test_that("occupancy budgets are exact and regional status uses the complete panel", {
  a <- array(0, c(4, 2, 2))
  a[1:2, 1, 1] <- 1
  a[2:4, 2, 1] <- 1
  a[1, 2, 2] <- 1
  a[4, 1, 2] <- 1
  q <- st2_occupancy(a, 1:3)$combination
  expect_equal(q$r1, q$r0 + q$gains - q$losses)
  expect_equal(q$gains + q$losses, q$expansion + q$contraction + 2 * q$redistribution)
  expect_true(q$pure_redistribution[1])
  expect_false(q$regional_appearance[2])
  full <- st2_occupancy(a)$combination
  expect_equal(full$expansion[1], 1)
  expect_equal(full$redistribution[1], 1)
})

test_that("Poisson-binomial joint probabilities agree with enumeration", {
  g <- c(.2, .7)
  l <- c(.1, .6)
  r <- st2_occupancy_probability(g, l)
  outcomes <- as.matrix(expand.grid(rep(list(0:1), 4)))
  probability <- apply(outcomes, 1, function(z) prod(ifelse(z == 1, c(g,l), 1-c(g,l))))
  a <- rowSums(outcomes[, 1:2])
  b <- rowSums(outcomes[, 3:4])
  expect_equal(r$expansion, sum(probability[a>b]))
  expect_equal(r$contraction, sum(probability[a<b]))
  expect_equal(r$pure_redistribution, sum(probability[a==b & a>0]))
  expect_equal(r$redistribution, sum(probability*pmin(a,b)))
  expect_equal(sum(r$joint_mass), 1)
  expect_equal(st2_poisson_binomial(numeric()), 1)
})
