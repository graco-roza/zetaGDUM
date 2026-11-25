source("docs/phase3/common.R")
source("docs/manuscript_validation/common.R")
set.seed(34912)

rho <- function(x, y) if (sd(x) > 0 && sd(y) > 0) cor(x, y, method = "spearman") else NA_real_
local_summary <- function(x, env, orders, label, source, reference = NULL) {
  n <- nrow(x)
  q <- x / rowSums(x)
  h <- sqrt(q)
  lcbd <- rowSums(sweep(h, 2L, colMeans(h))^2)
  lcbd <- lcbd / sum(lcbd)
  idx <- t(combn(n, 2))
  d <- suppressWarnings(ugdm_design(env, idx))
  r <- ugdm_response(x, idx, "relative")
  f <- ugdm_fit(d, r$y, r$trials)
  position <- drop(scale(env, d$environment_basis$center, d$environment_basis$scale) %*% f$lambda)
  local <- ugdm_local_distinctiveness(x, f, d)
  from_pair <- zetaGDUM:::.conditional_sample_mean(idx, r$y, rep(1, nrow(idx)), n)[, 1]
  rarity <- drop(q %*% (1 / pmax(colSums(x > 0), 1)))
  lczd <- zeta_lczd(x, orders)
  rows <- lapply(orders, function(k) {
    a <- lczd[lczd$order == k, ]
    raw <- zeta_lczd(x, k, "incidence")
    base <- if (is.null(reference)) rep(NA_real_, n) else
      reference$expected_dissimilarity[reference$order == k]
    data.frame(source, label, n, species = ncol(x), k, density = mean(x > 0),
      mean_relative_D = mean(a$expected_dissimilarity), relative_site_sd = sd(a$expected_dissimilarity),
      mean_raw_incidence_D = mean(raw$expected_dissimilarity), LCZD_sum = sum(a$LCZD),
      pair_identity_error = if (k == 2L) max(abs(from_pair - a$expected_dissimilarity)) else NA_real_,
      rho_LCBD = rho(a$expected_dissimilarity, lcbd),
      rho_mean_pair = rho(a$expected_dissimilarity, from_pair),
      rho_GDUM_position = rho(a$expected_dissimilarity, position),
      rho_GDUM_predicted = rho(a$expected_dissimilarity, local$predicted),
      rho_richness = rho(a$expected_dissimilarity, rowSums(x > 0)),
      rho_raw_incidence_richness = rho(raw$expected_dissimilarity, rowSums(x > 0)),
      rho_rare_abundance = rho(a$expected_dissimilarity, rarity),
      rho_totals = rho(a$expected_dissimilarity, rowSums(x)),
      rho_reference = if (is.null(reference)) NA_real_ else rho(a$expected_dissimilarity, base),
      max_change_reference = if (is.null(reference)) NA_real_ else max(abs(a$expected_dissimilarity - base)),
      pair_GDUM_D2 = f$D2, pair_GDUM_converged = f$converged)
  })
  do.call(rbind, rows)
}
sim <- simulate_metacommunity(72, 53, breadth = .65, seed = 4182)
x <- cbind(4L, sim$comm)
orders <- c(2L, 4L, 7L)
reference <- zeta_lczd(x, orders)
effort <- rep(1:6, length.out = nrow(x))
down <- matrix(rbinom(length(x), c(x), rep(seq(.12, .95, length.out = nrow(x)), ncol(x))), nrow(x))
down[, 1] <- pmax(1L, down[, 1])
rare <- diag(rep(0:4, length.out = nrow(x)))
relative_reweight <- sweep(x, 2L, seq(.5, 2, length.out = ncol(x)), "*")
cases <- list(base = x, count_scaling = x * effort, uneven_detection = down,
              rare_singletons = cbind(x, rare), taxon_reweight = relative_reweight)
out <- lapply(names(cases), function(label) {
  local_summary(cases[[label]], sim$env, orders, label, "simulation", reference)
})
write_result(do.call(rbind, out), "local_simulated.csv")

args <- commandArgs(TRUE)
if (length(args)) {
  specs <- read_manifest()
  out <- lapply(seq_len(nrow(specs)), function(j) {
    meta <- specs[j, ]
    dat <- load_dd2(meta, args[[1L]])
    local_summary(dat$x, dat$env, unlist(meta[c("low", "mid", "high")]), meta$community, "DD2")
  })
  write_result(do.call(rbind, out), "local_empirical.csv")
}

bench <- list()
for (n in c(72, 147, 241, 600)) for (density in c(.03, .2, .85)) {
  p <- 350
  k <- floor(.08 * n)
  x <- matrix(rexp(n * p) * rbinom(n * p, 1, density), n)
  x[, 1] <- 1
  for (type in c("incidence", "relative")) for (pairs in c(FALSE, TRUE)) {
    times <- numeric(12)
    for (i in seq_along(times)) {
      times[i] <- system.time(output <- conditional_zeta(x, k, type, pairs = pairs))[["elapsed"]]
    }
    bench[[length(bench) + 1L]] <- data.frame(n, p, k, density = mean(x > 0), type, pairs,
      seconds = median(times), input_kib = as.numeric(object.size(x)) / 1024,
      output_kib = as.numeric(object.size(output)) / 1024)
  }
}
write_result(do.call(rbind, bench), "conditional_cost.csv")
