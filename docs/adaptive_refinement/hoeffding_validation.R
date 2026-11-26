source("docs/phase5_refinement/common.R")
args <- commandArgs(TRUE)
B <- if (length(args)) as.integer(args[1L]) else 1000L
ns <- c(70L, 150L, 250L)
cases <- list(strong = c(1, 1.4), weak = c(.04, .08), cancellation = c(1, 1),
              spatial = c(1, 1))
map <- rbind(lower = c(1, 0), higher = c(0, 1), paired = c(-1, 1))
plot_data <- list()
for (n in ns) {
  set.seed(81540L + n)
  pair <- t(combn(n, 2))
  paired_products <- function(z) z[pair[, 1L]] * z[pair[, 2L]]
  b <- c(1, 2)
  runs <- lapply(names(cases), function(case) {
    a <- cases[[case]]
    result <- lapply(seq_len(B), function(r) {
      x <- rnorm(n)
      if (case == "spatial") {
        for (i in seq_len(n - 1L) + 1L) x[i] <- .65 * x[i - 1L] + sqrt(1 - .65^2) * x[i]
      }
      linear <- (x[pair[, 1L]] + x[pair[, 2L]]) / 2
      product <- paired_products(x)
      scores <- tcrossprod(linear, a) + tcrossprod(product, b)
      fit <- ugdm_pair_covariance(scores, pair, n)
      mapped <- function(V) diag(map %*% V %*% t(map))
      refvar <- if (case == "spatial") rep(NA_real_, 3L) else
        mapped(tcrossprod(a) / n + tcrossprod(b) / nrow(pair))
      point <- drop(map %*% colMeans(scores))
      # Compare the complete kernel with independent matched prefixes at a
      # controlled numerical budget. The kernel is the same pair average
      # inside a two- or higher-order set.
      gap <- cost <- rep(NA_real_, 3L)
      if (r <= min(B, 64L)) {
        high <- as.integer(round(.10 * n))
        paths <- sample_sets_coupled(matrix(1, n, 1L), c(2L, high), M = 3000L,
          exact_limit = 0L)
        sampled <- vapply(seq_along(paths), function(j) {
          idx <- paths[[j]]$idx
          z <- matrix(x[idx], nrow(idx))
          pair_average <- (rowSums(z)^2 - rowSums(z^2)) /
            (ncol(z) * (ncol(z) - 1L))
          mean(a[j] * rowMeans(z) + b[j] * pair_average)
        }, numeric(1))
        gap <- drop(map %*% sampled) - point
      }
      data.frame(n, case, replicate = r, contrast = rownames(map), point,
        first = mapped(fit$first), second = mapped(fit$second),
        overlap = mapped(fit$covariance), jackknife = mapped(fit$jackknife),
        exact = refvar, incomplete_error = gap)
    })
    do.call(rbind, result)
  })
  plot_data[[as.character(n)]] <- do.call(rbind, runs)
}
a <- do.call(rbind, plot_data)
out <- "docs/adaptive_refinement/results"
con <- gzfile(file.path(out, "hoeffding_replicates.csv.gz"), "wt", compression = 9)
write.csv(a, con, row.names = FALSE)
close(con)
q <- c(.50,.80,.90,.95)
summary <- do.call(rbind, lapply(split(a, a[c("n","case","contrast")]), function(g) {
  analytical <- g$exact[1L]
  ref <- if (is.na(analytical)) var(g$point) else analytical
  truth <- if (g$case[1L] != "spatial") 0 else {
    n <- g$n[1L]
    Epair <- sum((n - seq_len(n - 1L)) * .65^seq_len(n - 1L)) / choose(n, 2L)
    switch(g$contrast[1L], lower = Epair, higher = 2 * Epair, paired = Epair)
  }
  do.call(rbind, lapply(c("jackknife","overlap","reference"), function(m) {
    v <- if (m == "reference") rep(ref, nrow(g)) else g[[m]]
    covered <- vapply(q, function(level) {
      keep <- v > 0
      if (any(keep)) mean(abs(g$point[keep] - truth) <= qnorm((1 + level) / 2) * sqrt(v[keep])) else NA_real_
    }, numeric(1))
    data.frame(g[1L, c("n","case","contrast")], method = m, replicates = nrow(g),
      empirical_variance = var(g$point), reference = ref, mean_ratio = mean(v) / ref,
      covariance_RMSE = sqrt(mean((v - ref)^2)) / ref, invalid = mean(v <= 0),
      first_mean = mean(g$first) / ref, second_mean = mean(g$second) / ref,
      numerical_RMSE = sqrt(mean(g$incomplete_error^2, na.rm = TRUE) / ref),
      cover50 = covered[1L], cover80 = covered[2L], cover90 = covered[3L], cover95 = covered[4L])
  }))
}))
write_study(summary, out, "hoeffding_summary")
