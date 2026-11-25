source("docs/phase3/common.R")
Rcpp::sourceCpp("docs/manuscript_validation/calibration_helpers.cpp",
                cacheDir = Sys.getenv("CALIBRATION_CPP_CACHE", tempdir()))

args <- commandArgs(TRUE)
signal <- if (length(args)) match.arg(args[1L], c("null", "weak", "strong")) else "null"
reps <- if (length(args) > 1L) as.integer(args[2L]) else 48L
stopifnot(reps > 1)
M <- if (length(args) > 2L) as.integer(args[3L]) else 1600L
stopifnot(M >= 1000L)
B_boot <- 79L
overlap_pairs <- 180L
sites <- read.csv("docs/manuscript_validation/communities.csv")
anchor <- data.frame(gradient = seq(-3, 3, length.out = 41L))
basis <- ugdm_design(anchor, iid_sets(nrow(anchor), 2, 100))
effects <- rbind(spread = c(0, 1, 1, 1, 0), position = c(0, 0, 0, 0, 1))
strength <- c(null = 0, weak = .20, strong = 1)[[signal]]
coefs <- c(-.3, strength * c(.32, .35, .33), .035 * strength)

coef_fit <- function(X, y, start = rep(0, ncol(X)), keep = rep(TRUE, length(y))) {
  calibration_fit(X, y, unname(start), keep, rep(1, length(y)))
}

calibrate <- function(n, k, rep_id, seed) {
  set.seed(seed)
  env <- data.frame(gradient = rnorm(n))
  site_noise <- tanh(rnorm(n))
  input <- function(idx) {
    d <- ugdm_design(env, idx, ref = basis)
    X <- cbind(1, d$Xdir, d$Xpos)
    mu <- plogis(drop(X %*% coefs))
    y <- mu + .85 * pmin(mu, 1 - mu) * rowMeans(matrix(site_noise[idx], nrow(idx)))
    list(X = X, y = y)
  }
  m <- n %/% k
  groups <- M %/% m
  d <- sample_sets_matched(n, k, groups, m)
  a <- input(d$idx)
  point <- coef_fit(a$X, a$y)
  fit <- drop(effects %*% point)
  truth <- drop(effects %*% coefs)
  mu <- plogis(drop(a$X %*% point))
  H <- crossprod(a$X * sqrt(mu * (1 - mu))) / nrow(a$X)
  score <- a$X * (a$y - mu)
  contrast <- function(cov) diag(effects %*% cov %*% t(effects))
  normal <- function(v) sapply(v, function(vi) if (is.finite(vi) && vi > 0) sqrt(vi) else NA_real_)
  result <- list()

  measure <- system.time({
    estimate <- u_matched_covariance(score, d, H)
    v <- contrast(estimate$sandwich_incomplete)
  })[["elapsed"]]
  result$matched <- list(se = normal(v), psd = estimate$min_eigenvalue[["incomplete"]] >= 0,
                         seconds = measure, scores = nrow(a$X))

  measure <- system.time({
    kernel <- function(idx) {
      b <- input(idx)
      b$X * (b$y - plogis(drop(b$X %*% point)))
    }
    estimate <- u_overlap_covariance(n, k, kernel, overlap_pairs)
    cov <- (1 - 1 / groups) * estimate$complete +
      estimate$contrast[[k + 1L]] / nrow(a$X)
    extra <- zetaGDUM:::.sandwich_covariance(cov, H)
    v <- contrast(extra)
  })[["elapsed"]]
  result$overlap <- list(se = normal(v), psd = min(eigen(cov, symmetric = TRUE)$values) >= 0,
                         seconds = measure, scores = 2 * k * overlap_pairs)

  measure <- system.time({
    score_mean <- colMeans(score)
    included <- zetaGDUM:::.inclusion_sums(d$idx, cbind(1, score), n) / nrow(a$X)
    mass <- included[, 1L]
    pseudo <- (n - 1) * (included[, -1L] - outer(mass, score_mean)) / (1 - mass)
    influence <- sweep(pseudo, 2L, colMeans(pseudo)) %*% solve(H, t(effects))
    se <- sqrt(colSums(influence^2) / (n * (n - 1)))
  })[["elapsed"]]
  result$multiplier <- list(se = se, psd = TRUE, seconds = measure, scores = 0)

  measure <- system.time({
    deletion <- t(vapply(seq_len(n), function(i) {
      keep <- rowSums(d$idx == i) == 0
      drop(effects %*% coef_fit(a$X, a$y, point, keep))
    }, numeric(2)))
    se <- sqrt((n - 1) / n * colSums(sweep(deletion, 2L, colMeans(deletion))^2))
  })[["elapsed"]]
  result$jackknife <- list(se = se, psd = TRUE, seconds = measure, scores = n * nrow(a$X))

  measure <- system.time({
    omega <- site_weights(n, B_boot)
    bb <- t(vapply(seq_len(B_boot), function(j) {
      b <- input(calibration_cp(omega[, j], k, nrow(a$X)))
      drop(effects %*% coef_fit(b$X, b$y, point))
    }, numeric(2)))
    bounds <- apply(bb, 2L, quantile, c(.025, .975), na.rm = TRUE)
  })[["elapsed"]]
  result$bayesian <- list(se = apply(bb, 2, sd), bounds = bounds, psd = TRUE,
                         seconds = measure, scores = B_boot * nrow(a$X))

  rows <- lapply(names(result), function(method) {
    z <- result[[method]]
    q <- if (method %in% c("jackknife", "multiplier")) qt(.975, n - 1L) else qnorm(.975)
    bounds <- if (is.null(z$bounds)) rbind(fit - q * z$se, fit + q * z$se) else z$bounds
    data.frame(n, k, ratio = k / n, rep_id, signal, method, term = rownames(effects),
               point = unname(fit), truth = unname(truth), se = unname(z$se),
               lo = unname(bounds[1, ]), hi = unname(bounds[2, ]),
               seconds = z$seconds, kernel_rows = z$scores, psd = z$psd,
               budget = nrow(a$X), min_omission = nrow(a$X) * (1 - max(mass)),
               condition = kappa(H, exact = TRUE), design_kib = as.numeric(object.size(a)) / 1024)
  })
  do.call(rbind, rows)
}

out <- list()
for (j in seq_len(nrow(sites))) {
  n <- sites$sites[j]
  if (length(args) > 3L && !n %in% as.integer(strsplit(args[4L], ",", fixed = TRUE)[[1]])) next
  k <- sites$high[j]
  message(signal, " n=", n, " k=", k)
  for (r in seq_len(reps)) {
    out[[length(out) + 1L]] <- calibrate(n, k, r, 533000L + j * 1000L + r)
  }
}
write_result(do.call(rbind, out), paste0("calibration_", signal, if (M != 1600L) paste0("_M", M) else "", ".csv"))
