source("docs/phase5_refinement/common.R")
args <- commandArgs(TRUE)
mode <- match.arg(args[1L], c("regular", "boundary"))
n <- as.integer(args[2L])
case <- args[3L]
B <- if (length(args) >= 4L) as.integer(args[4L]) else 200L
M <- if (length(args) >= 5L) as.integer(args[5L]) else 3000L
orders <- c(2L, as.integer(round(.10 * n)))
n_draws <- 600L
cuts <- c(.50, .80, .90, .95)
target <- switch(case, null = c(-.3, 0, 0, -.3, 0, 0),
  weak = c(-.3, .25, .12, -.3, .15, .25),
  spatial = c(-.3, .25, .12, -.3, .15, .25),
  near = c(-.3, .10 * sqrt(70 / n), .12, -.3, .07 * sqrt(70 / n), .12),
  interior = c(-.3, .52, .12, -.3, .42, .12),
  stop("Unknown case"))
truth <- unname(drop(map6 %*% target))
if (!mode %in% c("regular","boundary") || is.na(B) || B < 1L) stop("Invalid run.")
sets <- function(M) sample_sets_coupled(matrix(1, n, 1), orders, M, exact_limit = 0L)
scenarios <- c("null", "weak", "spatial", "near", "interior")
records <- lapply(seq_len(B), function(r) {
  set.seed(49209L + n * 911L + match(case, scenarios) * 2101L + r * 31L)
  z <- rnorm(n)
  u <- rnorm(n)
  if (case == "spatial") for (i in seq_len(n - 1L) + 1L) {
    u[i] <- .60 * u[i - 1L] + sqrt(1 - .6^2) * u[i]
    z[i] <- .60 * z[i - 1L] + sqrt(1 - .6^2) * z[i]
  }
  u <- tanh(u)
  dist <- abs(outer(z, z, "-"))
  build <- function(j, idx) {
    k <- orders[j]
    spread <- internal("meanpair_dist_cpp")(dist, idx)
    position <- rowSums(matrix(z[idx], nrow(idx)))
    X <- cbind(intercept = 1, spread = spread, position = position / k)
    mu <- plogis(drop(X %*% target[(j - 1L) * 3L + seq_len(3L)]))
    y <- mu + .85 * pmin(mu, 1 - mu) * rowMeans(matrix(u[idx], nrow(idx)))
    design <- list(idx = idx, n = n, k = k, Xdir = cbind(g_linear = spread),
      Xpos = cbind(g = position), vars = "g", spread_width = 1L)
    list(design = design, y = y)
  }
  pool <- function(M) {
    sampled <- sets(M)
    lapply(seq_along(sampled), function(j) build(j, sampled[[j]]$idx))
  }
  p <- pool(M)
  pilot <- pool(max(5000L, M))
  pools <- replicate(3L, pool(M), simplify = FALSE)
  score <- function(a, theta, method = "inclusion") do.call(cbind, lapply(seq_along(a), function(j)
    ugdm_pool_scores(list(a[[j]]), theta[(j - 1L) * 3L + seq_len(3L)], method)))
  fit <- ugdm_bayes_gmm(p, pilot, moment_covariance = function(theta) {
    all <- lapply(pools, score, theta = theta)
    v <- cov(Reduce("+", all) / length(all))
    list(covariance = v, raw = all)
  }, prior_sd = 2, draws = n_draws, project = mode == "boundary")
  V <- fit$covariance_diagnostics$raw
  # pool and average are reported alongside the fully cross-pool estimators.
  covs <- list(
    single_ratio = cov(score(pools[[1]], fit$pilot_estimate, "centered")),
    pooled = fit$omega)
  if (mode == "regular") {
    covs$average <- ugdm_score_covariance(V, "average")$covariance
    covs$cross <- ugdm_score_covariance(V, "cross")$covariance
    covs$diagonal <- ugdm_score_covariance(V, "cross", "diagonal", alpha = .15)$covariance
    covs$block <- ugdm_score_covariance(V, "cross", "block", alpha = .15,
      blocks = rep(c("intercept", "spread", "position"), 2))$covariance
    covs$floor <- ugdm_score_covariance(V, "cross", "eigen", eigen_floor = .02)$covariance
    covs$psd <- ugdm_score_covariance(V, "cross", "psd")$covariance
  }
  info <- fit$derivative
  estimate <- fit$estimate
  rho <- function(V) V[3L, 6L] / sqrt(V[3L, 3L] * V[6L, 6L])
  summaries <- lapply(names(covs), function(method) {
    omega <- covs[[method]]
    spec <- eigen(omega / tcrossprod(sqrt(diag(fit$omega))), symmetric = TRUE)$values
    status <- all(is.finite(spec)) && min(spec) > max(spec) * 1e-9
    out <- data.frame(method, effect = rownames(map6), truth, level = rep(cuts, each = 6L),
       valid = status, covariance_min = min(spec), covariance_max = max(spec),
       position_correlation = suppressWarnings(rho(omega)))
    out$mean <- out$width <- out$lo <- out$hi <- out$p_positive <- out$p_practical <- NA_real_
    if (!status) return(out)
    Q <- n * crossprod(info, solve(omega, info))
    posterior <- solve(Q + diag(.25, length(estimate)))
    mu <- drop(map6 %*% (posterior %*% (Q %*% estimate)))
    sd <- sqrt(diag(map6 %*% posterior %*% t(map6)))
    half <- qnorm((1 + out$level) / 2) * rep(sd, times = length(cuts))
    out$mean <- rep(mu, times = length(cuts))
    out$lo <- out$mean - half
    out$hi <- out$mean + half
    out$width <- 2 * half
    out$p_positive <- rep(pnorm(mu / sd), times = length(cuts))
    out$p_practical <- rep(pnorm((mu - .1) / sd), times = length(cuts))
    out
  })
  if (mode == "boundary") {
    projected <- fit$projected %*% t(map6)
    observed_face <- ugdm_gmm_tangent(fit, c(2L, 5L), draws = n_draws)
    oracle_face <- ugdm_gmm_tangent(fit, c(2L, 5L), active =
      c(2L, 5L)[target[c(2L, 5L)] == 0], draws = n_draws)
    samples <- list(projected = projected,
      tangent_plugin = observed_face$increments %*% t(map6),
      tangent_oracle = oracle_face$increments %*% t(map6))
    point <- drop(map6 %*% observed_face$estimate)
    for (method in names(samples)) {
      X <- samples[[method]]
      add <- summaries[[1L]]
      add$method <- method
      add$valid <- TRUE
      add$mean <- if (method == "projected") colMeans(X) else point
      q <- t(vapply(seq_len(nrow(add)), function(i) {
        j <- 1L + (i - 1L) %% ncol(X)
        bounds <- quantile(X[, j], c((1 - add$level[i]) / 2, (1 + add$level[i]) / 2),
          names = FALSE)
        if (method == "projected") bounds else point[j] - rev(bounds)
      }, numeric(2)))
      add$lo <- q[, 1L]
      add$hi <- q[, 2L]
      add$width <- diff(t(q))[1L, ]
      add$p_positive <- if (method == "projected") colMeans(X > 0) else NA_real_
      add$p_practical <- if (method == "projected") colMeans(X > .1) else NA_real_
      summaries[[length(summaries) + 1L]] <- add
    }
  }
  out <- do.call(rbind, summaries)
  cbind(n, case, replicate = r, out)
})
data <- do.call(rbind, records)
out <- "docs/adaptive_refinement/results"
key <- paste("bayes", mode, n, case, sep = "_")
if (M != 3000L) key <- paste0(key, "_M", M)
con <- gzfile(file.path(out, paste0(key, ".csv.gz")), "wt", compression = 9)
write.csv(data, con, row.names = FALSE)
close(con)
summary <- do.call(rbind, lapply(split(data, data[c("method", "effect", "level")]), function(g) {
  keep <- g$valid
  val <- function(x) if (any(keep)) mean(x[keep], na.rm = TRUE) else NA_real_
  data.frame(g[1L, c("n","case","method","effect","level")], total = nrow(g),
    returned = sum(keep), coverage = val(g$lo <= g$truth + 1e-12 & g$hi >= g$truth - 1e-12),
    bias = val(g$mean - g$truth), width = val(g$width), positive = val(g$p_positive >= .95), negative = val(g$p_positive <= .05),
    practical = val(g$p_practical >= .95), mean_cor = mean(g$position_correlation, na.rm = TRUE),
    negative_spectrum = mean(g$covariance_min < 0))
}))
write_study(summary, out, paste0(key, "_summary"))
