source("docs/phase5_refinement/common.R")
set.seed(98031)
n <- 150L
orders <- c(6L, 14L)
z <- (rnorm(n) + stats::rexp(n, 2)) / 1.1
env <- data.frame(climate = z)
ref <- ugdm_design(env, cbind(seq_len(n - 1L), 2:n))
u <- tanh(rnorm(n))
path <- function(M) sample_sets_coupled(matrix(1, n, 1L), orders, M, exact_limit = 0)
raw_matrix <- function(d) {
  v <- ugdm_interaction_design(d, env, "nonlinear", "contextual")
  cbind(intercept = 1, spread = rowSums(v$base), position = v$pos[, 1],
    nonlinear = v$pos[, 2], pair = rowSums(v$pair_midpoint), context = rowSums(v$contextual))
}
fixed <- path(14000L)
scale <- vapply(seq_along(orders), function(j) {
  X <- raw_matrix(ugdm_design(env, fixed[[j]]$idx, ref = ref))
  c(1, apply(X[, -1L], 2L, sd))
}, numeric(6))
target <- cbind(c(-2, .24, .10, -.10, .28, .26), c(-2, .13, .27, -.14, .25, .38))
pool <- function(M) {
  sampled <- path(M)
  lapply(seq_along(orders), function(j) {
    d <- ugdm_design(env, sampled[[j]]$idx, ref = ref)
    X <- sweep(raw_matrix(d), 2L, scale[, j], "/")
    mu <- plogis(drop(X %*% target[, j]))
    y <- mu + .55 * pmin(mu, 1 - mu) * rowMeans(matrix(u[d$idx], nrow(d$idx)))
    list(design = d, y = y, X = X)
  })
}
a <- pool(8000L)
b <- pool(12000L)
raw <- replicate(4L, pool(12000L), simplify = FALSE)
design_matrix <- function(d) d$X
omega <- function(theta, method) {
  J <- lapply(raw, function(p) {
    do.call(cbind, lapply(seq_along(p), function(j) {
      d <- p[[j]]
      X <- d$X
      score <- internal("site_score_kernel_cpp")(t(d$design$idx), t(X), d$y,
        rep(1 / nrow(X), nrow(X)), theta[(j - 1L) * 6L + seq_len(6L)], n)
      full <- score[n + 1L, ]
      cond <- score[seq_len(n), ] * n / orders[j]
      sweep(sweep(cond, 2L, colMeans(cond)) * orders[j] * (n - 1) / (n - orders[j]),
        2L, full, "+")
    }))
  })
  if (method == "pooled") ugdm_score_covariance(J, "pooled") else
    ugdm_score_covariance(J, "cross", "eigen", eigen_floor = .02)
}
fits <- list(
  gaussian_pooled = ugdm_bayes_gmm(a, b, design_matrix = design_matrix,
    moment_covariance = function(x) omega(x, "pooled"), draws = 6000L, seed = 3010),
  gaussian_cross = ugdm_bayes_gmm(a, b, design_matrix = design_matrix,
    moment_covariance = function(x) omega(x, "cross"), draws = 6000L, seed = 3010),
  nonlinear_cross = ugdm_bayes_gmm(a, b, design_matrix = design_matrix,
    moment_covariance = function(x) omega(x, "cross"), integration = "mcmc",
    draws = 5500L, warmup = 1800L, seed = 7023))
evaluate <- function(draws) {
  change <- draws[, 7:12, drop = FALSE] - draws[, 1:6, drop = FALSE]
  cbind(position_amplitude_increase = sqrt(rowSums(draws[, 9:10, drop = FALSE]^2)) -
      sqrt(rowSums(draws[, 3:4, drop = FALSE]^2)), declining_spread = -change[, 2L],
    nonlinear_position = -draws[, 10L], pair_modification = draws[, 11L],
    contextual_modification = draws[, 12L], contextual_increase = change[, 6L])
}
true <- evaluate(matrix(as.vector(target), nrow = 1))
tab <- lapply(names(fits), function(m) {
  fit <- fits[[m]]
  x <- evaluate(fit$draws)
  q <- t(apply(x, 2L, quantile, c(.025,.975), names = FALSE))
  moments <- if (is.null(fit$rhat)) NULL else internal(".jel_chain_diagnostics")(
    lapply(split(seq_len(nrow(x)), fit$group), function(i) x[i, , drop = FALSE]))
  data.frame(method = m, effect = colnames(x), truth = drop(true), mean = colMeans(x),
    lo = q[, 1L], hi = q[, 2L], p_positive = colMeans(x > 0), p_practical = colMeans(x > .10),
    raw_min = min(fit$covariance_diagnostics$spectrum),
    adjusted_min = min(fit$covariance_diagnostics$adjusted_spectrum),
    min_ess = if (is.null(moments)) NA_real_ else min(fit$ess, moments$ess),
    max_rhat = if (is.null(moments)) NA_real_ else max(fit$rhat, moments$rhat),
    elapsed = fit$seconds)
})
write_study(do.call(rbind, tab), "docs/adaptive_refinement/results", "functional_posterior")
