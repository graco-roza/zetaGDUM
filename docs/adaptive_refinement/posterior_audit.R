source("docs/phase5_refinement/common.R")
args <- commandArgs(TRUE)
signal <- match.arg(args[1L], c("null","weak"))
old <- read.csv(gzfile("docs/phase5_refinement/results/posterior.csv.gz"))
extra <- read.csv("docs/interaction_refinement/results/bayes_additional_jel.csv")
j <- subset(rbind(old, extra), mode == "long_reference" & stable)
j <- j[j$signal == signal, ]
rows <- lapply(sort(unique(j$r)), function(r) {
  sim <- j_layout(70L, r, signal)
  a <- sim$pool(900L)
  pilot <- sim$pool(15000L)
  independent <- replicate(4L, sim$pool(10000L), simplify = FALSE)
  co <- function(theta, method) {
    scores <- lapply(independent, function(p) do.call(cbind, lapply(seq_along(p), function(i)
      ugdm_pool_scores(list(p[[i]]), theta[(i - 1L) * 3L + seq_len(3L)], "inclusion"))))
    blocks <- rep(c("intercept", "spread", "position"), 2)
    if (method == "pooled") ugdm_score_covariance(scores, "pooled") else
      ugdm_score_covariance(scores, "cross", "eigen", blocks = blocks, eigen_floor = .02)
  }
  gauss <- ugdm_bayes_gmm(a, pilot, moment_covariance = function(theta) co(theta, "pooled"),
    draws = 6000L, seed = 5013 + r)
  reg <- ugdm_bayes_gmm(a, pilot, moment_covariance = function(theta) co(theta, "cross"),
    draws = 6000L, seed = 5013 + r)
  ref <- ugdm_bayes_gmm(a, pilot, moment_covariance = function(theta) co(theta, "cross"),
    integration = "mcmc", draws = 5000L, warmup = 1500L, seed = 34021 + r)
  z <- j[j$r == r, ]
  z <- z[match(rownames(map6), z$effect), ]
  compare <- function(f, method) {
    x <- f$draws %*% t(map6)
    diagnostic <- if (is.null(f$rhat)) NULL else internal(".jel_chain_diagnostics")(
      lapply(split(seq_len(nrow(x)), f$group), function(i) x[i, , drop = FALSE]))
    q <- t(apply(x, 2L, stats::quantile, probs = c(.025,.975), names = FALSE))
    data.frame(signal, r, method, effect = rownames(map6), mean = colMeans(x),
      p_positive = colMeans(x > 0), p_practical = colMeans(x > .10), lo = q[, 1L], hi = q[, 2L],
      p_decline = colMeans(x < -.10),
      raw_spectrum_min = min(f$covariance_diagnostics$spectrum),
      adjusted_spectrum_min = min(f$covariance_diagnostics$adjusted_spectrum),
      ess = if (is.null(diagnostic)) NA_real_ else min(f$ess, diagnostic$ess),
      rhat = if (is.null(diagnostic)) NA_real_ else max(f$rhat, diagnostic$rhat),
      seconds = f$seconds)
  }
  template <- compare(ref, "nonlinear_cross")
  saved <- template
  saved$method <- "jel"
  saved$mean <- z$mean
  saved$p_positive <- z$p
  saved$p_practical <- z$practical
  saved$lo <- z$l95
  saved$hi <- z$h95
  saved[c("p_decline","raw_spectrum_min","adjusted_spectrum_min","ess","rhat","seconds")] <- NA_real_
  rbind(compare(gauss, "gaussian_pooled"), compare(reg, "gaussian_cross"), template, saved)
})
write_study(do.call(rbind, rows), "docs/adaptive_refinement/results", paste0("posterior_audit_", signal))
