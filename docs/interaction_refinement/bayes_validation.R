source("docs/phase5_refinement/common.R")
args <- commandArgs(TRUE)
mode <- match.arg(args[1L], c("calibration", "joint", "nonlinear"))
signal <- match.arg(args[2L], c("null", "weak"))
reps <- if (length(args) >= 3L) as.integer(args[3L]) else 16L
rows <- list()
grid <- if (mode == "joint") 900L else 5000L
reference <- NULL
if (mode == "joint") {
  reference <- read.csv(gzfile("docs/phase5_refinement/results/posterior.csv.gz"))
  additional <- "docs/interaction_refinement/results/bayes_additional_jel.csv"
  if (file.exists(additional)) reference <- rbind(reference, read.csv(additional))
  reference <- subset(reference, mode == "long_reference" & stable)
}
selected <- if (mode == "joint") sort(unique(reference$r[reference$signal == signal])) else seq_len(reps)
for (r in selected) {
  sim <- j_layout(70L, r, signal)
  a <- sim$pool(grid)
  b <- if (mode == "joint") a else sim$pool(15000L)
  g <- ugdm_bayes_gmm(a, b, draws = 5000L, seed = 6782 + r)
  if (r <= 4L && mode == "calibration") {
    projected <- ugdm_bayes_gmm(a, b, draws = 5000L, project = TRUE, seed = 6782 + r)
  }
  J <- do.call(cbind, lapply(seq_along(b), function(j)
    ugdm_pool_scores(list(b[[j]]), g$pilot_estimate[(j - 1L) * 3L + seq_len(3L)], "centered")))
  I <- sweep(J, 2L, colMeans(J)) %*% solve(g$derivative)
  n <- nrow(J)
  del <- sweep(-I, 2L, (n - 1) * g$estimate, "+") / (n - 1)
  sys <- function(f) {
    before <- proc.time()[["elapsed"]]
    out <- f()
    list(draws = out$draws, seconds = proc.time()[["elapsed"]] - before)
  }
  fit <- list(pseudovalue = sys(function()
    ugdm_pseudovalue_bootstrap(g$estimate, del, draws = 5000L, seed = r + 712L)),
    gaussian = g)
  if (r <= 4L && mode == "calibration")
    fit$projected <- list(draws = projected$projected, seconds = projected$seconds)
  if (mode %in% c("joint", "nonlinear")) {
    fit$nonlinear <- ugdm_bayes_gmm(a, b, integration = "mcmc", draws = 4000L,
      warmup = 1200L, seed = 8120 + r)

  }
  for (method in names(fit)) {
    f <- fit[[method]]
    contrast <- f$draws %*% t(map6)
    w <- if (is.null(f$weights)) rep(1 / nrow(contrast), nrow(contrast)) else f$weights
    diag <- if (is.null(f$rhat)) list(rhat = NA_real_, ess = NA_real_) else
      internal(".jel_chain_diagnostics")(lapply(split(seq_len(nrow(contrast)), f$group),
      function(i) contrast[i, , drop = FALSE]))
    relevant <- function(x, y, fn) if (all(is.na(c(x, y)))) NA_real_ else fn(c(x, y), na.rm = TRUE)
    rhat <- relevant(f$rhat, diag$rhat, max)
    ess <- relevant(f$ess, diag$ess, min)
    for (j in seq_len(nrow(map6))) {
      v <- contrast[, j]
      q <- qweight(v, w, c(.025, .25, .75, .975))
      rows[[length(rows) + 1L]] <- data.frame(signal, r, grid, method, effect = rownames(map6)[j],
        truth = drop(map6 %*% sim$truth)[j], mean = sum(v * w),
        p_positive = sum(w[v > 0]), p_practical = sum(w[v > .1]),
        l95 = q[1L], l50 = q[2L], u50 = q[3L], u95 = q[4L],
        ess, rhat, seconds = f$seconds, covariance_condition = g$covariance_condition)
    }
  }
  if (mode == "joint") {
    target_r <- r
    ref <- subset(reference, r == target_r)
    ref <- ref[ref$signal == signal, ]
    stopifnot(nrow(ref) == 6L)
    rows[[length(rows) + 1L]] <- data.frame(signal, r, grid, method = "jel", effect = ref$effect,
      truth = ref$truth, mean = ref$mean, p_positive = ref$p, p_practical = ref$practical,
      l95 = ref$l95, l50 = ref$l50, u50 = ref$h50, u95 = ref$h95,
      ess = NA_real_, rhat = NA_real_, seconds = NA_real_, covariance_condition = NA_real_)
  }
}
write_study(do.call(rbind, rows), "docs/interaction_refinement/results",
  paste("bayes", mode, signal, sep = "_"))
