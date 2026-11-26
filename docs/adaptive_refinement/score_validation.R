# From the repository root after installing zetaGDUM:
# Rscript docs/adaptive_refinement/score_validation.R
library(zetaGDUM)
out <- "docs/adaptive_refinement/results"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
methods <- c("raw", "inclusion", "both")
trials <- 16L
n_pool <- 4L
study <- function(n, type) {
  set.seed(n * 1609 + (type == "relative") * 189)
  orders <- if (n == 24L) c(2L, 4L) else if (n == 70L) c(2L, 8L) else c(3L, 12L)
  spatial <- sort(runif(n))
  z <- drop(scale(drop(t(chol(exp(-abs(outer(spatial, spatial, "-")) / .19))) %*%
                         rnorm(n)) + 1.2 * spatial))
  env <- data.frame(temperature = z)
  centres <- seq(-2.3, 2.3, length.out = 35L)
  breadth <- runif(length(centres), if (n >= 150L) .20 else .40, if (n >= 150L) .55 else 1.1)
  E <- exp(-.5 * sweep(outer(z, centres, "-"), 2L, breadth, "/")^2)
  x <- cbind(1L, matrix(rpois(length(E), 3 * E), n))
  # One universally occupied column ensures well-defined relative
  # abundance. When sparse, many sets share only that background species.
  model <- function(a) cbind(intercept = 1, spread = a$design$Xlinear[, 1L],
                            position = a$design$Xpos[, 1L] / a$design$k)
  ref <- ugdm_design(env, cbind(seq_len(n - 1L), seq.int(2L, n)))
  draws <- function(M, exhaust = FALSE) {
    sample_sets_coupled(x, orders, M, exact_limit = if (exhaust) M else 0L)
  }
  sample <- draws(if (n == 24L) 15000L else 80000L, exhaust = n == 24L)
  design <- ugdm_path_data(x, env, sample, ref = ref)
  theta <- lapply(design, function(a) suppressWarnings(glm.fit(model(a), a[[type]],
    family = binomial())$coefficients))
  H <- Map(function(a, b) {
    X <- model(a)
    mu <- plogis(drop(X %*% b))
    crossprod(X * sqrt(mu * (1 - mu))) / nrow(X)
  }, design, theta)
  contrast <- matrix(0, 6L, 6L, dimnames = list(
    c("lower_spread", "lower_position", "higher_spread", "higher_position",
      "paired_spread", "paired_position"), NULL))
  contrast[cbind(1:4, c(2, 3, 5, 6))] <- 1
  contrast[5, ] <- contrast[3, ] - contrast[1, ]
  contrast[6, ] <- contrast[4, ] - contrast[2, ]
  influence <- contrast
  for (j in 1:2) {
    ii <- (j - 1L) * 3L + 1:3
    influence[, ii] <- contrast[, ii] %*% solve(H[[j]])
  }
  plan <- function(paths = NULL, controls = character()) {
    ugdm_global_score_plan(x, env, orders, paths, theta = theta, type = type,
      controls = controls, ref = ref, design_matrix = model)
  }
  plain <- plan()
  centre <- function(a) sweep(a, 2L, colMeans(a))
  baseline <- function(path) centre(ugdm_global_site_scores(plain, path, theta)$pseudovalues)
  A <- baseline(sample)
  B <- if (n == 24L) A else baseline(draws(80000L))
  truth <- (A + B) / 2
  C <- cov(truth)
  V <- truth %*% t(influence)
  truth_variance <- diag(cov(V))
  reference_gap <- sqrt(colMeans(((A - B) %*% t(influence))^2) / 2)
  ecological_SDs <- sqrt(diag(C))
  cross_mask <- outer(rep(1:2, each = 3L), rep(1:2, each = 3L), "!=")
  finite_sharing <- vapply(design, function(a) mean(a$shared > 1), numeric(1L))
  rm(sample, design)
  build <- function(method, M) {
    started <- proc.time()[["user.self"]]
    cv <- switch(method, raw = character(), inclusion = "inclusion",
                 both = c("inclusion", "response"))
    pilot_m <- if (length(cv)) M %/% 5L else 0L
    p <- if (length(cv)) plan(draws(pilot_m), cv) else plain
    production <- (M - pilot_m) %/% n_pool
    evaluated <- replicate(n_pool, ugdm_global_site_scores(p, draws(production), theta),
                           simplify = FALSE)
    time <- proc.time()[["user.self"]] - started
    pools <- lapply(evaluated, "[[", "pseudovalues")
    lapply(c("raw", "average", "pooled", "cross", "cross_psd", "cross_block"), function(m) {
      covar <- if (startsWith(m, "cross_")) {
        ugdm_score_covariance(pools, regularization = if (m == "cross_psd") "psd" else "block",
          blocks = rep(c("intercept", "spread", "position"), 2L), contrasts = influence)
      } else ugdm_score_covariance(pools, method = m, contrasts = influence)
      Z <- if (m == "raw") pools[[1L]] else Reduce("+", pools) / n_pool
      err <- centre(Z) %*% t(influence) - V
      error <- (covar$covariance - C) / tcrossprod(ecological_SDs)
      variances <- diag(covar$contrasts$regularized)
      data.frame(covariance_method = m, effect = rownames(contrast),
        evaluated = 2L * (pilot_m + n_pool * production), pilot = 2L * pilot_m,
        time, site_RMSE = sqrt(colMeans(err^2)),
        site_bias = colMeans(err),
        cov_error = sqrt(mean(error^2)), cross_error = sqrt(mean(error[cross_mask]^2)),
        min_spectrum = min(covar$spectrum), min_adjusted = min(covar$adjusted_spectrum),
        adjustment = covar$adjustment,
        variance_ratio = variances / truth_variance, variance_ref = truth_variance,
        reference_gap = reference_gap,
        covariance_numerical_se = if (is.null(covar$contrast_numerical_se)) NA_real_ else
          diag(covar$contrast_numerical_se))
    })
  }
  sizes <- c(1800L, 7200L)
  timing <- vapply(methods, function(method) {
    vapply(sizes, function(M) median(replicate(3L, build(method, M)[[1L]]$time[1L])), numeric(1))
  }, numeric(2L))
  target_cpu <- max(timing[2L, ])
  trials_rows <- list()
  for (method in methods) {
    slope <- max(.001 / diff(sizes), diff(timing[, method]) / diff(sizes))
    cpu_budget <- as.integer(min(200000, max(sizes[1L], sizes[2L] + (target_cpu - timing[2L, method]) / slope)))
    for (budget in c("count", "cpu")) {
      M <- if (budget == "count") 7200L else cpu_budget
      message(n, " ", type, " ", method, " ", budget, " ", M)
      for (iteration in seq_len(trials)) {
        result <- do.call(rbind, build(method, M))
        trials_rows[[length(trials_rows) + 1L]] <- cbind(n, lower = orders[1L], higher = orders[2L],
          type, higher_nonbackground_sharing = unname(finite_sharing[2L]),
          iteration, method, budget, result)
      }
    }
  }
  do.call(rbind, trials_rows)
}
raw <- do.call(rbind, list(study(24L, "relative"), study(70L, "relative"),
                          study(70L, "incidence"), study(150L, "relative")))
connection <- gzfile(file.path(out, "score_repeats.csv.gz"), "wt", compression = 9)
write.csv(raw, connection, row.names = FALSE)
close(connection)
key <- c("n", "type", "lower", "higher", "method", "budget", "covariance_method", "effect")
summarise <- function(a) {
  data.frame(median_time = median(a$time), median_evaluated = median(a$evaluated),
    higher_nonbackground_sharing = a$higher_nonbackground_sharing[1L],
    site_RMSE = sqrt(mean(a$site_RMSE^2)), reference_gap = a$reference_gap[1L],
    covariance_RMSE = sqrt(mean(a$cov_error^2)),
    paired_covariance_RMSE = sqrt(mean(a$cross_error^2)),
    indefinite_fraction = mean(a$min_spectrum < -1e-12),
    contrast_variance_ratio = median(a$variance_ratio),
    median_numerical_contrast_cov_se = median(a$covariance_numerical_se, na.rm = TRUE))
}
tab <- do.call(rbind, lapply(split(seq_len(nrow(raw)), raw[key], drop = TRUE), function(i) {
  cbind(raw[i[1L], key], summarise(raw[i, ]))
}))
write.csv(tab, file.path(out, "score_summary.csv"), row.names = FALSE)
