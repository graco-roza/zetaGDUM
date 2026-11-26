# From the repository root after installing zetaGDUM:
# Rscript docs/adaptive_refinement/covariance_validation.R
library(zetaGDUM)
set.seed(768039)
R <- 4L
reps <- 240L
out <- "docs/adaptive_refinement/results"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
effects <- c("intercept", "spread", "position")
blocks <- rep(effects, 2L)
contrast <- rbind(
  lower_spread = c(0, 1, 0, 0, 0, 0),
  higher_spread = c(0, 0, 0, 0, 1, 0),
  paired_spread = c(0, -1, 0, 0, 1, 0),
  paired_position = c(0, 0, -1, 0, 0, 1)
)
methods <- c("raw", "average", "pooled", "cross", "cross_psd", "cross_eigen",
             "cross_diag", "cross_block")
rows <- list()
for (n in c(70L, 150L, 250L)) for (scenario in c("ordinary", "cancellation")) {
  within <- matrix(c(.6, .08, .12, .08, .24, .10, .12, .10, .43), 3L)
  rho <- if (scenario == "ordinary") .65 else .99
  S <- kronecker(matrix(c(1, rho, rho, 1), 2L), within)
  noise <- kronecker(matrix(c(1, .35, .35, 1.4), 2L), .35 * within)
  sd_true <- sqrt(diag(contrast %*% S %*% t(contrast)) / n)
  truth <- c(-.4, .30, .1, -.25, .25, .17)
  for (r in seq_len(reps)) {
    ecological <- matrix(rnorm(n * 6L), n) %*% chol(S)
    point <- truth + colMeans(ecological)
    reference <- cov(ecological)
    pools <- replicate(R, ecological + matrix(rnorm(n * 6L), n) %*% chol(noise),
                       simplify = FALSE)
    estimates <- list()
    for (method in methods) {
      if (!startsWith(method, "cross_")) {
        v <- ugdm_score_covariance(pools, method = method, contrasts = contrast)
      } else {
        kind <- switch(method, cross_psd = "psd", cross_eigen = "eigen",
                       cross_diag = "diagonal", cross_block = "block")
        v <- ugdm_score_covariance(pools, regularization = kind, blocks = blocks,
                                   eigen_floor = .03, contrasts = contrast)
      }
      A <- v$covariance
      se_diag <- diag(v$contrasts$regularized) / n
      bad <- min(v$adjusted_spectrum) <= .Machine$double.eps
      standardized_error <- (A - reference) / tcrossprod(sqrt(diag(S)))
      population_error <- (A - S) / tcrossprod(sqrt(diag(S)))
      tol <- c(.25, .25, -.04, .06)
      avg <- drop(contrast %*% point)
      se <- sqrt(pmax(0, se_diag))
      estimate_probability <- ifelse(se > 0, pnorm((avg - tol) / se), NA_real_)
      truth_probability <- pnorm((avg - tol) / sd_true)
      for (j in seq_len(nrow(contrast))) {
        differences <- rownames(contrast)[j]
        rows[[length(rows) + 1L]] <- data.frame(
          n, scenario, replicate = r, method, effect = differences,
          finite_cov_loss = sqrt(mean(standardized_error^2)),
          population_cov_loss = sqrt(mean(population_error^2)),
          max_reference_eigen = max(v$spectrum), min_reference_eigen = min(v$spectrum),
          min_adjusted_eigen = min(v$adjusted_spectrum), indefinite = min(v$spectrum) < -1e-12,
          singular_adjusted = bad, adjustment = v$adjustment, alpha = if (is.null(v$alpha)) NA else v$alpha,
          variance_ratio = se_diag[j] / sd_true[j]^2,
          finite_variance_ratio = se_diag[j] / drop(contrast[j, ] %*% reference %*% contrast[j, ]) * n,
          covariance_numerical_se = if (is.null(v$contrast_numerical_se)) NA else v$contrast_numerical_se[j, j],
          absolute_probability_error = abs(estimate_probability[j] - truth_probability[j])
        )
      }
    }
  }
}
raw <- do.call(rbind, rows)
connection <- gzfile(file.path(out, "covariance_repeats.csv.gz"), "wt", compression = 9)
write.csv(raw, connection, row.names = FALSE)
close(connection)
key <- c("n", "scenario", "method", "effect")
stats <- function(a) data.frame(
  covariance_RMSE_finite = sqrt(mean(a$finite_cov_loss^2)),
  covariance_RMSE_population = sqrt(mean(a$population_cov_loss^2)),
  cross_indefinite_fraction = mean(a$indefinite),
  singular_adjustment_fraction = mean(a$singular_adjusted),
  median_min_eigen = median(a$min_reference_eigen),
  median_alpha = median(a$alpha, na.rm = TRUE),
  median_adjustment = median(a$adjustment),
  median_numerical_contrast_cov_se = median(a$covariance_numerical_se, na.rm = TRUE),
  median_variance_ratio = median(a$variance_ratio),
  probability_available = mean(is.finite(a$absolute_probability_error)),
  probability_MAE = mean(a$absolute_probability_error, na.rm = TRUE))
out_table <- do.call(rbind, lapply(split(seq_len(nrow(raw)), raw[key], drop = TRUE), function(i) {
  cbind(raw[i[1L], key], stats(raw[i, ]))
}))
write.csv(out_table, file.path(out, "covariance_summary.csv"), row.names = FALSE)
