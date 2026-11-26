library(zetaGDUM)
out <- "docs/adaptive_refinement/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(82531)
n <- 70L
results <- lapply(c("uniform", "skewed", "symmetric_sets"), function(scenario) {
  z <- switch(scenario, uniform = sort(runif(n, -2, 2)),
    skewed = sort(rexp(n)), symmetric_sets = seq(-2, 2, length.out = n))
  env <- data.frame(climate = z)
  do.call(rbind, lapply(c(2L, 6L, 12L), function(k) {
    idx <- if (scenario == "symmetric_sets") {
      first <- if (k == 2L) matrix(sample.int(n / 2L, 2500L, replace = TRUE), ncol = 1L) else
        sample_sets(n / 2L, k / 2L, M = 2500L)
      cbind(first, n + 1L - first)
    } else sample_sets(n, k, M = 2500L)
    d <- ugdm_interaction_design(ugdm_design(env, idx), env, "nonlinear", "contextual")
    audit <- ugdm_contextual_diagnostics(d)
    e <- colMeans(d$contextual^2)^.5
    contrast <- rowSums(d$contextual)
    data.frame(scenario, n, k, rank_additive = audit$rank[["additive"]],
      rank_midpoint = audit$rank[["midpoint"]], rank_contextual = audit$rank[["context"]],
      free_columns = ncol(cbind(1, d$dir, d$pos)), context_ramp_rms = mean(contrast^2)^.5,
      context_basis_rms = max(e), max_context_residual = if (all(is.na(audit$context_residual_fraction)))
        NA_real_ else max(audit$context_residual_fraction, na.rm = TRUE))
  }))
})
write.csv(do.call(rbind, results), file.path(out, "context_geometry.csv"), row.names = FALSE)

quadratic <- function(b, P) drop(crossprod(b, P %*% b))
basis_penalty <- getFromNamespace(".context_penalty_basis", "zetaGDUM")
values <- lapply(c(1, 3, 7), function(width) {
  z <- seq(0, width, length.out = 30L)
  B <- splines2::iSpline(z, df = 3L, degree = 2L)
  P <- basis_penalty(B)
  do.call(rbind, lapply(c("linear", "zero_endpoint", "tilted"), function(pattern) {
    b <- switch(pattern, linear = c(1, 1, 1), zero_endpoint = c(1, -2, 1),
      tilted = c(.5, 1.3, 2.2))
    native <- quadratic(b, P$derivative)
    exact <- integrate(function(x) drop(predict(B, x, derivs = 1L) %*% b)^2, 0, width)$value
    data.frame(pattern = pattern, range = width, coefficient_penalty = sum(b^2), derivative_native = native,
      derivative_unit_range = width * native, max_quadrature_error = abs(native - exact),
      projected_linear_amplitude = sum(b * P$null_slope))
  }))
})
z <- do.call(rbind, values)
write.csv(z, file.path(out, "functional_penalty_units.csv"), row.names = FALSE)
