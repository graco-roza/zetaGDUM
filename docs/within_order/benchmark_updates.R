# Compare the same Rademacher score perturbations under three coalition solvers.
library(zetaGDUM)
n <- 70L
extended <- identical(commandArgs(TRUE), "--extended")
draws <- if (extended) 49L else 19L
M <- if (extended) 1600L else 1000L
jobs <- expand.grid(seed = if (extended) c(24L, 51L) else 24L, order = c(4L, 12L),
                    type = c("total", "relative"), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
results <- lapply(seq_len(nrow(jobs)), function(i) {
  started <- proc.time()[["elapsed"]]
  job <- jobs[i, ]
  sim <- simulate_metacommunity(n, n_species = 65L, spatial = 0, breadth = .85, seed = job$seed)
  x <- sim$comm
  x[rowSums(x) == 0, 1L] <- 1
  env <- cbind(sim$env, noise = rnorm(n))
  k <- job$order
  base <- ugdm_design(env, matrix(seq_len(k), nrow = 1L))
  make <- function(M) {
    idx <- sample_sets(n, k, M)
    d <- ugdm_interaction_design(ugdm_design(env, idx, ref = base), env, "nonlinear", "constant")
    list(design = d, y = ugdm_response(x, idx, job$type, working_weights = "constant")$y)
  }
  a <- make(M)
  eval <- make(1800)
  pools <- list(make(M), make(M))
  fit <- ugdm_fit_interactions(a$design, a$y, ridge = .002, distance = FALSE)
  shared_preparation <- proc.time()[["elapsed"]] - started
  methods <- c("one_step", "hybrid", "refit")
  runs <- lapply(methods, function(method) ugdm_within_order_resample(
    fit, a$design, a$y, evaluation = eval$design, evaluation_y = eval$y,
    method = "multiplier", score_update = method, score_pools = pools, draws = draws,
    seed = 419L, hybrid_tolerance = .0003))
  reference <- runs[[3L]]
  do.call(rbind, Map(function(z, method) {
    keep <- z$valid & reference$valid
    compare <- function(matrix) {
      lhs <- z[[matrix]][keep[z$valid], , drop = FALSE]
      rhs <- reference[[matrix]][keep[reference$valid], , drop = FALSE]
      sqrt(mean((lhs - rhs)^2))
    }
    data.frame(job, method, n = n, M = M, draws = draws, valid = sum(z$valid),
      paired = sum(keep), importance_RMSE = compare("predictor_draws"),
      contrast_RMSE = compare("difference_draws"), coalition_refits = sum(z$exact_score_refits, na.rm = TRUE),
      elapsed = z$elapsed, preparation_seconds = shared_preparation,
      total_elapsed = z$elapsed + shared_preparation)
  }, runs, methods))
})
filename <- if (extended) "benchmark_updates.csv" else "benchmark_updates_final.csv"
write.csv(do.call(rbind, results), file.path("docs/within_order/results", filename), row.names = FALSE)
