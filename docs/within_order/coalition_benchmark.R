library(zetaGDUM)
directory <- "docs/within_order/results"
dir.create(directory, showWarnings = FALSE, recursive = TRUE)
methods <- c("exact", "sample_40", "sample_160", "sample_640")
bench <- lapply(c(3L, 5L), function(p) {
  set.seed(4712L + p)
  n <- 70L
  env <- as.data.frame(matrix(rnorm(n * p), n))
  train <- seq_len(48L)
  holdout <- setdiff(seq_len(n), train)
  draws <- function(indices) matrix(indices[sample_sets(length(indices), 4L, M = 500L)], ncol = 4L)
  a <- ugdm_design(env, draws(train), basis_sites = train)
  b <- ugdm_design(env, draws(holdout), ref = a)
  da <- ugdm_interaction_design(a, env, "nonlinear")
  db <- ugdm_interaction_design(b, env, "nonlinear")
  site_effect <- rnorm(n, sd = .12)
  response <- function(d) plogis(-.3 + d$dir[, 2L] + .4 * d$dir[, 5L] +
    .8 * d$pos[, 1L] + .4 * d$pos[, 5L] +
    rowMeans(matrix(site_effect[d$idx], ncol = d$order)))
  ya <- response(da)
  yb <- response(db)
  f <- ugdm_fit_interactions(da, ya, distance = FALSE, ridge = .01)
  reference <- NULL
  do.call(rbind, lapply(methods, function(method) {
    rows <- lapply(seq_len(if (method == "exact") 1L else 5L), function(r) {
      m <- if (method == "exact") 40L else as.integer(sub("sample_", "", method))
      started <- proc.time()[["elapsed"]]
      q <- ugdm_within_order_importance(f, da, ya, evaluation = db,
        evaluation_y = yb, exact_max = if (method == "exact") 8L else 0L,
        permutations = m, seed = 1839L + r)
      elapsed <- proc.time()[["elapsed"]] - started
      cells <- q$hierarchy
      if (method == "exact") reference <<- cells$absolute
      data.frame(p, method, replicate = r, elapsed,
        coalitions = q$diagnostics$coalition_count, cells,
        reference = reference, error = cells$absolute - reference)
    })
    do.call(rbind, rows)
  }))
})
result <- do.call(rbind, bench)
write.csv(result, file.path(directory, "coalition_benchmark.csv"), row.names = FALSE)
