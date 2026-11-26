source("docs/phase5_refinement/common.R")
args <- commandArgs(TRUE)
B <- if (length(args)) as.integer(args[1L]) else 12L
n <- 70L
orders <- c(2L, 8L)
ans <- lapply(seq_len(B), function(r) {
  set.seed(8735L + 157L * r)
  env <- data.frame(temperature = rnorm(n))
  optima <- seq(-2.8, 2.8, length.out = 42L)
  breadth <- if (r %% 2L) .50 else .95
  rate <- exp(-.5 * (outer(env$temperature, optima, "-") / breadth)^2)
  x <- matrix(rpois(length(rate), 4 * rate + .10), n)
  x[, 1L] <- x[, 1L] + 1L
  ref <- ugdm_design(env, cbind(seq_len(n - 1L), 2:n))
  pool <- function(M) {
    d <- ugdm_path_data(x, env, sample_sets_coupled(x, orders, M, "uniform", exact_limit = 0),
      ref = ref)
    lapply(d, function(a) list(design = a$design, y = a$relative))
  }
  a <- pool(5000L)
  b <- pool(14000L)
  base <- ugdm_bayes_gmm(a, b, draws = 1500L, seed = 766L + r)
  results <- lapply(c(FALSE, TRUE), function(cv) {
    start <- proc.time()[["elapsed"]]
    plan <- ugdm_conditioned_score_plan(x, env, orders, draws = if (cv) 155L else 200L,
      theta = base$pilot_estimate, controls = if (cv) c("response", "predictors") else character(),
      pilot_draws = 45L, ref = ref, seed = 8531L + r)
    fit <- ugdm_bayes_gmm(a, b, site_scores = function(theta)
      ugdm_conditioned_scores(plan, theta)$pseudovalues, draws = 1500L, seed = 766L + r)
    list(covariance = fit$sandwich_covariance, omega = fit$omega,
      seconds = proc.time()[["elapsed"]] - start)
  })
  base_method <- list(covariance = base$sandwich_covariance, omega = base$omega, seconds = base$seconds)
  methods <- c(list(global = base_method), stats::setNames(results, c("conditional", "conditional_cv")))
  nm <- names(base$estimate)
  nc <- length(nm) / 2L
  K <- cbind(-diag(nc), diag(nc))
  K <- rbind(diag(length(nm)), K)
  data <- lapply(names(methods), function(method) {
    out <- methods[[method]]
    data.frame(replicate = r, breadth, method, contrast = c(nm, paste0("change.", nm[seq_len(nc)])),
      ecological_se = sqrt(diag(K %*% out$covariance %*% t(K))),
      covariance_difference = norm(out$omega - base$omega, "F") / norm(base$omega, "F"),
      seconds = out$seconds)
  })
  do.call(rbind, data)
})
write_study(do.call(rbind, ans), "docs/interaction_refinement/results", "bayes_conditioned")
