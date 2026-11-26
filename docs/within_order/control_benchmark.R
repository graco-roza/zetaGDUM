source("docs/within_order/simulation_design.R")
args <- commandArgs(TRUE)
B <- if (length(args)) as.integer(args[1L]) else 32L
fit_data <- expand.grid(order = c(2L, 8L), response = c("incidence", "relative"))
rows <- lapply(seq_len(nrow(fit_data)), function(j) {
  k <- fit_data$order[j]
  type <- as.character(fit_data$response[j])
  set.seed(3341L + k)
  sites <- within_sites(70L)
  point_data <- within_data(sites, k, 3000L, type, within_reference())
  fit <- within_fit(point_data)
  # The path integration API requires a reference whose site count matches
  # this community. point_data$base stores the external scales on these sites.
  ref <- point_data$base
  layout <- function(a) {
    d <- ugdm_interaction_design(a$design, sites$env, "nonlinear", "constant")
    cbind(intercept = 1, d$dir, d$pos)
  }
  paths <- function(M) sample_sets_coupled(sites$counts, k, M, exact_limit = 0L)
  exact <- function(M) sample_sets_coupled(sites$counts, k, M)
  plan_plain <- ugdm_global_score_plan(sites$counts, sites$env, k, type = type,
    controls = character(), ref = ref, design_matrix = layout)
  evaluate <- function(plan, pool) ugdm_global_site_scores(plan, pool, fit$coef)$pseudovalues
  reference1 <- evaluate(plan_plain, exact(70000L))
  reference2 <- evaluate(plan_plain, exact(70000L))
  reference <- (reference1 + reference2) / 2
  cost <- function(expr) {
    start <- proc.time()
    v <- force(expr)
    delta <- (proc.time() - start)[1:2]
    list(value = v, seconds = sum(delta))
  }
  use_cv <- function(M) {
    p <- ugdm_global_score_plan(sites$counts, sites$env, k, paths(M %/% 3L),
      fit$coef, type = type, ref = ref, design_matrix = layout)
    evaluate(p, paths(M - M %/% 3L))
  }
  use_plain <- function(M) evaluate(plan_plain, paths(M))
  counts <- c(1800L, 5400L)
  time_at <- function(fun) vapply(counts, function(M)
    median(replicate(5L, cost(fun(M))$seconds)), numeric(1L))
  times <- time_at(use_plain)
  budget <- median(replicate(5L, cost(use_cv(max(counts)))$seconds))
  slope <- diff(times) / diff(counts)
  nplain <- if (slope > 0) round(counts[1L] + (budget - times[1L]) / slope) else max(counts)
  nplain <- max(900L, min(18000L, nplain))
  scale <- sqrt(mean(reference^2))
  do.call(rbind, lapply(seq_len(B), function(r) {
    fits <- list(plain_equal_sets = cost(use_plain(max(counts))),
                 plain_cpu_target = cost(use_plain(nplain)),
                 controls = cost(use_cv(max(counts))))
    do.call(rbind, lapply(names(fits), function(m) {
      x <- fits[[m]]
      data.frame(order = k, response = type, replicate = r, method = m,
        seconds = x$seconds, evaluations = if (m == "plain_cpu_target") nplain else max(counts),
        reference_gap = sqrt(mean((reference1 - reference2)^2)) / scale,
        relative_score_error = sqrt(mean((x$value - reference)^2)) / scale,
        relative_covariance_error = norm(cov(x$value)-cov(reference),"F") / norm(cov(reference),"F"))
    }))
  }))
})
within_write(do.call(rbind, rows), "control_benchmark.csv")
