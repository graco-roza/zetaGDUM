library(zetaGDUM)

internal <- function(name) getFromNamespace(name, "zetaGDUM")
args <- commandArgs(TRUE)
site_count <- if (length(args)) as.integer(args[1L]) else 24L
replicates <- if (length(args) > 1L) as.integer(args[2L]) else 24L
cpu_only <- length(args) > 2L && args[3L] == "cpu_only"
cpu_replicates <- if (length(args) > 3L) as.integer(args[4L]) else 16L
k <- if (site_count <= 24L) c(2L, 4L) else c(2L, round(.11 * site_count))
n <- site_count
set.seed(n + 8149)
sites <- sort(runif(n))
gradient <- drop(t(chol(exp(-abs(outer(sites, sites, "-")) / .22))) %*% rnorm(n))
env <- data.frame(temperature = drop(scale(gradient + 1.6 * sites)))
center <- runif(34L, -2.5, 2.5)
width <- runif(length(center), .40, 1.7)
expected <- exp(-.5 * sweep(outer(env$temperature, center, "-"), 2L, width, "/")^2)
community <- matrix(rpois(length(expected), 3.4 * expected), n)
community <- cbind(community, rpois(n, 2) + 1L)
r <- ugdm_design(env, cbind(seq_len(n - 1L), seq.int(2L, n)))
build <- function(samples) ugdm_path_data(community, env, samples, ref = r)
design <- function(d) cbind(1, d$Xdir, d$Xpos / d$k)
center <- function(a) sweep(a, 2L, colMeans(a))
layout <- c(5L, 5L)
contrast <- rbind(lower_spread = c(0, 1, 1, 1, 0, rep(0, 5)),
  lower_position = c(0, 0, 0, 0, 1, rep(0, 5)))
contrast <- rbind(contrast, higher_spread = c(rep(0, 5), contrast[1L, 1:5]),
  higher_position = c(rep(0, 5), contrast[2L, 1:5]))
contrast <- rbind(contrast, paired_spread = contrast[3L, ] - contrast[1L, ],
                  paired_position = contrast[4L, ] - contrast[2L, ])
reference_sample <- function(M) sample_sets_coupled(community, k, M, exact_limit = 15000L)
first <- build(reference_sample(if (n <= 24L) 15000L else 250000L))
theta <- lapply(first, function(a) suppressWarnings(glm.fit(design(a$design), a$relative,
  family = binomial())$coefficients))
projection <- lapply(seq_along(first), function(j) {
  X <- design(first[[j]]$design)
  u <- plogis(drop(X %*% theta[[j]]))
  solve(crossprod(X * sqrt(u * (1 - u))) / nrow(X))
})
influence <- contrast
for (j in 1:2) {
  p <- (j - 1L) * 5L + 1:5
  influence[, p] <- contrast[, p] %*% projection[[j]]
}
global_scores <- function(pool, method) {
  out <- lapply(seq_along(pool), function(j) {
    a <- pool[[j]]
    s <- lapply(a, function(p) list(design = p$design, y = p$relative,
      sampling_weights = p$sampling_weights))
    P <- internal(".prepare_score_pools")(s)
    V <- internal(".pool_scores")(P, theta[[j]], method)
    G <- Reduce("+", lapply(P, function(p) internal(".site_score")(p, theta[[j]], "centered")$full)) /
      length(P)
    list(V = V, G = G)
  })
  list(V = do.call(cbind, lapply(out, "[[", "V")), G = unlist(lapply(out, "[[", "G")))
}
by_order <- function(...) {
  pools <- list(...)
  lapply(seq_along(k), function(j) lapply(pools, "[[", j))
}
Vfirst <- center(global_scores(by_order(first), "centered")$V)
if (n > 24L) {
  other <- build(reference_sample(250000L))
  Vother <- center(global_scores(by_order(other), "centered")$V)
  reference <- (Vfirst + Vother) / 2
  reference_gap <- sqrt(colMeans(((Vfirst - Vother) %*% t(influence))^2) / 2)
  rm(other, Vother)
} else {
  reference <- Vfirst
  reference_gap <- rep(0, nrow(influence))
}
truth <- reference %*% t(influence)
rm(first)
covariance <- function(a) c(cov(a[, 1L], a[, 3L]), cov(a[, 2L], a[, 4L]))
cross_truth <- covariance(truth)

methods <- c("deletion", "centered", "inclusion", "pooled", "focal_independent",
             "focal_matched", "response", "predictors", "adaptive")
pilot <- 30L
focal <- c("focal_independent", "focal_matched", "response", "predictors", "adaptive")
cpu <- function(fn) {
  start <- proc.time()[["user.self"]]
  out <- fn()
  out$seconds <- proc.time()[["user.self"]] - start
  out
}
simulate <- function(method, M) cpu(function() {
  if (method %in% focal) {
    controls <- switch(method, response = "response", predictors = c("response", "predictors"),
                       adaptive = c("response", "predictors"), character())
    m <- M - if (length(controls)) pilot else 0L
    plan <- ugdm_conditioned_score_plan(community, env, k, m, theta, controls = controls,
      pilot_draws = pilot, adaptive = method == "adaptive", minimum_draws = 12L,
      allocation_contrast = influence[5:6, , drop = FALSE], ref = r,
      coupled = method != "focal_independent")
    v <- ugdm_conditioned_scores(plan, theta, influence)
    return(list(V = v$contrast, G = v$score,
                cost = v$evaluations + v$pilot_evaluations, min_draws = min(v$counts)))
  }
  draw <- function(m) build(sample_sets_coupled(community, k, m * n, exact_limit = 0L))
  pools <- if (method == "pooled") by_order(draw(M %/% 2L), draw(M - M %/% 2L)) else by_order(draw(M))
  s <- global_scores(pools, if (method == "pooled") "conditional" else method)
  list(V = s$V %*% t(influence), G = s$G, cost = n * M * length(k), min_draws = NA_real_)
})

measure <- function(budget, method, target, i) {
  v <- simulate(method, target)
  shifted <- center(v$V)
  scale <- sqrt(diag(cov(truth))[1:2] * diag(cov(truth))[3:4])
  error <- sqrt(colMeans((shifted - truth)^2))
  columns <- seq_len(nrow(contrast))
  data.frame(n, n_species = ncol(community), lower = k[1], higher = k[2],
    budget, method, replicate = i, effect = rownames(contrast), evaluations = v$cost,
    seconds = v$seconds, min_draws = v$min_draws,
    score_point_shift = drop(influence %*% v$G),
    pseudo_RMSE = error, reference_gap, ecological_SD = apply(truth, 2L, sd),
    cross_covariance_error = c(rep(NA_real_, 4L), (covariance(shifted) - cross_truth) / scale),
    row.names = columns)
}
root <- "docs/interaction_refinement/results"
raw_file <- file.path(root, paste0("score_", n, "_repeats.csv.gz"))
if (cpu_only) {
  prior <- subset(read.csv(raw_file), budget != "cpu")
} else {
  set.seed(site_count + 26031)
  output <- list()
  for (M in c(90L, 240L)) for (method in methods) {
    message(n, " equal evaluations ", M, " ", method)
    for (i in seq_len(replicates)) output[[length(output) + 1L]] <- measure(
      paste0("evaluation_", M), method, M, i)
  }
  prior <- do.call(rbind, output)
}
set.seed(site_count + 53127)
timing <- vapply(methods, function(method) {
  vapply(c(90L, 240L), function(M) {
    median(replicate(3L, simulate(method, M)$seconds))
  }, numeric(1))
}, numeric(2))
target_time <- max(timing[2L, ])
output <- list()
for (method in methods) {
  low <- timing[1L, method]
  high <- timing[2L, method]
  slope <- max(.001 / 150, (high - low) / 150)
  M <- max(90L, round(240L + (target_time - high) / slope))
  message(n, " equal CPU ", M, " ", method)
  for (i in seq_len(cpu_replicates)) output[[length(output) + 1L]] <- measure("cpu", method, M, i)
}
all <- rbind(prior, do.call(rbind, output))
out <- gzfile(raw_file, "wt", compression = 9L)
write.csv(all, out, row.names = FALSE)
close(out)
summarise <- function(d) {
  a <- d[1L, c("n", "lower", "higher", "budget", "method", "effect")]
  cbind(a, min_draws = if (all(is.na(d$min_draws))) NA else median(d$min_draws),
    median_seconds = median(d$seconds), evaluations = median(d$evaluations),
    point_RMSE = sqrt(mean(d$score_point_shift^2)),
    site_RMSE = sqrt(mean(d$pseudo_RMSE^2)), reference_gap = median(d$reference_gap),
    ecological_SD = median(d$ecological_SD),
    cross_covariance_RMSE = sqrt(mean(d$cross_covariance_error^2)))
}
tab <- do.call(rbind, lapply(split(all, interaction(all$budget, all$method, all$effect)), summarise))
write.csv(tab, file.path(root, paste0("score_", n, "_summary.csv")), row.names = FALSE, na = "")
