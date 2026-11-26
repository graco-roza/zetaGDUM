# Run from the repository root using the installed development package.
source("docs/adaptive_refinement/ecology_common.R")
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Supply the DD2 data directory and a community manifest row.")
meta <- read_manifest()[as.integer(args[2L]), , drop = FALSE]
data <- load_dd2(meta, args[1L])
n <- nrow(data$x)
M <- if (length(args) >= 3L) as.integer(args[3L]) else 2400L
stopifnot(is.finite(M), M > 0L)
orders <- sort(unique(c(2L, 4L, 8L, 12L, meta$mid, meta$high)))
fold <- study_blocks(data$dist, partitions = 1L)[[1L]]
orders <- orders[orders <= min(table(fold))]
set.seed(55009L + n)
env <- data$env
dist <- data$dist
ridge <- .001
new_fit <- function(d, y) ugdm_fit_interactions(d, y, ridge = ridge)
make <- function(index, ref) study_d(env, index, ref, "nonlinear", "constant", dist)
response <- function(index, kind) {
  ugdm_response(data$x, index, if (kind == "incidence") "total" else "relative",
    working_weights = "constant")$y
}
allocate <- function(fit, d, y, ev = NULL, yy = NULL) {
  params <- list(fit = fit, design = d, y = y, evaluation = ev, evaluation_y = yy,
    permutations = 48L)
  # Pair orderings between the two numerical pools so permutation error is
  # not included in their difference. Restore the sampling RNG stream.
  previous <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(assign(".Random.seed", previous, envir = .GlobalEnv))
  params$seed <- 5003L
  if ("failure" %in% names(formals(ugdm_within_order_importance))) params$failure <- "warn"
  do.call(ugdm_within_order_importance, params)
}
fit_info <- function(fit, d) {
  layout <- internal(".interaction_layout")(d, fit$has_distance)
  mu <- ugdm_predict(fit, d)
  X <- layout$X * sqrt(mu * (1 - mu))
  scale <- sqrt(colSums(X^2))
  scale[scale <= .Machine$double.eps] <- 1
  spectrum <- svd(sweep(X, 2L, scale, "/"), nu = 0L, nv = 0L)$d
  rank <- sum(spectrum > max(spectrum) * sqrt(.Machine$double.eps))
  condition <- if (rank == ncol(X)) (max(spectrum) / min(spectrum))^2 else Inf
  data.frame(converged = fit$converged, information_rank = rank,
    parameters = ncol(X), condition, active = sum(fit$coef[layout$dir] <= 1e-8))
}
describe <- function(k, kind, group, contribution) {
  data.frame(meta[c("id", "community")], order = k, response = kind,
    group, term = names(contribution$absolute), absolute = as.numeric(contribution$absolute),
    relative = as.numeric(contribution$relative), total = contribution$total, row.names = NULL)
}
all_sites <- seq_len(n)
model_rows <- point_rows <- fold_rows <- list()
for (k in orders) {
  message(meta$community, " at order ", k)
  train1 <- study_idx(all_sites, k, M)
  train2 <- study_idx(all_sites, k, M)
  assess <- study_idx(all_sites, k, max(2L * M, 6000L))
  base <- study_reference(env, train1, dist)
  designs <- lapply(list(train1 = train1, train2 = train2, assess = assess), make, ref = base)
  sharing <- mean(internal("shared_species_cpp")(1L * (data$x > 0), assess) > 0)
  variation <- c(apply(designs$assess$base, 2L, stats::sd),
    apply(designs$assess$pos, 2L, stats::sd))
  spatial <- lapply(sort(unique(fold)), function(label) {
    train <- which(fold != label)
    test <- which(fold == label)
    tr <- study_idx(train, k, max(600L, M %/% 2L))
    te <- study_idx(test, k, max(600L, M %/% 2L))
    ref <- study_reference(env, tr, dist, train)
    list(label = label, tr = tr, te = te, train = make(tr, ref), test = make(te, ref))
  })
  for (kind in c("incidence", "relative")) {
    yy <- lapply(list(train1 = train1, train2 = train2, assess = assess), response, kind = kind)
    fit <- new_fit(designs$train1, yy$train1)
    second <- new_fit(designs$train2, yy$train2)
    values <- allocate(fit, designs$train1, yy$train1)
    common1 <- allocate(fit, designs$train1, yy$train1, designs$assess, yy$assess)
    common2 <- allocate(second, designs$train2, yy$train2, designs$assess, yy$assess)
    point <- rbind(describe(k, kind, "predictor", values$predictors),
                   describe(k, kind, "component", values$components))
    point$allocation_mc_se <- c(
      if (is.null(values$predictors$mc_se)) rep(NA_real_, nrow(point) - length(values$components$absolute)) else values$predictors$mc_se,
      if (is.null(values$components$mc_se)) rep(NA_real_, length(values$components$absolute)) else values$components$mc_se)
    point$common_assessment_absolute <- c(common1$predictors$absolute, common1$components$absolute)
    point$second_pool_absolute <- c(common2$predictors$absolute, common2$components$absolute)
    point$pool_difference <- point$common_assessment_absolute - point$second_pool_absolute
    point_rows[[length(point_rows) + 1L]] <- point
    row <- data.frame(meta[c("id", "community")], sites = n, order = k, response = kind,
      fit_sets = M, assessment_sets = nrow(assess), sd = stats::sd(yy$assess), sharing,
      smallest_design_sd = min(variation), median_design_sd = stats::median(variation),
      fit_info(fit, designs$train1), second_converged = second$converged, D2 = values$predictors$total,
      same_sites_assessment_D2 = common1$predictors$total,
      largest_pool_gap = if (any(is.finite(point$pool_difference)))
        max(abs(point$pool_difference), na.rm = TRUE) else NA_real_)
    spatial_results <- lapply(spatial, function(v) {
      u <- response(v$tr, kind)
      t <- response(v$te, kind)
      f <- new_fit(v$train, u)
      null <- loss(t, rep(mean(u), length(t)), rep(1 / length(t), length(t)))
      dev <- loss(t, ugdm_predict(f, v$test), rep(1 / length(t), length(t)))
      data.frame(meta[c("id", "community")], response = kind, order = k, block = v$label,
        n_training_sites = sum(fold != v$label), n_test_sites = sum(fold == v$label),
        training_sets = nrow(v$tr), test_sets = nrow(v$te), converged = f$converged,
        null, deviance = dev, heldout_D2 = if (f$converged && null > 0) 1 - dev / null else NA_real_)
    })
    f <- do.call(rbind, spatial_results)
    ok <- f$converged & is.finite(f$heldout_D2)
    row$spatial_D2 <- if (all(ok)) 1 - sum(f$deviance) / sum(f$null) else NA_real_
    row$spatial_positive_blocks <- sum(ok & f$heldout_D2 > 0)
    row$spatial_valid_blocks <- sum(ok)
    model_rows[[length(model_rows) + 1L]] <- row
    fold_rows[[length(fold_rows) + 1L]] <- f
  }
}
folder <- "docs/within_order/results"
write_study(do.call(rbind, point_rows), folder, paste0("empirical_importance_", meta$id))
write_study(do.call(rbind, model_rows), folder, paste0("empirical_diagnostics_", meta$id))
write_study(do.call(rbind, fold_rows), folder, paste0("empirical_spatial_", meta$id))
