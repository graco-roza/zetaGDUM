# Base-R summaries from fixed replicate-level validation outputs.
root <- "docs/adaptive_refinement/results"
read_many <- function(pattern) {
  data <- do.call(rbind, lapply(list.files(root, pattern, full.names = TRUE), read.csv))
  names(data)[names(data) == "form"] <- "model"
  data
}
use <- function(x) x[is.finite(x)]
safe <- function(x, f) if (length(use(x))) f(use(x)) else NA_real_
save <- function(value, id) write.csv(value, file.path(root, paste0(id,".csv")), row.names = FALSE)
group <- function(d, fields, calc) {
  ff <- do.call(interaction, c(unname(as.list(d[fields])), list(drop = TRUE, lex.order = TRUE)))
  do.call(rbind, lapply(split(seq_len(nrow(d)), ff), function(ii) {
    a <- d[ii, , drop = FALSE]
    cbind(a[1L, fields, drop = FALSE], calc(a), row.names = NULL)
  }))
}
err <- function(a, column, f) safe(a[[column]] - a[[paste0("true_", column)]], f)
rmse <- function(v) sqrt(mean(v^2))

sim <- read_many("^context_(linear|redundant|positive|negative|collinear|spatial|sparse)[.]csv$")
stopifnot(length(unique(sim$scenario)) == 7L, nrow(sim) == 3360L)
point_stats <- function(a) data.frame(replicates = nrow(a),
  converged = mean(a$convergence), structural_rank = safe(a$rank, median),
  fit_parameters = safe(a$parameters, median),
  raw_context_novelty = safe(a$context_incremental_rank, median),
  residual_norm_fraction = safe(a$context_orthogonal_norm, median),
  holdout_RMSE = safe(a$test_RMSE, mean), conditional_RMSE = safe(a$reference_RMSE, mean),
  support_gap = safe(a$support_gap, median), sharing = safe(a$sharing, median),
  response_spread_bias = err(a, "response_spread", mean),
  response_spread_SD = err(a, "response_spread", stats::sd),
  response_spread_RMSE = err(a, "response_spread", rmse),
  response_position_bias = err(a, "response_position", mean),
  response_position_SD = err(a, "response_position", stats::sd),
  response_position_RMSE = err(a, "response_position", rmse),
  response_joint_bias = err(a, "response_joint", mean),
  response_joint_SD = err(a, "response_joint", stats::sd),
  response_joint_RMSE = err(a, "response_joint", rmse),
  link_joint_RMSE = err(a, "link_joint", rmse),
  n_context = length(use(a$fitted_context)), context_bias = err(a, "fitted_context", mean),
  context_RMSE = err(a, "fitted_context", rmse),
  true_context = safe(a$true_context, mean))
# The fitted-context target has a deliberately explicit name in the raw files.
sim$true_fitted_context <- sim$true_context
save(group(sim, c("scenario", "n", "order", "model"), point_stats), "ecology_context_summary")

pairs <- do.call(rbind, lapply(split(sim, sim$scenario), function(a) {
  k <- sort(unique(a$order))
  do.call(rbind, lapply(seq_len(length(k) - 1L), function(i) {
    available <- c("scenario", "n", "replicate", "model")
    lower <- a[a$order == k[i], ]
    higher <- a[a$order == k[i + 1L], ]
    same <- merge(lower, higher, by = available, suffixes = c(".low", ".high"))
    out <- same[available]
    out$low <- k[i]
    out$high <- k[i + 1L]
    for (v in c("response_spread", "response_position", "response_joint", "link_joint")) {
      out[[v]] <- same[[paste0(v, ".high")]] - same[[paste0(v, ".low")]]
      out[[paste0("true_", v)]] <- same[[paste0("true_",v,".high")]] -
        same[[paste0("true_",v,".low")]]
    }
    out
  }))
}))
save(group(pairs, c("scenario", "n", "low", "high", "model"), function(a) data.frame(
  response_spread_RMSE = err(a, "response_spread", rmse),
  response_position_RMSE = err(a, "response_position", rmse),
  response_joint_RMSE = err(a, "response_joint", rmse),
  link_joint_RMSE = err(a, "link_joint", rmse))), "ecology_context_paired_summary")

a <- sim[sim$order == 6L & sim$scenario %in% c("redundant", "positive", "negative") &
  sim$model %in% c("D","D_weak","D_derivative","D_strong"), ]
a$error <- a$response_joint - a$true_response_joint
groups <- c("D", "D_weak", "D_derivative", "D_strong")
p <- c("#5c707c","#bb573d","#28718a")
C <- c("redundant", "positive", "negative")
plot_file <- file.path(root, "adaptive_interaction_path.png")
png(plot_file, 1280, 590, res = 140)
par(mfrow = c(1, 2), mar = c(4.5, 4.6, 2, 1))
yr <- range(a$error)
plot(c(1, 4), yr, type = "n", xaxt = "n", xlab = "Derivative penalty",
     ylab = "Signed response-contrast error")
axis(1, 1:4, c("0", "1e-6", ".005", ".15"))
abline(h = 0, col = "grey78")
for (j in seq_along(C)) {
  d <- a[a$scenario == C[j], ]
  vals <- lapply(groups, function(m) d$error[d$model == m])
  average <- vapply(vals, mean, 0.0)
  limits <- sapply(vals, quantile, c(.25, .75), names = FALSE)
  x <- seq_along(vals)
  lines(x,average, col = p[j], lwd = 2, type = "b", pch = 16)
  arrows(x, limits[1L, ], x, limits[2L, ], angle = 90, code = 3, length = .04, col = p[j])
}
y_high <- max(vapply(split(a$error, interaction(a$scenario, a$model)), rmse, 0.0))
plot(c(1, 4), c(0, y_high * 1.12), type = "n", xaxt = "n", xlab = "Derivative penalty",
     ylab = "Response-contrast RMSE")
axis(1, 1:4, c("0", "1e-6", ".005", ".15"))
for (j in seq_along(C)) {
  d <- a[a$scenario == C[j], ]
  values <- vapply(groups, function(m) rmse(d$error[d$model == m]), 0.0)
  lines(seq_along(values), values, col = p[j], lwd = 2, type = "b", pch = 16)
}
legend("topright", c("Redundant turnover", "Positive context", "Negative context"),
       col = p, lwd = 2, bty = "n", cex = .86)
invisible(dev.off())

sel <- read_many("^selection_(linear|redundant|positive|spatial)[.]csv$")
if (NROW(sel)) {
  sel$contrast_RMSE <- sqrt((sel$spread_error^2 + sel$position_error^2 + sel$joint_error^2) / 3)
  choose <- function(which) group(sel[sel$candidate == sel[[which]], ], "scenario", function(a)
    data.frame(policy = which, replicates = nrow(a), median_penalty_ID = median(a$candidate),
      mean_contrast_RMSE = mean(a$contrast_RMSE), median_contrast_RMSE = median(a$contrast_RMSE),
      mean_agreement = mean(a$number_agreeing), oracle_agreement = mean(a$candidate == a$oracle)))
  save(do.call(rbind, lapply(c("rule", "loss_only", "oracle"), choose)), "ecology_selection_summary")
}
cv <- read_many("^ecology_cv_.*[.]csv$")
if (NROW(cv)) {
  cv_summary <- group(cv, c("id", "community", "response", "candidate"), function(a) data.frame(
    selected = unique(a$selected) == a$candidate[1L], mean_loss = mean(a$validation_loss),
    spread_SD = safe(a$spread, sd), position_SD = safe(a$position, sd),
    joint_SD = safe(a$joint, sd), convergence = mean(a$converged)))
  save(cv_summary, "ecology_cv_summary")
}
mods <- read_many("^ecology_models_.*[.]csv$")
if (NROW(mods)) save(mods, "ecology_empirical_summary")
sites <- read_many("^ecology_sites_.*[.]csv$")
if (NROW(sites)) {
  cor_ok <- function(a, b) {
    good <- is.finite(a) & is.finite(b)
    if (sum(good) <= 2L || stats::sd(a[good]) < 1e-12 || stats::sd(b[good]) < 1e-12) return(NA_real_)
    stats::cor(a[good], b[good], method = "spearman")
  }
  keep <- c("observed", "predicted", "spread", "position", "modification")
  summary <- group(sites, c("id", "community", "response", "order", "model"), function(a) data.frame(
    sites = nrow(a), median_retention = median(a$characteristic_retention),
    retention_SD = sd(a$characteristic_retention), median_order_survival = median(a$retention),
    attribution_error = max(abs(a$predicted - (a$base + a$spread + a$position + a$modification))),
    as.list(vapply(a[keep], function(b) cor_ok(a$characteristic_retention, b), 0.0))))
  save(summary, "ecology_site_summary")
}
