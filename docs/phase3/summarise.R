source("docs/phase3/common.R")
options(width = 155)
results <- function(file) read.csv(file.path("docs/phase3/results", file), na.strings = "")
table_by <- function(data, key, measure) {
  groups <- split(data, interaction(data[key], drop = TRUE))
  rownames_out <- lapply(groups, function(g) cbind(g[1L, key, drop = FALSE], as.data.frame(as.list(measure(g)))))
  result <- do.call(rbind, rownames_out)
  print(result, row.names = FALSE, digits = 4)
}
rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))
fitting <- function(g) c(
  fallback = mean(g$method != g$used),
  spread = rmse(g$spread_error), position = rmse(g$position_error),
  importance = rmse(g$importance_error), shapley = rmse(g$phi_error),
  rank_agreement = mean(g$same_importance_rank, na.rm = TRUE),
  prediction = rmse(g$prediction_rmse), local_shapley = rmse(g$local_shapley_rmse),
  loss_difference = mean(g$loss_difference), bias_spread = mean(g$spread_error),
  bias_position = mean(g$position_error), minimum_ESS = mean(g$min_delete_ess),
  site_cv = mean(g$inclusion_cv), pair_cv = mean(g$pair_cv),
  sharing = mean(g$sharing), nonzero = mean(g$nonzero), median_seconds = median(g$elapsed),
  design_KiB = median(g$design_kib), moment_seconds = median(g$moments_seconds))

exact <- results("fitting_exact.csv")
scale_results <- results("fitting_scale.csv")
cat("\nExhaustive comparison\n")
table_by(exact, c("M", "type", "method"), fitting)
cat("\nPaired change error (order five minus order three)\n")
keys <- c("type", "M", "method", "repeat_id")
paired <- merge(exact[exact$k == 3, ], exact[exact$k == 5, ], by = keys, suffixes = c("_lo", "_hi"))
table_by(paired, keys[1:3], function(g) c(spread = rmse(g$spread_error_hi - g$spread_error_lo),
  position = rmse(g$position_error_hi - g$position_error_lo), importance = rmse(g$importance_error_hi - g$importance_error_lo),
  shapley = rmse(g$phi_error_hi - g$phi_error_lo)))
cat("\nManuscript dimensions with a 30,000-set reference\n")
table_by(scale_results, c("type", "method"), fitting)
cat("\nManuscript dimensions individually\n")
table_by(scale_results, c("n", "type", "method"), fitting)

cat("\nFixed kernel variance estimates as ratios to truth\n")
fix <- results("fixed_kernel.csv")
table_by(fix, c("n", "k"), function(g) c(
  empirical_complete = var(g$observed_complete) / mean(g$truth_complete),
  matched = mean(g$complete_matched / g$truth_complete),
  overlap = mean(g$complete_overlap / g$truth_complete),
  matched_incomplete = mean(g$incomplete_matched / g$truth_incomplete),
  overlap_incomplete = mean(g$incomplete_overlap / g$truth_incomplete),
  sd_matched = sd(g$complete_matched / g$truth_complete),
  sd_overlap = sd(g$complete_overlap / g$truth_complete),
  negative_overlap = mean(g$complete_overlap < 0)))

calibration <- function(g) {
  valid <- is.finite(g$lo) & is.finite(g$hi)
  coverage <- g$lo[valid] <= g$truth[valid] & g$hi[valid] >= g$truth[valid]
  reject <- g$lo[valid] > 0 | g$hi[valid] < 0
  c(valid = sum(valid), total = nrow(g), coverage = mean(coverage),
    rejection = mean(reject), available_rejection = sum(reject) / nrow(g),
    bias = mean(g$point - g$truth), sd = sd(g$point), mean_se = mean(g$se, na.rm = TRUE),
    psd = mean(g$psd), median_seconds = median(g$seconds), median_kernel_rows = median(g$kernel_rows))
}
cal <- do.call(rbind, lapply(c("null", "weak", "strong"), function(s) results(paste0("calibration_", s, ".csv"))))
cat("\nInterval calibration at approximately 1,600 sampled sets\n")
table_by(cal, c("signal", "term", "method"), calibration)
cat("\nCalibration by manuscript site count and high order\n")
table_by(cal, c("n", "k", "signal", "term", "method"), calibration)
check_files <- list.files("docs/phase3/results", "^calibration_.*_M[0-9]+\\.csv$")
if (length(check_files)) {
  cat("\nLarger-budget sensitivity\n")
  hi <- do.call(rbind, lapply(check_files, results))
  table_by(hi, c("n", "k", "signal", "term", "method"), calibration)
}
