out <- "docs/interaction_refinement/results"
recovery <- read.csv(gzfile(file.path(out, "model_recovery.csv.gz")))
precision <- read.csv(file.path(out, "model_pool_precision.csv"))
summarise <- function(data, keys, calc) {
  rows <- split(seq_len(nrow(data)), interaction(data[keys], drop = TRUE))
  result <- do.call(rbind, lapply(rows, function(i) {
    cbind(data[i[1L], keys, drop = FALSE], as.data.frame(as.list(calc(data[i, , drop = FALSE]))))
  }))
  rownames(result) <- NULL
  result
}
avg <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
mid <- function(x) if (any(!is.na(x))) median(x, na.rm = TRUE) else NA_real_
rmse <- function(x) sqrt(avg(x^2))
metrics <- function(x) {
  threshold <- .05
  c(n = nrow(x), convergence = mean(x$converged), rank_deficient = mean(x$rank < x$parameters),
    median_condition = mid(x$condition), median_active = mid(x$active),
    median_support_gap = mid(x$support_gap), heldout_response_RMSE = avg(x$pred_RMSE),
    conditional_response_RMSE = avg(x$truth_RMSE), spread_SD_error = rmse(x$spread_SD - x$true_spread_SD),
    position_SD_error = rmse(x$position_SD - x$true_position_SD),
    link_RMSE = rmse(x$link - x$true_link), response_RMSE = rmse(x$response - x$true_response),
    explicit_bias = avg(x$explicit - x$true_explicit), explicit_RMSE = rmse(x$explicit - x$true_explicit),
    true_explicit = avg(x$true_explicit), explicit_exceeds_threshold = mean(abs(x$explicit) > threshold),
    explicit_sign_recovery = avg(ifelse(abs(x$true_explicit) > threshold,
      sign(x$explicit) == sign(x$true_explicit), NA_real_)),
    nonlinear_position_residual = mid(x$position_residual), interaction_residual = mid(x$interaction_residual),
    interaction_D2 = mid(x$interaction_D2), position_D2 = mid(x$position_D2),
    spread_D2 = mid(x$spread_D2), full_D2 = mid(x$D2))
}
result <- summarise(recovery, c("scenario", "k", "model"), metrics)
write.csv(result, file.path(out, "model_recovery_summary.csv"), row.names = FALSE)
selected <- subset(recovery, model == "site_selected")
if (nrow(selected)) {
  selection <- summarise(selected, c("scenario", "k", "penalty"), function(x) {
    c(chosen = nrow(x), average_omission_agreement = mean(x$selected_folds) / 5)
  })
  write.csv(selection, file.path(out, "model_smoothing_selection.csv"), row.names = FALSE)
}
metrics_precision <- function(x) c(link_SD = sd(x$link), response_SD = sd(x$response),
                                   explicit_SD = sd(x$explicit))
a <- summarise(precision, c("scenario", "M", "r", "k"), metrics_precision)
ks <- sort(unique(precision$k))
for (j in seq_len(length(ks) - 1L)) {
  key <- c("scenario", "M", "r", "p")
  b <- merge(precision[precision$k == ks[j], ], precision[precision$k == ks[j+1L], ], by = key)
  for (name in c("link", "response", "explicit")) b[[name]] <- b[[paste0(name, ".y")]] - b[[paste0(name, ".x")]]
  b <- summarise(b, c("scenario", "M", "r"), metrics_precision)
  b$k <- paste(ks[j], ks[j+1L], sep = " to ")
  a <- rbind(a, b)
}
precision_summary <- summarise(a, c("scenario", "M", "k"), function(x) {
  c(n = nrow(x), median_link_SD = median(x$link_SD), median_response_SD = median(x$response_SD),
    median_explicit_SD = median(x$explicit_SD))
})
write.csv(precision_summary, file.path(out, "model_pool_precision_summary.csv"), row.names = FALSE)

curves <- read.csv(file.path(out, "model_curve_recovery.csv"))
curve_summary <- summarise(curves, c("scenario", "k", "model"), function(x) {
  c(n = nrow(x), directional_RMSE = median(x$directional_RMSE),
    position_RMSE = median(x$position_RMSE), endpoint_RMSE = median(x$endpoint_RMSE))
})
write.csv(curve_summary, file.path(out, "model_curve_recovery_summary.csv"), row.names = FALSE)

png(file.path(out, "interaction_recovery.png"), width = 1480, height = 670, res = 130)
par(mfrow = c(1,3), mar = c(4.5, 4.5, 2.7, .8), mgp = c(2.75, .8, 0))
names <- c("curved_turnover", "positive", "negative")
colors <- c("#395E82", "#AC5127", "#237158")
groups <- lapply(names, function(case) {
  a <- subset(recovery, scenario == case & model %in% c("varying_spread", "joint", "joint_penalty"))
  boxplot(a$explicit - a$true_explicit ~ interaction(a$model, a$k), plot = FALSE)
})
limits <- range(unlist(lapply(groups, function(a) a$stats)), finite = TRUE) * 1.06
for (i in seq_along(names)) {
  a <- subset(recovery, scenario == names[i] & model %in% c("varying_spread", "joint", "joint_penalty"))
  group <- interaction(a$model, a$k)
  boxplot(a$explicit - a$true_explicit ~ group, names = rep(c("2", "6", "12"), each = 3L),
          col = rep(colors, 3L), ylab = "Link-scale explicit interaction error", xlab = "Zeta order",
          main = c("Curved turnover only", "Positive interaction", "Negative interaction")[i],
          outline = FALSE, ylim = limits)
  abline(h = 0, lty = 2, col = "#405366")
  if (i == 1L) legend("topright", c("Both extensions", "Both with penalties", "Varying spread only"),
                         fill = colors, bty = "n", cex = .75)
}
invisible(dev.off())
