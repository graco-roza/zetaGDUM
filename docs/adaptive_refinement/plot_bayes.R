folder <- "docs/adaptive_refinement/results"
data <- read.csv(file.path(folder, "bayes_regular_70_weak_summary.csv"))
methods <- c(pooled = "Pooled", average = "Average", block = "Block shrinkage", floor = "Eigenvalue floor")
cols <- c(pooled = "#256384", average = "#956429", block = "#347b63", floor = "#8e526d")
effects <- c("lower_spread", "higher_spread", "paired_spread", "paired_position")
labels <- c("Lower-order spread", "Higher-order spread", "Paired spread change", "Paired position change")
png(file.path(folder, "bayes_coverage.png"), width = 1350, height = 960, res = 145)
par(mfrow = c(2L, 2L), mar = c(4, 4.4, 3, 1), mgp = c(2.5, .7, 0), las = 1)
for (i in seq_along(effects)) {
  d <- subset(data, effect == effects[i] & method %in% names(methods))
  plot(NA_real_, xlim = c(.50, .96), ylim = c(.40, 1.025), xaxt = "n",
       xlab = "Nominal interval coverage", ylab = "Observed coverage", main = labels[i])
  axis(1L, at = c(.5,.8,.9,.95), labels = c("50%", "80%", "90%", "95%"))
  abline(0, 1, col = "#757575", lty = 2)
  for (j in seq_along(methods)) {
    a <- d[d$method == names(methods)[j], ]
    a <- a[order(a$level), ]
    lines(a$level, a$coverage, type = "b", pch = 14 + j, lwd = 2, col = cols[j])
  }
  if (i == 1L) legend("bottomright", unname(methods), col = cols, lty = 1,
                       pch = seq_along(methods) + 14, bty = "n", cex = .77)
}
invisible(dev.off())
