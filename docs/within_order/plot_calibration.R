source("docs/within_order/simulation_design.R")
a <- read.csv(file.path(within_dir, "calibration_summary.csv"))
a <- subset(a, scenario == "baseline" & B == 200L & draw_budget == 79L)
methods <- c("bayesian_bootstrap", "gaussian_gmm", "multiplier", "site_bootstrap")
colors <- c("#543F72", "#226B8A", "#A45131", "#39765C")
method_names <- c("Direct Bayesian bootstrap", "Gaussian GMM", "Score multiplier", "Site bootstrap")
effects <- c("strong","weak","irrelevant","strong - weak")
pdf(file.path(within_dir, "calibration_coverage.pdf"), width = 9, height = 6)
par(mfrow = c(2L, 2L), mar = c(4.5, 4.1, 2.2, .8), mgp = c(2.6, .65, 0))
stopifnot(nrow(a) > 0L)
for (j in seq_along(effects)) {
  d <- a[a$estimand == effects[j], ]
  if (nrow(d)) {
    keys <- unique(d[c("n", "order", "response")])
    keys <- keys[order(keys$order, keys$n, keys$response), ]
    x <- match(paste(d$n, d$order, d$response), paste(keys$n, keys$order, keys$response))
    xlabel <- paste0("k", keys$order, ", ", keys$n, "\n", substr(keys$response, 1L, 3L))
  } else {
    keys <- data.frame(n=0)
    xlabel <- ""
    x <- numeric()
  }
  plot(NA_real_, xlim = c(.5, nrow(keys) + .5), ylim = c(0, 1.02),
       xaxt = "n", xlab = "Order, sites and response", ylab = "Empirical coverage",
       main = effects[j])
  axis(1, at = seq_len(nrow(keys)), labels = xlabel, cex.axis = .77)
  abline(h = c(.8,.95), lty = 2L, col = "grey60")
  for (m in seq_along(methods)) {
    use <- d$method == methods[m]
    at <- x[use] + (m - 2.5) * .11
    points(at, d$coverage80[use], pch = 1L, col = colors[m])
    points(at, d$coverage95[use], pch = 16L, col = colors[m])
  }
  if (j == 1L) legend("bottomleft", method_names[methods %in% a$method],
    col = colors[methods %in% a$method], pch = 16L, bty = "n", cex = .67)
}
invisible(dev.off())
