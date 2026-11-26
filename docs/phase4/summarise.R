source("docs/phase4/common.R")
root <- "docs/phase4/results"
stack <- function(pattern, filename) {
  paths <- list.files("docs/phase4/runs", pattern, full.names = TRUE)
  output <- file.path(root, filename)
  if (!length(paths)) return(read.csv(output))
  x <- do.call(rbind, lapply(paths, read.csv, na.strings = c("", "NA")))
  write.csv(x, output, row.names = FALSE)
  x
}
grouped <- function(x, keys, fn) {
  ids <- split(seq_len(nrow(x)), interaction(x[keys], drop = TRUE))
  a <- lapply(ids, function(i) cbind(x[i[1L], keys, drop = FALSE], as.data.frame(as.list(fn(x[i, ])))))
  y <- do.call(rbind, a)
  y[do.call(order, y[keys]), , drop = FALSE]
}
profile <- stack("^profile_[NS]", "empirical_profiles.csv")
precision <- stack("^precision_[NS]", "empirical_precision.csv")
boot <- stack("^bootstrap_", "empirical_bootstrap.csv")
invisible(stack("^pilot_[NS]", "empirical_pilots.csv"))
invisible(stack("^composition_", "simulated_communities.csv"))
invisible(stack("^profile_sim", "simulated_profiles.csv"))
invisible(stack("^precision_sim", "simulated_precision.csv"))
invisible(stack("^pilot_sim", "simulated_pilots.csv"))
invisible(stack("^sensitivity_[NS]", "sensitivity_profiles.csv"))
invisible(stack("^sensitivity_precision_", "sensitivity_precision.csv"))
both <- merge(precision, boot, by = c("community", "component", "n", "predictor", "order", "preceding", "effect"))
both$ratio <- both$bootstrap_se / both$jackknife_1
both$contrast <- ifelse(is.na(both$preceding), "level", "change")
tab <- grouped(both, c("component", "effect", "contrast"), function(d) {
  val <- d$ratio[is.finite(d$ratio)]
  c(n = length(val), median_ratio = median(val), median_JK = median(d$jackknife_1, na.rm = TRUE),
    median_BB = median(d$bootstrap_se, na.rm = TRUE))
})
write.csv(tab, file.path(root, "summary_intervals.csv"), row.names = FALSE)

metric <- function(d) {
  cpu <- median(d$cpu)
  error <- d$point - d$reference
  use <- is.finite(error) & as.logical(d$regular)
  rmse <- sqrt(mean(error[use]^2))
  c(valid = sum(use), RMSE = rmse, bias = mean(error[use]), CPU = cpu,
    RMSE_sqrt_CPU = rmse * sqrt(cpu), sharing = mean(d$sharing), shared_draws = mean(d$shared_draws),
    min_exclusion_ESS = median(d$exclusion_ess), fallback = mean(d$fallback),
    mean_mass_draws = mean(d$mass_draws), mean_proposal_count = mean(d$proposals))
}
comp <- stack("^comparison_", "sampler_replicates.csv")
keys <- c("case", "n", "k", "component", "M", "sampler", "effect")
if (!is.null(comp)) {
  tab <- grouped(comp, keys, metric)
  write.csv(tab, file.path(root, "summary_sampler.csv"), row.names = FALSE)
  same <- c("case", "n", "component", "M", "sampler", "r", "effect")
  high <- ave(comp$k, comp$case, FUN = max) == comp$k
  paired <- merge(comp[high, ], comp[!high, ], by = same, suffixes = c("", "_lo"))
  paired$point <- paired$point - paired$point_lo
  paired$reference <- paired$reference - paired$reference_lo
  paired$regular <- as.logical(paired$regular) & as.logical(paired$regular_lo)
  paired$cpu <- paired$cpu + paired$cpu_lo
  paired$sharing <- paired$shared_draws <- paired$exclusion_ess <- NA_real_
  paired$fallback <- as.logical(paired$fallback) | as.logical(paired$fallback_lo)
  paired$mass_draws <- paired$mass_draws + paired$mass_draws_lo
  paired$proposals <- paired$proposals + paired$proposals_lo
  write.csv(grouped(paired, keys, metric), file.path(root, "summary_paired.csv"), row.names = FALSE)
}

p_png <- function(filename, width, height) png(file.path(root, filename), width, height, res = 145)
manifest <- read_manifest()
colors <- c(total = "#914b10", relative = "#087e84")
p_png("empirical_profiles.png", 1450, 1720)
par(mfrow = c(5, 2), mar = c(2.25, 4.6, 2.3, 1), oma = c(2.8, 0.3, 3, 0.3),
    cex.axis = .83, mgp = c(2.5, .7, 0))
for (i in seq_len(nrow(manifest))) {
  m <- manifest[i, ]
  a <- profile[profile$community == m$community & profile$order < m$sites, ]
  u <- a[a$component == "total", ]
  plot(u$order, u$sharing_estimate, type = "n", ylim = c(0, 1), xlab = "",
       ylab = "Probability of sharing", main = m$community, xaxt = "n")
  axis(1, u$order)
  polygon(c(u$order, rev(u$order)), c(u$sharing_lower, rev(u$sharing_upper)),
          col = "#cad8ec", border = NA)
  abline(v = m$high, lty = 3, col = "#6b6b6b")
  lines(u$order, u$sharing_estimate, type = "b", pch = 16, col = "#325584")
  if (i == 1L) legend("bottomleft", c("Estimated", "Analytical bounds"), pch = c(16, 15),
                       col = c("#325584", "#cad8ec"), bty = "n", cex = .8)
  floor <- 1e-5
  plot(u$order, pmax(u$response_sd, floor), ylim = c(floor, .25), log = "y", type = "n",
       ylab = "Dissimilarity SD", xlab = "", main = if (i == 1L) "Observed dissimilarity" else "", xaxt = "n")
  axis(1, u$order)
  abline(h = .005, lty = 2, col = "#626262")
  abline(v = m$high, lty = 3, col = "#6b6b6b")
  for (component in names(colors)) {
    g <- a[a$component == component, ]
    lines(g$order, pmax(g$response_sd, floor), type = "b",
          pch = ifelse(g$pass, 19, 1), col = colors[component], lwd = 1.2)
  }
  if (i == 1L) legend("bottomleft", c("Normalised incidence", "Relative abundance"), pch = 1,
                       col = colors, bty = "n", cex = .8)
}
mtext("Sharing and response variation across manuscript communities", line = 1.25,
      side = 3, outer = TRUE, font = 2, cex = 1.35)
mtext("Zeta order", line = .15, side = 1, outer = TRUE, cex = .9)
mtext("Filled response points passed. Dashed lines show the SD tolerance. Vertical lines mark the highest manuscript orders.",
      line = 1.35, side = 1, outer = TRUE, cex = .77)
invisible(dev.off())

verified <- read.csv(file.path(root, "verification_profile.csv"))
v <- read.csv(file.path(root, "verification_precision.csv"))
p <- read.csv(file.path(root, "verification_pilot.csv"))
v <- v[v$effect %in% c("spread", "position") & (v$order == 3 | !is.na(v$preceding)), ]
p_png("budget_verification.png", 1300, 925)
par(mfrow = c(2, 2), mar = c(4, 4.5, 2.6, .9), oma = c(.5, .7, 3, .2), mgp = c(2.6, .7, 0))
for (i in seq_len(nrow(v))) {
  a <- v[i, ]
  use <- p$order == a$order & p$effect == a$effect &
    (is.na(p$preceding) == is.na(a$preceding))
  j <- p[use, ]
  x <- sort(unique(c(j$budget, a$candidate_budget)))
  target <- .30 * a$ecological_se
  required <- verified$candidate_budget[match(a$order, verified$order)]
  g <- seq(min(x), max(x), length.out = 100)
  label <- if (is.na(a$preceding)) "Order 3" else "Order 5 minus 3"
  plot(j$budget, sqrt(j$mc_variance), log = "xy", main = paste(label, a$effect, sep = " · "),
       ylim = range(c(sqrt(j$mc_variance), a$confirmation_mc_upper, target)),
       xlab = "Combinations", ylab = "Monte Carlo SD (logit contrast)", pch = 16, col = "#325584")
  lines(g, sqrt(a$A / g), lty = 3, col = "#325584")
  abline(h = target, lty = 2, col = "#333333")
  points(required, a$confirmation_mc_upper, pch = 18, cex = 1.6,
         col = if (a$confirmation_pass) "#087e84" else "#ab4a22")
}
mtext("Independent verification of the predicted combination budget", side = 3,
      outer = TRUE, line = 1.3, font = 2, cex = 1.3)
invisible(dev.off())
