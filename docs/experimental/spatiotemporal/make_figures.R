# Run after run_simulations.R from the repository root.
root <- "docs/experimental/spatiotemporal"
source(file.path(root, "load.R"))
st_load(root)
output <- Sys.getenv("ST_OUTPUT", root)
dir.create(file.path(output, "figures"), recursive = TRUE, showWarnings = FALSE)
read_result <- function(name) utils::read.csv(file.path(output, "results", paste0(name, ".csv")))
draw <- function(name, width, height, expr) {
  grDevices::svg(file.path(output, "figures", paste0(name, ".svg")), width, height,
                 family = "sans", pointsize = 12)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4.3, 4.4, 3, 1), mgp = c(2.5, 0.7, 0), las = 1, bty = "l",
                col.axis = "#2B3540", col.lab = "#192934", col.main = "#192934")
  force(expr)
}
ink <- .st_colours
surface <- read_result("order_surface")
draw("01-decline-retention", 13.2, 4.1, {
  par(mfrow = c(1, 3), mar = c(3.8, 4.1, 3.9, 0.8))
  plot_st_surface(surface, "zeta", "Mean shared richness", main = "A   Joint decline")
  plot_st_surface(surface, "spatial_retention", "Adjacent-order ratio", limits = c(0, 1),
                  main = "B   Spatial retention")
  plot_st_surface(surface, "temporal_retention", "Adjacent-order ratio", limits = c(0, 1),
                  main = "C   Temporal retention")
})
draw("02-order-information", 13.2, 4.1, {
  par(mfrow = c(1, 3), mar = c(3.8, 4.1, 3.9, 0.8))
  plot_st_surface(surface, "zero_fraction", "Fraction of rectangles with zero zeta", limits = c(0, 1),
                  main = "A   Zero sharing")
  plot_st_surface(surface, "effective_species", "Inverse concentration of contributions", digits = 1,
                  main = "B   Effective species")
  plot_st_surface(surface, "relative_position_information", "Relative to order (1, 1)", limits = c(0, 1),
                  main = "C   Position information")
})

scenarios <- read_result("scenario_zeta")
draw("03-directional-sharing", 12, 4, {
  par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3, 0.9))
  labels <- c(gains_only = "A   Gains only", losses_only = "B   Losses only", redistribution = "C   Redistribution")
  for (process in names(labels)) {
    y <- scenarios[scenarios$scenario == process, ]
    z <- aggregate(cbind(gains, losses, change) ~ spatial_order, y, mean)
    plot(z$spatial_order, z$gains, type = "n", ylim = range(z[c("gains", "losses", "change")], 0),
         main = labels[process], xlab = "Spatial order", ylab = "Mean shared-species events")
    abline(h = 0, col = "#C5CCD0")
    points(z$spatial_order, z$gains, type = "b", pch = 16, col = ink["blue"], lwd = 2)
    points(z$spatial_order, z$losses, type = "b", pch = 1, col = ink["orange"], lty = 2, lwd = 2)
    points(z$spatial_order, z$change, type = "b", pch = 3, col = ink["grey"], lty = 3, lwd = 2)
    if (process == "gains_only") legend("topright", c("Gains", "Losses", "Net change"), pch = c(16, 1, 3),
                                        lty = 1:3, lwd = 2, col = ink[c("blue","orange","grey")],
                                        bty = "n", cex = 0.82)
  }
})

draw("04-occupancy-abundance", 11, 8.2, {
  par(mfrow = c(2, 2), mar = c(4.4, 4.6, 2.8, 1.0), oma = c(1.4, 0, 0, 0))
  labels <- c(expansion = "A   Occupancy expansion", contraction = "B   Occupancy contraction",
              abundance_only = "C   Positive abundance change", combined = "D   Expansion, redistribution and abundance")
  for (process in names(labels)) {
    y <- read_result(paste0("species_", process))
    x <- y$occupancy_change
    a <- y$persistent_conditional_change
    mark <- ifelse(y$redistribution, 1, 16)
    plot(x, a, type = "n", xlim = c(-30, 30),
         ylim = c(-26, 26), main = labels[process],
         ylab = "Mean count change at persistent sites", xlab = "Net change in occupied sites")
    abline(h = 0, v = 0, col = "#C5CCD0")
    points(x, a, pch = mark, cex = 0.85, lwd = 1.2, col = ifelse(mark == 1, ink["orange"], ink["blue"]))
    if (process == "expansion") legend("topleft", c("Occupancy events in one direction or neither",
                                                    "Both newly occupied and lost sites"),
                                       pch = c(16, 1), col = ink[c("blue","orange")], bty = "n", cex = 0.78)
  }
  mtext("Points represent species in one replicate. Coincident points can overlap.", side = 1,
        outer = TRUE, line = 0.1, cex = 0.77)
})

recovery <- read_result("parameter_recovery")
draw("05-effect-recovery", 9.5, 5.9, {
  keys <- c("gains_only/all/gain/position_change.gradient" = "Gains / environmental change",
            "losses_only/low_persistence/loss/position_change.gradient" = "Losses / environmental change",
            "interaction/all/gain/position_by_change.gradient" = "Gains / initial position × change",
            "uncommon_effect/target/gain/position_change.gradient" = "Gains / uncommon species",
            "widespread_effect/target/gain/position_change.gradient" = "Gains / widespread species",
            "persistent_effect/target/gain/position_change.gradient" = "Gains / persistent species",
            "abundance_only/all/positive_abundance/position_change.gradient" = "Positive abundance / environmental change",
            "null/all/gain/position_change.gradient" = "Null / gains",
            "detection_failure/all/gain/position_change.gradient" = "Imperfect detection / apparent gains")
  key <- with(recovery, paste(scenario, pool, event, term, sep = "/"))
  averages <- t(vapply(names(keys), function(k) {
    err <- recovery$error[key == k & is.finite(recovery$error)]
    c(bias = mean(err), margin = 1.96 * sd(err) / sqrt(length(err)))
  }, numeric(2)))
  par(mar = c(4.4, 19.5, 3.3, 1.3))
  y <- rev(seq_along(keys))
  x <- averages[, "bias"]
  margin <- averages[, "margin"]
  plot(x, y, type = "n", xlim = range(x - margin, x + margin, 0), yaxt = "n",
       xlab = "Mean estimated coefficient minus known coefficient",
       ylab = "", main = "Recovery error across simulated panels")
  axis(2, y, unname(keys), las = 1, tick = FALSE, cex.axis = 0.83)
  abline(v = 0, col = "#687986", lty = 2)
  segments(x - margin, y, x + margin, y, col = ink["grey"], lwd = 1.7)
  points(x, y, pch = ifelse(y == 1, 1, 16), cex = 1.05,
         col = ifelse(y == 1, ink["orange"], ink["blue"]))
  mtext("Intervals show Monte Carlo uncertainty in the mean", side = 3, line = 0.25, cex = 0.76)
})

draw("06-environmental-functions", 12.5, 8.5, {
  par(mfrow = c(2, 2))
  s <- simulate_st_gdum("spatial_filtering", n_sites = 44, n_species = 64, seed = 14)
  fit <- st_spatial_gdum(s$latent, s$env)[[1]]
  basis <- fit$design$environment_basis$splines[[1]]
  bounds <- attr(basis, "Boundary.knots")
  spread <- seq(0, 2 * min(abs(bounds)), length.out = 101)
  # Matching opposite centered positions holds summed position at zero.
  difference <- predict(basis, spread / 2) - predict(basis, -spread / 2)
  y <- plogis(fit$fit$intercept + difference %*% fit$fit$beta[1:3])
  plot(spread, y, type = "l", lwd = 2.4, col = ink["blue"], ylim = c(0, 1),
       xlab = "Environmental separation (SD)", ylab = "Simpson dissimilarity",
       main = "A   Spatial spread")
  s <- simulate_st_gdum(c("interaction", "losses_only"), n_sites = 50, n_species = 72, seed = 190)
  d <- st_design(s$latent, s$env)
  term <- c("position_initial.gradient", "position_change.gradient", "position_by_change.gradient")
  f <- st_gdum(s$latent, d, term, loss_terms = term[1:2])
  delta <- seq(-1.3, 1.3, length.out = 121)
  grid <- function(position) {
    z <- d
    z$X <- cbind(position, delta, position * delta)
    colnames(z$X) <- term
    z
  }
  plot_st_gdum(f, grid(-1), delta, ylim = c(0, 0.75), xlab = "Environmental change (SD)",
               main = "B   Gains by initial position")
  plot_st_gdum(f, grid(1), delta, lines_only = TRUE, col = ink["orange"], lty = 2)
  legend("topleft", c("Initial position −1 SD", "Initial position +1 SD"), bty = "n",
         col = ink[c("blue", "orange")], lty = 1:2, lwd = 2, cex = 0.85)
  plot_st_gdum(f, grid(0), delta, event = "loss", ylim = c(0, 0.75), xlab = "Environmental change (SD)",
               main = "C   Losses at mean initial position")
  a <- simulate_st_gdum("abundance_only", abundance = "gamma", seed = 90)
  da <- st_design(a$latent, a$env)
  fa <- st_abundance_gdum(a$latent, da, term[2])
  mu <- exp(fa$fit$coefficients[1] + fa$fit$coefficients[2] * delta)
  plot(delta, mu, type = "l", lwd = 2.4, col = ink["green"], ylim = c(0, max(mu)),
       xlab = "Environmental change (SD)", ylab = "Mean abundance given presence",
       main = "D   Conditional abundance")
})

draw("07-detection-and-prediction", 10.8, 4.1, {
  par(mfrow = c(1, 2), mar = c(4.3, 4.4, 3, 1.1))
  a <- scenarios[scenarios$scenario == "detection_failure", ]
  z <- aggregate(cbind(gains, losses) ~ observation + spatial_order, a, mean)
  obs <- z[z$observation == "observed", ]
  ref <- z[z$observation == "matched_latent", ]
  plot(obs$spatial_order, obs$gains, type = "b", pch = 16, col = ink["blue"], lwd = 2,
       ylim = range(0, z$gains, z$losses), xlab = "Spatial order", ylab = "Mean shared-species events",
       main = "A   Static occupancy with imperfect detection")
  lines(obs$spatial_order, obs$losses, type = "b", pch = 1, col = ink["orange"], lty = 2, lwd = 2)
  lines(ref$spatial_order, ref$gains, col = ink["grey"], lty = 3, lwd = 2)
  legend("topright", c("Apparent gains", "Apparent losses", "True gains and losses"),
         col = ink[c("blue", "orange", "grey")], lty = 1:3, lwd = 2, bty = "n", cex = 0.85)
  h <- read_result("site_holdout")
  h <- aggregate(gain_log_loss ~ replicate + spatial_order + architecture, h, mean)
  wide <- reshape(h, idvar = c("replicate", "spatial_order"), timevar = "architecture", direction = "wide")
  wide$improvement <- wide$gain_log_loss.additive - wide$gain_log_loss.interaction
  boxplot(improvement ~ spatial_order, wide, col = "#D4E7F0", border = ink["blue"],
          xlab = "Spatial order", ylab = "Improvement in held-site gain log loss",
          main = "B   Predictive value of the interaction")
  abline(h = 0, lty = 2, col = ink["grey"])
})

draw("08-overlap", 10.5, 4.4, {
  par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3, 1))
  a <- read_result("overlap_dependence")
  a$ratio <- a$site_se / a$independent_se
  a$group <- factor(paste("Order", a$spatial_order, a$event),
                    levels = c("Order 1 gain", "Order 1 loss", "Order 2 gain", "Order 2 loss"))
  labels <- c("1 Gain", "1 Loss", "2 Gain", "2 Loss")
  boxplot(ratio ~ group, a, names = labels, col = rep(c("#D4E7F0", "#F5D9CD"), 2),
          border = ink["grey"], xlab = "Spatial order and event", ylab = "Site SE / independent-trial SE",
          main = "A   Uncertainty from overlapping sets")
  abline(h = 1, lty = 2, col = ink["grey"])
  a$independent <- abs(a$estimate) > 1.96 * a$independent_se
  a$site <- abs(a$estimate) > 1.96 * a$site_se
  rate <- t(aggregate(cbind(independent, site) ~ group, a, mean)[, -1])
  barplot(rate, beside = TRUE, col = c("#9CC9DE", "#F1B797"), border = ink["grey"],
          names.arg = labels, ylim = c(0, 0.55), xlab = "Spatial order and event",
          ylab = "Fraction of null slopes rejected", main = "B   False attribution under a clean null")
  abline(h = 0.05, lty = 2, col = ink["grey"])
  legend("topleft", c("Independent-trial SE", "Delete-site SE"), bty = "n",
         fill = c("#9CC9DE", "#F1B797"), border = ink["grey"], cex = 0.82)
})
