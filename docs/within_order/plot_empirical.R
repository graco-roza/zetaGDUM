# Draw the descriptive figures or render intervals after calibrated inference.
root <- "docs/within_order/results"
manifest <- read.csv("docs/manuscript_validation/communities.csv")
read_profile <- function(kind) {
  files <- file.path(root, paste0("empirical_", kind, "_", manifest$id, ".csv"))
  if (!all(file.exists(files))) stop("Complete all empirical profiles first.")
  do.call(rbind, lapply(files, read.csv))
}
point <- read_profile("importance")
diagnostic <- read_profile("diagnostics")
colors <- c("#116681", "#ba4c29", "#567e44", "#76509c", "#a76b05", "#272d32",
            "#467eb0", "#9e5276", "#006b62")
types <- c(16, 17, 15, 18, 8, 12, 3, 2, 4)
style <- function() {
  graphics::par(mfrow = c(nrow(manifest), 2), mar = c(3.6, 4.4, 2.8, 1.2),
                mgp = c(2.1, .62, 0), cex = .83, las = 1, family = "sans")
}
order_axis <- function(k) {
  graphics::axis(1L, at = sort(unique(k)))
}
annotation <- function(message) {
  xy <- graphics::par("usr")
  graphics::text(mean(xy[1:2]), mean(xy[3:4]), message, col = "#47545c")
}
limits <- function(x) {
  x <- range(c(0, x[is.finite(x)]))
  if (diff(x) <= .Machine$double.eps) c(-1, 1) else
    c(x[1] - diff(x) * .08, x[2] + diff(x) * .08)
}
plot_importance <- function(group, measure, title, y_label) {
  grDevices::cairo_pdf(file.path(root, paste0("empirical_", title, ".pdf")),
    width = 12, height = 15, onefile = TRUE)
  on.exit(grDevices::dev.off())
  style()
  for (id in manifest$id) {
    rows <- point[point$id == id & point$group == group, ]
    labels <- unique(rows$term)
    yrange <- limits(rows[[measure]])
    for (kind in c("incidence", "relative")) {
      d <- rows[rows$response == kind, ]
      r <- diagnostic[diagnostic$id == id & diagnostic$response == kind, ]
      has <- is.finite(r$spatial_D2) & r$spatial_D2 > 0
      graphics::plot(range(r$order), yrange, type = "n", xlab = "Zeta order",
        ylab = y_label, xaxt = "n", main = paste(unique(d$community),
          if (kind == "relative") "Relative abundance" else "Incidence", sep = " | "))
      order_axis(r$order)
      graphics::abline(h = 0, col = "#aaa", lty = "dotted")
      if (!any(is.finite(d[[measure]]))) annotation("No finite allocation at these orders")
      for (j in seq_along(labels)) {
        e <- d[d$term == labels[j], ]
        e <- e[order(e$order), ]
        col <- rep(colors, length.out = length(labels))[j]
        pch <- rep(types, length.out = length(labels))[j]
        graphics::lines(e$order, e[[measure]], type = "b", lwd = 1.5, col = col, pch = pch)
      }
      graphics::rug(r$order[has], col = "#58616a", ticksize = .025, lwd = 1.5)
      if (kind == "incidence") graphics::legend("topright", legend = labels,
        col = rep(colors, length.out = length(labels)), pch = rep(types, length.out = length(labels)),
        horiz = FALSE, ncol = if (length(labels) > 3L) 2L else 1L, cex = .62, bty = "n")
    }
  }
}
plot_importance("predictor", "absolute", "absolute", "Absolute explained deviance")
plot_importance("predictor", "relative", "relative", "Proportion of explained deviance")
plot_importance("component", "absolute", "components", "Absolute explained deviance")
plot_composition <- function() {
  grDevices::cairo_pdf(file.path(root, "empirical_composition.pdf"), width = 12, height = 15)
  on.exit(grDevices::dev.off())
  style()
  for (id in manifest$id) for (kind in c("incidence", "relative")) {
    d <- diagnostic[diagnostic$id == id & diagnostic$response == kind, ]
    d <- d[order(d$order), ]
    graphics::plot(d$order, d$sd^2, type = "b", col = colors[1], pch = 16, lwd = 1.6,
      ylim = c(0, max(diagnostic$sd[diagnostic$id == id])^2 * 1.10),
      xlab = "Zeta order", ylab = "Variance in dissimilarity", xaxt = "n",
      main = paste(unique(d$community),
        if (kind == "relative") "Relative abundance" else "Incidence", sep = " | "))
    order_axis(d$order)
    has <- is.finite(d$spatial_D2) & d$spatial_D2 > 0
    graphics::rug(d$order[has], col = "#58616a", ticksize = .025, lwd = 1.5)
  }
}
plot_composition()
