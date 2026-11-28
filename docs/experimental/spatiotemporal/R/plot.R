.st_colours <- c(blue = "#0072B2", orange = "#D55E00", green = "#00866B", grey = "#4F5C68")

plot_st_surface <- function(surface, field = "zeta", legend = field, digits = 2L,
                            limits = NULL, main = NULL) {
  k <- sort(unique(surface$spatial_order))
  m <- sort(unique(surface$temporal_order))
  z <- matrix(NA_real_, length(k), length(m))
  z[cbind(match(surface$spatial_order, k), match(surface$temporal_order, m))] <- surface[[field]]
  if (is.null(limits)) {
    finite <- z[is.finite(z)]
    limits <- if (length(finite)) range(finite) else c(0, 1)
  }
  if (diff(limits) == 0) limits <- limits + c(0, 1)
  pal <- grDevices::colorRampPalette(c("#F0F7FC", "#80B5CE", "#00527F"))(100)
  graphics::plot(NA, xlim = c(0.5, length(k) + 0.5), ylim = c(0.5, length(m) + 0.5),
                 xaxs = "i", yaxs = "i", axes = FALSE, xlab = "Spatial order", ylab = "Temporal order",
                 main = main)
  graphics::axis(1, seq_along(k), k, tick = FALSE)
  graphics::axis(2, seq_along(m), m, tick = FALSE, las = 1)
  for (i in seq_along(k)) {
    for (j in seq_along(m)) {
      v <- z[i, j]
      scaled <- if (is.finite(v)) max(0, min(1, (v - limits[1]) / diff(limits))) else NA_real_
      colour <- if (is.na(scaled)) "#F1F0ED" else pal[1L + floor(99 * scaled)]
      graphics::rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, col = colour, border = "white")
      graphics::text(i, j, if (is.finite(v)) formatC(v, digits = digits, format = "f") else "–",
                     cex = 0.84, col = if (is.finite(scaled) && scaled > 0.67) "white" else "#152734")
    }
  }
  graphics::mtext(legend, side = 3, line = 0.3, cex = 0.77)
  invisible(z)
}

# Conditional effect on explicitly supplied counterfactual design rows.
plot_st_gdum <- function(object, design = object$design, x, event = c("gain", "loss"),
                         lines_only = FALSE, col = .st_colours["blue"], lty = 1L,
                         xlab = "Environmental contrast", ...) {
  event <- match.arg(event)
  X <- .st_model_matrix(design, object$terms[[event]])
  if (!is.numeric(x) || length(x) != nrow(X) || any(!is.finite(x))) stop("Supply a value for each design row.")
  y <- .st_predict_risk(object$fit[[event]], X)
  i <- order(x)
  if (lines_only) {
    graphics::lines(x[i], y[i], col = col, lty = lty, lwd = 2, ...)
  } else {
    graphics::plot(x[i], y[i], type = "l", col = col, lty = lty, lwd = 2,
                   xlab = xlab, ylab = paste("Conditional", event, "probability"), ...)
  }
  invisible(data.frame(x = x, probability = y))
}
