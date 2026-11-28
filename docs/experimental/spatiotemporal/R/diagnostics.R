# Exact means, approximate rectangle support above budget, and working
# information for one common environmental position contrast across orders.
st_order_diagnostic <- function(x, spatial_max, temporal_max, environment = NULL,
                                budget = 1000L, time_limit = 1e5) {
  p <- .st_panel(x)
  d <- dim(p$x)
  budget <- .st_one(budget, .Machine$integer.max, "budget")
  tab <- st_zeta_retention(st_zeta_decline(p, spatial_max, temporal_max, time_limit))
  if (!is.null(environment)) {
    if (!is.matrix(environment) || !is.numeric(environment) ||
        !identical(dim(environment), d[1:2]) || any(!is.finite(environment)) ||
        stats::sd(as.vector(environment)) <= 0) stop("Provide one varying, aligned sites x times contrast.")
    environment <- matrix(scale(as.vector(environment)), d[1], d[2])
  }
  columns <- c("rectangles_evaluated", "zero_fraction", "zero_mc_se", "design_rank",
               "minimum_eigenvalue", "condition", "position_information", "site_jackknife_se")
  diagnostic <- t(vapply(seq_len(nrow(tab)), function(r) {
    k <- tab$spatial_order[r]
    m <- tab$temporal_order[r]
    exhaustive <- tab$rectangles[r] <= budget
    if (exhaustive) {
      a <- .st_indices(d[1], k, budget)
      b <- .st_indices(d[2], m, budget)
      sites <- a[rep(seq_len(nrow(a)), times = nrow(b)), , drop = FALSE]
      times <- b[rep(seq_len(nrow(b)), each = nrow(a)), , drop = FALSE]
    } else {
      sites <- t(replicate(budget, sample.int(d[1], k)))
      times <- t(replicate(budget, sample.int(d[2], m)))
      sites <- matrix(sites, budget, k)
      times <- matrix(times, budget, m)
    }
    y <- vapply(seq_len(nrow(sites)), function(i) {
      sum(colSums(p$x[sites[i, ], times[i, ], , drop = FALSE] > 0, dims = 2L) == k * m)
    }, numeric(1))
    zero <- mean(y == 0)
    out <- setNames(c(length(y), zero, if (exhaustive) 0 else sqrt(zero * (1 - zero) / length(y)),
                      rep(NA_real_, 5L)), columns)
    if (is.null(environment)) return(out)
    contrast <- vapply(seq_len(nrow(sites)), function(i) {
      mean(environment[sites[i, ], times[i, ], drop = FALSE])
    }, numeric(1))
    X <- cbind(intercept = 1, position = contrast)
    opportunities <- rep(d[3], length(y))
    fit <- .st_fit_risk(X, y, opportunities, ridge = 0.1)
    out["design_rank"] <- fit$rank
    if (fit$status == "boundary") {
      out["minimum_eigenvalue"] <- out["position_information"] <- 0
    }
    if (fit$status != "estimated") return(out)
    mu <- .st_predict_risk(fit, X)
    w <- mu * (1 - mu)
    out["minimum_eigenvalue"] <- min(fit$curvature)
    out["condition"] <- max(fit$curvature) / min(fit$curvature)
    out["position_information"] <- mean(w * (contrast - stats::weighted.mean(contrast, w))^2)
    if (exhaustive && k < d[1] && d[1] >= 3L) {
      slopes <- vapply(seq_len(d[1]), function(i) {
        keep <- rowSums(sites == i) == 0
        .st_fit_risk(X[keep, , drop = FALSE], y[keep], opportunities[keep], 0.1)$coefficients[2]
      }, numeric(1))
      out["site_jackknife_se"] <- sqrt((d[1] - 1) / d[1] * sum((slopes - mean(slopes))^2))
    }
    out
  }, setNames(numeric(length(columns)), columns)))
  out <- cbind(tab, as.data.frame(diagnostic))
  ref <- out$position_information[out$spatial_order == 1 & out$temporal_order == 1]
  out$relative_position_information <- .st_ratio(out$position_information, ref)
  out
}
