#' Pre-fit compositional and environmental information by zeta order
#'
#' Profiles the finite community under uniform sampling of distinct sites.
#' Relative information divides each response variance by its order-two
#' variance. It can exceed one, rise with order, or be undefined. Exact zeta
#' moments and means use [zeta_order_support()]. McGeoch et al. (2019,
#' doi:10.1002/ecy.2832) established incidence-based zeta variance.
#'
#' Incidence is one minus shared richness divided by mean set richness.
#' Relative abundance is one minus the sum of within-set minima of species
#' proportions. Its occupancy component sums mean proportions for species
#' missing from at least one site. Its heterogeneity component sums mean
#' minus minimum proportions for species occurring at every site. The two
#' components add to relative-abundance dissimilarity. Their variances also
#' depend on their covariance. These are descriptive partitions of the response.
#'
#' Environmental position is a set mean. Linear and nonlinear position
#' variances are computed with the finite-population correction. Spread is
#' the mean pairwise difference, with exact variance from edge inclusion
#' probabilities. Joint rank and conditioning use centred, unit-scaled
#' design columns. Geographic splines use a common distance range across orders.
#' Raw predictor-specific variances and the nonlinear basis variances are
#' returned in the `environment_detail` attribute.
#'
#' The pre-fit flags identify provable degeneracy and otherwise indicate
#' candidates for model checking. They do not establish ecological inference
#' or set a transferable threshold. Use [ugdm_within_order_diagnostic()] for
#' post-fit confirmation with independent prediction and calibrated intervals.
#' Numerical standard errors quantify combination integration on fixed sites.
#'
#' @param community Numeric sites by species matrix. Every site must have
#'   a positive finite total for the abundance diagnostics.
#' @param environment Optional numeric sites by predictors data frame or matrix
#'   in community order. Predictors with fewer than four distinct values
#'   are reported as unsupported by the spline design.
#' @param coordinates Optional two-column site coordinates, used with
#'   environmental predictors for the geographic distance term.
#' @param orders Distinct orders from two to the number of sites. The default
#'   stops at the smaller of 30 and the number of sites.
#' @param response One of `"relative"`, `"incidence"`, `"occupancy"`,
#'   `"heterogeneity"`, or `"all"`.
#' @param method `"auto"` uses analytical moments and enumerates small
#'   populations when other quantities require integration. `"sample"` uses
#'   independent pools for non-analytical quantities. `"exact"` requires
#'   enumeration of quantities without an analytical result, within `exact_limit`.
#' @param combinations Number of IID uniform combinations per numerical pool.
#' @param pools Number of independent numerical pools, at least two.
#' @param exact_limit Maximum number of combinations to enumerate at an order.
#' @param threshold_limit Maximum number of nonzero abundance thresholds for
#'   exact covariance of minimum-abundance indicators. Zero disables that
#'   calculation. Its memory use is quadratic in the retained thresholds.
#' @param seed Optional random seed.
#' @param lonlat Interpret coordinates as longitude and latitude in degrees.
#' @return An order data frame of class `zeta_order_information`. Multiple
#'   responses give one row per response and order with shared ecological and
#'   environmental columns. The `candidate_orders` attribute contains the
#'   pre-fit candidate orders by response. Exact numerical errors are zero.
#'   Sharing bounds and the sharing interval are both reported.
#' @export
zeta_order_information <- function(community, environment = NULL, coordinates = NULL,
    orders = 2:30, response = "relative", method = "auto", combinations = 1000L,
    pools = 4L, exact_limit = 5000L, seed = NULL, lonlat = FALSE, threshold_limit = 256L) {
  types <- c("incidence", "relative", "occupancy", "heterogeneity")
  response <- match.arg(response, c(types, "all"))
  method <- match.arg(method, c("auto", "sample", "exact"))
  selected <- if (response == "all") types else response
  if (anyNA(community)) stop("Community values must be finite.")
  x <- .conditional_community(community)
  q <- .conditional_community(x, relative = TRUE)
  n <- nrow(x)
  if (n < 2L) stop("At least two sites are required.")
  if (missing(orders)) orders <- seq.int(2L, min(30L, n))
  if (!length(orders) || anyDuplicated(orders)) stop("Provide distinct orders.")
  orders <- vapply(orders, .zeta_order_integer, integer(1), n = n, minimum = 2L)
  combinations <- .set_design_integer(combinations, "combinations", 2L)
  pools <- .set_design_integer(pools, "pools", 2L)
  exact_limit <- .set_design_integer(exact_limit, "exact_limit", 1L)
  threshold_limit <- .set_design_integer(threshold_limit, "threshold_limit", 0L)
  if (!is.null(seed)) set.seed(seed)
  grid <- sort(unique(c(2L, orders)))
  ecological <- zeta_order_support(x, seq.int(2L, max(grid)))
  breakpoint <- zeta_break(data.frame(k = ecological$order, zeta = ecological$expected_richness,
                                     retention = ecological$incidence_retention))
  env <- .order_environment(environment, coordinates, n, lonlat)
  presence <- (x > 0) * 1L
  richness <- rowSums(presence)
  occupancy <- colSums(presence)
  constant_q <- all(q == rep(q[1L, ], each = n))
  constant_positive <- all(vapply(seq_len(ncol(q)), function(s) {
    v <- q[presence[, s] > 0L, s]
    length(v) <= 1L || all(v == v[1L])
  }, logical(1L)))
  threshold <- if (!constant_q && any(selected != "incidence")) {
    .order_minimum_thresholds(q, threshold_limit)
  }
  rows <- vector("list", length(grid))
  detail <- vector("list", length(grid))

  for (j in seq_along(grid)) {
    k <- grid[j]
    supporting <- ecological[k - 1L, ]
    moments <- .order_response_exact(q, richness, occupancy, constant_q, constant_positive,
                                     threshold, k, supporting)
    exact_sharing <- supporting$sharing_lower == supporting$sharing_upper
    # Endpoints with one possible set are cheap even if all quantities are analytical.
    can_enumerate <- lchoose(n, k) <= log(exact_limit) + 1e-12
    need_pool <- !exact_sharing || anyNA(moments$variance[selected]) || (!is.null(env$basis) && k < n)
    enumerate <- can_enumerate && method != "sample" && need_pool
    if (need_pool && method == "exact" && !enumerate) {
      stop("Order ", k, " requires enumeration above exact_limit.")
    }
    indices <- if (!need_pool) list() else if (enumerate) {
      list(t(utils::combn(n, k)))
    } else lapply(seq_len(pools), function(i) .iid_uniform_sets(n, k, combinations))
    take <- lapply(indices, .order_information_pool, q = q, presence = presence,
                   richness = richness, env = env, types = selected)
    summaries <- lapply(take, "[[", "response")
    variances <- moments$variance[selected]
    averages <- moments$mean[selected]
    source <- ifelse(is.finite(variances), "analytical", if (enumerate) "enumerated" else "sampled")
    standard_error <- stats::setNames(rep(0, length(selected)), selected)
    if (length(summaries)) {
      missing <- is.na(variances)
      if (any(missing)) {
        v <- do.call(cbind, lapply(summaries, function(z) apply(z, 2L, .order_variance, exact = enumerate)))
        variances[missing] <- rowMeans(v)[missing]
        if (!enumerate) standard_error[missing] <- (apply(v, 1L, stats::sd) / sqrt(pools))[missing]
      }
      missing_mean <- is.na(averages)
      if (any(missing_mean)) {
        a <- do.call(cbind, lapply(summaries, colMeans))
        averages[missing_mean] <- rowMeans(a)[missing_mean]
      }
    }
    sharing <- .order_sharing(take, supporting, exact_sharing, enumerate)
    geometry <- .order_geometry(env, k, n, take)
    detail[[j]] <- geometry$detail
    rows[[j]] <- data.frame(order = k, response = selected,
      mean_zeta = ifelse(selected == "incidence", supporting$expected_richness, supporting$abundance_zeta),
      retention_rate = ifelse(selected == "incidence", supporting$incidence_retention, supporting$abundance_retention),
      incidence_zeta = supporting$expected_richness, abundance_zeta = supporting$abundance_zeta,
      incidence_retention = supporting$incidence_retention, abundance_retention = supporting$abundance_retention,
      sharing_probability = sharing[["estimate"]], sharing_se = sharing[["se"]],
      sharing_interval_lower = sharing[["lo"]], sharing_interval_upper = sharing[["hi"]],
      sharing_lower = supporting$sharing_lower, sharing_upper = supporting$sharing_upper,
      effective_incidence_species = supporting$effective_incidence_species,
      effective_abundance_species = supporting$effective_abundance_species,
      effective_species = ifelse(selected == "incidence", supporting$effective_incidence_species,
                                 supporting$effective_abundance_species),
      shared_richness_variance = supporting$richness_variance, response_mean = unname(averages),
      response_variance = unname(variances), variance_se = unname(standard_error),
      variance_method = source, geometry$summary, integration_method = if (!length(take)) {
        "analytical"
      } else if (enumerate) "enumerated" else "sampled",
      evaluated_sets = sum(lengths(lapply(take, "[[", "shared"))),
      numerical_pools = if (enumerate) 0L else length(take), row.names = NULL)
  }
  out <- do.call(rbind, rows)
  baseline <- out[out$order == 2L, ]
  b <- match(out$response, baseline$response)
  denom <- baseline$response_variance[b]
  denom_se <- baseline$variance_se[b]
  out$relative_information <- out$relative_information_se <- NA_real_
  positive <- denom > 0
  out$relative_information[positive] <- out$response_variance[positive] / denom[positive]
  out$relative_information_se[positive] <- sqrt((out$variance_se[positive] / denom[positive])^2 +
    (out$relative_information[positive] * denom_se[positive] / denom[positive])^2)
  out$relative_information_se[positive & out$order == 2L] <- 0
  out$pairwise_variance <- denom
  out$pairwise_variance_method <- baseline$variance_method[b]
  out$retention_breakpoint <- breakpoint
  structural_zero <- out$response_variance == 0 & out$variance_method != "sampled"
  exact_rank_deficient <- out$integration_method == "enumerated" &
    !is.na(out$design_rank) & out$design_rank < out$design_columns
  out$candidate <- if (is.null(env)) NA else
    !structural_zero & !exact_rank_deficient & out$unsupported_predictors == 0L
  out$recommended_status <- if (is.null(env)) "ecology_only" else ifelse(
    out$candidate, ifelse(out$response_variance > 0, "candidate", "verify_variation"), "unsupported")
  out <- out[order(match(out$order, orders), match(out$response, selected)), ]
  out <- out[out$order %in% orders, ]
  rownames(out) <- NULL
  structure(out, class = c("zeta_order_information", "data.frame"),
    environment_detail = if (is.null(env)) NULL else do.call(rbind, detail[grid %in% orders]),
    candidate_orders = if (is.null(env)) NULL else split(out$order[out$candidate], out$response[out$candidate]))
}

.order_variance <- function(x, exact = FALSE) {
  if (diff(range(x)) <= 8 * .Machine$double.eps * max(abs(x))) return(0)
  stats::var(x) * if (exact) (length(x) - 1) / length(x) else 1
}

.order_response_exact <- function(q, richness, occupancy, constant_q, constant_positive,
                                  threshold, k, support) {
  relative_mean <- max(0, 1 - support$abundance_zeta)
  used <- occupancy >= k
  shared_mean <- sum(colMeans(q)[used] *
    exp(lchoose(occupancy[used] - 1L, k - 1L) - lchoose(nrow(q) - 1L, k - 1L)))
  occupancy_mean <- max(0, 1 - shared_mean)
  averages <- c(incidence = NA_real_, relative = relative_mean, occupancy = occupancy_mean,
                heterogeneity = max(0, shared_mean - support$abundance_zeta))
  variances <- stats::setNames(rep(NA_real_, 4L), names(averages))
  if (all(richness == richness[1L])) {
    averages["incidence"] <- 1 - support$expected_richness / richness[1L]
    richness_variance <- support$richness_variance
    # Near a cancellation in E(Z^2) - E(Z)^2, integrate the centred response
    # unless we can establish constancy from species occupancies alone.
    constant_shared <- all(occupancy < k | occupancy == nrow(q))
    if (constant_shared || richness_variance > 64 * .Machine$double.eps * support$richness_second_moment) {
      variances["incidence"] <- richness_variance / richness[1L]^2
    }
  }
  if (support$sharing_upper == 0) {
    averages["incidence"] <- 1
    variances[] <- 0
  }
  if (!is.null(threshold) && is.na(variances["relative"])) {
    p <- exp(lchoose(0:nrow(q), k) - lchoose(nrow(q), k))
    a <- p[threshold$size + 1L]
    joint <- matrix(p[threshold$intersection + 1L], nrow = length(a))
    terms <- (joint - tcrossprod(a)) * threshold$increment_product
    magnitude <- sum(abs(terms))
    value <- sum(terms)
    # Near complete cancellation, use a centred variance of actual responses.
    if (magnitude == 0 || value > 64 * .Machine$double.eps * magnitude) {
      variances["relative"] <- value
    }
  }
  if (constant_positive) {
    variances["heterogeneity"] <- 0
    variances["occupancy"] <- variances["relative"]
  }
  if (all(occupancy %in% c(0, nrow(q)))) {
    variances["occupancy"] <- 0
    variances["heterogeneity"] <- variances["relative"]
  }
  if (constant_q) variances[] <- 0
  if (k == nrow(q)) {
    averages["incidence"] <- 1 - support$expected_richness / mean(richness)
    variances[] <- 0
  }
  list(mean = averages, variance = variances)
}

.order_minimum_thresholds <- function(q, limit) {
  if (limit == 0L) return(NULL)
  levels <- lapply(seq_len(ncol(q)), function(j) {
    values <- sort(q[q[, j] > 0, j])
    if (length(values) < 2L) numeric() else sort(unique(values[-length(values)]))
  })
  sizes <- lengths(levels)
  if (!sum(sizes) || sum(sizes) > limit) return(NULL)
  above <- do.call(cbind, lapply(which(sizes > 0L), function(j) {
    outer(q[, j], levels[[j]], ">=")
  }))
  increments <- unlist(lapply(levels, function(v) diff(c(0, v))), use.names = FALSE)
  list(size = colSums(above), intersection = crossprod(above),
       increment_product = tcrossprod(increments))
}

.order_information_pool <- function(idx, q, presence, richness, env, types) {
  shared <- shared_species_cpp(presence, idx)
  out <- matrix(NA_real_, nrow(idx), length(types), dimnames = list(NULL, types))
  if ("incidence" %in% types) {
    out[, "incidence"] <- 1 - shared / rowMeans(matrix(richness[idx], nrow(idx)))
  }
  if (any(c("occupancy", "heterogeneity") %in% types)) {
    parts <- relative_partition_cpp(q, idx)
    if ("occupancy" %in% types) out[, "occupancy"] <- parts[, 1L]
    if ("heterogeneity" %in% types) out[, "heterogeneity"] <- parts[, 2L]
    if ("relative" %in% types) out[, "relative"] <- rowSums(parts)
  } else if ("relative" %in% types) {
    out[, "relative"] <- pmax(0, 1 - min_abundance_cpp(q, idx))
  }
  X <- NULL
  if (!is.null(env$basis)) {
    basis <- ugdm_design(env$eligible, idx, env$distance, ref = env$basis)
    positions <- vapply(seq_len(ncol(env$position_sites)), function(j) {
      rowMeans(matrix(env$position_sites[idx, j], nrow(idx)))
    }, numeric(nrow(idx)))
    X <- cbind(basis$Xdir, matrix(positions, nrow(idx)), basis$Xdist)
  }
  list(response = out, shared = shared > 0, design = X)
}

.order_sharing <- function(take, support, exact, enumerate) {
  if (exact) {
    p <- support$sharing_lower
    return(c(estimate = p, se = 0, lo = p, hi = p))
  }
  shared <- unlist(lapply(take, "[[", "shared"), use.names = FALSE)
  p <- mean(shared)
  if (enumerate) return(c(estimate = p, se = 0, lo = p, hi = p))
  N <- length(shared)
  z <- stats::qnorm(.975)
  center <- (p + z^2 / (2 * N)) / (1 + z^2 / N)
  half <- z * sqrt(p * (1 - p) / N + z^2 / (4 * N^2)) / (1 + z^2 / N)
  c(estimate = p, se = sqrt(p * (1 - p) / N), lo = max(0, center - half), hi = min(1, center + half))
}

.order_edge_moments <- function(edge_mean, edge_ss, row_sums, n) {
  c(square = max(0, edge_ss), overlap = sum((row_sums - (n - 1) * edge_mean)^2))
}

.order_abs_moments <- function(x) {
  x <- sort(x - mean(x))
  n <- length(x)
  i <- seq_len(n)
  prefix <- cumsum(x)
  row_sums <- i * x - prefix + prefix[n] - prefix - (n - i) * x
  edge_mean <- mean(row_sums) / (n - 1)
  .order_edge_moments(edge_mean, n * sum(x^2) - choose(n, 2) * edge_mean^2,
                      row_sums, n)
}

.order_spread_variance <- function(moments, k, n) {
  if (k == n) return(rep(0, ncol(moments)))
  p2 <- k * (k - 1) / (n * (n - 1))
  if (n == 3L) return(p2 * moments["square", ])
  p3 <- p2 * (k - 2) / (n - 2)
  only_edge <- p2 * (n - k) * (n - k - 1) / ((n - 2) * (n - 3))
  shared_vertex <- p3 * (n - k) / (n - 3)
  (only_edge * moments["square", ] + shared_vertex * moments["overlap", ]) / choose(k, 2)^2
}

.order_environment <- function(environment, coordinates, n, lonlat) {
  if (is.null(environment)) return(NULL)
  x <- as.matrix(as.data.frame(environment))
  if (!ncol(x) || nrow(x) != n || !is.numeric(x) || any(!is.finite(x)) ||
      anyDuplicated(colnames(x))) stop("Provide a finite numeric environment with unique predictor names and matching sites.")
  eligible <- vapply(seq_len(ncol(x)), function(j) length(unique(x[, j])) >= 4L, logical(1L))
  distance <- if (is.null(coordinates)) NULL else {
    xy <- as.matrix(coordinates)
    if (!is.numeric(xy) || !identical(dim(xy), c(n, 2L)) || any(!is.finite(xy))) {
      stop("Coordinates must be a finite numeric matrix with two columns and matching sites.")
    }
    coords_dist(xy, lonlat)
  }
  raw_variance <- apply(x, 2L, stats::var)
  scales <- sqrt(raw_variance)
  z <- scale(x, scale = ifelse(scales > 0, scales, 1))
  out <- list(eligible = x[, eligible, drop = FALSE], distance = distance,
    raw_variance = raw_variance, linear_variance = ifelse(scales > 0, 1, 0),
    raw_spread = apply(x, 2L, .order_abs_moments),
    spread = apply(z, 2L, .order_abs_moments), unsupported = sum(!eligible))
  if (!is.null(distance)) {
    edges <- distance[lower.tri(distance)]
    m <- mean(edges)
    out$geographic <- matrix(.order_edge_moments(m, sum((edges - m)^2), rowSums(distance), n),
                             ncol = 1L, dimnames = list(c("square", "overlap"), NULL))
  }
  if (any(eligible)) {
    out$basis <- ugdm_design(out$eligible, matrix(c(1L, 2L), 1L), distance,
                            dist_range = if (is.null(distance)) NULL else .pair_distance_range(distance))
    e <- .environment_for_moments(out$eligible, out$basis)
    out$position_sites <- do.call(cbind, lapply(seq_len(ncol(e$z)), function(j) {
      bounds <- attr(out$basis$environment_basis$splines[[j]], "Boundary.knots")
      t <- (e$z[, j] - mean(bounds)) / (diff(bounds) / 2)
      cbind(e$z[, j], t^2, t^3)
    }))
    out$spline_spread <- apply(e$spline, 2L, .order_abs_moments)
  }
  out
}

.order_geometry <- function(env, k, n, take) {
  columns <- c("position_variance", "nonlinear_position_variance", "spread_variance",
               "spline_spread_variance", "geographic_spread_variance", "design_rank",
               "design_columns", "design_condition", "unsupported_predictors")
  blank <- as.list(stats::setNames(rep(NA_real_, length(columns)), columns))
  if (is.null(env)) return(list(summary = blank, detail = NULL))
  factor <- (1 - k / n) / k
  position <- env$linear_variance * factor
  spread <- .order_spread_variance(env$spread, k, n)
  geographic <- if (is.null(env$geographic)) NA_real_ else .order_spread_variance(env$geographic, k, n)
  detail <- data.frame(order = k, predictor = names(env$raw_variance),
    position_variance = position, raw_position_variance = env$raw_variance * factor,
    spread_variance = spread, raw_spread_variance = .order_spread_variance(env$raw_spread, k, n),
    quadratic_variance = NA_real_, cubic_variance = NA_real_,
    spread_basis1_variance = NA_real_, spread_basis2_variance = NA_real_, spread_basis3_variance = NA_real_)
  included <- match(colnames(env$eligible), detail$predictor)
  nonlinear <- spline <- numeric()
  blank$design_columns <- 6L * nrow(detail) + 3L * !is.null(env$distance)
  blank$design_rank <- 0L
  blank$design_condition <- Inf
  if (length(included)) {
    variance <- t(matrix(apply(env$position_sites, 2L, stats::var) * factor, nrow = 3L))
    nonlinear <- variance[, 2:3, drop = FALSE]
    spline <- t(matrix(.order_spread_variance(env$spline_spread, k, n), nrow = 3L))
    detail[included, c("quadratic_variance", "cubic_variance")] <- nonlinear
    detail[included, c("spread_basis1_variance", "spread_basis2_variance", "spread_basis3_variance")] <- spline
    X <- do.call(rbind, lapply(take, "[[", "design"))
    if (length(X) && nrow(X) > 1L) {
      centered <- scale(X, scale = FALSE)
      scale <- sqrt(colSums(centered^2))
      active <- scale > sqrt(nrow(X)) * 16 * .Machine$double.eps * apply(abs(X), 2L, max)
      if (any(active)) {
        d <- svd(sweep(centered[, active, drop = FALSE], 2L, scale[active], "/"), nu = 0, nv = 0)$d
        blank$design_rank <- sum(d > max(d) * max(dim(X)) * .Machine$double.eps)
        if (blank$design_rank == blank$design_columns) blank$design_condition <- max(d) / min(d)
      }
    }
  }
  blank$position_variance <- mean(position)
  blank$nonlinear_position_variance <- if (length(nonlinear)) mean(nonlinear) else NA_real_
  blank$spread_variance <- mean(spread)
  blank$spline_spread_variance <- if (length(spline)) mean(spline) else NA_real_
  blank$geographic_spread_variance <- geographic
  blank$unsupported_predictors <- env$unsupported
  list(summary = blank, detail = detail)
}

#' Plot the zeta order information profile
#' @param x A profile from [zeta_order_information()].
#' @param ... Reserved for future use.
#' @return The input, invisibly.
#' @export
plot.zeta_order_information <- function(x, ...) {
  has_env <- any(is.finite(x$position_variance))
  old <- graphics::par(mfrow = if (has_env) c(3L, 3L) else c(2L, 3L),
                       mar = c(3.8, 3.8, 2.1, .8), mgp = c(2.4, .7, 0))
  on.exit(graphics::par(old))
  profile <- x[order(x$order, x$response), ]
  baseline <- profile[!duplicated(profile$order), ]
  palette <- c("#226B8A", "#9A4D24", "#58734D", "#75568B", "#505050")
  panel <- function(order, columns, title, ylab) {
    y <- as.matrix(columns)
    limits <- range(y[is.finite(y)], 0)
    if (diff(limits) == 0) limits <- limits + c(0, 1)
    colors <- palette[seq_len(ncol(y))]
    graphics::matplot(order, y, type = "o", lty = 1, pch = seq_len(ncol(y)),
                      cex = .55, col = colors, ylim = limits, xlab = "Zeta order",
                      ylab = ylab, main = title, cex.main = .95)
    if (ncol(y) > 1L) graphics::legend("topright", colnames(y), bty = "n",
      cex = .63, col = colors, lty = 1, pch = seq_len(ncol(y)))
  }
  orders <- unique(profile$order)
  response <- unique(profile$response)
  by_response <- function(metric) {
    cols <- lapply(response, function(type) profile[[metric]][profile$response == type])
    do.call(cbind, stats::setNames(cols, response))
  }
  panel(baseline$order, cbind(Incidence = baseline$incidence_zeta,
    Abundance = baseline$abundance_zeta), "Mean shared composition", "Zeta")
  panel(baseline$order, cbind(Incidence = baseline$incidence_retention,
    Abundance = baseline$abundance_retention, Sharing = baseline$sharing_probability),
    "Retention and sharing", "Fraction")
  panel(baseline$order, cbind(Incidence = baseline$effective_incidence_species,
    Abundance = baseline$effective_abundance_species), "Contributing species", "Effective number")
  panel(orders, by_response("response_variance"), "Compositional variation", "Response variance")
  panel(orders, by_response("relative_information"), "Relative information", "Variance / pairwise variance")
  if (has_env) for (j in seq_along(response)) {
    candidate <- profile[profile$response == response[j] & profile$candidate, "order"]
    if (length(candidate)) graphics::rug(candidate, col = palette[j], ticksize = .02 * j)
  }
  panel(orders, by_response("variance_se"), "Numerical uncertainty", "Variance estimate SE")
  if (has_env) {
    panel(baseline$order, cbind(Position = baseline$position_variance,
      Spread = baseline$spread_variance, `Nonlinear position` = baseline$nonlinear_position_variance,
      `Spline spread` = baseline$spline_spread_variance), "Environmental information", "Variance")
    panel(baseline$order, cbind(Observed = baseline$design_rank, Full = baseline$design_columns),
      "Joint support", "Centred design rank")
    condition <- log10(baseline$design_condition)
    deficient <- is.infinite(condition)
    condition[deficient] <- max(1, condition[!deficient & is.finite(condition)]) + 1
    panel(baseline$order, condition, "Conditioning", "Log10 scaled condition")
    if (any(deficient)) graphics::points(baseline$order[deficient], condition[deficient], pch = 4)
  }
  invisible(x)
}
