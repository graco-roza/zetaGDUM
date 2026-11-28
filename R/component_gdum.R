#' Relative-abundance components for multisite GDUM
#'
#' Occupancy sums mean site proportions for taxa absent from any member of
#' the set. Heterogeneity sums mean minus minimum proportions among taxa
#' present throughout. Both use the original site total as denominator.
#' @param community Nonnegative sites-by-species matrix with positive row totals.
#' @param sets Integer matrix of distinct site indices, one set per row.
#' @return A response matrix with `total`, `occupancy` and `heterogeneity`
#'   columns and the original set indices.
#' @export
ugdm_component_response <- function(community, sets) {
  if (anyNA(community)) stop("Supply complete community observations.")
  q <- .conditional_community(community, relative = TRUE)
  idx <- .species_indices(sets, nrow(q))
  storage.mode(idx) <- "integer"
  parts <- relative_partition_cpp(q, idx)
  parts <- parts / pmax(1, rowSums(parts))
  structure(cbind(total = pmin(1, rowSums(parts)), occupancy = parts[, 1L],
    heterogeneity = parts[, 2L]), sets = idx, class = c("ugdm_component_response", "matrix", "array"))
}

#' Design for component-specific multisite GDUM
#'
#' Uses the existing nonlinear spread and position bases. The additive
#' I-splines can already produce different effects for the same raw spread
#' at different positions; position terms are an addition to, not a
#' replacement of, the spread curves.
#' @param environment Numeric environmental data, one row per site.
#' @param sets Site index matrix as in [ugdm_component_response()].
#' @param model Spread alone, or spread with additive nonlinear position.
#'   Defaults to the reference model if supplied.
#' @param coordinates Optional site coordinates for a geographic distance term.
#' @param reference Optional training design or component fit. Reuses its basis.
#' @param lonlat Use great-circle distance for longitude and latitude.
#' @return A component design containing training coordinates and set geometry.
#' @export
ugdm_component_design <- function(environment, sets,
    model = c("additive", "spread_only"), coordinates = NULL,
    reference = NULL, lonlat = FALSE) {
  if (!is.null(reference) && !inherits(reference, c("ugdm_component_fit", "ugdm_component_design"))) {
    stop("Supply a component training design or fit as reference.")
  }
  if (missing(model) && !is.null(reference)) model <- reference$model
  model <- match.arg(model)
  idx <- .species_indices(sets, nrow(environment))
  storage.mode(idx) <- "integer"
  env <- as.data.frame(environment)
  if (!ncol(env) || anyDuplicated(names(env))) stop("Provide uniquely named environmental predictors.")
  distance <- if (!is.null(coordinates)) {
    xy <- as.matrix(coordinates)
    if (!identical(dim(xy), c(nrow(env), 2L)) || !is.numeric(xy) || any(!is.finite(xy))) {
      stop("Supply two finite coordinate columns in site order.")
    }
    coords_dist(xy, lonlat)
  }
  if (!is.null(reference) && xor(is.null(distance), is.null(reference$reference$dist_basis))) {
    stop("Training and evaluation must use the same geographic terms.")
  }
  base <- ugdm_design(env, idx, dist_matrix = distance, ref = reference$reference,
    dist_range = if (!is.null(distance)) .pair_distance_range(distance))
  design <- ugdm_position_design(base, env, "nonlinear")
  design$model <- model
  design$reference <- base[c("environment_basis", "dist_basis", "dist_range")]
  values <- as.matrix(env)
  design$raw_position <- vapply(seq_len(ncol(values)), function(j) {
    rowMeans(matrix(values[idx, j], nrow(idx)))
  }, numeric(nrow(idx)))
  design$raw_spread <- do.call(cbind, lapply(seq_len(ncol(values)), function(j) {
    z <- values[, j] - mean(values[, j])
    mean_pair_isp_cpp(matrix(z, ncol = 1L), z, idx)
  }))
  dim(design$raw_position) <- c(nrow(idx), ncol(values))
  colnames(design$raw_position) <- colnames(design$raw_spread) <- names(env)
  class(design) <- c("ugdm_component_design", class(design))
  design
}

.component_problem <- function(design, ridge) {
  a <- .position_layout(design, !is.null(design$geo))
  groups <- lapply(seq_along(design$vars), function(j) {
    b <- 3L * (j - 1L) + seq_len(3L)
    direction <- list(spread = a$dir[b])
    if (design$model != "spread_only") direction$position <- a$pos[b]
    direction
  })
  names(groups) <- design$vars
  P <- .position_penalties(design, a, c(ridge = ridge))
  if (a$distance) {
    label <- make.unique(c(names(groups), "geographic_distance"))[length(groups) + 1L]
    groups[[label]] <- list(distance = utils::tail(a$dir, ncol(design$geo)))
  }
  list(X = a$X, P = P, groups = groups)
}

.component_transform <- function(problem, keep) {
  diag(ncol(problem$X))[, keep, drop = FALSE]
}

.component_fits <- function(X, P, y, w, architecture, bounded, start = NULL) {
  total <- y[, "total"]
  means <- stats::setNames(colSums(y * (w / sum(w))), colnames(y))
  solve <- function(y, weights, bound, name) {
    initial <- if (is.null(start)) NULL else start[[name]]$coef
    a <- .variant_solution(X, y, weights, P, bound, initial)
    if (!a$converged && !is.null(initial)) a <- .variant_solution(X, y, weights, P, bound)
    a$constant_value <- if (isTRUE(a$constant)) y[which(weights > 0)[1L]] else NA_real_
    a
  }
  out <- list(total = solve(total, w, bounded, "total"))
  if (architecture == "independent") {
    out$occupancy <- solve(y[, "occupancy"], w, integer(), "occupancy")
    out$heterogeneity <- solve(y[, "heterogeneity"], w, integer(), "heterogeneity")
  } else {
    proportion <- y[, "heterogeneity"] / ifelse(total > 0, total, 1)
    w2 <- w * total
    out$share <- if (any(w2 > 0)) solve(pmin(1, proportion), w2, integer(), "share") else
      list(coef = rep(0, ncol(X)), constant_value = .5, converged = TRUE, constant = TRUE,
           projected_gradient = 0)
  }
  attr(out, "share_identified") <- architecture == "coherent" && means["total"] > 0
  out
}

.component_values <- function(fits, X, architecture) {
  fitted <- lapply(fits, function(f) {
    if (is.finite(f$constant_value)) rep(f$constant_value, nrow(X)) else
      stats::plogis(drop(X %*% f$coef))
  })
  total <- fitted$total
  heterogeneity <- if (architecture == "coherent") total * fitted$share else fitted$heterogeneity
  occupancy <- if (architecture == "coherent") total - heterogeneity else fitted$occupancy
  data.frame(shared = 1 - total, total = total, occupancy = occupancy,
    heterogeneity = heterogeneity, closure_error = occupancy + heterogeneity - total)
}

#' Fit independent or coherent component multisite GDUM
#'
#' The coherent model fits total dissimilarity and its heterogeneity share.
#' The share fit uses total as a working weight, so rows with zero total
#' contribute only to the total fit. It estimates the ratio of expected
#' heterogeneity to expected total. Predictions preserve compositional closure.
#' The independent architecture fits all three marginals separately.
#' Spread monotonicity is optional for total and unrestricted for the
#' components or conditional share. The same ridge penalty applies to the
#' spread-only and additive models, so they are comparable.
#' @param design Component design from [ugdm_component_design()].
#' @param response Aligned component response from [ugdm_component_response()].
#' @param architecture Fit independent marginal regressions or a coherent model.
#' @param ridge Nonnegative main-effect ridge penalty.
#' @param monotone_total Constrain the total model's directional coefficients.
#' @param sampling_weights Optional uniform-target importance weights.
#' @return A fit for [ugdm_predict_components()] and [ugdm_component_importance()].
#' @export
ugdm_fit_components <- function(design, response, architecture = c("coherent", "independent"),
    ridge = .001, monotone_total = TRUE, sampling_weights = NULL) {
  architecture <- match.arg(architecture)
  .component_match_response(design, response)
  if (!is.logical(monotone_total) || length(monotone_total) != 1L || is.na(monotone_total) ||
      !is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) stop("Invalid penalties or monotonicity.")
  problem <- .component_problem(design, ridge)
  keep <- sort(c(1L, unlist(problem$groups, use.names = FALSE)))
  A <- .component_transform(problem, keep)
  w <- .model_weights(nrow(response), sampling_weights = sampling_weights)
  w <- w / max(w)
  X <- problem$X %*% A
  direction <- unlist(lapply(problem$groups, function(g) unname(g[setdiff(names(g), "position")])),
                     use.names = FALSE)
  bounds <- if (monotone_total) which(keep %in% direction) else integer()
  fits <- .component_fits(X, crossprod(A, problem$P %*% A), response, w, architecture, bounds)
  structure(list(model = design$model, architecture = architecture, fits = fits,
    keep = keep, transform = A, ridge = ridge,
    direction = if (monotone_total) direction else integer(), weights = w,
    reference = design$reference, order = design$order, vars = design$vars,
    sets = design$idx, share_trials = if (architecture == "coherent") response[, "total"],
    converged = all(vapply(fits, function(f) f$converged, logical(1L))),
    share_identified = attr(fits, "share_identified"), support = .component_support_reference(design, w)),
    class = "ugdm_component_fit")
}

.component_match_response <- function(design, response) {
  if (!inherits(design, "ugdm_component_design") ||
      !inherits(response, "ugdm_component_response") ||
      !identical(attr(response, "sets"), design$idx) ||
      !identical(colnames(response), c("total", "occupancy", "heterogeneity")) ||
      any(!is.finite(response) | response < 0 | response > 1) ||
      any(abs(response[, "total"] - response[, "occupancy"] - response[, "heterogeneity"]) > 1e-12)) {
    stop("Supply valid component responses in the design's set order.")
  }
}

.component_prediction_matrix <- function(fit, design) {
  if (!inherits(fit, "ugdm_component_fit") || !inherits(design, "ugdm_component_design") ||
      design$model != fit$model || design$order != fit$order || !identical(design$vars, fit$vars) ||
      !identical(design$reference, fit$reference)) {
    stop("Rebuild the evaluation design with the fitted training reference and order.")
  }
  .position_layout(design, !is.null(design$geo))$X %*% fit$transform
}

#' Predict the component dissimilarities
#' @param fit Fitted component model.
#' @param design Matching new design built with the training reference.
#' @param mask Replace predictions outside observed spread-position support
#'   with missing values.
#' @param ... Support settings for [ugdm_component_support()].
#' @return A data frame with predicted shared composition, total dissimilarity,
#'   occupancy, heterogeneity and the signed component closure error.
#' @export
ugdm_predict_components <- function(fit, design, mask = FALSE, ...) {
  out <- .component_values(fit$fits, .component_prediction_matrix(fit, design), fit$architecture)
  if (isTRUE(mask)) out[!ugdm_component_support(fit, design, ...)$supported, ] <- NA_real_
  out
}
