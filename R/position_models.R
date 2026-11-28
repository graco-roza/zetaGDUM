#' Environmental position for additive MS-GDUM
#'
#' Extends the original MS-GDUM design with environmental position, either
#' the mean standardised predictor (linear position) or the mean of a
#' site-level cubic spline (a single polynomial segment with linear
#' continuation outside the training bounds; nonlinear position). Spread
#' columns are the position-independent I-spline spread columns of
#' [ugdm_design()]. Coefficients of the linear and nonlinear position terms
#' are unrestricted. Zero nonlinear position coefficients recover the
#' original model. The I-spline spread terms of a single predictor can
#' already produce different effects for the same raw spread at different
#' positions, so position terms should be interpreted as an addition to,
#' not a replacement of, the spread curves.
#'
#' @param design Original [ugdm_design()] or a design from [ugdm_path_data()].
#' @param env Numeric site predictors in the same row order as `design`.
#' @param position Linear or nonlinear local environmental position.
#' @return Design for [ugdm_fit_position()].
#' @export
ugdm_position_design <- function(design, env, position = c("linear", "nonlinear")) {
  position <- match.arg(position)
  if (!is.null(design$spread_width)) stop("Supply the original environmental I-spline design.")
  details <- .environment_for_moments(env, design)
  z <- details$z
  k <- design$k
  nr <- nrow(design$idx)
  p <- ncol(z)
  mean_z <- design$Xpos / k
  bounds <- vapply(design$environment_basis$splines, function(b) {
    attr(b, "Boundary.knots")
  }, numeric(2L))
  half <- (bounds[2L, ] - bounds[1L, ]) / 2
  mid <- colMeans(bounds)
  raw_u <- sweep(sweep(mean_z, 2L, bounds[1L, ]), 2L, 2 * half, "/")
  means <- function(x) rowMeans(matrix(x[design$idx], nr, k))
  square_mean <- vapply(seq_len(p), function(j) means(z[, j]^2), numeric(nr))
  variance <- matrix(square_mean, nr, p) - mean_z^2
  pos <- mean_z
  if (position == "nonlinear") {
    blocks <- lapply(seq_len(p), function(j) {
      t <- (z[, j] - mid[j]) / half[j]
      clipped <- pmax(-1, pmin(1, t))
      outside <- t - clipped
      block <- cbind(mean_z[, j],
        means(clipped^2 + 2 * clipped * outside),
        means(clipped^3 + 3 * clipped^2 * outside))
      colnames(block) <- paste0(design$vars[j], c("", "_quadratic", "_cubic"))
      block
    })
    pos <- do.call(cbind, blocks)
  }
  structure(list(dir = design$Xdir, pos = pos, base = design$Xdir, geo = design$Xdist,
    curvature = ugdm_curvature_penalty(design), position_curvature = 1 / half^3,
    position = position, outside = raw_u < 0 | raw_u > 1,
    mean_position = mean_z, position_variance = pmax(variance, 0), order = k,
    basis = design$environment_basis, dist_basis = design$dist_basis, dist_range = design$dist_range,
    idx = design$idx, n = design$n, vars = design$vars), class = "ugdm_position_design")
}

.position_layout <- function(design, distance) {
  has_geo <- distance && !is.null(design$geo)
  direction <- if (has_geo) cbind(design$dir, design$geo) else design$dir
  nd <- ncol(direction)
  list(X = cbind(intercept = 1, direction, design$pos), nd = nd,
       dir = seq_len(nd) + 1L, pos = seq_len(ncol(design$pos)) + 1L + nd,
       distance = has_geo)
}

# Null-space slope of the integrated curvature basis (used only by the optional
# linear-component penalties spread_null and position_null).
.spread_penalty_basis <- function(basis) {
  bounds <- attr(basis, "Boundary.knots")
  cuts <- sort(unique(c(bounds, attr(basis, "knots"))))
  width <- diff(bounds)
  half <- diff(cuts) / 2
  middle <- utils::head(cuts, -1L) + half
  nodes <- rep(middle, 3L) + rep(c(-sqrt(3 / 5), 0, sqrt(3 / 5)), each = length(half)) * half
  weights <- rep(c(5 / 9, 8 / 9, 5 / 9), each = length(half)) * half
  deriv <- stats::predict(basis, nodes, derivs = 1L)
  values <- stats::predict(basis, nodes)
  slope <- 12 * drop(crossprod(values, weights * ((nodes - bounds[1L]) / width - .5))) / width
  list(derivative = crossprod(deriv * sqrt(weights)), null_slope = slope, width = width)
}

.position_penalties <- function(d, layout, strengths, penalty_scale = "native") {
  strength <- function(name) if (name %in% names(strengths)) strengths[[name]] else 0
  ridge <- strength("ridge")
  P <- diag(c(0, rep(ridge, ncol(layout$X) - 1L)))
  curvature <- strength("spread_smooth")
  position <- strength("position_smooth")
  spread_null <- strength("spread_null")
  position_null <- strength("position_null")
  for (j in seq_along(d$vars)) {
    at <- 1L + 3L * (j - 1L) + seq_len(3L)
    basis <- d$basis$splines[[j]]
    range <- diff(attr(basis, "Boundary.knots"))
    rescale <- if (penalty_scale == "range") range else 1
    M <- if (spread_null > 0) .spread_penalty_basis(basis)
    C <- curvature * d$curvature[[j]] * rescale^3
    if (spread_null > 0) C <- C + spread_null * tcrossprod(M$null_slope * rescale / range)
    P[at, at] <- P[at, at] + C
    at <- layout$pos
    if (d$position == "nonlinear") {
      at <- at[3L * (j - 1L) + seq_len(3L)]
      P[at, at] <- P[at, at] + position * rescale^3 *
        d$position_curvature[j] * diag(c(0, 8, 24))
      slope <- c(1, 0, 1.2 / range)
    } else {
      at <- at[j]
      slope <- 1
    }
    P[at, at] <- P[at, at] + position_null * tcrossprod(slope * rescale)
  }
  P
}

.position_fit_core <- function(X, y, w, P, nonnegative, start = NULL) {
  fit <- .variant_solution(X, y, w, P, nonnegative, start)
  fit$D2 <- if (!fit$constant && fit$null_deviance > 1e-15) {
    1 - fit$deviance / fit$null_deviance
  } else NA_real_
  fit
}

#' Fit an additive MS-GDUM with environmental position
#'
#' The directional (spread) and geographic coefficients are constrained
#' nonnegative; position coefficients are unrestricted. Three independent
#' penalties act on directional curvature, position curvature and slope
#' magnitude (ridge). With `spread_only = TRUE` the same design is fitted
#' with the position columns removed (the spread-only MS-GDM comparator), so
#' that bases, links, constraints and penalties are identical to the position
#' model.
#'
#' @param design Design from [ugdm_position_design()].
#' @param y,trials,sampling_weights,distance As in [ugdm_fit_variant()].
#' @param spread_smooth,position_smooth,ridge Nonnegative penalty strengths
#'   in the average working loss.
#' @param penalty_scale Use native standardised predictor units or express
#'   curvature penalties on a common unit range.
#' @param spread_null,position_null Optional penalties on the linear component
#'   obtained by function-space projection, separate from curvature.
#' @param start Optional previous fit or numeric vector in the same basis and
#'   coefficient layout for warm-starting the optimiser.
#' @param spread_only Fit without environmental position terms.
#' @return A fit for [ugdm_predict()].
#' @export
ugdm_fit_position <- function(design, y, trials = rep(1, length(y)),
                              sampling_weights = NULL, distance = TRUE,
                              spread_smooth = 0, position_smooth = 0, ridge = 0,
                              penalty_scale = c("native", "range"),
                              spread_null = 0, position_null = 0, start = NULL, spread_only = FALSE) {
  penalty_scale <- match.arg(penalty_scale)
  if (!inherits(design, "ugdm_position_design")) stop("Supply a position design.")
  if (length(y) != nrow(design$idx) || any(!is.finite(y) | y < 0 | y > 1)) {
    stop("Supply one finite fractional response per design row.")
  }
  supplied <- list(spread_smooth = spread_smooth, position_smooth = position_smooth,
                   ridge = ridge, spread_null = spread_null, position_null = position_null)
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  if (!all(vapply(supplied, valid, TRUE))) {
    stop("Supply finite nonnegative penalty strengths.")
  }
  strengths <- vapply(supplied, function(x) unname(x), numeric(1L))
  w <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  a <- .position_layout(design, distance)
  P <- .position_penalties(design, a, strengths, penalty_scale)
  if (spread_only) {
    if (!is.null(start)) stop("Warm starts are not available for the spread-only model.")
    keep <- c(1L, a$dir)
    X <- a$X[, keep, drop = FALSE]
    fit <- .position_fit_core(X, y, w, P[keep, keep, drop = FALSE], a$dir)
    fit$coef <- stats::setNames(fit$coef, colnames(X))
    return(structure(c(fit, list(strengths = strengths, penalty = P[keep, keep, drop = FALSE], vars = design$vars,
      penalty_scale = penalty_scale, nd = a$nd, has_distance = a$distance, basis = design$basis, order = design$order,
      dist_basis = if (a$distance) attributes(design$dist_basis)[c("degree", "knots", "Boundary.knots", "intercept")],
      dist_range = if (a$distance) design$dist_range, position = design$position, spread_only = TRUE)),
      class = "ugdm_position_fit"))
  }
  if (inherits(start, "ugdm_position_fit")) {
    if (!identical(start$basis, design$basis) || !identical(names(start$coef), colnames(a$X)) ||
        (a$distance && (!identical(start$dist_range, design$dist_range) ||
         !identical(start$dist_basis, attributes(design$dist_basis)[names(start$dist_basis)])))) {
      stop("Warm start must use the same basis and coefficient layout.")
    }
    start <- start$coef
  }
  if (!is.null(start) && (!is.numeric(start) || length(start) != ncol(a$X) ||
      any(!is.finite(start)))) stop("Supply one finite initial value per coefficient.")
  fit <- .position_fit_core(a$X, y, w, P, a$dir, start)
  fit$coef <- stats::setNames(fit$coef, colnames(a$X))
  structure(c(fit, list(strengths = strengths, penalty = P, vars = design$vars,
    penalty_scale = penalty_scale,
    nd = a$nd, has_distance = a$distance, basis = design$basis, order = design$order,
    dist_basis = if (a$distance) attributes(design$dist_basis)[c("degree", "knots", "Boundary.knots", "intercept")],
    dist_range = if (a$distance) design$dist_range,
    position = design$position)), class = "ugdm_position_fit")
}

.position_predict <- function(fit, design, type) {
  if (!inherits(design, "ugdm_position_design") ||
      !identical(design$position, fit$position) ||
      !identical(design$basis, fit$basis) || (fit$has_distance &&
      (!identical(design$dist_range, fit$dist_range) ||
       !identical(attributes(design$dist_basis)[names(fit$dist_basis)], fit$dist_basis)))) {
    stop("Supply a matching position design built with the training reference.")
  }
  X <- .position_layout(design, fit$has_distance)$X
  if (isTRUE(fit$spread_only)) X <- X[, c(1L, .position_layout(design, fit$has_distance)$dir), drop = FALSE]
  if (!identical(colnames(X), names(fit$coef))) stop("Prediction columns differ from the fitted model.")
  eta <- drop(X %*% fit$coef)
  if (type == "link") eta else stats::plogis(eta)
}

#' Position-model contrasts on the link and dissimilarity scales
#'
#' @param fit A fit from [ugdm_fit_position()].
#' @param design Evaluation design using the fitting reference basis.
#' @param weights Contrast matrix with one column per environmental set.
#'   For a two-by-two comparison, specify e.g. `c(1,-1,-1,1)`. The
#'   configurations and their common support should be declared beforehand.
#' @return Link and response contrasts, a row per supplied contrast.
#' @export
ugdm_position_contrast <- function(fit, design, weights) {
  W <- if (is.null(dim(weights))) matrix(weights, 1L) else as.matrix(weights)
  if (!is.numeric(W) || ncol(W) != nrow(design$idx) || any(!is.finite(W))) {
    stop("Supply finite contrast weights with one column per design row.")
  }
  eta <- ugdm_predict(fit, design, "link")
  data.frame(link = drop(W %*% eta), response = drop(W %*% stats::plogis(eta)))
}

#' Identifiability and attribution for additive MS-GDUM with position
#'
#' Scaled information is the unpenalised fractional-logit working
#' information, so a stabilising penalty cannot hide a weak direction.
#' Nonlinear position columns are assessed after projection off the spread
#' and linear position columns. The spread and position allocations of the
#' deviance explained are obtained from the matched restricted fits (spread
#' only, position only) under the same penalty.
#'
#' @param fit,design Matching fit and design.
#' @param y,trials,sampling_weights As in [ugdm_fit_position()].
#' @param contrasts Optional finite matrix of link-scale coefficient contrasts
#'   for reporting projection into weak working-information directions.
#' @return Rank, scaled condition, active constraints, predictor support and
#'   spread and position deviances explained.
#' @export
ugdm_position_diagnostics <- function(fit, design, y, trials = rep(1, length(y)),
                                      sampling_weights = NULL, contrasts = NULL) {
  w <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  a <- .position_layout(design, fit$has_distance)
  mu <- ugdm_predict(fit, design)
  v <- sqrt(w * mu * (1 - mu) / sum(w))
  Xw <- a$X * v
  scale <- sqrt(colSums(Xw^2))
  scale[scale < 1e-14] <- 1
  spec <- eigen(crossprod(sweep(Xw, 2L, scale, "/")), symmetric = TRUE)
  ev <- spec$values
  rank <- sum(ev > max(ev) * ncol(a$X) * .Machine$double.eps)
  contrast <- NULL
  if (!is.null(contrasts)) {
    L <- as.matrix(contrasts)
    if (ncol(L) != ncol(a$X) || any(!is.finite(L))) stop("Supply finite contrasts matching the coefficients.")
    projection <- sweep(L, 2L, scale, "/") %*% spec$vectors
    total <- rowSums(projection^2)
    weak <- ev <= max(ev) * 1e-8
    contrast <- data.frame(weak_direction_fraction = ifelse(total > 0,
      rowSums(projection[, weak, drop = FALSE]^2) / total, 0),
      inverse_curvature = rowSums(sweep(projection[, !weak, drop = FALSE]^2, 2L, ev[!weak], "/")))
  }
  P <- .position_penalties(design, a, fit$strengths,
    if (is.null(fit$penalty_scale)) "native" else fit$penalty_scale)
  restricted <- function(idx) {
    f <- .position_fit_core(a$X[, idx, drop = FALSE], y, w, P[idx, idx, drop = FALSE],
                            which(idx %in% a$dir))
    if (f$converged) f$D2 else NA_real_
  }
  both <- fit$D2
  spread_only <- restricted(c(1L, a$dir))
  pos_only <- restricted(c(1L, a$pos))
  variation <- function(X, Z) {
    Z <- Z * v
    norm <- colSums(Z^2)
    R <- qr.resid(qr(X * v), Z)
    ifelse(norm > 1e-24, sqrt(colSums(R^2) / norm), NA_real_)
  }
  nonlinear <- if (design$position == "nonlinear") a$pos[rep(c(FALSE, TRUE, TRUE), length(design$vars))] else integer()
  pos_resid <- if (length(nonlinear)) {
    stats::setNames(variation(a$X[, -nonlinear, drop = FALSE], a$X[, nonlinear, drop = FALSE]),
                    colnames(a$X)[nonlinear])
  } else NULL
  p <- design$mean_position
  sd <- function(X) apply(X, 2L, stats::sd)
  list(rank = rank, parameters = ncol(a$X),
    condition = if (rank == ncol(a$X)) max(ev) / min(ev) else Inf,
    active = names(fit$coef)[a$dir[fit$coef[a$dir] <= 1e-8]],
    projected_gradient = fit$projected_gradient, converged = fit$converged,
    contrast_information = contrast,
    nonlinear_position_residual = pos_resid,
    support = data.frame(variable = design$vars, mean_min = apply(p, 2L, min),
      mean_max = apply(p, 2L, max), mean_sd = sd(p),
      within_set_variance_sd = sd(design$position_variance),
      mean_variance_correlation = vapply(seq_along(design$vars), function(j) {
        spread <- design$position_variance[, j]
        if (length(spread) < 2L || min(stats::sd(spread), stats::sd(p[, j])) == 0) return(NA_real_)
        stats::cor(p[, j], spread)
      }, numeric(1L)),
      fraction_extrapolating_interpolation = colMeans(design$outside)),
    deviance_allocation = c(spread = (both + spread_only - pos_only) / 2,
      position = (both + pos_only - spread_only) / 2),
    full_D2 = fit$D2, spread_only_D2 = spread_only, position_only_D2 = pos_only)
}

#' Select additive MS-GDUM penalties by holding out groups of sites
#'
#' Each fold's environmental scales, knot bounds and distance scaling come
#' exclusively from the training sites. Training and held-out sets use
#' disjoint sites and identical set pools across candidate penalties. Folds
#' contribute equally to the validation loss. Use spatial blocks to assess
#' transfer across space. The fold-omission choices diagnose sensitivity
#' of the selected strengths to the validation blocks.
#'
#' @param env Numeric site predictors.
#' @param folds One block label per site, or a list or matrix of repeated
#'   partitions with one label per site per partition. Each partition needs
#'   at least two blocks.
#' @param order Order of the GDUM response.
#' @param response Function accepting an index matrix with original site
#'   numbers. Return a numeric response or a list with `y`, optionally
#'   `trials` and `keep`, as in [ugdm_response()].
#' @param penalties Data frame with candidate rows and columns
#'   `spread_smooth`, `position_smooth` and `ridge`. Optionally include
#'   `spread_null` and `position_null`.
#' @param position As in [ugdm_position_design()].
#' @param dist_matrix Optional site distance matrix.
#' @param M Combinations per training and per held-out fold.
#' @param seed Optional seed for combination sampling.
#' @param penalty_scale As in [ugdm_fit_position()].
#' @param warm_start Reuse the preceding converged fit along the regularisation path.
#' @return Candidate scores, fold scores, fold-omission choices and the
#'   index of the selected row. Ties follow the supplied candidate order.
#' @export
ugdm_select_penalties <- function(env, folds, order, response, penalties,
                                  position = "nonlinear",
                                  dist_matrix = NULL, M = 1000L, seed = NULL,
                                  penalty_scale = c("native", "range"),
                                  warm_start = TRUE) {
  penalty_scale <- match.arg(penalty_scale)
  fields <- c("spread_smooth", "position_smooth", "ridge")
  if (!is.data.frame(penalties) || !nrow(penalties) ||
      !all(fields %in% names(penalties)) ||
      !all(vapply(penalties[fields], is.numeric, TRUE)) ||
      any(!is.finite(as.matrix(penalties[fields])) | as.matrix(penalties[fields]) < 0)) {
    stop("Provide a candidate data frame with three finite nonnegative strengths.")
  }
  fields <- c(fields, intersect(c("spread_null", "position_null"), names(penalties)))
  if (any(!vapply(penalties[fields], is.numeric, TRUE)) ||
      any(!is.finite(as.matrix(penalties[fields])) | as.matrix(penalties[fields]) < 0)) {
    stop("Supply finite nonnegative penalties.")
  }
  if (is.matrix(folds)) folds <- lapply(seq_len(ncol(folds)), function(j) folds[, j])
  if (!is.list(folds)) folds <- list(folds)
  valid_partition <- function(f) length(f) == nrow(env) && !anyNA(f) && length(unique(f)) >= 2L
  if (!length(folds) || !all(vapply(folds, valid_partition, TRUE))) {
    stop("Supply a block label for every site and at least two blocks per partition.")
  }
  if (!is.function(response)) stop("Supply a set-level response function.")
  k <- .set_design_integer(order, "order", 2L)
  M <- .set_design_integer(M, "M", 1L)
  groups <- lapply(folds, function(cluster) split(seq_len(nrow(env)), cluster, drop = TRUE))
  blocks <- if (length(groups) == 1L) groups[[1L]] else
    do.call(c, stats::setNames(groups, paste0("partition", seq_along(groups))))
  partition <- rep(seq_along(groups), lengths(groups))
  if (any(lengths(blocks) < k) || any(nrow(env) - lengths(blocks) < max(k, 4L))) {
    stop("Every fold needs at least order sites, and each training fold needs at least four sites.")
  }
  if (!is.null(dist_matrix)) {
    if (!identical(dim(dist_matrix), rep(nrow(env), 2L)) || !is.numeric(dist_matrix) ||
        any(!is.finite(dist_matrix)) || !isSymmetric(dist_matrix)) {
      stop("Supply a finite symmetric site distance matrix.")
    }
  }
  if (!is.null(seed)) set.seed(seed)
  sample <- function(sites) matrix(sites[sample_sets(length(sites), k, M)], ncol = k)
  # The traversal changes computation only. Ties are still resolved in
  # caller-supplied order, and failed warm fits are retried from scratch.
  penalty_norm <- vapply(penalties[fields], max, numeric(1L))
  penalty_norm[penalty_norm == 0] <- 1
  path <- order(rowSums(sweep(as.matrix(penalties[fields]), 2L, penalty_norm, "/")),
                decreasing = TRUE)
  scores <- lapply(blocks, function(test) {
    train <- setdiff(seq_len(nrow(env)), test)
    distance <- if (is.null(dist_matrix)) NULL else .pair_distance_range(dist_matrix[train, train, drop = FALSE])
    define <- function(idx, ref = NULL) {
      values <- response(idx)
      if (!is.list(values)) values <- list(y = values)
      keep <- if (is.null(values$keep)) rep(TRUE, nrow(idx)) else values$keep
      w <- if (is.null(values$trials)) rep(1, nrow(idx)) else values$trials
      if (!is.logical(keep) || length(keep) != nrow(idx) || anyNA(keep) || !any(keep) ||
          length(values$y) != nrow(idx) || length(w) != nrow(idx) ||
          any(!is.finite(values$y[keep]) | values$y[keep] < 0 | values$y[keep] > 1)) {
        stop("The response function supplied no valid rows or invalid fractional responses.")
      }
      design <- ugdm_design(env, idx[keep, , drop = FALSE], dist_matrix, ref = ref,
                            basis_sites = train, dist_range = distance)
      list(design = ugdm_position_design(design, env, position),
           ref = design, y = values$y[keep], w = .model_weights(sum(keep), w[keep]))
    }
    a <- define(sample(train))
    b <- define(sample(test), a$ref)
    if (!any(b$w > 0)) stop("Each validation fold needs at least one positive working weight.")
    out <- numeric(nrow(penalties)) + Inf
    current <- NULL
    for (j in path) {
      args <- c(list(design = a$design, y = a$y, trials = a$w, penalty_scale = penalty_scale),
                as.list(penalties[j, fields, drop = FALSE]))
      f <- do.call(ugdm_fit_position, c(args, list(start = current)))
      if (!f$converged && !is.null(current)) f <- do.call(ugdm_fit_position, args)
      if (!f$converged) next
      if (warm_start) current <- f
      out[j] <- .binom_dev(b$y, ugdm_predict(f, b$design), b$w / sum(b$w))
    }
    out
  })
  fold_loss <- do.call(rbind, scores)
  rownames(fold_loss) <- names(blocks)
  average <- function(mat, rows = seq_along(blocks)) {
    split_rows <- split(rows, partition[rows])
    colMeans(do.call(rbind, lapply(split_rows, function(i) colMeans(mat[i, , drop = FALSE]))))
  }
  mean_loss <- average(fold_loss)
  best <- function(p) if (any(is.finite(p))) which.min(replace(p, !is.finite(p), Inf)) else NA_integer_
  choice <- best(mean_loss)
  alternative <- vapply(seq_along(blocks), function(i) {
    rows <- setdiff(seq_along(blocks), i)
    best(average(fold_loss, rows))
  }, integer(1L))
  list(penalties = penalties, mean_loss = mean_loss, fold_loss = fold_loss,
       partition = partition, selected = choice,
       selected_without_fold = stats::setNames(alternative, names(blocks)))
}
