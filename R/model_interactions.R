#' Environmental position and position-dependent spread for MS-GDUM
#'
#' Extends the original MS-GDUM design with the mean of a site-level cubic
#' spline (a single polynomial segment with linear continuation outside
#' training bounds) and linear interpolation between monotone spread curves.
#' The interpolation coordinate is based on mean standardised position,
#' clipped to the original training range. Coefficients of the linear and
#' nonlinear position terms are unrestricted. Equal spread endpoints and
#' zero nonlinear position coefficients recover the original model.
#'
#' @param design Original [ugdm_design()] or a design from [ugdm_path_data()].
#' @param env Numeric site predictors in the same row order as `design`.
#' @param position Linear or nonlinear local environmental position.
#' @param spread Position-independent, set-position-dependent or contextual
#'   spread. The contextual option combines pair-midpoint and set-mean
#'   endpoint interpolation with nonnegative coefficients. It uses no
#'   product of the two interpolation coordinates.
#' @return Design for [ugdm_fit_interactions()].
#' @export
ugdm_interaction_design <- function(design, env, position = c("linear", "nonlinear"),
                                    spread = c("constant", "varying", "contextual")) {
  position <- match.arg(position)
  spread <- match.arg(spread)
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
  u <- pmin(pmax(raw_u, 0), 1)
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
  base <- design$Xdir
  dir <- base
  if (spread == "varying") {
    dir <- do.call(cbind, lapply(seq_len(p), function(j) {
      block <- base[, 3L * (j - 1L) + seq_len(3L), drop = FALSE]
      out <- cbind(block * (1 - u[, j]), block * u[, j])
      colnames(out) <- paste0(rep(colnames(block), 2L), rep(c("_L", "_U"), each = 3L))
      out
    }))
  }
  context <- NULL
  if (spread == "contextual") {
    context <- .context_geometry(z, details$spline, design, mean_z, bounds, u)
    dir <- context$dir
  }
  structure(list(dir = dir, pos = pos, base = base, geo = design$Xdist,
    contextual = context$context, pair_midpoint = context$pair, bounded_context = context$effective,
    curvature = ugdm_curvature_penalty(design), position_curvature = 1 / half^3,
    position = position, spread = spread, u = u, outside = raw_u < 0 | raw_u > 1,
    mean_position = mean_z, position_variance = pmax(variance, 0), order = k,
    basis = design$environment_basis, dist_basis = design$dist_basis, dist_range = design$dist_range,
    idx = design$idx, vars = design$vars), class = "ugdm_interaction_design")
}

.interaction_layout <- function(design, distance) {
  has_geo <- distance && !is.null(design$geo)
  direction <- if (has_geo) cbind(design$dir, design$geo) else design$dir
  nd <- ncol(direction)
  list(X = cbind(intercept = 1, direction, design$pos), nd = nd,
       dir = seq_len(nd) + 1L, pos = seq_len(ncol(design$pos)) + 1L + nd,
       distance = has_geo)
}

.interaction_penalties <- function(d, layout, strengths,
                                   interaction_penalty = "coefficients", penalty_scale = "native") {
  strength <- function(name) if (name %in% names(strengths)) strengths[[name]] else 0
  ridge <- strength("ridge")
  P <- diag(c(0, rep(ridge, ncol(layout$X) - 1L)))
  curvature <- strength("spread_smooth")
  position <- strength("position_smooth")
  interaction <- strength("interaction_shrink")
  spread_null <- strength("spread_null")
  position_null <- strength("position_null")
  I <- diag(3L)
  for (j in seq_along(d$vars)) {
    width <- switch(d$spread, constant = 3L, varying = 6L, contextual = 12L)
    at <- 1L + width * (j - 1L) + seq_len(width)
    basis <- d$basis$splines[[j]]
    range <- diff(attr(basis, "Boundary.knots"))
    rescale <- if (penalty_scale == "range") range else 1
    M <- if (interaction_penalty == "derivative" || spread_null > 0) .context_penalty_basis(basis)
    C <- d$curvature[[j]] * rescale^3
    C <- curvature * C
    if (spread_null > 0) C <- C + spread_null * tcrossprod(M$null_slope * rescale / range)
    D <- if (interaction_penalty == "derivative") M$derivative * rescale else I
    evaluate <- function(A, G) crossprod(A, G %*% A)
    if (d$spread == "constant") {
      P[at, at] <- P[at, at] + C
    } else if (d$spread == "varying") {
      delta <- cbind(-I, I)
      A <- matrix(0, 6L, 6L)
      A[1:3, 1:3] <- A[4:6, 4:6] <- C / 2
      P[at, at] <- P[at, at] + A + interaction * evaluate(delta, D)
    } else {
      zero <- matrix(0, 3L, 3L)
      low_pair <- cbind(I, zero, zero, zero)
      high_pair <- cbind(zero, I, zero, zero)
      low_set <- cbind(zero, zero, I, zero)
      high_set <- cbind(zero, zero, zero, I)
      corners <- list(low_pair + low_set, high_pair + low_set, low_pair + high_set, high_pair + high_set)
      A <- Reduce("+", lapply(corners, evaluate, G = C)) / 4
      delta_pair <- high_pair - low_pair
      delta_set <- high_set - low_set
      # Use physical coordinates so the penalty is invariant to moving
      # a common function between the pair and set endpoint representation.
      R <- evaluate(Reduce("+", corners) / 4, I) +
        evaluate(delta_pair + delta_set, I) + evaluate(delta_set, I)
      P[at, at] <- A + ridge * R + interaction * (evaluate(delta_pair, D) + evaluate(delta_set, D))
    }
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

.interaction_fit_core <- function(X, y, w, P, nonnegative, start = NULL) {
  fit <- .variant_solution(X, y, w, P, nonnegative, start)
  fit$D2 <- if (!fit$constant && fit$null_deviance > 1e-15) {
    1 - fit$deviance / fit$null_deviance
  } else NA_real_
  fit
}

#' Fit an extended multisite GDUM
#'
#' The directional endpoints and geographic coefficients are constrained
#' nonnegative. Four independent penalties act on directional curvature,
#' position curvature, endpoint coefficient differences and slope magnitude.
#' The directional curvature penalty averages the two endpoint penalties
#' under varying spread. With equal endpoints its curvature is exactly
#' that of the corresponding position-independent spread.
#'
#' @param design Design from [ugdm_interaction_design()].
#' @param y,trials,sampling_weights,distance As in [ugdm_fit_variant()].
#' @param spread_smooth,position_smooth,interaction_shrink,ridge Nonnegative
#'   penalty strengths in the average working loss.
#' @param interaction_penalty Penalise endpoint coefficient differences
#'   or the integrated squared derivative difference on the environmental range.
#' @param penalty_scale Use native standardised predictor units or express
#'   derivative and curvature penalties on a common unit range.
#' @param spread_null,position_null Optional penalties on the linear component
#'   obtained by function-space projection, separate from curvature.
#' @param start Optional previous fit or numeric vector in the same basis and
#'   coefficient layout for warm-starting the optimiser.
#' @return An extended GDUM fit for [ugdm_predict()].
#' @export
ugdm_fit_interactions <- function(design, y, trials = rep(1, length(y)),
                                  sampling_weights = NULL, distance = TRUE,
                                  spread_smooth = 0, position_smooth = 0,
                                  interaction_shrink = 0, ridge = 0,
                                  interaction_penalty = c("coefficients", "derivative"),
                                  penalty_scale = c("native", "range"),
                                  spread_null = 0, position_null = 0, start = NULL) {
  interaction_penalty <- match.arg(interaction_penalty)
  penalty_scale <- match.arg(penalty_scale)
  if (!inherits(design, "ugdm_interaction_design")) stop("Supply an interaction design.")
  if (length(y) != nrow(design$idx) || any(!is.finite(y) | y < 0 | y > 1)) {
    stop("Supply one finite fractional response per design row.")
  }
  supplied <- list(spread_smooth = spread_smooth, position_smooth = position_smooth,
                   interaction_shrink = interaction_shrink, ridge = ridge,
                   spread_null = spread_null, position_null = position_null)
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  if (!all(vapply(supplied, valid, TRUE))) {
    stop("Supply finite nonnegative penalty strengths.")
  }
  strengths <- vapply(supplied, function(x) unname(x), numeric(1L))
  w <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  a <- .interaction_layout(design, distance)
  P <- .interaction_penalties(design, a, strengths, interaction_penalty, penalty_scale)
  if (inherits(start, "ugdm_interaction_fit")) {
    if (!identical(start$basis, design$basis) || !identical(names(start$coef), colnames(a$X)) ||
        (a$distance && (!identical(start$dist_range, design$dist_range) ||
         !identical(start$dist_basis, attributes(design$dist_basis)[names(start$dist_basis)])))) {
      stop("Warm start must use the same basis and coefficient layout.")
    }
    start <- start$coef
  }
  if (!is.null(start) && (!is.numeric(start) || length(start) != ncol(a$X) ||
      any(!is.finite(start)))) stop("Supply one finite initial value per coefficient.")
  fit <- .interaction_fit_core(a$X, y, w, P, a$dir, start)
  fit$coef <- stats::setNames(fit$coef, colnames(a$X))
  structure(c(fit, list(strengths = strengths, penalty = P, vars = design$vars,
    interaction_penalty = interaction_penalty, penalty_scale = penalty_scale,
    nd = a$nd, has_distance = a$distance, basis = design$basis, order = design$order,
    dist_basis = if (a$distance) attributes(design$dist_basis)[c("degree", "knots", "Boundary.knots", "intercept")],
    dist_range = if (a$distance) design$dist_range,
    position = design$position, spread = design$spread)), class = "ugdm_interaction_fit")
}

.interaction_predict <- function(fit, design, type) {
  if (!inherits(design, "ugdm_interaction_design") ||
      !identical(design$spread, fit$spread) || !identical(design$position, fit$position) ||
      !identical(design$basis, fit$basis) || (fit$has_distance &&
      (!identical(design$dist_range, fit$dist_range) ||
       !identical(attributes(design$dist_basis)[names(fit$dist_basis)], fit$dist_basis)))) {
    stop("Supply a matching interaction design built with the training reference.")
  }
  X <- .interaction_layout(design, fit$has_distance)$X
  if (!identical(colnames(X), names(fit$coef))) stop("Prediction columns differ from the fitted model.")
  eta <- drop(X %*% fit$coef)
  if (type == "link") eta else stats::plogis(eta)
}

#' Extended GDUM contrasts on the link and dissimilarity scales
#'
#' @param fit An extended fit from [ugdm_fit_interactions()].
#' @param design Evaluation design using the fitting reference basis.
#' @param weights Contrast matrix with one column per environmental set.
#'   For a two-by-two comparison, specify e.g. `c(1,-1,-1,1)`. The
#'   configurations and their common support should be declared beforehand.
#' @return Link and response contrasts, a row per supplied contrast.
#' @export
ugdm_interaction_contrast <- function(fit, design, weights) {
  W <- if (is.null(dim(weights))) matrix(weights, 1L) else as.matrix(weights)
  if (!is.numeric(W) || ncol(W) != nrow(design$idx) || any(!is.finite(W))) {
    stop("Supply finite contrast weights with one column per design row.")
  }
  eta <- ugdm_predict(fit, design, "link")
  data.frame(link = drop(W %*% eta), response = drop(W %*% stats::plogis(eta)))
}

#' Identifiability and hierarchical attribution for extended MS-GDUM
#'
#' Scaled information is the unpenalised fractional-logit working
#' information, so a stabilising penalty cannot hide a weak direction.
#' Position-dependent columns are assessed after projection off the spread
#' and position main effects. The base spread and position allocations
#' sum to the fitted additive deviance explained. Interaction attribution is
#' the difference between the full model and additive fit on the same rows.
#'
#' @param fit,design Matching extended fit and design.
#' @param y,trials,sampling_weights As in [ugdm_fit_interactions()].
#' @param contrasts Optional finite matrix of link-scale coefficient contrasts
#'   for reporting projection into weak working-information directions.
#' @return Rank, scaled condition, active constraints, predictor support and
#'   hierarchical spread, position and interaction deviances explained.
#' @export
ugdm_interaction_diagnostics <- function(fit, design, y, trials = rep(1, length(y)),
                                         sampling_weights = NULL, contrasts = NULL) {
  w <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  a <- .interaction_layout(design, fit$has_distance)
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
  additive <- design
  additive$dir <- additive$base
  additive$spread <- "constant"
  b <- .interaction_layout(additive, fit$has_distance)
  P <- .interaction_penalties(additive, b, fit$strengths,
    if (is.null(fit$interaction_penalty)) "coefficients" else fit$interaction_penalty,
    if (is.null(fit$penalty_scale)) "native" else fit$penalty_scale)
  # Equal endpoint coefficients reproduce the original design and curvature.
  # Their total ridge penalty has twice the environmental spread strength.
  if (design$spread == "varying") {
    ii <- 1L + seq_len(ncol(design$base))
    diag(P)[ii] <- diag(P)[ii] + fit$strengths[["ridge"]]
  }
  restricted <- function(idx) {
    f <- .interaction_fit_core(b$X[, idx, drop = FALSE], y, w, P[idx, idx, drop = FALSE],
                               which(idx %in% b$dir))
    if (f$converged) f$D2 else NA_real_
  }
  both <- if (design$spread == "constant") fit$D2 else restricted(seq_len(ncol(b$X)))
  spread_only <- restricted(c(1L, b$dir))
  pos_only <- restricted(c(1L, b$pos))
  variation <- function(X, Z) {
    Z <- Z * v
    norm <- colSums(Z^2)
    R <- qr.resid(qr(X * v), Z)
    ifelse(norm > 1e-24, sqrt(colSums(R^2) / norm), NA_real_)
  }
  nonlinear <- if (design$position == "nonlinear") b$pos[rep(c(FALSE, TRUE, TRUE), length(design$vars))] else integer()
  pos_resid <- if (length(nonlinear)) {
    stats::setNames(variation(b$X[, -nonlinear, drop = FALSE], b$X[, nonlinear, drop = FALSE]),
                    colnames(b$X)[nonlinear])
  } else NULL
  int_resid <- if (design$spread == "varying") {
    up_cols <- as.vector(outer(4:6, 6L * (seq_along(design$vars) - 1L), "+"))
    Z <- design$dir[, up_cols, drop = FALSE] - design$base / 2
    stats::setNames(variation(b$X, Z), colnames(Z))
  } else NULL
  p <- design$mean_position
  sd <- function(X) apply(X, 2L, stats::sd)
  list(rank = rank, parameters = ncol(a$X),
    condition = if (rank == ncol(a$X)) max(ev) / min(ev) else Inf,
    active = names(fit$coef)[a$dir[fit$coef[a$dir] <= 1e-8]],
    projected_gradient = fit$projected_gradient, converged = fit$converged,
    contrast_information = contrast,
    nonlinear_position_residual = pos_resid, interaction_residual = int_resid,
    contextual = if (design$spread == "contextual") ugdm_contextual_diagnostics(design, v^2),
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
      position = (both + pos_only - spread_only) / 2, interaction = fit$D2 - both),
    additive_D2 = both, full_D2 = fit$D2)
}

#' Select extended GDUM penalties by holding out groups of sites
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
#'   `spread_smooth`, `position_smooth`, `interaction_shrink` and `ridge`.
#'   Optionally include `spread_null` and `position_null`.
#' @param position,spread As in [ugdm_interaction_design()].
#' @param dist_matrix Optional site distance matrix.
#' @param M Combinations per training and per held-out fold.
#' @param seed Optional seed for combination sampling.
#' @param interaction_penalty,penalty_scale As in [ugdm_fit_interactions()].
#' @param contrast_loss Optional function of fit, held-out interaction
#'   design, response and working weights returning one nonnegative loss.
#' @param selection Selection criterion. Use contrast loss alone or require
#'   candidates to lie within a predeclared absolute held-out deviance
#'   tolerance of the minimum before comparing contrasts.
#' @param prediction_tolerance Held-out deviance tolerance for
#'   `selection = "contrast_within_prediction"`.
#' @param warm_start Reuse the preceding converged fit along the regularisation path.
#' @return Candidate scores, fold scores, fold-omission choices and the
#'   index of the selected row. Ties follow the supplied candidate order.
#' @export
ugdm_select_interactions <- function(env, folds, order, response, penalties,
                                     position = "nonlinear", spread = "varying",
                                     dist_matrix = NULL, M = 1000L, seed = NULL,
                                     interaction_penalty = c("coefficients", "derivative"),
                                     penalty_scale = c("native", "range"),
                                     contrast_loss = NULL,
                                     selection = c("predictive", "contrast", "contrast_within_prediction"),
                                     prediction_tolerance = 0, warm_start = TRUE) {
  interaction_penalty <- match.arg(interaction_penalty)
  penalty_scale <- match.arg(penalty_scale)
  selection <- match.arg(selection)
  fields <- c("spread_smooth", "position_smooth", "interaction_shrink", "ridge")
  if (!is.data.frame(penalties) || !nrow(penalties) ||
      !all(fields %in% names(penalties)) ||
      !all(vapply(penalties[fields], is.numeric, TRUE)) ||
      any(!is.finite(as.matrix(penalties[fields])) | as.matrix(penalties[fields]) < 0)) {
    stop("Provide a candidate data frame with four finite nonnegative strengths.")
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
  if ((!is.null(contrast_loss) && !is.function(contrast_loss)) ||
      (selection != "predictive" && is.null(contrast_loss))) {
    stop("Supply a contrast-loss function when selecting a contrast criterion.")
  }
  if (!is.numeric(prediction_tolerance) || length(prediction_tolerance) != 1L ||
      !is.finite(prediction_tolerance) || prediction_tolerance < 0) {
    stop("Supply a finite nonnegative prediction tolerance.")
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
      list(design = ugdm_interaction_design(design, env, position, spread),
           ref = design, y = values$y[keep], w = .model_weights(sum(keep), w[keep]))
    }
    a <- define(sample(train))
    b <- define(sample(test), a$ref)
    if (!any(b$w > 0)) stop("Each validation fold needs at least one positive working weight.")
    out <- matrix(Inf, nrow(penalties), 2L)
    current <- NULL
    for (j in path) {
      args <- c(list(design = a$design, y = a$y, trials = a$w,
                     interaction_penalty = interaction_penalty, penalty_scale = penalty_scale),
                as.list(penalties[j, fields, drop = FALSE]))
      f <- do.call(ugdm_fit_interactions, c(args, list(start = current)))
      if (!f$converged && !is.null(current)) f <- do.call(ugdm_fit_interactions, args)
      if (!f$converged) next
      if (warm_start) current <- f
      out[j, 1L] <- .binom_dev(b$y, ugdm_predict(f, b$design), b$w / sum(b$w))
      if (!is.null(contrast_loss)) {
        value <- contrast_loss(f, b$design, b$y, b$w)
        if (is.numeric(value) && length(value) == 1L && !is.na(value) && value >= 0) {
          out[j, 2L] <- value
        } else stop("contrast_loss must return one nonnegative number.")
      }
    }
    out
  })
  stack <- function(j) do.call(rbind, lapply(scores, function(s) t(s[, j, drop = FALSE])))
  fold_loss <- stack(1L)
  rownames(fold_loss) <- names(blocks)
  fold_contrast <- if (!is.null(contrast_loss)) {
    value <- stack(2L)
    rownames(value) <- names(blocks)
    value
  }
  average <- function(mat, rows = seq_along(blocks)) {
    split_rows <- split(rows, partition[rows])
    colMeans(do.call(rbind, lapply(split_rows, function(i) colMeans(mat[i, , drop = FALSE]))))
  }
  mean_loss <- average(fold_loss)
  mean_contrast <- if (is.null(contrast_loss)) NULL else average(fold_contrast)
  best <- function(p, c) {
    score <- if (selection == "predictive") p else c
    if (selection == "contrast_within_prediction") score[p > min(p) + prediction_tolerance] <- Inf
    if (any(is.finite(p) & is.finite(score))) which.min(replace(score, !is.finite(p), Inf)) else NA_integer_
  }
  choice <- best(mean_loss, mean_contrast)
  alternative <- vapply(seq_along(blocks), function(i) {
    rows <- setdiff(seq_along(blocks), i)
    best(average(fold_loss, rows), if (!is.null(fold_contrast)) average(fold_contrast, rows))
  }, integer(1L))
  list(penalties = penalties, mean_loss = mean_loss, fold_loss = fold_loss,
       mean_contrast_loss = mean_contrast, fold_contrast_loss = fold_contrast,
       partition = partition, selected = choice, selection = selection,
       selected_without_fold = stats::setNames(alternative, names(blocks)))
}
