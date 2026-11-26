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
#' @param spread Position-independent or position-dependent directional spread.
#' @return Design for [ugdm_fit_interactions()].
#' @export
ugdm_interaction_design <- function(design, env, position = c("linear", "nonlinear"),
                                    spread = c("constant", "varying")) {
  position <- match.arg(position)
  spread <- match.arg(spread)
  if (!is.null(design$spread_width)) stop("Supply the original environmental I-spline design.")
  z <- .environment_for_moments(env, design)$z
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
  structure(list(dir = dir, pos = pos, base = base, geo = design$Xdist,
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

.interaction_penalties <- function(d, layout, strengths) {
  P <- diag(c(0, rep(strengths[["ridge"]], ncol(layout$X) - 1L)))
  for (j in seq_along(d$vars)) {
    width <- if (d$spread == "varying") 6L else 3L
    at <- 1L + width * (j - 1L) + seq_len(width)
    C <- d$curvature[[j]]
    if (d$spread == "constant") {
      P[at, at] <- P[at, at] + strengths[["spread_smooth"]] * C
    } else {
      interaction <- cbind(-diag(3L), diag(3L))
      A <- matrix(0, 6L, 6L)
      A[1:3, 1:3] <- A[4:6, 4:6] <- C / 2
      P[at, at] <- P[at, at] + strengths[["spread_smooth"]] * A +
        strengths[["interaction_shrink"]] * crossprod(interaction)
    }
    if (d$position == "nonlinear") {
      at <- layout$pos[3L * (j - 1L) + 2:3]
      P[at, at] <- P[at, at] + strengths[["position_smooth"]] *
        d$position_curvature[j] * diag(c(8, 24))
    }
  }
  P
}

.interaction_fit_core <- function(X, y, w, P, nonnegative) {
  fit <- .variant_solution(X, y, w, P, nonnegative)
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
#' @return An extended GDUM fit for [ugdm_predict()].
#' @export
ugdm_fit_interactions <- function(design, y, trials = rep(1, length(y)),
                                  sampling_weights = NULL, distance = TRUE,
                                  spread_smooth = 0, position_smooth = 0,
                                  interaction_shrink = 0, ridge = 0) {
  if (!inherits(design, "ugdm_interaction_design")) stop("Supply an interaction design.")
  if (length(y) != nrow(design$idx) || any(!is.finite(y) | y < 0 | y > 1)) {
    stop("Supply one finite fractional response per design row.")
  }
  supplied <- list(spread_smooth = spread_smooth, position_smooth = position_smooth,
                   interaction_shrink = interaction_shrink, ridge = ridge)
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  if (!all(vapply(supplied, valid, TRUE))) {
    stop("Supply four finite nonnegative penalty strengths.")
  }
  strengths <- vapply(supplied, function(x) unname(x), numeric(1L))
  w <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  a <- .interaction_layout(design, distance)
  P <- .interaction_penalties(design, a, strengths)
  fit <- .interaction_fit_core(a$X, y, w, P, a$dir)
  fit$coef <- stats::setNames(fit$coef, colnames(a$X))
  structure(c(fit, list(strengths = strengths, penalty = P, vars = design$vars,
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
  P <- .interaction_penalties(additive, b, fit$strengths)
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
#' @param folds One block label per site. At least two blocks are required.
#' @param order Order of the GDUM response.
#' @param response Function accepting an index matrix with original site
#'   numbers. Return a numeric response or a list with `y`, optionally
#'   `trials` and `keep`, as in [ugdm_response()].
#' @param penalties Data frame with candidate rows and columns
#'   `spread_smooth`, `position_smooth`, `interaction_shrink` and `ridge`.
#' @param position,spread As in [ugdm_interaction_design()].
#' @param dist_matrix Optional site distance matrix.
#' @param M Combinations per training and per held-out fold.
#' @param seed Optional seed for combination sampling.
#' @return Candidate scores, fold scores, fold-omission choices and the
#'   index of the selected row. Ties follow the supplied candidate order.
#' @export
ugdm_select_interactions <- function(env, folds, order, response, penalties,
                                     position = "nonlinear", spread = "varying",
                                     dist_matrix = NULL, M = 1000L, seed = NULL) {
  fields <- c("spread_smooth", "position_smooth", "interaction_shrink", "ridge")
  if (!is.data.frame(penalties) || !nrow(penalties) ||
      !all(fields %in% names(penalties)) ||
      !all(vapply(penalties[fields], is.numeric, TRUE)) ||
      any(!is.finite(as.matrix(penalties[fields])) | as.matrix(penalties[fields]) < 0)) {
    stop("Provide a candidate data frame with four finite nonnegative strengths.")
  }
  if (length(folds) != nrow(env) || anyNA(folds) || length(unique(folds)) < 2L) {
    stop("Supply a block label for every site and at least two blocks.")
  }
  if (!is.function(response)) stop("Supply a set-level response function.")
  k <- .set_design_integer(order, "order", 2L)
  M <- .set_design_integer(M, "M", 1L)
  blocks <- split(seq_len(nrow(env)), folds, drop = TRUE)
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
  fold_loss <- t(vapply(blocks, function(test) {
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
    vapply(seq_len(nrow(penalties)), function(j) {
      f <- do.call(ugdm_fit_interactions, c(list(design = a$design, y = a$y, trials = a$w),
                                         as.list(penalties[j, fields, drop = FALSE])))
      if (f$converged) .binom_dev(b$y, ugdm_predict(f, b$design), b$w / sum(b$w)) else Inf
    }, numeric(1L))
  }, numeric(nrow(penalties))))
  if (nrow(penalties) == 1L) fold_loss <- t(fold_loss)
  mean_loss <- colMeans(fold_loss)
  best <- function(x) if (any(is.finite(x))) which.min(x) else NA_integer_
  choice <- best(mean_loss)
  alternative <- vapply(seq_along(blocks), function(i) {
    best(colMeans(fold_loss[-i, , drop = FALSE]))
  }, integer(1L))
  list(penalties = penalties, mean_loss = mean_loss, fold_loss = fold_loss,
       selected = choice, selected_without_fold = stats::setNames(alternative, names(blocks)))
}
