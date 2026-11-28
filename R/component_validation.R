.component_support_reference <- function(design, w) {
  keep <- w > 0
  list(position = design$raw_position[keep, , drop = FALSE],
       spread = design$raw_spread[keep, , drop = FALSE],
       idx = design$idx[keep, , drop = FALSE])
}

.component_cells <- function(pos, spread, bins, pos_range = range(pos), spread_range = range(spread)) {
  index <- function(x, r) {
    outside <- x < r[1L] | x > r[2L]
    i <- if (diff(r) == 0) rep(1L, length(x)) else
      pmin(bins, floor(bins * (x - r[1L]) / diff(r)) + 1L)
    i[outside] <- NA_integer_
    as.integer(i)
  }
  i <- index(pos, pos_range)
  j <- index(spread, spread_range)
  data.frame(position_bin = i, spread_bin = j, cell = (i - 1L) * bins + j)
}

.component_support_cells <- function(reference, predictor, bins) {
  pos <- reference$position[, predictor]
  spread <- reference$spread[, predictor]
  cells <- .component_cells(pos, spread, bins)
  idx <- t(vapply(seq_len(nrow(reference$idx)), function(i) sort.int(reference$idx[i, ]),
                  integer(ncol(reference$idx))))
  keep <- !duplicated(idx)
  rows <- split(which(keep), cells$cell[keep])
  sets <- sites <- integer(bins^2)
  used <- as.integer(names(rows))
  sets[used] <- lengths(rows)
  sites[used] <- vapply(rows, function(i) length(unique(as.vector(idx[i, , drop = FALSE]))),
                        integer(1L))
  list(sets = sets, sites = sites, position_range = range(pos), spread_range = range(spread))
}

#' Observed spread-position support for component GDUM predictions
#'
#' A rectangular joint grid is formed separately for each predictor. Training
#' cells must contain the declared number of distinct site sets and sites.
#' This is a display and extrapolation diagnostic. It does not establish
#' overlap of all predictors jointly or turn reused sites into independent data.
#' @param fit,design As in [ugdm_predict_components()].
#' @param bins Number of equal-width bins along each axis.
#' @param min_sets,min_sites Minimum distinct training combinations and sites per cell.
#' @return Per-row support and a list of per-predictor support tables.
#' @export
ugdm_component_support <- function(fit, design, bins = 8L, min_sets = 8L, min_sites = 12L) {
  .component_prediction_matrix(fit, design)
  bins <- .set_design_integer(bins, "bins", 1L, 100L)
  min_sets <- .set_design_integer(min_sets, "min_sets", 1L)
  min_sites <- .set_design_integer(min_sites, "min_sites", 1L)
  detail <- lapply(fit$vars, function(p) {
    train <- .component_support_cells(fit$support, p, bins)
    cell <- .component_cells(design$raw_position[, p], design$raw_spread[, p],
                             bins, train$position_range, train$spread_range)
    cell$sets <- train$sets[cell$cell]
    cell$sites <- train$sites[cell$cell]
    cell$supported <- !is.na(cell$cell) & cell$sets >= min_sets & cell$sites >= min_sites
    cell
  })
  names(detail) <- fit$vars
  list(supported = Reduce("&", lapply(detail, "[[", "supported")), by_predictor = detail)
}

#' Component response surface on observed site-set configurations
#'
#' Averages model predictions over supplied site sets in each raw spread and
#' position cell. Other predictors and within-set configurations follow their
#' distribution among those sets. The surface is a conditional predictive
#' display. Empty or unsupported cells have missing predictions.
#' @param fit,design As in [ugdm_predict_components()].
#' @param predictor Name of the displayed environmental predictor.
#' @param bins,min_sets,min_sites As in [ugdm_component_support()].
#' @return Grid with raw position and spread, support counts and predictions.
#' @export
ugdm_component_surface <- function(fit, design, predictor, bins = 8L,
                                   min_sets = 8L, min_sites = 12L) {
  if (length(predictor) != 1L || !predictor %in% fit$vars) stop("Choose one fitted predictor.")
  supported <- ugdm_component_support(fit, design, bins, min_sets, min_sites)
  cell <- supported$by_predictor[[predictor]]
  train <- .component_support_cells(fit$support, predictor, bins)
  centers <- function(r) r[1L] + (seq_len(bins) - .5) * diff(r) / bins
  out <- data.frame(position = rep(centers(train$position_range), each = bins),
    spread = rep(centers(train$spread_range), times = bins),
    training_sets = train$sets, training_sites = train$sites)
  pred <- ugdm_predict_components(fit, design)
  rows <- split(which(supported$supported), cell$cell[supported$supported])
  out$evaluated_sets <- integer(nrow(out))
  at <- as.integer(names(rows))
  out$evaluated_sets[at] <- lengths(rows)
  out$supported <- out$evaluated_sets > 0L
  for (name in names(pred)) {
    out[[name]] <- NA_real_
    out[[name]][at] <- vapply(rows, function(i) mean(pred[[name]][i]), numeric(1L))
  }
  out
}

#' Working-information diagnostics for component MS-GDUM
#'
#' Reports unpenalised scaled design and information rank, condition and
#' active constraints for each fitted part. For the coherent share, only rows
#' with positive total contribute.
#' @param fit,design As in [ugdm_predict_components()].
#' @param response Training responses.
#' @return A list of stage diagnostics and the identification status of the share.
#' @export
ugdm_component_diagnostics <- function(fit, design, response) {
  .component_match_response(design, response)
  X <- .component_prediction_matrix(fit, design)
  if (!identical(fit$sets, design$idx)) stop("Supply the original training design.")
  rank <- function(Z) {
    scale <- sqrt(colSums(Z^2))
    S <- sweep(Z, 2L, ifelse(scale > 0, scale, 1), "/")
    d <- svd(S, nu = 0, nv = 0)$d
    r <- sum(d > max(d) * max(dim(S)) * .Machine$double.eps)
    c(rank = r, condition = if (r == ncol(S)) (max(d) / min(d))^2 else Inf)
  }
  rows <- lapply(names(fit$fits), function(name) {
    a <- fit$fits[[name]]
    w <- fit$weights * if (name == "share") response[, "total"] else 1
    mu <- if (is.finite(a$constant_value)) rep(a$constant_value, nrow(X)) else plogis(drop(X %*% a$coef))
    design_info <- rank(X * sqrt(w))
    info <- rank(X * sqrt(w * mu * (1 - mu)))
    bounded <- if (name == "total") which(fit$keep %in% fit$direction) else integer()
    data.frame(stage = name, parameters = ncol(X), design_rank = design_info[["rank"]],
      information_rank = info[["rank"]], condition = info[["condition"]],
      active = sum(a$coef[bounded] <= 1e-8), converged = a$converged,
      projected_gradient = a$projected_gradient, nonzero_weight_sets = sum(w > 0))
  })
  list(stages = do.call(rbind, rows), share_identified = fit$share_identified)
}

.component_assessment <- function(response, predictions, means, weight) {
  components <- c("total", "occupancy", "heterogeneity")
  if (!length(weight) || sum(weight) == 0) {
    return(data.frame(component = components, RMSE = NA_real_, bias = NA_real_,
        calibration_slope = NA_real_, deviance = NA_real_, null_deviance = NA_real_, D2 = NA_real_, n = 0L))
  }
  w <- weight / sum(weight)
  do.call(rbind, lapply(components, function(p) {
    y <- response[, p]
    mu <- predictions[[p]]
    null <- .binom_dev(y, rep(means[[p]], length(y)), w)
    dev <- .binom_dev(y, mu, w)
    centered <- mu - sum(w * mu)
    variance <- sum(w * centered^2)
    slope <- if (variance > 1e-15) sum(w * centered * y) / variance else NA_real_
    data.frame(component = p, RMSE = sqrt(sum(w * (mu - y)^2)),
      bias = sum(w * (mu - y)), calibration_slope = slope, deviance = dev, null_deviance = null,
      D2 = if (null > 1e-10) 1 - dev / null else NA_real_, n = length(y))
  }))
}

#' Assess component predictions on independent site sets
#' @param fit,design As in [ugdm_predict_components()].
#' @param response Evaluation response.
#' @param training_response Response used to fit the model for null predictions.
#' @param mask Restrict assessment to rows supported for every predictor.
#' @param ... Support settings for [ugdm_component_support()].
#' @return Component RMSE, calibration bias and slope, deviance and explained deviance
#'   relative to the training component means. Closure error is reported separately.
#' @export
ugdm_assess_components <- function(fit, design, response, training_response, mask = FALSE, ...) {
  .component_match_response(design, response)
  if (!inherits(training_response, "ugdm_component_response") ||
      !identical(attr(training_response, "sets"), fit$sets)) stop("Supply the original training response.")
  pred <- ugdm_predict_components(fit, design, mask = mask, ...)
  keep <- is.finite(pred$total)
  means <- colSums(training_response * (fit$weights / sum(fit$weights)))
  metrics <- .component_assessment(response[keep, , drop = FALSE], pred[keep, ], means,
                                  rep(1, sum(keep)))
  list(metrics = metrics, evaluated_fraction = mean(keep),
       closure_RMSE = if (any(keep)) sqrt(mean(pred$closure_error[keep]^2)) else NA_real_)
}

#' Matched independent spatial validation of component MS-GDUM
#'
#' Sets contain only training sites or only validation sites. Each fold
#' estimates its basis exclusively on training sites. All models share
#' sampled sets within a fold. Fold-level scores are descriptive, especially
#' with few spatial blocks.
#' @param community,environment,coordinates As in [ugdm_component_design()]
#'   and [ugdm_component_response()]. Coordinates are required to construct
#'   folds when fold labels are absent.
#' @param order Number of distinct sites in a set.
#' @param fold Optional site-level validation-block labels.
#' @param blocks Number of coordinate clusters when fold labels are absent.
#' @param models,architectures Component model and architecture choices.
#' @param combinations Training and held-out combinations per fold.
#' @param ridge,monotone_total As in [ugdm_fit_components()].
#' @param geographic Include a distance predictor in addition to spatial folds.
#' @param lonlat Use great-circle distance for geographic coordinates.
#' @param seed Optional random seed.
#' @return Fold-level predictive metrics, overlap, closure and convergence.
#' @export
ugdm_validate_components <- function(community, environment, coordinates = NULL, order = 2L,
    fold = NULL, blocks = 4L, models = c("spread_only", "additive"),
    architectures = c("independent", "coherent"), combinations = 1000L,
    ridge = .001, monotone_total = TRUE,
    geographic = FALSE, lonlat = FALSE, seed = NULL) {
  if (anyNA(community)) stop("Supply complete community observations.")
  x <- .conditional_community(community)
  env <- as.data.frame(environment)
  n <- nrow(x)
  if (nrow(env) != n) stop("Community and environment rows must agree.")
  order <- .zeta_order_integer(order, n, 2L)
  combinations <- .set_design_integer(combinations, "combinations", 2L)
  models <- unique(match.arg(models, c("spread_only", "additive"), several.ok = TRUE))
  architectures <- unique(match.arg(architectures, c("independent", "coherent"), several.ok = TRUE))
  if (is.null(fold) || geographic) {
    xy <- as.matrix(coordinates)
    if (!is.numeric(xy) || !identical(dim(xy), c(n, 2L)) || any(!is.finite(xy))) {
      stop("Supply two finite coordinate columns in site order.")
    }
  }
  if (is.null(fold)) {
    blocks <- .set_design_integer(blocks, "blocks", 2L, n)
    h <- stats::hclust(stats::as.dist(coords_dist(xy, lonlat)), method = "complete")
    fold <- stats::cutree(h, k = blocks)
  }
  if (length(fold) != n || anyNA(fold) || length(unique(fold)) < 2L) stop("Supply at least two complete fold labels.")
  counts <- table(fold)
  if (min(counts) < order || min(n - counts) < max(4L, order)) {
    stop("Every fold needs enough training and validation sites for the requested order.")
  }
  if (!is.null(seed)) set.seed(seed)
  result <- lapply(unique(fold), function(b) {
    train <- which(fold != b)
    test <- which(fold == b)
    i <- sample_sets(length(train), order, combinations)
    j <- sample_sets(length(test), order, combinations)
    y <- ugdm_component_response(x[train, , drop = FALSE], i)
    v <- ugdm_component_response(x[test, , drop = FALSE], j)
    coords <- function(z) if (geographic) xy[z, , drop = FALSE] else NULL
    do.call(rbind, lapply(models, function(m) {
      d <- ugdm_component_design(env[train, , drop = FALSE], i, m, coords(train), lonlat = lonlat)
      heldout <- ugdm_component_design(env[test, , drop = FALSE], j,
        coordinates = coords(test), reference = d, lonlat = lonlat)
      do.call(rbind, lapply(architectures, function(a) {
        t <- proc.time()[["elapsed"]]
        fit <- ugdm_fit_components(d, y, a, ridge, monotone_total)
        seconds <- proc.time()[["elapsed"]] - t
        out <- ugdm_assess_components(fit, heldout, v, y)
        support <- ugdm_assess_components(fit, heldout, v, y, mask = TRUE)
        data.frame(fold = b, order, model = m, architecture = a, sites_train = length(train),
          sites_test = length(test), converged = fit$converged, out$metrics,
          supported_D2 = support$metrics$D2, supported_fraction = support$evaluated_fraction,
          closure_RMSE = out$closure_RMSE, fit_seconds = seconds)
      }))
    }))
  })
  do.call(rbind, result)
}
