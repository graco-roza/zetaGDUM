.within_resample_summary <- function(predictor, threshold) {
  interval <- function(x) {
    if (!length(x)) return(rep(NA_real_, 5L))
    stats::quantile(x, c(.5, .1, .9, .025, .975), names = FALSE)
  }
  quantiles <- function(a) {
    out <- t(vapply(seq_len(ncol(a)), function(j) interval(a[, j]), numeric(5L)))
    colnames(out) <- c("median", "lo80", "hi80", "lo95", "hi95")
    out
  }
  lead <- rep(NA_real_, ncol(predictor))
  if (nrow(predictor)) {
    top <- predictor == apply(predictor, 1L, max)
    lead <- colMeans(top / rowSums(top))
  }
  predictors <- data.frame(predictor = colnames(predictor), quantiles(predictor),
                           p_leading = lead, row.names = NULL)
  g <- ncol(predictor)
  if (g < 2L) return(list(predictors = predictors, differences = NULL, difference_draws = NULL))
  pairs <- utils::combn(g, 2L)
  delta <- predictor[, pairs[1L, ], drop = FALSE] - predictor[, pairs[2L, ], drop = FALSE]
  left <- colnames(predictor)[pairs[1L, ]]
  right <- colnames(predictor)[pairs[2L, ]]
  colnames(delta) <- paste(left, right, sep = " - ")
  probability <- function(condition) if (nrow(delta)) colMeans(condition) else rep(NA_real_, ncol(delta))
  differences <- data.frame(predictor = left, comparator = right, quantiles(delta),
    p_greater = probability(delta > 0), p_less = probability(delta < 0),
    p_practical = probability(delta > threshold),
    p_practical_reverse = probability(delta < -threshold), row.names = NULL)
  list(predictors = predictors, differences = differences, difference_draws = delta)
}

.within_score_inputs <- function(p, design, n) {
  idx <- .species_indices(design$idx, n)
  effective <- p$weights / sum(p$weights)
  mass <- score_inclusion_cpp(idx, matrix(1, nrow(idx), 1L), effective * 0,
    effective, n)[seq_len(n), 1L]
  if (any(mass >= 1 - 1e-12)) stop("A site deletion has no support in a score pool.")
  list(X = p$X, y = p$y, idx = t(idx), effective = effective, mass = mass)
}

.within_centered_score <- function(p, keep, theta, n) {
  inc <- site_score_kernel_cpp(p$idx, t(p$X[, keep, drop = FALSE]), p$y, p$effective, theta, n)
  deleted <- sweep(-inc[seq_len(n), , drop = FALSE], 2L, inc[n + 1L, ], "+") / (1 - p$mass)
  sweep(deleted, 2L, colMeans(deleted)) * (1 - n)
}

.within_score_state <- function(p, training, keep, theta, bounded, pools, n, groups, score_control) {
  X <- p$X[, keep, drop = FALSE]
  ww <- training$effective
  P <- p$penalty[keep, keep, drop = FALSE]
  mu <- stats::plogis(drop(X %*% theta))
  g <- drop(crossprod(X, ww * (p$y - mu)) - P %*% theta)
  H <- crossprod(X * sqrt(ww * mu * (1 - mu))) + P
  scale <- sqrt(diag(H))
  if (any(!is.finite(scale) | scale <= 0)) stop("A coalition has a zero information direction.")
  scaled <- H / tcrossprod(scale)
  chol_H <- chol(scaled)
  inv <- chol2inv(chol_H) / tcrossprod(scale)
  H <- crossprod(chol_H) * tcrossprod(scale)
  raw <- if (is.null(score_control)) {
    Reduce("+", lapply(pools, .within_centered_score, keep = keep, theta = theta, n = n)) / length(pools)
  } else as.matrix(score_control(keep, theta))
  if (!is.numeric(raw) || !identical(dim(raw), c(n, length(keep))) || any(!is.finite(raw))) {
    stop("The score function must return one finite row per site and column per coalition coefficient.")
  }
  raw <- sweep(raw, 2L, colMeans(raw))
  by_group <- rowsum(raw, groups, reorder = TRUE)
  G <- nrow(by_group)
  if (G < 2L) stop("Score resampling needs two or more independent groups.")
  by_group <- by_group * sqrt(G / (G - 1)) / n
  Omega <- crossprod(by_group)
  units <- sqrt(diag(Omega))
  units[units == 0] <- 1
  ev <- eigen(Omega / tcrossprod(units), symmetric = TRUE, only.values = TRUE)$values
  rank <- sum(ev > max(ev) * max(1L, length(ev)) * .Machine$double.eps)
  condition <- if (rank == length(ev) && min(ev) > 0) max(ev) / min(ev) else Inf
  project <- if (length(bounded)) .gmm_project(H, bounded) else identity
  center <- colSums(X * ww)
  center[1L] <- 0
  optimizer_scale <- sqrt(colSums(sweep(X, 2L, center)^2 * ww))
  optimizer_scale[optimizer_scale < 1e-12] <- 1
  list(beta = theta, g = g, inv = inv, project = project,
    scores = by_group, covariance_condition = condition, covariance_rank = rank,
    information_condition = kappa(scaled, exact = TRUE), X = X,
    P = P, bounded = bounded, optimizer_center = center, optimizer_scale = optimizer_scale)
}

#' Within-order coalition resampling for extended MS-GDUM
#'
#' Each draw shares the same site weights, combination sample, evaluation rows
#' and permutation orders across coalitions. Exact Bayesian bootstrap draws
#' sample with probability proportional to products of Dirichlet site weights.
#' The site bootstrap samples site copies with replacement. Score methods
#' use a common random multiplier for every coalition, with optional block
#' sums, and project each quadratic update in the coalition loss metric.
#' At interior fits the Gaussian score perturbation is the flat-prior
#' Gaussian GMM approximation for square estimating equations. At constraints
#' it is a local approximation to the perturbed constrained loss.
#'
#' @param fit,design,y,trials,sampling_weights,evaluation,evaluation_y,evaluation_trials,evaluation_sampling_weights
#'   As in \code{\link{ugdm_within_order_importance}}. Score methods evaluate
#'   on fixed sets, using training sets when no separate evaluation is supplied.
#' @param method Direct Bayesian bootstrap, Gaussian GMM approximation,
#'   Rademacher multiplier approximation, or full site resampling.
#' @param draws Number of resampling draws.
#' @param M Number of sets per direct resampling draw.
#' @param build Function of a set index matrix and original site IDs
#'   (function(idx, sites)). Returns design, y, and
#'   optionally trials, and must keep the original fitting basis.
#'   A duplicated original site ID represents an independent resampled copy.
#' @param blocks Optional block labels in original site order. Common
#'   weights act within blocks for Bayesian bootstrap and score resampling.
#' @param threshold Nonnegative minimum practical difference in absolute
#'   explained-deviance units.
#' @param seed Optional seed. Per-draw seeds make serial and forked runs match.
#' @param permutations,exact_max,min_D2 As in \code{\link{ugdm_within_order_importance}}.
#' @param workers Number of fork workers on supporting platforms.
#' @param score_pools Optional list of independent score pools. Each element
#'   supplies design, y, optional trials and
#'   sampling_weights. Pool pseudovalues are averaged before covariance
#'   estimation. The default uses the training sets.
#' @param score_control Optional function of (keep, theta) returning an
#'   aligned n-by-length(keep) matrix of site scores from a separate design,
#'   such as independently piloted analytical control variates. Use score
#'   units of the average working loss, excluding the deterministic penalty.
#' @param score_update Use a projected quadratic update, solve the exact
#'   perturbed fitting objective, or refit when a hybrid check finds a large
#'   projected score residual. These options apply to score resampling.
#' @param hybrid_tolerance Maximum absolute projected score residual in the
#'   scaled optimisation coordinates before a hybrid draw refits a coalition.
#' @param n Number of sites in older extended designs that do not store it.
#' @return Draw matrices, conditional importance intervals and ranking
#'   probabilities, fit failures, elapsed time and score diagnostics.
#' @export
ugdm_within_order_resample <- function(fit, design, y, trials = rep(1, length(y)),
    sampling_weights = NULL, evaluation = NULL, evaluation_y = NULL,
    evaluation_trials = NULL, evaluation_sampling_weights = NULL,
    method = c("bayesian_bootstrap", "gaussian_gmm", "multiplier", "site_bootstrap"),
    draws = 199L, M = nrow(design$idx), build = NULL, blocks = NULL,
    threshold = 0, seed = NULL, permutations = 1000L, exact_max = 8L,
    min_D2 = 1e-4, workers = 1L, score_pools = NULL, score_control = NULL,
    n = design$n, score_update = c("one_step", "refit", "hybrid"),
    hybrid_tolerance = .001) {
  started <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  score_update <- match.arg(score_update)
  if (length(hybrid_tolerance) != 1L || !is.finite(hybrid_tolerance) || hybrid_tolerance < 0) {
    stop("hybrid_tolerance must be a nonnegative finite scalar.")
  }
  draws <- .set_design_integer(draws, "draws", 1L)
  n <- .set_design_integer(n, "n", 3L)
  workers <- .set_design_integer(workers, "workers", 1L)
  k <- design$order
  if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
    stop("threshold must be a nonnegative finite scalar.")
  }
  if (k >= n || n < max(design$idx)) stop("Resampling requires valid original sites and order < n.")
  if (!is.null(seed)) set.seed(seed)
  p <- .within_importance_setup(fit, design, y, trials, sampling_weights,
    evaluation, evaluation_y, evaluation_trials, evaluation_sampling_weights)
  plan <- .within_importance_plan(p$groups, permutations, exact_max)
  if (is.null(blocks)) blocks <- seq_len(n)
  if (length(blocks) != n || anyNA(blocks)) stop("Provide one non-missing block label per site.")
  groups <- match(blocks, unique(blocks))
  G <- max(groups)
  if (G < 2L) stop("Provide at least two blocks.")
  exact <- method %in% c("bayesian_bootstrap", "site_bootstrap")
  if (exact && !is.function(build)) stop("Direct site resampling requires a build function.")
  M <- .set_design_integer(M, "M", 1L)
  if (!is.null(score_control) && !is.function(score_control)) stop("score_control must be a function.")
  if (method == "site_bootstrap" && G != n) stop("The full site bootstrap requires individual sites.")
  if (exact && score_update != "one_step") stop("Score updates apply only to score resampling.")

  training <- states <- NULL
  if (!exact) {
    training <- .within_score_inputs(p, design, n)
    pools <- if (is.null(score_pools)) list(training) else lapply(score_pools, function(a) {
      if (a$design$order != k || (!is.null(a$design$n) && a$design$n != n)) stop("Align score-pool sites and orders.")
      pp <- .within_importance_setup(fit, a$design, a$y,
        if (is.null(a$trials)) rep(1, length(a$y)) else a$trials, a$sampling_weights,
        evaluation, evaluation_y, evaluation_trials, evaluation_sampling_weights)
      .within_score_inputs(pp, a$design, n)
    })
    if (!length(pools)) stop("Supply at least one score pool.")
    states <- new.env(parent = emptyenv())
    store <- function(keep, start) {
      bounded <- which(keep %in% p$bounded)
      a <- .variant_solution(p$X[, keep, drop = FALSE], p$y, p$weights,
                             p$penalty[keep, keep, drop = FALSE], bounded, start)
      if (!a$converged) a <- .variant_solution(p$X[, keep, drop = FALSE], p$y, p$weights,
        p$penalty[keep, keep, drop = FALSE], bounded)
      if (a$converged) {
        state <- .within_score_state(p, training, keep, a$coef, bounded, pools, n, groups, score_control)
        assign(paste(keep, collapse = ","), state, envir = states)
      }
      a
    }
    .within_importance_game(p, permutations, exact_max, min_D2, fit_coalition = store, plan = plan)
  }
  seeds <- sample.int(.Machine$integer.max, draws)
  run <- function(i) {
    start_time <- proc.time()[["elapsed"]]
    set.seed(seeds[i])
    tryCatch({
      refits <- 0L
      if (exact) {
        sites <- if (method == "site_bootstrap") sample.int(n, n, replace = TRUE) else seq_len(n)
        weights <- if (method == "bayesian_bootstrap") site_weights(G, 1L)[groups] else rep(1, n)
        idx <- sample_sets_weighted(weights, k, M)
        a <- build(idx, sites)
        q <- .within_importance_setup(fit, a$design, a$y,
          if (is.null(a$trials)) rep(1, length(a$y)) else a$trials, a$sampling_weights,
          evaluation, evaluation_y, evaluation_trials, evaluation_sampling_weights)
        ans <- .within_importance_game(q, permutations, exact_max, min_D2, plan = plan)
      } else {
        e <- if (method == "gaussian_gmm") stats::rnorm(G) else sample(c(-1, 1), G, replace = TRUE)
        update <- function(keep, start) {
          z <- get(paste(keep, collapse = ","), envir = states, inherits = FALSE)
          change <- drop(crossprod(z$scores, e))
          coef <- z$project(z$beta + drop(z$inv %*% (z$g + change)))
          if (score_update == "one_step") return(list(coef = coef, converged = all(is.finite(coef))))
          if (score_update == "hybrid") {
            mu <- stats::plogis(drop(z$X %*% coef))
            gradient <- drop(crossprod(z$X, training$effective * (mu - p$y)) +
              z$P %*% coef - change)
            gradient <- (gradient - gradient[1L] * z$optimizer_center) / z$optimizer_scale
            active <- z$bounded[coef[z$bounded] <= 1e-8]
            gradient[active] <- pmin(gradient[active], 0)
            if (max(abs(gradient)) <= hybrid_tolerance) {
              return(list(coef = coef, converged = all(is.finite(coef))))
            }
          }
          refits <<- refits + 1L
          .variant_solution(z$X, p$y, p$weights, z$P, z$bounded, coef, score_shift = change)
        }
        ans <- .within_importance_game(p, permutations, exact_max, min_D2,
                                        fit_coalition = update, plan = plan)
      }
      if (!all(is.finite(ans$predictors$absolute))) {
        stop("Undefined predictor attribution for this evaluation sample.")
      }
      list(value = ans, refits = refits, elapsed = proc.time()[["elapsed"]] - start_time)
    }, error = function(e) list(error = conditionMessage(e),
      elapsed = proc.time()[["elapsed"]] - start_time))
  }
  result <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(draws), run, mc.cores = workers, mc.set.seed = FALSE)
  } else lapply(seq_len(draws), run)
  ok <- vapply(result, function(x) !is.null(x$value) && all(is.finite(x$value$predictors$absolute)),
    logical(1L))
  value <- lapply(result[ok], "[[", "value")
  extract <- function(group, scale, labels) {
    if (!length(value)) return(matrix(numeric(), 0, length(labels), dimnames = list(NULL, labels)))
    do.call(rbind, lapply(value, function(x) x[[group]][[scale]]))
  }
  predictor <- extract("predictors", "absolute", names(p$groups))
  component <- extract("components", "absolute", unique(unlist(lapply(p$groups, names))))
  relative <- extract("predictors", "relative", names(p$groups))
  hierarchy <- do.call(rbind, lapply(value, function(x) {
    h <- x$hierarchy
    stats::setNames(h$absolute, paste(h$predictor, h$component, sep = ":"))
  }))
  sums <- .within_resample_summary(predictor, threshold)
  diagnostics <- if (exact) NULL else do.call(rbind, lapply(ls(states), function(key) {
    z <- states[[key]]
    data.frame(coalition = key, covariance_condition = z$covariance_condition,
               covariance_rank = z$covariance_rank, coefficients = length(z$beta),
               information_condition = z$information_condition)
  }))
  c(sums, list(predictor_draws = predictor, component_draws = component,
    relative_draws = relative, hierarchy_draws = hierarchy,
    components = .within_resample_summary(component, threshold)$predictors,
    relative = .within_resample_summary(relative[stats::complete.cases(relative), , drop = FALSE],
                                        threshold)$predictors,
    valid = ok, errors = vapply(result, function(x) if (is.null(x$error)) "" else x$error, ""),
    draw_seconds = vapply(result, "[[", numeric(1L), "elapsed"),
    exact_score_refits = vapply(result, function(x) if (is.null(x$refits)) NA_integer_ else x$refits, 0L),
    elapsed = proc.time()[["elapsed"]] - started, method = method,
    score_diagnostics = diagnostics, n = n, groups = G, order = k))
}
