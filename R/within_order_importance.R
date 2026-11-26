#' Within-order hierarchical explained-deviance attribution
#'
#' Each predictor is one group of its directional and position terms.
#' Between predictors the allocation is Shapley. Within each predictor the
#' Owen allocation averages the order of its spread and position terms,
#' conditional on each other-predictor coalition. These component allocations
#' sum to the corresponding predictor contributions. Geographic distance,
#' when fitted, is a single-component group.
#'
#' Coalitions share the declared penalty, coefficient bounds and evaluation
#' sets. Held-out gains can be negative. The intercept-only reference is
#' estimated on the training sets, including for held-out evaluation.
#' @param fit Fitted [ugdm_fit_interactions()] model.
#' @param design Its training design.
#' @param y,trials,sampling_weights Training response and weights.
#' @param evaluation Optional held-out design using the training basis.
#' @param evaluation_y,evaluation_trials,evaluation_sampling_weights Held-out
#'   response and weights, required with evaluation.
#' @param permutations Number of random blocked predictor orderings when
#'   the number of predictors exceeds `exact_max`.
#' @param seed Optional permutation seed.
#' @param min_D2 Minimum positive total for relative shares.
#' @param exact_max Maximum group count for exact subset calculation, from
#'   zero through eight. Zero forces a permutation approximation.
#' @param failure Use `"error"` to stop or `"warn"` to return affected
#'   allocations as missing with coalition failure details.
#' @return Predictor Shapley and aggregated component Owen contributions.
#'   The hierarchy table gives one row per predictor-component pair.
#'   Monte Carlo standard errors are conditional on the fitted data and
#'   measure permutation error only. Diagnostics identify invalid coalition
#'   fits and undefined intercept-only deviance.
#' @export
ugdm_within_order_importance <- function(fit, design, y, trials = rep(1, length(y)),
    sampling_weights = NULL, evaluation = NULL, evaluation_y = NULL,
    evaluation_trials = NULL, evaluation_sampling_weights = NULL,
    permutations = 1000L, seed = NULL, min_D2 = 1e-4,
    exact_max = 8L, failure = c("error", "warn")) {
  problem <- .within_importance_setup(fit, design, y, trials, sampling_weights,
    evaluation, evaluation_y, evaluation_trials, evaluation_sampling_weights)
  if (!is.null(seed)) set.seed(seed)
  .within_importance_game(problem, permutations, exact_max, min_D2,
    failure = match.arg(failure))
}

.within_importance_setup <- function(fit, design, y, trials = NULL,
    sampling_weights = NULL, evaluation = NULL, evaluation_y = NULL,
    evaluation_trials = NULL, evaluation_sampling_weights = NULL) {
  if (!inherits(fit, "ugdm_interaction_fit") || !inherits(design, "ugdm_interaction_design")) {
    stop("Use a fitted extended MS-GDUM and matching training design.")
  }
  verify <- function(d) {
    if (!inherits(d, "ugdm_interaction_design") ||
        !identical(d$order, fit$order) || !identical(d$vars, fit$vars)) {
      stop("Each importance design must match the fitted predictors and zeta order.")
    }
    # Prediction also checks environmental and geographic bases and all column names.
    invisible(.interaction_predict(fit, d, "link"))
  }
  verify(design)
  if (!is.null(evaluation)) {
    if (is.null(evaluation_y)) stop("Supply the held-out response.")
    verify(evaluation)
  }
  train <- .interaction_layout(design, fit$has_distance)
  X <- train$X
  check_response <- function(z, n) {
    if (!is.numeric(z) || length(z) != n || any(!is.finite(z) | z < 0 | z > 1)) {
      stop("Supply one finite fractional response per design row.")
    }
  }
  check_response(y, nrow(X))
  test_X <- if (is.null(evaluation)) X else .interaction_layout(evaluation, fit$has_distance)$X
  test_y <- if (is.null(evaluation)) y else evaluation_y
  check_response(test_y, nrow(test_X))
  weights <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  test_weights <- if (is.null(evaluation)) weights else
    .model_weights(length(test_y), evaluation_trials, sampling_weights = evaluation_sampling_weights)

  penalty <- fit$penalty
  expected <- .interaction_penalties(design, train, fit$strengths,
    if (is.null(fit$interaction_penalty)) "coefficients" else fit$interaction_penalty,
    if (is.null(fit$penalty_scale)) "native" else fit$penalty_scale)
  if (!isTRUE(all.equal(unname(penalty), unname(expected), tolerance = 1e-10))) {
    stop("The fitted penalty must match the training design and its coefficient coordinates.")
  }
  span <- switch(design$spread, constant = 3L, varying = 6L, contextual = 12L)
  pos_span <- if (design$position == "linear") 1L else 3L
  groups <- lapply(seq_along(design$vars), function(j) {
    list(spread = train$dir[span * (j - 1L) + seq_len(span)],
         position = train$pos[pos_span * (j - 1L) + seq_len(pos_span)])
  })
  names(groups) <- design$vars
  if (fit$has_distance) {
    label <- if ("distance" %in% design$vars) "geographic_distance" else "distance"
    groups[[label]] <- list(distance = setdiff(train$dir, seq_len(ncol(design$dir)) + 1L))
  }
  list(X = X, y = y, weights = weights, test_X = test_X, test_y = test_y,
    test_weights = test_weights, bounded = train$dir, penalty = penalty, groups = groups,
    start = fit$coef, idx = design$idx, n = design$n, order = design$order)
}

.within_importance_plan <- function(groups, permutations, exact_max) {
  exact_max <- .set_design_integer(exact_max, "exact_max", 0L, 8L)
  groups <- lapply(groups, unname)
  blocks <- lengths(groups)
  group <- rep(seq_along(groups), blocks)
  cells <- unlist(groups, recursive = FALSE, use.names = FALSE)
  cell_ids <- split(seq_along(cells), group)
  group_count <- length(groups)
  exact <- group_count <= exact_max
  path <- function(ordering) {
    order <- unlist(lapply(cell_ids[ordering], function(i) i[sample.int(length(i))]),
      use.names = FALSE)
    list(before = lapply(seq_along(order), function(i) order[seq_len(i - 1L)]), enter = order)
  }
  if (!exact) {
    permutations <- .set_design_integer(permutations, "permutations", 2L)
    return(list(cells = cells, group = group, exact = FALSE,
      paths = lapply(seq_len(permutations), function(i) path(sample.int(group_count)))))
  }
  subset <- function(ids) {
    masks <- seq_len(2^length(ids)) - 1L
    lapply(masks, function(mask) ids[bitwAnd(mask, 2^(seq_along(ids) - 1L)) != 0L])
  }
  terms <- lapply(seq_along(cells), function(i) {
    same <- group[i]
    outer <- subset(setdiff(seq_along(groups), same))
    inner <- subset(setdiff(cell_ids[[same]], i))
    before <- unlist(lapply(outer, function(other) {
      ext <- unlist(cell_ids[other], use.names = FALSE)
      lapply(inner, function(partial) c(ext, partial))
    }), recursive = FALSE)
    outer_weight <- 1 / (group_count * choose(group_count - 1L, lengths(outer)))
    inner_weight <- 1 / (blocks[same] * choose(blocks[same] - 1L, lengths(inner)))
    list(before = before, weight = rep(outer_weight, each = length(inner)) *
      rep(inner_weight, times = length(outer)))
  })
  list(cells = cells, group = group, exact = TRUE, terms = terms)
}

.within_importance_game <- function(problem, permutations = 1000L, exact_max = 8L,
    min_D2 = 1e-4, fit_coalition = NULL, failure = "error", plan = NULL) {
  if (!is.numeric(min_D2) || length(min_D2) != 1L || !is.finite(min_D2) || min_D2 < 0) {
    stop("Supply a finite nonnegative minimum explained deviance.")
  }
  if (is.null(plan)) plan <- .within_importance_plan(problem$groups, permutations, exact_max)
  ww <- problem$test_weights / sum(problem$test_weights)
  null_mu <- stats::weighted.mean(problem$y, problem$weights)
  cache <- new.env(parent = emptyenv())
  errors <- list()
  if (is.null(fit_coalition)) {
    fit_coalition <- function(keep, start) {
      X <- problem$X[, keep, drop = FALSE]
      penalty <- problem$penalty[keep, keep, drop = FALSE]
      bounded <- which(keep %in% problem$bounded)
      solve <- function(initial) .variant_solution(X, problem$y, problem$weights,
        penalty, bounded, initial)
      f <- solve(start)
      if (!isTRUE(f$converged)) f <- solve(NULL)
      f
    }
  }
  key <- function(selected) if (length(selected)) paste(sort(selected), collapse = ",") else "empty"
  solve <- function(code, keep, start) {
    fitted <- tryCatch(fit_coalition(keep, start), error = function(e) list(message = conditionMessage(e)))
    if (isTRUE(fitted$converged) && is.numeric(fitted$coef) &&
        length(fitted$coef) == length(keep) && all(is.finite(fitted$coef))) return(fitted$coef)
    reason <- if (is.null(fitted$message)) "optimiser did not converge" else fitted$message
    errors[[code]] <<- list(columns = colnames(problem$X)[keep], reason = reason)
    if (failure == "error") stop("Coalition ", code, " failed: ", reason, call. = FALSE)
    NULL
  }
  intercept <- pmin(pmax(null_mu, 1e-12), 1 - 1e-12)
  beta_null <- solve("empty", 1L, stats::qlogis(intercept))
  null <- if (is.null(beta_null)) NA_real_ else
    .binom_dev(problem$test_y, rep(plogis(beta_null), length(ww)), ww)
  if (all(problem$test_y == problem$y[1L]) &&
      all(problem$y == problem$y[1L]) && !is.null(beta_null) &&
      abs(plogis(beta_null) - null_mu) < 1e-9) null <- 0
  cache$empty <- list(value = 0, keep = 1L, coef = beta_null)
  score <- function(selected, parent = integer()) {
    code <- key(selected)
    cached <- cache[[code]]
    if (!is.null(cached)) return(cached$value)
    keep <- sort(c(1L, unlist(plan$cells[selected], use.names = FALSE)))
    previous <- cache[[key(parent)]]
    start <- problem$start[keep]
    if (length(parent) && !is.null(previous$coef)) {
      start[] <- 0
      start[match(previous$keep, keep)] <- previous$coef
    }
    coefficient <- solve(code, keep, start)
    value <- if (is.null(coefficient)) NA_real_ else {
      mu <- plogis(drop(problem$test_X[, keep, drop = FALSE] %*% coefficient))
      1 - .binom_dev(problem$test_y, mu, ww) / null
    }
    cache[[code]] <- list(value = value, keep = keep, coef = coefficient)
    value
  }
  if (is.null(problem$start)) problem$start <- c(stats::qlogis(intercept),
    numeric(ncol(problem$X) - 1L))

  good_null <- is.finite(null) && null > 0
  cell_count <- length(plan$cells)
  margin <- function(before, enter) {
    first <- score(before)
    score(c(before, enter), before) - first
  }
  draws <- NULL
  if (!good_null) {
    estimate <- rep(NA_real_, cell_count)
  } else if (plan$exact) {
    estimate <- vapply(seq_along(plan$cells), function(j) {
      term <- plan$terms[[j]]
      sum(term$weight * vapply(term$before, margin, numeric(1L), enter = j))
    }, numeric(1L))
  } else {
    draws <- t(vapply(plan$paths, function(path) {
      out <- numeric(cell_count)
      out[path$enter] <- vapply(seq_along(path$enter), function(i) {
        margin(path$before[[i]], path$enter[i])
      }, numeric(1L))
      out
    }, numeric(cell_count)))
    dim(draws) <- c(length(plan$paths), cell_count)
    estimate <- colMeans(draws)
  }
  total <- if (good_null) score(seq_along(plan$cells)) else NA_real_
  if (length(errors) && failure == "warn") {
    warning(length(errors), " coalition fits failed. Affected allocations are missing.", call. = FALSE)
  }
  groups <- problem$groups
  predictor <- rep(names(groups), lengths(groups))
  component <- unlist(lapply(groups, names), use.names = FALSE)
  mc_se <- function(x) if (is.null(draws)) ifelse(is.finite(x), 0, NA_real_) else
    stats::sd(x) / sqrt(length(x))
  cells <- data.frame(predictor, component, absolute = estimate,
    relative = if (is.finite(total) && total > min_D2) estimate / total else NA_real_,
    mc_se = if (is.null(draws)) mc_se(estimate) else apply(draws, 2L, mc_se))
  allocation <- function(group) {
    index <- split(seq_along(group), factor(group, levels = unique(group)))
    absolute <- vapply(index, function(i) sum(estimate[i]), numeric(1L))
    se <- if (is.null(draws)) mc_se(absolute) else
      vapply(index, function(i) mc_se(rowSums(draws[, i, drop = FALSE])), numeric(1L))
    list(absolute = absolute, total = total, relative = if (is.finite(total) && total > min_D2)
      absolute / total else absolute * NA_real_, mc_se = se,
      permutations = if (plan$exact) 0L else length(plan$paths))
  }
  list(predictors = allocation(predictor), components = allocation(component),
    hierarchy = cells, diagnostics = list(exact = plan$exact, null_deviance = null,
      coalition_count = length(ls(cache)), failures = errors, valid_null = good_null))
}
