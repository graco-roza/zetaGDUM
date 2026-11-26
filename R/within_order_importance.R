#' Within-order hierarchical explained-deviance attribution
#'
#' Each ecological predictor is a group of its directional and position
#' basis columns. Spread, position and distance are also reported as
#' separate component groupings. Coalitions are refitted with the same
#' quadratic penalty, coefficient bounds and combination corrections as
#' the supplied extended MS-GDUM design. Negative held-out gains remain.
#' @param fit Fitted [ugdm_fit_interactions()] model.
#' @param design Its training design.
#' @param y,trials,sampling_weights Training response and weights.
#' @param evaluation Optional matching held-out design.
#' @param evaluation_y,evaluation_trials,evaluation_sampling_weights Held-out
#'   response and weights, required with evaluation.
#' @param permutations Number of random predictor orderings above eight groups.
#' @param seed Optional permutation seed.
#' @return Absolute predictor, component and joint Shapley allocations;
#'   relative shares only when full explained deviance exceeds `min_D2`.
#' @param min_D2 Minimum positive total for relative shares.
#' @export
ugdm_within_order_importance <- function(fit, design, y, trials = rep(1, length(y)),
    sampling_weights = NULL, evaluation = NULL, evaluation_y = NULL,
    evaluation_trials = NULL, evaluation_sampling_weights = NULL,
    permutations = 1000L, seed = NULL, min_D2 = 1e-4) {
  if (!inherits(fit, "ugdm_interaction_fit") || !inherits(design, "ugdm_interaction_design")) {
    stop("Use a fitted extended MS-GDUM and matching design.")
  }
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(evaluation) && (is.null(evaluation_y) ||
      !inherits(evaluation, "ugdm_interaction_design"))) {
    stop("Supply a matching evaluation design and response.")
  }
  train <- .interaction_layout(design, fit$has_distance)
  X <- train$X
  weights <- .model_weights(length(y), trials, sampling_weights = sampling_weights)
  test <- if (is.null(evaluation)) train else .interaction_layout(evaluation, fit$has_distance)
  yy <- if (is.null(evaluation)) y else evaluation_y
  ww <- if (is.null(evaluation)) weights else
    .model_weights(length(yy), evaluation_trials,
                   sampling_weights = evaluation_sampling_weights)
  if (ncol(test$X) != ncol(X) || !identical(colnames(test$X), colnames(X)) ||
      length(yy) != nrow(test$X) || any(!is.finite(yy) | yy < 0 | yy > 1)) {
    stop("Evaluation design and response differ from the fitted layout.")
  }
  span <- if (design$spread == "constant") 3L else
    if (design$spread == "varying") 6L else 12L
  position_span <- if (design$position == "linear") 1L else 3L
  environmental <- seq_along(design$vars)
  directions <- lapply(environmental, function(j)
    1L + span * (j - 1L) + seq_len(span))
  positions <- lapply(environmental, function(j)
    train$pos[position_span * (j - 1L) + seq_len(position_span)])
  distance <- if (fit$has_distance) setdiff(train$dir,
    unlist(directions, use.names = FALSE)) else integer()
  predictors <- Map(c, directions, positions)
  names(predictors) <- design$vars
  if (length(distance)) predictors$distance <- distance
  components <- list(spread = unlist(directions, use.names = FALSE),
                     position = unlist(positions, use.names = FALSE))
  if (length(distance)) components$distance <- distance
  null_mu <- sum(weights * y) / sum(weights)
  null <- .binom_dev(yy, rep(null_mu, length(yy)), ww / sum(ww))
  score <- function(indices) {
    if (!length(indices)) return(0)
    keep <- sort(unique(c(1L, indices)))
    bounded <- which(keep %in% train$dir)
    f <- .variant_solution(X[, keep, drop = FALSE], y, weights,
                           fit$penalty[keep, keep, drop = FALSE], bounded)
    if (!f$converged || !is.finite(null) || null <= 0) return(NA_real_)
    mu <- plogis(drop(test$X[, keep, drop = FALSE] %*% f$coef))
    1 - .binom_dev(yy, mu, ww / sum(ww)) / null
  }
  allocate <- function(groups) {
    g <- length(groups)
    cache <- new.env(parent = emptyenv())
    coalition <- function(selected) {
      key <- if (length(selected)) paste(sort(selected), collapse = ",") else "empty"
      if (!exists(key, envir = cache, inherits = FALSE)) {
        assign(key, score(unlist(groups[selected], use.names = FALSE)), envir = cache)
      }
      get(key, envir = cache, inherits = FALSE)
    }
    nperm <- if (g <= 8L) factorial(g) else
      .set_design_integer(permutations, "permutations")
    contribution <- numeric(g)
    for (r in seq_len(nperm)) {
      order <- if (g <= 8L) .importance_permutation(g, r) else sample.int(g)
      selected <- integer()
      previous <- coalition(selected)
      for (j in order) {
        selected <- c(selected, j)
        current <- coalition(selected)
        contribution[j] <- contribution[j] + current - previous
        previous <- current
      }
    }
    out <- stats::setNames(contribution / nperm, names(groups))
    total <- coalition(seq_len(g))
    list(absolute = out, total = total,
         relative = if (is.finite(total) && total > min_D2) out / total else
           stats::setNames(rep(NA_real_, g), names(groups)),
         permutations = nperm)
  }
  list(predictors = allocate(predictors), components = allocate(components))
}

.importance_permutation <- function(g, r) {
  available <- seq_len(g)
  out <- integer(g)
  index <- r - 1L
  for (j in seq_len(g)) {
    divisor <- factorial(g - j)
    at <- index %/% divisor + 1L
    out[j] <- available[at]
    available <- available[-at]
    index <- index %% divisor
  }
  out
}
