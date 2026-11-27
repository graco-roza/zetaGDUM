#' Within-order predictive attribution for component MS-GDUM
#'
#' Refits both stages or all independent marginals for each coalition. The
#' scoring function evaluates the final prediction of each response, including
#' the products from the coherent model. Predictor groups use Shapley weights.
#' Within a group, spread and position use Owen weights. Explicit interaction
#' enters after both main effects and receives the conditional gain above the
#' predictor's additive restriction. Equal endpoints enforce that restriction.
#' This hierarchy avoids double counting. Negative held-out allocations are
#' possible. Relative allocations require positive explained deviance.
#'
#' @param fit,design Training component fit and design.
#' @param response Training component responses.
#' @param evaluation,evaluation_response Optional independent design and responses.
#' @param resolution Full hierarchy or predictor-only coalitions.
#' @param permutations,exact_max,seed,min_D2 As in [ugdm_within_order_importance()].
#' @return Predictor allocations, a predictor-term hierarchy and failed coalition
#'   identifiers. Standard errors describe permutation sampling only.
#' @export
ugdm_component_importance <- function(fit, design, response, evaluation = NULL,
    evaluation_response = NULL, resolution = c("hierarchy", "predictor"),
    permutations = 200L, exact_max = 4L, seed = NULL, min_D2 = 1e-4) {
  resolution <- match.arg(resolution)
  .component_match_response(design, response)
  .component_prediction_matrix(fit, design)
  if (is.null(evaluation)) {
    evaluation <- design
    evaluation_response <- response
  }
  .component_match_response(evaluation, evaluation_response)
  .component_prediction_matrix(fit, evaluation)
  if (!is.numeric(min_D2) || length(min_D2) != 1L || !is.finite(min_D2) || min_D2 < 0) {
    stop("Supply a nonnegative finite minimum explained deviance.")
  }
  if (!is.null(seed)) set.seed(seed)
  problem <- .component_problem(design, fit$ridge, fit$interaction_shrink)
  Z <- .interaction_layout(evaluation, !is.null(evaluation$geo))$X
  groups <- problem$groups
  if (resolution == "predictor") {
    groups <- lapply(groups, function(g) list(predictor = unlist(g, use.names = FALSE)))
  }
  plan <- .within_importance_plan(groups, permutations, exact_max, interaction_last = TRUE)
  w <- fit$weights
  means <- colSums(response * (w / sum(w)))
  components <- names(means)
  n <- nrow(Z)
  eval_w <- rep(1 / n, n)
  loss <- function(pred) vapply(components, function(p) {
    .binom_dev(evaluation_response[, p], pred[[p]], eval_w)
  }, numeric(1L))
  null <- loss(lapply(means, rep, times = n))
  zero <- stats::setNames(rep(NA_real_, length(null)), components)
  zero[null > 1e-10] <- 0
  cache <- new.env(parent = emptyenv())
  cache$empty <- zero
  failures <- character()
  score <- function(select) {
    if (!length(select)) return(zero)
    key <- paste(sort(select), collapse = ",")
    stored <- cache[[key]]
    if (!is.null(stored)) return(stored)
    if (length(select) == length(plan$cells)) {
      pred <- .component_values(fit$fits, Z %*% fit$transform, fit$architecture)
      ok <- fit$converged
    } else {
      keep <- sort(c(1L, unlist(plan$cells[select], use.names = FALSE)))
      A <- .component_transform(problem, keep)
      candidate <- .component_fits(problem$X %*% A, crossprod(A, problem$P %*% A),
        response, w, fit$architecture, which(keep %in% fit$direction))
      ok <- all(vapply(candidate, function(x) x$converged, logical(1L)))
      pred <- .component_values(candidate, Z %*% A, fit$architecture)
    }
    out <- zero
    if (ok) out[null > 1e-10] <- (1 - loss(pred) / null)[null > 1e-10] else {
      failures <<- c(failures, key)
      out[] <- NA_real_
    }
    cache[[key]] <- out
    out
  }
  margin <- function(before, enter) score(c(before, enter)) - score(before)
  cells <- seq_along(plan$cells)
  evaluate <- function(g) do.call(cbind, lapply(cells, g))
  draws <- NULL
  if (plan$exact) {
    absolute <- evaluate(function(i) {
      term <- plan$terms[[i]]
      drop(vapply(term$before, margin, zero, enter = i) %*% term$weight)
    })
  } else {
    draws <- lapply(plan$paths, function(p) {
      evaluate(function(i) margin(p$before[[match(i, p$enter)]], i))
    })
    absolute <- Reduce("+", draws) / length(draws)
  }
  total <- score(cells)
  col_labels <- rep(names(groups), lengths(groups))
  term_labels <- unlist(lapply(groups, names), use.names = FALSE)
  compile <- function(index, name, term = NULL) {
    result <- rowSums(absolute[, index, drop = FALSE])
    se <- if (is.null(draws)) ifelse(is.finite(result), 0, NA_real_) else
      apply(vapply(draws, function(x) rowSums(x[, index, drop = FALSE]), zero), 1L, stats::sd) /
        sqrt(length(draws))
    out <- data.frame(response = components, predictor = name, absolute = result,
      relative = ifelse(is.finite(total) & total > min_D2, result / total, NA_real_),
      D2 = unname(total), mc_se = se, row.names = NULL)
    if (!is.null(term)) out$term <- term
    out
  }
  list(predictors = do.call(rbind, lapply(names(groups), function(p) {
    compile(which(col_labels == p), p)
  })), hierarchy = do.call(rbind, lapply(cells, function(i) {
    compile(i, col_labels[i], term_labels[i])
  })), diagnostics = list(exact = plan$exact, failed_coalitions = failures,
    evaluated_coalitions = length(ls(cache)), null_deviance = null))
}
