#' Assess design and precision at a single zeta order
#'
#' Supply caller-declared thresholds in the units of dissimilarity or
#' explained deviance. Independently repeated importance integrations
#' keep sites, model specification and ecological draw fixed. Widths
#' describe uncertainty across ecological draws. A successful flag
#' requires a method calibrated for the conditions of the analysis.
#' @param fit,design Matching fitted additive MS-GDUM or component model and its design.
#' @param y Dissimilarity for each set.
#' @param integration Sites-fixed repeat integrations by predictor matrix.
#' @param draws Ecological resampling draws by predictor matrix.
#' @param heldout_D2 Explained deviance in independent site holdouts.
#' @param sharing Optional logical indicators of shared composition.
#' @param sampling_weights Optional importance weights targeting the original
#'   set distribution. These are distinct from the model's working weights.
#' @param min_response_sd,min_heldout_D2,max_mc_se,max_width,max_condition
#'   Declared criteria for response variation, prediction, numerical SE,
#'   ecological 95 percent interval width and scaled information.
#' @param calibrated Whether the ecological uncertainty procedure has
#'   coverage evidence for the conditions of this fit.
#' @param component Response name for a component fit. Coherent occupancy and
#'   heterogeneity require information in both fitted stages.
#' @return Response, sharing, design and precision diagnostics with status.
#' @export
ugdm_within_order_diagnostic <- function(fit, design, y, integration, draws,
    heldout_D2, sharing = NULL, min_response_sd, min_heldout_D2, max_mc_se,
    max_width, max_condition, calibrated = FALSE, sampling_weights = NULL,
    component = c("total", "occupancy", "heterogeneity")) {
  compositional <- inherits(fit, "ugdm_component_fit")
  component <- match.arg(component)
  if (!compositional && (!inherits(fit, "ugdm_interaction_fit") ||
      !inherits(design, "ugdm_interaction_design") || design$spread != "constant")) {
    stop("Supply the additive extended MS-GDUM and matching design.")
  }
  if (length(y) != nrow(design$idx) || any(!is.finite(y) | y < 0 | y > 1) ||
      length(y) < 2L || (!is.null(sharing) && (!is.logical(sharing) ||
      length(sharing) != length(y) || anyNA(sharing)))) stop("Invalid response or sharing indicators.")
  limits <- c(min_response_sd, min_heldout_D2, max_mc_se, max_width, max_condition)
  if (length(limits) != 5L || any(!is.finite(limits)) ||
      any(limits[c(1L, 3L, 4L, 5L)] < 0) || length(calibrated) != 1L ||
      !is.logical(calibrated) || is.na(calibrated)) stop("Invalid diagnostic criteria.")
  if (!is.numeric(heldout_D2) || length(heldout_D2) != 1L) stop("Supply one held-out D2.")
  if (compositional && !identical(fit$sets, design$idx)) stop("Supply the original component design.")
  target_names <- if (compositional) names(.component_problem(design, fit$ridge, fit$interaction_shrink)$groups) else
    names(.within_importance_setup(fit, design, y)$groups)
  if (compositional && is.null(sampling_weights)) sampling_weights <- fit$weights
  weights <- .model_weights(length(y), sampling_weights = sampling_weights)
  weights <- weights / max(weights)
  weights <- weights / sum(weights)
  if (sum(weights > 0) < 2L) stop("At least two sets must have positive sampling weight.")
  weighted_sd <- function(values) {
    average <- sum(weights * values)
    sqrt(sum(weights * (values - average)^2) / (1 - sum(weights^2)))
  }
  contrast_values <- function(value) {
    if (is.null(value)) return(NULL)
    value <- as.matrix(value)
    if (!is.numeric(value) || !identical(colnames(value), target_names)) {
      stop("Supply draws in predictor order with matching column names.")
    }
    if (length(target_names) == 1L) return(value)
    pairs <- utils::combn(seq_along(target_names), 2L)
    differences <- value[, pairs[1L, ], drop = FALSE] - value[, pairs[2L, ], drop = FALSE]
    colnames(differences) <- paste(target_names[pairs[1L, ]], target_names[pairs[2L, ]], sep = " - ")
    cbind(value, differences)
  }
  numeric_draws <- contrast_values(integration)
  ecological_draws <- contrast_values(draws)
  mc_se <- if (!is.null(numeric_draws) && nrow(numeric_draws) >= 3L &&
      all(is.finite(numeric_draws))) apply(numeric_draws, 2L, stats::sd) else NA_real_
  width <- if (!is.null(ecological_draws) && nrow(ecological_draws) >= 20L &&
      all(is.finite(ecological_draws))) apply(ecological_draws, 2L, function(value)
        diff(stats::quantile(value, c(.025, .975), names = FALSE))) else NA_real_
  X <- if (compositional) .component_prediction_matrix(fit, design) else
    .interaction_layout(design, fit$has_distance)$X
  working_information <- function(mu, working) {
    Z <- X * sqrt(working * mu * (1 - mu))
    scale <- sqrt(colSums(Z^2))
    scale[scale == 0] <- 1
    eigenvalues <- eigen(crossprod(sweep(Z, 2L, scale, "/")), symmetric = TRUE,
                         only.values = TRUE)$values
    rank <- sum(eigenvalues > max(eigenvalues) * ncol(X) * .Machine$double.eps)
    c(rank = rank, condition = if (rank == ncol(X)) max(eigenvalues) / min(eigenvalues) else Inf)
  }
  information <- if (compositional) {
    stages <- if (component != "total" && fit$architecture == "coherent") c("total", "share") else component
    values <- vapply(stages, function(stage) {
      f <- fit$fits[[stage]]
      mu <- if (is.finite(f$constant_value)) rep(f$constant_value, nrow(X)) else plogis(drop(X %*% f$coef))
      working_information(mu, weights * if (stage == "share") fit$share_trials else 1)
    }, numeric(2L))
    c(rank = min(values[1L, ]), condition = max(values[2L, ]))
  } else working_information(ugdm_predict(fit, design), weights)
  rank <- information[["rank"]]
  condition <- information[["condition"]]
  sd <- weighted_sd(y)
  variation <- sd >= min_response_sd && sd > 0
  predict <- fit$converged && is.finite(heldout_D2) && heldout_D2 > min_heldout_D2
  supported <- variation && predict && calibrated && condition <= max_condition &&
    all(is.finite(mc_se) & mc_se <= max_mc_se) &&
    all(is.finite(width) & width <= max_width)
  status <- if (!variation) "insufficient_response_variation" else
    if (!predict) "prediction_unsupported" else
      if (supported) "predictor_inference_supported" else "prediction_only"
  list(status = status, response_variance = sd^2,
       sharing = if (is.null(sharing)) NA_real_ else sum(weights * sharing),
       predictor_sd = apply(X[, -1L, drop = FALSE], 2L, weighted_sd),
       effective_rank = rank, columns = ncol(X), condition = condition,
       sampling_ESS = 1 / sum(weights^2), mc_se = mc_se,
       interval_width95 = width, heldout_D2 = heldout_D2)
}
