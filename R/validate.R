#' Predictions of a multi-site GDUM
#'
#' Fitted dissimilarity (or linear predictor) of the sets of a design, as in the observed versus predicted
#' plots of generalised dissimilarity models.
#'
#' @param fit An original or extended GDUM fit from \code{\link{ugdm_fit}} or
#'   \code{\link{ugdm_fit_interactions}}.
#' @param design A matching original or extended design. Reuse the training
#'   reference when constructing evaluation sets.
#' @param type `"response"` (default; expected dissimilarity) or `"link"` (linear predictor, the "predicted ecological distance").
#' @return A numeric vector with one value per set.
#' @export
ugdm_predict <- function(fit, design, type = c("response", "link")) {
  type <- match.arg(type)
  if (inherits(fit, "ugdm_interaction_fit")) return(.interaction_predict(fit, design, type))
  dm <- .design_matrix(design, fit$vars, fit$has_distance)
  eta <- fit$intercept + drop(dm$Xdir %*% fit$beta) + drop(dm$Xpos %*% fit$lambda)
  if (type == "link") eta else stats::plogis(eta)
}

#' Evaluate a multi-site GDUM on held-out sites
#'
#' Sets of sites share sites, so the fit to the sets used in a model overstates how well it predicts new data. This
#' function repeatedly splits the sites into training and test sites, fits the model to sets of training sites, and
#' predicts sets made only of test sites, as in the evaluation of multi-site generalised dissimilarity models
#' (a 70/30 split of the sites). Predictor standardisation and spline
#' bounds are estimated from training sites. Held-out spread bases are
#' clamped at training bounds, while position predictors retain their
#' training-scaled values.
#'
#' @param comm Sites x species matrix.
#' @param env Data frame of sites x predictors.
#' @param coords Optional site coordinates.
#' @param lonlat `TRUE` if `coords` are longitude and latitude.
#' @param order Order (number of sites per set).
#' @param component Response, see [ugdm_response()].
#' @param test_frac Fraction of sites held out.
#' @param M Number of sets for training, and for testing if `M_test` is omitted.
#' @param nrep Number of random splits.
#' @param seed Random seed.
#' @param keep_predictions If `TRUE`, the observed and predicted values of the test sets of the first split are returned.
#' @param sampling Training-set sampler, either `"conventional"` or
#'   `"species"` (see [sample_sets_species()]).
#' @param alpha Targeted mixture probability for species-informed sampling.
#' @param working_weights Response weighting mode from [ugdm_response()].
#' @param test_sampling Held-out sampler. Use `"conventional"` to compare
#'   fits on identical held-out sets for a given split and seed.
#' @param M_test Number of held-out sets, defaulting to `M`.
#' @return A list with `summary` (one row per split: `D2_train`, `D2_test` (deviance explained in the test sets relative
#'   to the mean of the training sets), `r2_train` and `r2_test` (squared Pearson correlation of observed and predicted
#'   dissimilarity), `slope_test`, the slope of observed on predicted
#'   dissimilarity through the origin, root mean squared prediction error and
#'   sampling diagnostics) and, if requested, `predictions`. Test deviance and
#'   RMSE use working and sampling weights. Correlation and slope retain their
#'   original set-level target and use sampling weights alone.
#' @export
#' @examples
#' sim <- simulate_metacommunity(n_sites = 80, n_species = 60, seed = 1)
#' v <- ugdm_validate(sim$comm, sim$env, sim$coords, order = 3, M = 2000, nrep = 3)
#' v$summary
ugdm_validate <- function(comm, env, coords = NULL, lonlat = FALSE, order = 2, component = "turnover",
                          test_frac = 0.3, M = 1e4, nrep = 20, seed = 1, keep_predictions = FALSE,
                          sampling = c("conventional", "species"), alpha = 0.5,
                          working_weights = c("default", "constant"), test_sampling = sampling, M_test = M) {
  sampling <- match.arg(sampling)
  test_sampling <- match.arg(test_sampling, c("conventional", "species"))
  working_weights <- match.arg(working_weights)
  .check_alpha(alpha)
  comm <- .as_abundance(comm)
  P <- (comm > 0) * 1L
  n <- nrow(comm)
  env <- as.data.frame(env)
  if (nrow(env) != n) stop("`comm` and `env` must have the same number of sites.")
  dm <- if (!is.null(coords)) coords_dist(coords, lonlat) else NULL
  rows <- vector("list", nrep)
  preds <- NULL
  for (r in seq_len(nrep)) {
    set.seed(seed + r)
    test <- sort(sample.int(n, round(test_frac * n)))
    train <- setdiff(seq_len(n), test)
    if (length(test) < order || length(train) < order) stop("Too few sites for this order and test fraction.")
    evaluate_sets <- function(sites, sampler, budget, draw_seed) {
      sets <- .model_sets(P[sites, , drop = FALSE], order, budget, sampler, alpha, seed = draw_seed)
      idx <- matrix(sites[sets$idx], ncol = order)
      response <- ugdm_response(comm, idx, component, working_weights = working_weights)
      keep <- response$keep
      if (!any(keep)) stop("No sets have a defined response in split ", r, ".")
      sw <- sets$sampling_weights[keep]
      list(idx = idx[keep, , drop = FALSE], y = response$y[keep], trials = response$trials[keep],
           sw = sw, diagnostics = .sampling_diagnostics(sets, keep, sw))
    }
    tr <- evaluate_sets(train, sampling, M, NULL)
    te <- evaluate_sets(test, test_sampling, M_test, seed + nrep + r)
    dr <- if (!is.null(dm)) .pair_distance_range(dm[train, train, drop = FALSE]) else NULL
    dtr <- suppressWarnings(ugdm_design(env, tr$idx, dm, dist_range = dr, basis_sites = train))
    dte <- suppressWarnings(ugdm_design(env, te$idx, dm, ref = dtr))
    fit <- ugdm_fit(dtr, tr$y, tr$trials, sampling_weights = tr$sw)
    ptr <- ugdm_predict(fit, dtr)
    pte <- ugdm_predict(fit, dte)
    wtr <- .model_weights(length(tr$y), tr$trials, sampling_weights = tr$sw)
    wte <- .model_weights(length(te$y), te$trials, sampling_weights = te$sw)
    p0 <- stats::weighted.mean(tr$y, wtr)
    null_dev <- .binom_dev(te$y, rep(p0, length(pte)), wte)
    d2_test <- if (null_dev > 0 && !all(te$y[wte > 0] == p0)) 1 - .binom_dev(te$y, pte, wte) / null_dev else NA_real_
    rows[[r]] <- data.frame(split = r, n_train_sets = length(tr$y), n_test_sets = length(te$y),
                            D2_train = fit$D2, D2_test = d2_test,
                            r2_train = .weighted_r2(tr$y, ptr, tr$sw), r2_test = .weighted_r2(te$y, pte, te$sw),
                            slope_test = sum(te$sw * te$y * pte) / sum(te$sw * pte^2),
                            rmse_test = sqrt(stats::weighted.mean((te$y - pte)^2, wte)),
                            ess_train = unname(tr$diagnostics["ess"]), ess_test = unname(te$diagnostics["ess"]),
                            shared_train = unname(tr$diagnostics["n_shared"]), shared_test = unname(te$diagnostics["n_shared"]))
    if (keep_predictions && r == 1) {
      preds <- data.frame(observed = te$y, predicted = pte, link = ugdm_predict(fit, dte, "link"),
                          trials = te$trials, sampling_weights = te$sw)
    }
  }
  list(summary = do.call(rbind, rows), predictions = preds)
}

.weighted_r2 <- function(y, predicted, w) {
  if (length(unique(y[w > 0])) < 2L || length(unique(predicted[w > 0])) < 2L) return(NA_real_)
  .wcov(y, predicted, w)^2 / (.wcov(y, y, w) * .wcov(predicted, predicted, w))
}
