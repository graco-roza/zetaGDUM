#' Compare observed set configurations in a fitted MS-GDUM
#'
#' Reports changes on the link and dissimilarity scales, plus the actual
#' predictor changes and number of replaced sites. Position columns are
#' expressed on the fixed site scale as set means. Use near-zero position
#' changes to diagnose contrasts intended to vary spread at fixed position.
#' @param fit,design A fitted model and its matching design.
#' @param from,to Equal-length vectors of design row indices.
#' @return A list of effect changes and the underlying design differences.
#' @export
ugdm_configuration_effects <- function(fit, design, from, to) {
  if (!length(from) || length(from) != length(to) ||
      any(!is.finite(c(from, to)) | c(from, to) < 1 | c(from, to) > nrow(design$idx) |
          c(from, to) != floor(c(from, to)))) stop("Supply matched valid design row indices.")
  x <- .design_matrix(design, fit$vars, fit$has_distance)
  S <- x$Xdir[to, , drop = FALSE] - x$Xdir[from, , drop = FALSE]
  P <- x$Xpos[to, , drop = FALSE] - x$Xpos[from, , drop = FALSE]
  replace <- vapply(seq_along(from), function(i) {
    sum(!design$idx[to[i], ] %in% design$idx[from[i], ])
  }, integer(1))
  effects <- cbind(spread = drop(S %*% fit$beta), position = drop(P %*% fit$lambda))
  eta <- ugdm_predict(fit, design, "link")[from]
  end <- stats::plogis(eta + rowSums(effects))
  start <- stats::plogis(eta)
  sp <- .5 * (stats::plogis(eta + effects[, 1L]) - start +
                end - stats::plogis(eta + effects[, 2L]))
  list(effects = data.frame(from, to, sites_replaced = replace, spread_link = effects[, 1L],
    position_link = effects[, 2L], dissimilarity_from = start, dissimilarity_to = end,
    spread_response = sp, position_response = end - start - sp),
    spread_design_change = S, mean_position_change = P / design$k)
}

#' Pilot and confirm a combination budget for declared MS-GDUM contrasts
#'
#' The callback generates a fresh fit or jointly fitted set of orders.
#' Samplers are compared using the worst variance-times-cost among the
#' declared, tolerance-standardised contrasts. A noisy inverse-size trend
#' triggers extra pilot replication. Confirmation uses independent draws.
#' @param fit_once Function of budget and sampler returning a named numeric
#'   vector of predeclared contrasts. Include differences in the vector
#'   when using a paired design. Time all relevant work inside the callback.
#' @param samplers Names understood by the callback.
#' @param pilots At least three distinct positive combination counts.
#' @param tolerance Positive absolute numerical SD floor and fallback for
#'   every contrast, supplied as a named numeric vector.
#' @param ecological_se Optional named ecological standard errors.
#' @param calibrated Optional named logical flags for usable ecological SEs.
#' @param fraction Fraction of the calibrated ecological SE used when it
#'   exceeds the corresponding absolute floor.
#' @param nrep Number of independent fits per pilot count and confirmation.
#' @param max_budget Computational ceiling. Infeasible predictions are retained.
#' @param seed Optional random seed.
#' @return Sampler choice, pilot diagnostics, contrast-level predicted counts,
#'   and independent confirmation. Numerical confirmation is separate from
#'   ecological or statistical certification.
#' @export
ugdm_contrast_budget <- function(fit_once, samplers, pilots, tolerance,
                                    ecological_se = NULL, calibrated = NULL, fraction = .2,
                                    nrep = 8L, max_budget = 50000L, seed = NULL) {
  nrep <- .set_design_integer(nrep, "nrep", 4L)
  max_budget <- .set_design_integer(max_budget, "max_budget")
  if (!length(samplers) || anyNA(samplers) || anyDuplicated(samplers)) stop("Supply distinct samplers.")
  pilots <- sort(unique(vapply(pilots, .set_design_integer, integer(1), name = "pilot budget")))
  terms <- names(tolerance)
  if (length(pilots) < 3L || max(pilots) > max_budget || !length(terms) ||
      anyNA(terms) || any(!nzchar(terms)) || anyDuplicated(terms) ||
      any(!is.finite(tolerance) | tolerance <= 0) || length(fraction) != 1L ||
      !is.finite(fraction) || fraction <= 0) stop("Supply valid pilot counts and named tolerances.")
  if (!is.null(seed)) set.seed(seed)
  target <- tolerance
  if (!is.null(ecological_se) && !is.null(calibrated)) {
    se <- ecological_se[terms]
    good <- calibrated[terms] %in% TRUE & is.finite(se) & se > 0
    target[good] <- pmax(tolerance[good], fraction * se[good])
  }
  fit <- function(M, method, reps) {
    seconds <- numeric(reps)
    q <- vapply(seq_len(reps), function(j) {
      tick <- sum(proc.time()[1:2])
      v <- fit_once(M, method)
      seconds[j] <<- sum(proc.time()[1:2]) - tick
      if (!all(terms %in% names(v))) stop("fit_once must return every named contrast.")
      unname(v[terms])
    }, numeric(length(terms)))
    if (is.null(dim(q))) q <- matrix(q, nrow = length(terms))
    list(values = q, seconds = seconds)
  }
  slope <- function(s) {
    if (any(!is.finite(s) | s <= 0)) return(NA_real_)
    unname(stats::coef(stats::lm(log(s) ~ log(pilots)))[2L])
  }
  result <- lapply(samplers, function(method) {
    runs <- lapply(pilots, fit, method = method, reps = nrep)
    vars <- function(runs) vapply(runs, function(a) apply(a$values, 1L, stats::var),
                                  numeric(length(terms)))
    v <- matrix(vars(runs), length(terms))
    slopes <- apply(v, 1L, slope)
    noisy <- any(!is.finite(slopes) | abs(slopes + 1) > .6)
    if (noisy) {
      more <- lapply(pilots, fit, method = method, reps = nrep)
      runs <- Map(function(a, b) list(values = cbind(a$values, b$values),
                                      seconds = c(a$seconds, b$seconds)), runs, more)
      v <- matrix(vars(runs), length(terms))
      slopes <- apply(v, 1L, slope)
    }
    A <- apply(sweep(v, 2L, pilots, "*"), 1L, stats::median)
    cpu <- vapply(runs, function(a) max(stats::median(a$seconds), .001), numeric(1))
    cost <- stats::median(cpu / pilots)
    list(A = A, slope = slopes, cost = if (all(is.finite(A))) max(A / target^2) * cost else Inf,
      pilot = data.frame(sampler = method, budget = rep(pilots, each = length(terms)),
        effect = rep(terms, length(pilots)), MC_variance = c(v), CPU = rep(cpu, each = length(terms)),
        replicates = if (noisy) 2L * nrep else nrep))
  })
  efficiency <- vapply(result, function(x) x$cost, numeric(1))
  j <- which.min(efficiency)
  if (!length(j) || !is.finite(efficiency[j])) return(list(sampler = NA_character_,
    pilots = do.call(rbind, lapply(result, function(x) x$pilot)), reason = "pilot fits unavailable"))
  best <- result[[j]]
  adjust <- (nrep - 1) / stats::qchisq(.1, nrep - 1)
  predicted <- ceiling(best$A / target^2)
  proposed <- pmax(min(pilots), ceiling(predicted * adjust))
  tab <- data.frame(effect = terms, tolerance = unname(target), A = best$A,
    slope = best$slope, predicted = predicted, candidate = proposed,
    verified_at = NA_real_, confirmed_sd = NA_real_, confirmed_upper = NA_real_, passed = FALSE)
  for (M in sort(unique(proposed[is.finite(proposed) & proposed <= max_budget]))) {
    eval <- fit(M, samplers[j], nrep)
    at <- which(proposed == M)
    sd <- apply(eval$values, 1L, stats::sd)[at]
    tab$verified_at[at] <- M
    tab$confirmed_sd[at] <- sd
    tab$confirmed_upper[at] <- sd * sqrt(adjust)
    tab$passed[at] <- is.finite(sd) & sd * sqrt(adjust) <= target[at]
  }
  list(sampler = samplers[j], cost = stats::setNames(efficiency, samplers),
       pilots = do.call(rbind, lapply(result, function(x) x$pilot)), contrasts = tab)
}

#' Classify ecological, statistical and numerical support separately
#' @param order,label Order and effect identifiers.
#' @param ecological,estimable,numerical Logical vectors or missing values
#'   for separately declared ecological support, statistical precision,
#'   and independent numerical confirmation.
#' @return Per-contrast status and highest passing order per axis and effect.
#' @export
ugdm_order_status <- function(order, label, ecological, estimable, numerical) {
  n <- length(order)
  if (!n || anyNA(order) || anyNA(label) ||
      any(lengths(list(label, ecological, estimable, numerical)) != n) ||
      !all(vapply(list(ecological, estimable, numerical), is.logical, TRUE))) {
    stop("Supply aligned orders, effect identifiers and logical diagnostics.")
  }
  status <- ifelse(ecological %in% FALSE, "ecologically uninformative",
    ifelse(is.na(ecological), "ecological support unresolved",
      ifelse(!estimable %in% TRUE, "statistical precision unresolved",
        ifelse(!numerical %in% TRUE, "computationally unresolved", "supported"))))
  x <- data.frame(order, effect = label, ecological, estimable, numerical, status)
  hi <- function(k, q) if (any(q %in% TRUE)) max(k[q %in% TRUE]) else NA_real_
  top <- lapply(split(x, x$effect), function(a) {
    c(ecological = hi(a$order, a$ecological),
      statistical = hi(a$order, a$ecological & a$estimable),
      computational = hi(a$order, a$ecological & a$estimable & a$numerical))
  })
  list(profile = x, highest = do.call(rbind, top))
}
