.precision_inputs <- function(comm, env, orders, M, nrep, coords = NULL, lonlat = FALSE) {
  comm <- .as_abundance(comm)
  env <- as.data.frame(env)
  if (nrow(comm) != nrow(env) || nrow(comm) < 4L) stop("Supply matching community and environment rows for at least four sites.")
  if (!length(orders) || anyNA(orders) || any(orders != floor(orders) | orders < 2 | orders > nrow(comm)) ||
      anyDuplicated(orders)) stop("Supply distinct integer orders between 2 and the number of sites.")
  if (!length(M) || any(!is.finite(M) | M < 1 | M != floor(M))) stop("Budgets must be positive integers.")
  if (length(nrep) != 1L || !is.finite(nrep) || nrep < 2 || nrep != floor(nrep)) {
    stop("At least two independent sampling replicates are required.")
  }
  if (!is.null(coords) && nrow(coords) != nrow(comm)) stop("Coordinates must match community rows.")
  dm <- if (!is.null(coords)) coords_dist(coords, lonlat) else NULL
  list(comm = comm, env = env, P = (comm > 0) * 1L,
       dm = dm, dist_range = if (!is.null(dm)) .pair_distance_range(dm) else NULL)
}

.precision_run <- function(input, k, M, nrep, alpha, component, working_weights,
                           free_sign, groups, seed) {
  scores <- vector("list", nrep)
  fits <- vector("list", nrep)
  seconds <- numeric(nrep)
  for (i in seq_len(nrep)) {
    t0 <- sum(proc.time()[1:2])
    s <- .inference_sets(input$P, k, M, alpha, seed = if (!is.null(seed)) seed + i)
    f <- .inference_fit(input$comm, input$env, s, component, working_weights, free_sign, groups,
                        dist_matrix = input$dm, dist_range = input$dist_range)
    seconds[i] <- sum(proc.time()[1:2]) - t0
    fits[[i]] <- f$effects
    info <- f$information
    d <- .sampling_diagnostics(s, f$keep, f$weights)
    scores[[i]] <- data.frame(ess = unname(d["ess"]), shared = mean(s$shared > 0),
      sharing_estimate = mean(s$sampling_weights * (s$shared > 0)), exact = s$exact,
      active = if (free_sign) 0L else info$active, alpha_effective = s$alpha,
      rank = info$rank, parameters = info$parameters, condition = info$condition,
      min_eigenvalue = info$min_eigenvalue, converged = f$fit$converged,
      D2 = f$fit$D2, seconds = seconds[i])
  }
  list(effects = do.call(rbind, fits), diagnostics = do.call(rbind, scores))
}

#' Compare mixture fractions using an independent fitting pilot
#'
#' The criterion is the sum of the sampling variances of the spread heights
#' and position coefficients, multiplied by median process CPU time per
#' fit. The criterion uses common spline bounds and the same one-SD position
#' scale for all candidates. Use a fresh seed for the final analysis.
#'
#' @param comm,env Sites by species and sites by predictor data.
#' @param coords,lonlat Optional coordinates and longitude/latitude flag as
#'   in [zeta_orders()]. Distance uses common site-pair bounds.
#' @param order Number of sites per set.
#' @param fractions Candidate targeted mixture probabilities.
#' @param M Number of sets per pilot fit.
#' @param nrep Independent pilot fits per candidate.
#' @param component,working_weights As in [ugdm_response()].
#' @param free_sign As in [ugdm_fit()].
#' @param seed Pilot seed.
#' @return Selected `alpha`, and candidate variance, CPU time, cost and
#'   fit diagnostics. Unidentified or nonconverged candidates are ineligible.
#'   The conventional sampler is used for exhaustive pairs or if the
#'   pilot cannot estimate a positive variance for any candidate.
#' @export
ugdm_sampling_pilot <- function(comm, env, order, fractions = c(0, .1, .25, .5, .75, .9),
                                M = 1000L, nrep = 8L, component = "turnover",
                                working_weights = "default", free_sign = FALSE, seed = 1L,
                                coords = NULL, lonlat = FALSE) {
  if (length(order) != 1L || length(M) != 1L) stop("Supply one order and one pilot budget.")
  p <- .precision_inputs(comm, env, order, M, nrep, coords, lonlat)
  if (!length(fractions) || anyDuplicated(fractions)) stop("Supply distinct candidate fractions.")
  for (alpha in fractions) .check_alpha(alpha)
  if (!0 %in% fractions) stop("Include zero as the conventional control.")
  if (choose(nrow(p$P), order) <= M && order == 2L) {
    return(list(alpha = 0, candidates = NULL, reason = "exhaustive pairs"))
  }
  invisible(.precision_run(p, order, min(M, 100L), 1L, 0, component, working_weights,
                           free_sign, NULL, seed))
  rows <- lapply(seq_along(fractions), function(j) {
    alpha <- fractions[j]
    a <- .precision_run(p, order, M, nrep, alpha, component, working_weights,
                        free_sign, NULL, seed + 10000L * j)
    d <- a$diagnostics
    effects <- a$effects[, grepl("^(spread|position)\\.", colnames(a$effects)), drop = FALSE]
    variance <- sum(apply(effects, 2L, stats::var))
    sec <- stats::median(d$seconds)
    usable <- all(d$converged & d$rank == d$parameters & is.finite(d$D2)) &&
      is.finite(variance) && variance > 0 && sec > 0
    data.frame(alpha = alpha, variance = variance, seconds = sec,
               cost = if (usable) variance * sec else Inf,
               median_ess = stats::median(d$ess), median_shared_fraction = stats::median(d$shared),
               active_fit_fraction = mean(d$active > 0), effective_alpha = stats::median(d$alpha_effective))
  })
  results <- do.call(rbind, rows)
  good <- is.finite(results$cost)
  if (all(results$effective_alpha == 0)) {
    return(list(alpha = 0, candidates = results, reason = "no effective targeted branch"))
  }
  list(alpha = if (any(good)) results$alpha[which.min(results$cost)] else 0,
       candidates = results, reason = if (any(good)) NULL else "no usable variance estimate")
}

.ecological_reference <- function(p, k, M, alpha, component, working_weights,
                                  free_sign, groups, exact_limit) {
  sets <- .inference_sets(p$P, k, M, alpha, exact_limit)
  f <- .inference_fit(p$comm, p$env, sets, component, working_weights, free_sign, groups,
                      dist_matrix = p$dm, dist_range = p$dist_range)
  out <- matrix(NA_real_, nrow(p$P), length(f$effects), dimnames = list(NULL, names(f$effects)))
  ess <- numeric(nrow(out))
  if (k < nrow(out)) {
    excluded <- split(rep(seq_len(nrow(f$design$idx)), k), as.vector(f$design$idx))
    for (i in seq_len(nrow(out))) {
      keep <- rep(TRUE, length(f$weights))
      keep[excluded[[as.character(i)]]] <- FALSE
      w <- f$weights[keep]
      ess[i] <- .kish(w)
      if (sum(w > 0) < f$information$parameters) next
      d <- .subset_order_design(f$design, keep)
      y <- f$y[keep]
      tr <- f$trials[keep]
      fit <- ugdm_fit(d, y, tr, free_sign = free_sign, sampling_weights = w)
      info <- .order_information(d, fit, .model_weights(length(y), tr, sampling_weights = w))
      if (info$rank < info$parameters || !fit$converged || !is.finite(fit$D2)) next
      imp <- ugdm_importance(fit, d, trials = tr, sampling_weights = w)
      shp <- if (!is.null(groups)) ugdm_shapley(d, y, tr, groups, free_sign = free_sign, sampling_weights = w)
      out[i, ] <- c(.order_effects(fit), imp, shp)
    }
  }
  list(delete = out, min_delete_ess = min(ess), exact = sets$exact)
}

#' Diagnose Monte Carlo precision across zeta orders
#'
#' Compares independently drawn set samples at each budget to paired
#' delete-one-site ecological scales. Adjacent-order changes use the same
#' site deletions at both orders. The smallest passing common budget uses
#' fresh independent replicates for confirmation. A pass measures relative
#' numerical precision and does not establish inferential calibration.
#'
#' @param comm,env Sites by species and sites by predictor data.
#' @param coords,lonlat Optional coordinates and longitude/latitude flag as
#'   in [zeta_orders()]. Distance uses common site-pair bounds.
#' @param orders Increasingly ordered site-set sizes.
#' @param budgets Numbers of combinations to evaluate per order.
#' @param nrep Number of independent draws of the whole set sample at each budget.
#' @param alpha A fixed targeted mixture probability, or `"adaptive"` to run
#'   [ugdm_sampling_pilot()] independently by order.
#' @param pilot_M,pilot_nrep Sets and replicates for each adaptive-pilot candidate.
#' @param reference_M Number of draws for the independent ecological reference.
#' @param exact_limit Enumerate the ecological reference when this many sets suffice.
#' @param mc_fraction Tolerated ratio of one-fit Monte Carlo SD to the paired
#'   deletion SD. Must be positive. This is a user accuracy target.
#' @param groups Optional groups for Shapley estimates as in [ugdm_shapley()].
#' @param component,working_weights As in [ugdm_response()].
#' @param free_sign As in [ugdm_fit()].
#' @param seed Random seed.
#' @param assess Optional names of effects to use in the stopping comparison.
#'   The default assesses every returned spread, position, importance, and
#'   requested Shapley effect, including adjacent-order changes.
#' @param stop_after_confirmation If `TRUE`, skip larger budgets once the
#'   stopping criterion is independently confirmed. The default evaluates
#'   the entire grid for comparison.
#' @return A list with effect-level precision curves, order-level diagnostics,
#'   pilot results, the ecological retention breakpoint and the first
#'   confirmed budget (or `NA`). No universal informative order is inferred.
#' @export
zeta_order_precision <- function(comm, env, orders, budgets = c(500, 1000, 2000, 4000, 10000, 20000, 50000),
                                 nrep = 8L, alpha = 0, pilot_M = 1000L, pilot_nrep = 8L,
                                 reference_M = 50000L, exact_limit = 50000L, mc_fraction = .2,
                                 groups = NULL, component = "turnover", working_weights = "default",
                                 free_sign = FALSE, seed = 1L, assess = NULL, stop_after_confirmation = FALSE,
                                 coords = NULL, lonlat = FALSE) {
  p <- .precision_inputs(comm, env, orders, c(budgets, reference_M), nrep, coords, lonlat)
  orders <- sort(orders)
  budgets <- sort(unique(budgets))
  if (length(mc_fraction) != 1L || !is.finite(mc_fraction) || mc_fraction <= 0) {
    stop("`mc_fraction` must be positive.")
  }
  if (length(exact_limit) != 1L || !is.finite(exact_limit) || exact_limit < 0) {
    stop("`exact_limit` must be non-negative.")
  }
  adaptive <- identical(alpha, "adaptive")
  if (!adaptive) .check_alpha(alpha)
  pilots <- vector("list", length(orders))
  fractions <- numeric(length(orders))
  reference <- vector("list", length(orders))
  for (j in seq_along(orders)) {
    k <- orders[j]
    if (adaptive) {
      pilots[[j]] <- ugdm_sampling_pilot(p$comm, p$env, k, M = pilot_M, nrep = pilot_nrep,
        component = component, working_weights = working_weights, free_sign = free_sign, seed = seed + j * 100000L, coords = coords, lonlat = lonlat)
    }
    fractions[j] <- if (adaptive) pilots[[j]]$alpha else alpha
    set.seed(seed + j * 100000L + 80000L)
    reference[[j]] <- .ecological_reference(p, k, reference_M, fractions[j], component,
                                           working_weights, free_sign, groups, exact_limit)
  }
  make_table <- function(runs, M, confirmed) {
    n <- length(orders)
    out <- vector("list", n)
    for (j in seq_len(n)) {
      results <- function(effects, delete, contrast, regular) {
        mc <- apply(effects, 2L, stats::sd)
        eco <- .jackknife_sd(delete)
        data.frame(contrast = contrast, effect = names(mc), budget = M, confirmation = confirmed,
                   mc_sd = unname(mc), ecological_sd = unname(eco), ratio = unname(mc / eco),
                   pass = regular & is.finite(mc) & is.finite(eco) & eco > 0 & mc <= mc_fraction * eco)
      }
      regular <- function(r) with(r$diagnostics, all(converged & rank == parameters & is.finite(D2)))
      a <- runs[[j]]
      out[[j]] <- results(a$effects, reference[[j]]$delete, paste0("k", orders[j]), regular(a))
      if (j > 1L) {
        b <- runs[[j - 1L]]
        out[[j]] <- rbind(out[[j]], results(a$effects - b$effects, reference[[j]]$delete - reference[[j - 1L]]$delete,
                             paste0("k", orders[j], " - k", orders[j - 1L]), regular(a) && regular(b)))
      }
    }
    do.call(rbind, out)
  }
  run_budget <- function(i, confirmation) lapply(seq_along(orders), function(j) {
    .precision_run(p, orders[j], budgets[i], nrep, fractions[j], component, working_weights,
                   free_sign, groups, seed + 2000000L + 10000L * i + 500L * j + 250L * confirmation)
  })
  passes <- function(table) {
    selected <- if (is.null(assess)) table else table[table$effect %in% assess, , drop = FALSE]
    if (!is.null(assess) && !all(assess %in% selected$effect)) stop("Unknown effect in `assess`.")
    all(selected$pass)
  }
  precision <- diagnostics <- list()
  first <- NA_integer_
  for (i in seq_along(budgets)) {
    runs <- run_budget(i, FALSE)
    table <- make_table(runs, budgets[i], FALSE)
    precision[[length(precision) + 1L]] <- table
    if (is.na(first) && passes(table)) {
      confirm <- make_table(run_budget(i, TRUE), budgets[i], TRUE)
      precision[[length(precision) + 1L]] <- confirm
      if (passes(confirm)) first <- budgets[i]
    }
    for (j in seq_along(orders)) {
      k <- orders[j]
      share <- .species_zeta_presence(p$comm, k)
      d <- runs[[j]]$diagnostics
      diagnostics[[length(diagnostics) + 1L]] <- cbind(
        order = k, k_over_n = k / nrow(p$comm), budget = budgets[i], alpha = fractions[j],
        replicate = seq_len(nrep), expected_shared_species = sum(share),
        sharing_lower = max(c(0, share)), sharing_upper = min(1, sum(share)),
        min_delete_ess = reference[[j]]$min_delete_ess, exact_reference = reference[[j]]$exact, d)
    }
    if (isTRUE(stop_after_confirmation) && !is.na(first)) break
  }
  list(precision = do.call(rbind, precision), diagnostics = do.call(rbind, diagnostics),
       confirmed_budget = first, pilots = stats::setNames(pilots, orders),
       retention_break = zeta_break(zeta_curve(p$comm, kmax = nrow(p$comm))),
       mc_fraction = mc_fraction)
}
