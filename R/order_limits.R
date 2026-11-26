.limit_draw <- function(p, k, M, sampler, alpha, exact_limit, mass_draws) {
  n <- nrow(p$P)
  exact <- lchoose(n, k) <= log(min(M, exact_limit)) + 1e-12
  if (exact) {
    idx <- t(utils::combn(n, k))
    sets <- list(idx = idx, sampling_weights = rep(1, nrow(idx)))
  } else {
    sets <- switch(sampler,
      species = .sample_sets_mixture(p$P, k, M, alpha),
      stratified = sample_sets_stratified(p$comm, k, M, mass_draws = max(mass_draws, M),
                                         exact_limit = exact_limit),
      conventional = list(idx = .iid_uniform_sets(n, k, M), sampling_weights = rep(1, M)),
      site = list(idx = sample_sets_balanced(n, k, M), sampling_weights = rep(1, M)))
  }
  if (is.null(sets$shared)) sets$shared <- shared_species_cpp(p$P, sets$idx)
  sets$exact <- exact
  sets
}

.limit_effect_spec <- function(env, design, predictor) {
  j <- match(predictor, design$vars)
  b <- design$environment_basis
  q <- (stats::quantile(env[[j]], c(.25, .75)) - b$center[j]) / b$scale[j]
  list(j = j, spread = .1 * colSums(diff(stats::predict(b$splines[[j]], q))),
       position = unname(diff(q)))
}

.limit_effects <- function(f, d, y, w, spec, free_sign, shapley) {
  effects <- c(spread = sum(f$beta[3L * (spec$j - 1L) + 1:3] * spec$spread),
               position = unname(f$lambda[spec$j] * spec$position))
  if (shapley) {
    x <- cbind(d$Xdir, d$Xdist)
    spread <- .fit_logit_nnls(if (free_sign) cbind(-x, x) else -x, y, w)
    position <- .fit_logit_nnls(cbind(d$Xpos, -d$Xpos), y, w)
    valid <- spread$converged && position$converged && all(is.finite(c(f$D2, spread$D2, position$D2))) &&
      max(spread$D2, position$D2) <= f$D2 + 1e-6
    phi <- if (valid) .5 * (f$D2 + position$D2 - spread$D2) else NA_real_
    effects <- c(effects, shapley_spread = f$D2 - phi, shapley_position = phi)
  }
  effects
}

.limit_fit <- function(p, sets, component, working_weights, free_sign, spec, shapley,
                       delete = FALSE) {
  y <- ugdm_response(p$comm, sets$idx, component, working_weights = working_weights)
  keep <- y$keep & sets$sampling_weights > 0
  if (!any(keep)) {
    terms <- c("spread", "position", if (shapley) c("shapley_spread", "shapley_position"))
    effects <- stats::setNames(rep(NA_real_, length(terms)), terms)
    de <- matrix(NA_real_, nrow(p$P), length(terms), dimnames = list(NULL, terms))
    return(list(effects = effects, delete = de, regular = FALSE,
      diagnostics = c(response_mean = NA_real_, response_sd = 0, working_response_sd = 0,
        fitted_sd = 0, response_defined = 0, shared_draws = sum(sets$shared > 0),
        sharing_estimate = mean(sets$sampling_weights * (sets$shared > 0)), ess = 0,
        minimum_delete_ess = 0, min_predictor_sd = NA_real_, joint_support = NA_real_,
        rank = 0, parameters = 1 + 4 * ncol(p$env) + 3 * !is.null(p$dm),
        condition = Inf, min_eigenvalue = NA_real_, active = NA_real_,
        converged = FALSE, exact = sets$exact, sampler_fallback = !is.null(sets$fallback), D2 = NA_real_)))
  }
  w <- sets$sampling_weights[keep]
  y <- list(y = y$y[keep], trials = y$trials[keep])
  d <- ugdm_design(p$env, sets$idx[keep, , drop = FALSE], dist_matrix = p$dm,
                   dist_range = p$dist_range, ref = p$basis)
  f <- ugdm_fit(d, y$y, y$trials, free_sign = free_sign, sampling_weights = w)
  ww <- w * y$trials
  info <- .order_information(d, f, ww)
  v <- function(x, w) {
    w <- w / sum(w)
    sqrt(sum(w * (x - sum(w * x))^2))
  }
  X <- info$X[, -1L, drop = FALSE]
  p_w <- ww / sum(ww)
  center <- sweep(X, 2L, colSums(X * p_w))
  co <- crossprod(center * sqrt(p_w))
  scale <- sqrt(pmax(diag(co), 0))
  support <- if (all(scale > 0)) .minimum_eigenvalue(co / outer(scale, scale)) else 0
  effects <- .limit_effects(f, d, y$y, ww, spec, free_sign, shapley)
  omissions <- .inclusion_sums(d$idx, cbind(p_w, p_w^2), d$n)
  omit1 <- 1 - omissions[, 1L]
  omit2 <- sum(p_w^2) - omissions[, 2L]
  ess <- ifelse(omit1 > 0 & omit2 > 0 & d$k < d$n, omit1^2 / omit2, 0)
  regular <- f$converged && is.finite(f$D2) && info$rank == info$parameters
  de <- matrix(NA_real_, d$n, length(effects), dimnames = list(NULL, names(effects)))
  if (delete && regular && d$k < d$n) {
    rows <- split(rep(seq_len(nrow(d$idx)), d$k), as.vector(d$idx))
    for (i in seq_len(d$n)) {
      use <- rep(TRUE, length(ww))
      use[rows[[as.character(i)]]] <- FALSE
      if (sum(use) < info$parameters) next
      dd <- .subset_order_design(d, use)
      ff <- ugdm_fit(dd, y$y[use], y$trials[use], free_sign = free_sign, sampling_weights = w[use])
      ii <- .order_information(dd, ff, ww[use])
      if (ff$converged && ii$rank == ii$parameters && is.finite(ff$D2)) {
        de[i, ] <- .limit_effects(ff, dd, y$y[use], ww[use], spec, free_sign, shapley)
      }
    }
  }
  list(effects = effects, delete = de, regular = regular,
       diagnostics = c(response_mean = stats::weighted.mean(y$y, w),
         response_sd = v(y$y, w), working_response_sd = v(y$y, ww),
         fitted_sd = v(info$mu, w), response_defined = mean(keep),
         shared_draws = sum(sets$shared > 0),
         sharing_estimate = mean(sets$sampling_weights * (sets$shared > 0)),
         ess = .kish(w), minimum_delete_ess = min(ess),
         min_predictor_sd = min(scale), joint_support = max(0, support),
         rank = info$rank, parameters = info$parameters, condition = info$condition,
         min_eigenvalue = info$min_eigenvalue, active = if (free_sign) 0 else info$active,
         converged = f$converged, exact = sets$exact,
         sampler_fallback = !is.null(sets$fallback), D2 = f$D2))
}

.limit_eco <- function(a, b, noise_fraction) {
  e1 <- .jackknife_sd(a)
  e2 <- .jackknife_sd(b)
  e <- pmin(e1, e2)
  noise <- .jackknife_sd(a - b) / sqrt(2)
  ratio <- noise / e
  list(se = e, first = e1, second = e2, noise_ratio = ratio,
       reliable = is.finite(e) & e > 0 & is.finite(ratio) & ratio <= noise_fraction)
}

#' Profile MS-GDUM order feasibility and verify a sampling budget
#'
#' The ecological gate combines predeclared tolerances with two independently
#' sampled, matched site deletion analyses. Candidate budgets use inverse-size
#' pilot variance estimates followed by fresh confirmation draws. All fits
#' share a fixed environmental and geographic basis.
#'
#' @param comm,env Matching sites by species and sites by environmental predictor tables.
#' @param orders Site-set sizes to assess. Order differences use the previous
#'   supplied order, with the same deletions in both fits.
#' @param predictor Name of the environmental predictor for the declared
#'   spread and position contrasts. The full supplied model is fitted.
#' @param max_ecological_se Named positive tolerances for spread and position
#'   on the logit scale. Applied to both levels and changes.
#' @param min_response_sd Minimum meaningful SD of the original dissimilarity.
#' @param sampler Independent uniform, species mixture, disjoint site blocks,
#'   or zero/positive-sharing stratification.
#' @param alpha Targeted species mixture fraction.
#' @param pilot_budgets At least two distinct set counts for estimating the
#'   coefficient of the inverse sample size variance relationship.
#' @param reference_M Set count for each of two independent ecological references.
#' @param max_budget Maximum count for independent candidate confirmation.
#' @param nrep Independent draws per pilot and confirmation. At least four.
#' @param epsilon Maximum Monte Carlo SD as a fraction of ecological SE.
#' @param ecological_noise Maximum ratio of between-pool deletion noise to
#'   the smaller matched jackknife SE.
#' @param max_condition Maximum allowable working Hessian condition number.
#' @param min_delete_ess Minimum exclusion ESS using combined fitting weights.
#' @param min_sharing Minimum positive-sharing probability, certified by the
#'   analytical lower bound. Zero disables this optional ecological gate.
#' @param stability Maximum ratio between the largest and smallest estimated
#'   coefficients of inverse-size Monte Carlo variance.
#' @param mc_confidence One-sided normal-theory confidence level for the Monte
#'   Carlo SD used at confirmation. Set to zero to use the observed SD.
#' @param require_change Require a precision pass on the preceding order
#'   contrast as well as estimates at the current order.
#' @param shapley Compute spread and position deviance Shapley summaries.
#'   They are evaluated separately and never gate the budget.
#' @param exact_limit Maximum population size permitting enumeration.
#' @param mass_draws Minimum fresh stratum probability pilot draws per
#'   stratified sample. The pilot uses at least as many draws as the fit.
#' @param component,working_weights As in [ugdm_response()].
#' @param free_sign As in [ugdm_fit()].
#' @param dist_matrix Optional site-by-site geographic distance matrix.
#' @param seed Optional seed covering all independent set draws.
#' @return Exact ecological support and a grid of model diagnostics and
#'   decisions, with pilot curves, matched ecological precision and fresh
#'   confirmation results. The retention breakpoint and maximum passing
#'   tested order are returned separately. The latter can be missing.
#' @export
ugdm_order_limits <- function(comm, env, orders, predictor, max_ecological_se, min_response_sd,
                              sampler = c("conventional", "species", "site", "stratified"),
                              alpha = .5, pilot_budgets = c(1000L, 4000L), reference_M = 10000L,
                              max_budget = 50000L, nrep = 10L, epsilon = .25, ecological_noise = .35,
                              max_condition = 1e7, min_delete_ess = 100, min_sharing = 0,
                              stability = 4, mc_confidence = .9, require_change = TRUE,
                              shapley = FALSE, exact_limit = 50000L, mass_draws = 10000L,
                              component = "total", working_weights = "constant", free_sign = FALSE,
                              dist_matrix = NULL, seed = NULL) {
  sampler <- match.arg(sampler)
  p <- .precision_inputs(comm, env, orders, c(pilot_budgets, reference_M, max_budget), nrep)
  n <- nrow(p$P)
  orders <- as.integer(sort(orders))
  pilot_budgets <- sort(unique(pilot_budgets))
  effects <- c("spread", "position")
  if (length(predictor) != 1L || !is.character(predictor) || !predictor %in% names(p$env) ||
      !is.numeric(max_ecological_se) || !all(effects %in% names(max_ecological_se)) ||
      any(!is.finite(max_ecological_se[effects]) | max_ecological_se[effects] <= 0)) {
    stop("Supply one predictor and positive named spread and position ecological SE tolerances.")
  }
  max_ecological_se <- max_ecological_se[effects]
  scalar <- function(x, lower, upper, strict = FALSE) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x <= upper &&
      (if (strict) x > lower else x >= lower)
  }
  if (length(pilot_budgets) < 2L || any(pilot_budgets < 2) || nrep < 4L ||
      !scalar(reference_M, 2, .Machine$integer.max) ||
      !scalar(max_budget, max(pilot_budgets), .Machine$integer.max) ||
      !scalar(min_response_sd, 0, 1, TRUE) || !scalar(min_sharing, 0, 1) ||
      !scalar(epsilon, 0, 1, TRUE) || !scalar(ecological_noise, 0, 1, TRUE) ||
      !scalar(max_condition, 1, Inf) || !scalar(min_delete_ess, 1, Inf) ||
      !scalar(stability, 1, Inf) || !scalar(mc_confidence, 0, 1 - .Machine$double.eps)) {
    stop("Invalid budgets, precision tolerances or diagnostic limits.")
  }
  exact_limit <- .set_design_integer(exact_limit, "exact_limit", 0L)
  .check_alpha(alpha)
  if (!is.null(dist_matrix)) {
    if (!is.matrix(dist_matrix) || !identical(dim(dist_matrix), c(n, n)) ||
        any(!is.finite(dist_matrix)) || any(dist_matrix < 0) || !isSymmetric(dist_matrix)) {
      stop("Supply a finite symmetric site-by-site distance matrix.")
    }
    p$dm <- dist_matrix
    p$dist_range <- .pair_distance_range(dist_matrix)
  }
  if (!is.null(seed)) set.seed(seed)
  p$basis <- ugdm_design(p$env, cbind(seq_len(n - 1L), seq.int(2L, n)),
                         dist_matrix = p$dm, dist_range = p$dist_range)
  spec <- .limit_effect_spec(p$env, p$basis, predictor)
  support <- zeta_order_support(p$comm, orders, abundance = all(rowSums(p$comm) > 0))
  fit <- function(k, M, delete = FALSE) {
    t0 <- proc.time()[["elapsed"]]
    s <- .limit_draw(p, k, M, sampler, alpha, exact_limit, mass_draws)
    f <- .limit_fit(p, s, component, working_weights, free_sign, spec, shapley, delete)
    f$diagnostics <- c(f$diagnostics, elapsed = proc.time()[["elapsed"]] - t0)
    f
  }
  reference <- lapply(orders, function(k) {
    a <- fit(k, reference_M, TRUE)
    b <- if (isTRUE(as.logical(a$diagnostics[["exact"]]))) a else fit(k, reference_M, TRUE)
    list(a = a, b = b)
  })
  runs <- function(k, M) {
    fits <- lapply(seq_len(nrep), function(i) fit(k, M))
    list(effect = do.call(rbind, lapply(fits, function(f) f$effects)),
         valid = all(vapply(fits, function(f) {
           f$regular && f$diagnostics[["condition"]] <= max_condition
         }, TRUE)))
  }
  pilots <- lapply(pilot_budgets, function(M) lapply(orders, runs, M = M))
  confidence_factor <- if (mc_confidence > 0) {
    sqrt((nrep - 1) / stats::qchisq(1 - mc_confidence, nrep - 1))
  } else 1
  keys <- summary <- verification <- list()
  detail <- function(j, change) {
    out <- function(run) {
      a <- run[[j]]
      if (!change) return(a)
      b <- run[[j - 1L]]
      list(effect = a$effect - b$effect, valid = a$valid && b$valid)
    }
    a <- reference[[j]]
    da <- a$a$delete
    db <- a$b$delete
    if (change) {
      da <- da - reference[[j - 1L]]$a$delete
      db <- db - reference[[j - 1L]]$b$delete
    }
    point <- a$a$effects
    if (change) point <- point - reference[[j - 1L]]$a$effects
    eco <- .limit_eco(da, db, ecological_noise)
    width <- stats::qt(.975, n - 1L) * eco$first
    var <- do.call(rbind, lapply(pilots, function(run) apply(out(run)$effect, 2L, stats::var)))
    eligible <- vapply(pilots, function(run) out(run)$valid, TRUE)
    A <- var * pilot_budgets
    hi <- apply(A, 2L, max)
    lo <- apply(A, 2L, min)
    stable <- (hi == 0 & lo == 0) | (lo > 0 & hi / lo <= stability)
    cap <- round(if (change) max(choose(n, orders[c(j - 1L, j)])) else choose(n, orders[j]))
    exact_cap <- is.finite(cap) && cap <= min(exact_limit, max_budget)
    reliable <- eco$reliable & is.finite(hi) & hi > 0 & !is.na(stable) & stable & all(eligible)
    needed <- ifelse(reliable, ceiling(hi / (epsilon * eco$se)^2), NA_real_)
    planned <- ceiling(needed * confidence_factor^2)
    budget <- if (exact_cap) pmin(planned, cap, na.rm = TRUE) else planned
    budget <- ifelse(eco$reliable, budget, NA_real_)
    data <- data.frame(order = orders[j], preceding = if (change) orders[j - 1L] else NA_integer_,
        effect = names(hi), estimate = unname(point), jackknife_low = unname(point - width),
        jackknife_high = unname(point + width), ecological_se = unname(eco$se), jackknife_1 = unname(eco$first),
        jackknife_2 = unname(eco$second), deletion_noise_ratio = unname(eco$noise_ratio),
        ecology_reliable = unname(eco$reliable), A = unname(hi), scaling_ratio = unname(hi / lo),
        scaling_supported = !is.na(stable) & stable & hi > 0,
        predicted_budget = unname(needed), candidate_budget = unname(pmax(if (exact_cap) 2 else
          min(pilot_budgets), budget)))
    list(data = data, pilot = do.call(rbind, lapply(seq_along(pilot_budgets), function(i) {
      data.frame(order = orders[j], preceding = if (change) orders[j - 1L] else NA_integer_,
        effect = colnames(var), budget = pilot_budgets[i], mc_variance = unname(var[i, ]),
        fits_regular = eligible[i])
    })))
  }
  for (j in seq_along(orders)) {
    k <- orders[j]
    a <- reference[[j]]$a
    b <- reference[[j]]$b
    di <- a$diagnostics
    selected <- detail(j, FALSE)
    keys[[length(keys) + 1L]] <- selected$pilot
    if (j > 1L) {
      change <- detail(j, TRUE)
      keys[[length(keys) + 1L]] <- change$pilot
      selected$data <- rbind(selected$data, change$data)
    }
    tab <- selected$data
    check <- tab$effect %in% effects & (require_change | is.na(tab$preceding))
    tab$acceptable_ecological_se <- tab$ecology_reliable &
      tab$ecological_se <= unname(max_ecological_se[tab$effect])
    proposed <- max(tab$candidate_budget[check])
    required <- if (require_change && j > 1L) c(j - 1L, j) else j
    reference_fits <- unlist(reference[required], recursive = FALSE, use.names = FALSE)
    information_ok <- all(vapply(reference_fits, function(f) {
      f$regular && f$diagnostics[["condition"]] <= max_condition
    }, TRUE))
    deletion_ok <- all(vapply(reference_fits, function(f) {
      f$diagnostics[["minimum_delete_ess"]] >= min_delete_ess
    }, TRUE))
    gate <- c(response = min(di[["response_sd"]], b$diagnostics[["response_sd"]]) >= min_response_sd,
      sharing = support$sharing_lower[j] >= min_sharing,
      information = information_ok,
      deletion = deletion_ok,
      ecology = all(tab$acceptable_ecological_se[check]),
      feasible_budget = is.finite(proposed) && proposed <= max_budget)
    tab$confirmation_mc_sd <- tab$confirmation_mc_upper <- NA_real_
    tab$confirmation_pass <- NA
    if (all(gate)) {
      fresh <- runs(k, proposed)
      sd <- apply(fresh$effect, 2L, stats::sd)
      observed <- unname(sd[tab$effect])
      valid <- rep(fresh$valid, nrow(tab))
      has_change <- !is.na(tab$preceding)
      if (any(has_change)) {
        prev <- runs(orders[j - 1L], proposed)
        delta <- apply(fresh$effect - prev$effect, 2L, stats::sd)
        observed[has_change] <- unname(delta[tab$effect[has_change]])
        valid[has_change] <- valid[has_change] & prev$valid
      }
      upper <- observed * confidence_factor
      tab$confirmation_mc_sd <- observed
      tab$confirmation_mc_upper <- upper
      tab$confirmation_pass <- valid & is.finite(upper) & upper <= epsilon * tab$ecological_se
    }
    passed <- all(gate) && isTRUE(all(tab$confirmation_pass[check]))
    reason <- if (!all(gate)) paste(names(gate)[!gate], collapse = ", ") else
      if (!passed) "independent Monte Carlo confirmation" else "confirmed"
    summary[[j]] <- cbind(support[j, , drop = FALSE], as.data.frame(as.list(di)),
      candidate_budget = if (is.finite(proposed)) proposed else NA_real_,
      confirmed_budget = if (passed) proposed else NA_real_, pass = passed, reason = reason)
    verification[[j]] <- tab
  }
  profile <- do.call(rbind, summary)
  list(profile = profile, precision = do.call(rbind, verification), pilots = do.call(rbind, keys),
       effects = do.call(rbind, lapply(reference, function(r) r$a$effects)),
       deletes = lapply(reference, function(r) list(first = r$a$delete, second = r$b$delete)),
       retention_break = zeta_break(zeta_curve(p$comm, kmax = n)),
       highest_supported_order = if (any(profile$pass)) max(profile$order[profile$pass]) else NA_integer_)
}
