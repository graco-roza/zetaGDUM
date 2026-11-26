.gmm_blocks <- function(data, design_matrix) {
  if (!length(data)) stop("Provide at least one matching order.")
  if (any(vapply(data, function(x) !is.null(x$fit) && !isTRUE(x$fit$free_sign), TRUE))) {
    stop("Supply set data for regular interior inference. Use project explicitly for a boundary sensitivity.")
  }
  p <- lapply(data, function(x) {
    a <- .prepare_score_pools(list(x))[[1L]]
    if (!is.null(design_matrix)) {
      X <- as.matrix(design_matrix(x))
      if (nrow(X) != length(a$y) || !ncol(X) || !is.numeric(X) || any(!is.finite(X)) ||
          is.null(colnames(X)) || anyNA(colnames(X)) || anyDuplicated(colnames(X))) {
        stop("Custom model matrices require finite values and distinct coefficient names.")
      }
      a$X <- X
      a$design_by_row <- t(X)
    }
    a
  })
  n <- p[[1L]]$n
  if (anyDuplicated(vapply(p, function(x) x$k, numeric(1L)))) stop("Supply each order once.")
  if (any(vapply(p, function(x) x$n != n, TRUE))) stop("Use the same sites in every order.")
  p
}

.gmm_inverse <- function(A, label, tolerance) {
  scale <- sqrt(diag(A))
  if (any(!is.finite(scale) | scale <= 0)) stop(label, " has zero or invalid marginal information.")
  correlations <- A / tcrossprod(scale)
  values <- eigen(correlations, symmetric = TRUE, only.values = TRUE)$values
  if (min(values) <= max(values) * tolerance) stop(label, " is numerically rank deficient.")
  list(inverse = chol2inv(chol(correlations)) / tcrossprod(scale),
       condition = max(values) / min(values))
}

.gmm_make <- function(p, index, penalty) {
  point <- lapply(seq_along(p), function(j) {
    x <- p[[j]]
    K <- penalty[[j]]
    at <- x$effective > 0
    if (length(unique(x$y[at])) < 2L) stop("Responses must vary within each order.")
    fit <- suppressWarnings(stats::glm.fit(x$X, x$y, weights = x$effective,
      family = stats::binomial(), control = stats::glm.control(epsilon = 1e-10, maxit = 100L)))
    if (!fit$converged || fit$rank != ncol(x$X) || any(!is.finite(fit$coefficients))) {
      stop("The unrestricted score root is nonregular.")
    }
    if (!any(K != 0)) return(fit$coefficients)
    if (any(x$X[, 1L] != 1)) stop("Penalised model matrices must begin with an intercept.")
    result <- .variant_solution(x$X, x$y, x$effective, K, integer(), fit$coefficients)
    if (!result$converged) stop("The unrestricted penalised score root did not converge.")
    result$coef
  })
  list(point = unlist(point, use.names = FALSE), evaluate = function(theta, gradient = FALSE) {
    values <- lapply(seq_along(p), function(j) {
      a <- p[[j]]
      mu <- stats::plogis(drop(a$X %*% theta[index[[j]]]))
      K <- sum(a$effective) * penalty[[j]]
      if (gradient) crossprod(a$X * sqrt(a$effective * mu * (1 - mu))) + K else
        drop(crossprod(a$X, a$effective * (a$y - mu)) - K %*% theta[index[[j]]])
    })
    if (!gradient) return(unlist(values, use.names = FALSE))
    out <- matrix(0, length(theta), length(theta))
    for (j in seq_along(p)) out[index[[j]], index[[j]]] <- values[[j]]
    out
  })
}

.gmm_project <- function(precision, bounded) {
  free <- setdiff(seq_len(nrow(precision)), bounded)
  scale <- sqrt(diag(precision))
  metric <- precision / tcrossprod(scale)
  joint <- metric[free, bounded, drop = FALSE]
  map <- solve(metric[free, free, drop = FALSE], joint)
  profile <- metric[bounded, bounded, drop = FALSE] - crossprod(joint, map)
  U <- chol(profile)
  function(x) {
    if (!any(x[bounded] < 0)) return(x)
    z <- x * scale
    b <- nnls::nnls(U, drop(U %*% z[bounded]))
    if (b$mode != 1L) stop("Experimental loss-metric projection failed.")
    z[free] <- z[free] - drop(map %*% (b$x - z[bounded]))
    z[bounded] <- b$x
    z / scale
  }
}

#' Frozen-covariance Bayesian estimating equations for MS-GDUM
#'
#' Use regular unrestricted score roots and a normal prior on the canonical
#' coefficients. The frozen moment covariance is estimated at the unrestricted
#' pilot root from aligned site pseudovalues. The nonlinear kernel evaluates
#' the same fixed-covariance quadratic moment loss at every draw.
#' @param data List, one entry per order, with a fixed design, fractional y, and
#'   optional trials and sampling_weights. Order and site indexing must agree
#'   with pilot. Position coefficients are in mean units.
#' @param pilot Independent larger set data of the same layout for the moment
#'   covariance. The default uses data and does not check integration precision.
#' @param site_scores Optional function taking the named joint pilot coefficient
#'   vector and returning sites by coefficients score pseudovalues in that order.
#'   Supply aligned, independently checked moments. By default, use centred
#'   deletion scores from pilot.
#' @param design_matrix Optional function of one data entry returning a
#'   custom, named, unpenalised model matrix. Use the same coefficient order
#'   in site_scores. Supply the same predictor basis in data and pilot.
#'   For example, an interior sensitivity analysis can use a nonlinear
#'   position matrix with position-independent spread.
#' @param moment_covariance Optional function of the joint pilot coefficients
#'   returning a fixed covariance matrix or a list with a covariance field and
#'   diagnostics. This replaces site_scores and must use the same units.
#' @param penalty Optional list of positive-semidefinite quadratic penalty
#'   matrices in average working-loss units, one per order in the fitted
#'   coordinates. This changes the estimand. Its derivative enters the mean
#'   score and information, while centred site scores remain stochastic only.
#' @param prior_sd Positive normal prior SD, scalar or one per coefficient.
#'   Use Inf for prior-free Gaussian coordinates. Specify the prior separately
#'   from the penalty defining the estimand.
#' @param integration Gaussian linearisation or nonlinear MCMC.
#' @param draws Gaussian draws, or retained iterations per MCMC chain.
#' @param warmup,chains Warmup per chain and number of independent MCMC chains.
#' @param project Experimental projection onto nonnegative environmental spread
#'   coefficients in the local GMM loss metric. The unprojected draws remain
#'   available. Neither method supplies calibrated constrained probabilities.
#' @param bounded Optional indices of coefficients to project. Required for
#'   a custom design_matrix and otherwise defaults to environmental spread.
#' @param tolerance Relative eigenvalue threshold for scaled covariance and
#'   derivative diagnostics. Singular systems stop without covariance inflation.
#' @param seed Optional random seed.
#' @return Joint draws, score roots, Gaussian and sandwich covariance, frozen
#'   moment covariance, and numerical diagnostics. The Gaussian approximation
#'   includes the specified prior. The sandwich covariance omits it.
#' @export
ugdm_bayes_gmm <- function(data, pilot = data, site_scores = NULL, prior_sd = 2,
    integration = c("gaussian", "mcmc"), draws = 1000L, warmup = 700L,
    chains = 4L, project = FALSE, tolerance = 1e-9, seed = NULL,
    design_matrix = NULL, bounded = NULL, moment_covariance = NULL, penalty = NULL) {
  integration <- match.arg(integration)
  if (!is.null(seed)) set.seed(seed)
  started <- proc.time()[["elapsed"]]
  draws <- .set_design_integer(draws, "draws", 20L)
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1) {
    stop("tolerance must be strictly between zero and one.")
  }
  if (!is.null(design_matrix) && !is.function(design_matrix)) stop("design_matrix must be a function.")
  if (project && !is.null(design_matrix) && is.null(bounded)) {
    stop("Supply bounded explicitly for a custom projection.")
  }
  a <- .gmm_blocks(data, design_matrix)
  b <- .gmm_blocks(pilot, design_matrix)
  same <- length(a) == length(b) && all(vapply(seq_along(a), function(j)
    a[[j]]$n == b[[j]]$n && a[[j]]$k == b[[j]]$k &&
      identical(colnames(a[[j]]$X), colnames(b[[j]]$X)) &&
      identical(data[[j]]$design$environment_basis, pilot[[j]]$design$environment_basis) &&
      identical(data[[j]]$design$dist_range, pilot[[j]]$design$dist_range), TRUE))
  if (!same) stop("Global and pilot pools must share the sites, orders and fixed bases.")
  widths <- vapply(a, function(x) ncol(x$X), 1L)
  index <- split(seq_len(sum(widths)), rep(seq_along(a), widths))
  n <- a[[1L]]$n
  d <- sum(widths)
  if (n <= d) stop("More sites than joint coefficients are needed for a full-rank covariance.")
  nm <- make.unique(unlist(lapply(a, function(x) paste0("order", x$k, ".", colnames(x$X)))))
  if (length(prior_sd) == 1L) prior_sd <- rep(prior_sd, d)
  if (length(prior_sd) != d || anyNA(prior_sd) || any(prior_sd <= 0) ||
      (integration == "mcmc" && any(!is.finite(prior_sd)))) stop("Invalid prior SD for this integration.")
  if (is.null(penalty)) penalty <- lapply(widths, function(j) matrix(0, j, j))
  if (!is.list(penalty) || length(penalty) != length(widths)) stop("Provide one penalty matrix per order.")
  penalty <- Map(function(K, j) {
    if (!is.matrix(K) || !is.numeric(K) || !identical(dim(K), c(j, j)) ||
        any(!is.finite(K)) || !isSymmetric(unname(K))) stop("Invalid quadratic penalty matrix.")
    ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -max(1, max(ev)) * 1e-10) stop("Penalty matrices must be positive-semidefinite.")
    K
  }, penalty, widths)
  main <- .gmm_make(a, index, penalty)
  ref <- if (identical(a, b)) main else .gmm_make(b, index, penalty)
  theta <- stats::setNames(main$point, nm)
  pilot_theta <- stats::setNames(ref$point, nm)
  covariance_diagnostics <- NULL
  if (!is.null(moment_covariance)) {
    if (!is.null(site_scores) || !is.function(moment_covariance)) {
      stop("Supply one covariance function or a site-score function.")
    }
    result <- moment_covariance(pilot_theta)
    omega <- if (is.list(result)) result$covariance else result
    if (is.list(result)) covariance_diagnostics <- result[setdiff(names(result), "covariance")]
    if (!is.matrix(omega) || !is.numeric(omega) || !identical(dim(omega), c(d, d)) ||
        any(!is.finite(omega)) || !isSymmetric(unname(omega))) {
      stop("Provide a finite symmetric joint moment covariance.")
    }
  } else {
    if (is.null(site_scores)) {
      J <- do.call(cbind, lapply(seq_along(b), function(j)
        .site_score(b[[j]], pilot_theta[index[[j]]], "centered")$pseudovalues))
    } else {
      if (!is.function(site_scores)) stop("site_scores must be a function.")
      J <- as.matrix(site_scores(pilot_theta))
    }
    if (!identical(dim(J), c(n, d)) || !is.numeric(J) || any(!is.finite(J))) {
      stop("Supply finite site pseudovalues for every site and joint score coordinate.")
    }
    omega <- stats::cov(J)
  }
  dimnames(omega) <- list(nm, nm)
  om <- .gmm_inverse(omega, "Frozen moment covariance", tolerance)
  H <- main$evaluate(theta, gradient = TRUE)
  derivative <- .gmm_inverse(H, "Global moment derivative", tolerance)
  loss <- n * crossprod(H, om$inverse %*% H)
  prior_precision <- 1 / prior_sd^2
  precision <- loss + diag(prior_precision, d)
  covariance <- .gmm_inverse(precision, "Gaussian posterior precision", tolerance)$inverse
  mean <- stats::setNames(drop(covariance %*% (loss %*% theta)), nm)
  U <- chol(covariance)
  W <- t(U)
  sandwich <- derivative$inverse %*% omega %*% derivative$inverse / n
  convert <- function(z) {
    x <- sweep(z %*% U, 2L, mean, "+")
    colnames(x) <- nm
    x
  }
  diag <- list()
  group <- rep(1L, draws)
  if (integration == "gaussian") {
    posterior <- convert(matrix(stats::rnorm(draws * d), draws, d))
  } else {
    warmup <- .set_design_integer(warmup, "warmup")
    chains <- .set_design_integer(chains, "chains", 2L)
    density <- function(z) {
      value <- mean + drop(W %*% z)
      G <- main$evaluate(value)
      list(value = -n / 2 * sum(G * (om$inverse %*% G)) -
        sum((value / prior_sd)^2) / 2, status = "valid")
    }
    target <- list(log_density = density, d = d, point = mean, W = W)
    mode <- .jel_mode(target)
    run <- .jel_mcmc(target, mode, draws, warmup, chains)
    posterior <- convert(do.call(rbind, run$samples))
    group <- rep(seq_len(chains), each = draws)
    diag <- c(list(acceptance = run$acceptance, mode_convergence = mode$convergence,
      mode = drop(convert(matrix(mode$center, 1L))), clipped = mode$clipped),
      .jel_chain_diagnostics(lapply(run$samples, convert)))
  }
  projected <- NULL
  if (project) {
    if (is.null(bounded)) bounded <- unlist(lapply(seq_along(a), function(j) {
      width <- ncol(data[[j]]$design$Xdir)
      if (width > 0L) index[[j]][1L + seq_len(width)] else integer()
    }), use.names = FALSE)
    if (!length(bounded) || !is.numeric(bounded) || any(!is.finite(bounded) |
        bounded != floor(bounded) | bounded < 1 | bounded > d) ||
        anyDuplicated(bounded) || length(bounded) == d) {
      stop("Supply valid distinct projection indices and at least one free coefficient.")
    }
    projection <- .gmm_project(loss, bounded)
    projected <- t(vapply(seq_len(nrow(posterior)), function(i)
      projection(posterior[i, ]), numeric(d)))
    colnames(projected) <- nm
  }
  c(list(draws = posterior, projected = projected, group = group, estimate = theta,
    pilot_estimate = pilot_theta, gaussian_mean = mean, gaussian_covariance = covariance,
    sandwich_covariance = sandwich, omega = omega, derivative = H, loss_precision = loss,
    covariance_diagnostics = covariance_diagnostics,
    covariance_condition = om$condition, derivative_condition = derivative$condition,
    prior_sd = prior_sd, integration = integration, seconds = proc.time()[["elapsed"]] - started), diag)
}

#' Tangent-cone sampling approximation for bound-constrained MS-GDUM
#'
#' Approximate the constrained M-estimator by projecting its unrestricted
#' root in the estimating-equation derivative metric. Project centred
#' sandwich Gaussian increments onto the tangent cone for the specified
#' active face. These are sampling increments for basic intervals. They
#' differ from projected posterior draws and have an explicitly chosen face.
#' @param fit Result of [ugdm_bayes_gmm()].
#' @param bounded Indices of nonnegative coefficients.
#' @param active Active face, a subset of bounded. By default uses the face of
#'   the projected root. Supply integer(0) for an interior comparator.
#' @param draws Number of Gaussian increments.
#' @param seed Optional random seed.
#' @return Projected estimate, active indices, and centred sampling increments.
#' @export
ugdm_gmm_tangent <- function(fit, bounded, active = NULL, draws = 1000L, seed = NULL) {
  theta <- fit$estimate
  H <- fit$derivative
  V <- fit$sandwich_covariance
  d <- length(theta)
  valid <- function(i) is.numeric(i) && !anyDuplicated(i) &&
    all(is.finite(i) & i == floor(i) & i >= 1 & i <= d)
  if (!length(bounded) || !valid(bounded) || !is.matrix(H) || !is.matrix(V)) {
    stop("Supply a GMM fit and valid nonnegative coefficient indices.")
  }
  if (!is.null(seed)) set.seed(seed)
  draws <- .set_design_integer(draws, "draws", 20L)
  estimate <- .gmm_project(H, bounded)(theta)
  if (is.null(active)) active <- bounded[estimate[bounded] == 0]
  if (!valid(active) || !all(active %in% bounded)) stop("active must be a subset of bounded.")
  increments <- matrix(stats::rnorm(draws * d), draws, d) %*% chol(V)
  if (length(active)) {
    project <- .gmm_project(H, active)
    increments <- t(vapply(seq_len(draws), function(i) project(increments[i, ]), numeric(d)))
  }
  colnames(increments) <- names(theta)
  list(estimate = stats::setNames(estimate, names(theta)), active = active, increments = increments)
}
