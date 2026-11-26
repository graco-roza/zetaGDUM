.global_score_data <- function(x, env, paths, orders, type, weighting, controls,
                               dist_matrix, ref, design_matrix) {
  if (length(paths) != length(orders) ||
      !identical(vapply(paths, function(a) ncol(a$idx), integer(1), USE.NAMES = FALSE), orders)) {
    stop("Paths must follow all requested orders in order.")
  }
  path_data <- ugdm_path_data(x, env, paths, dist_matrix, ref)
  lapply(path_data, function(a) {
    d <- a$design
    if (is.null(design_matrix)) {
      X <- cbind(intercept = 1, d$Xdir, d$Xdist, d$Xpos / d$k)
      colnames(X) <- c("intercept", colnames(d$Xdir), colnames(d$Xdist), paste0(d$vars, ".position"))
    } else X <- as.matrix(design_matrix(a))
    y <- a[[type]]
    trial <- if (weighting == "richness") a$mean_richness else rep(1, length(y))
    w <- a$sampling_weights
    y[!is.finite(y) & trial == 0] <- 0
    if (!is.numeric(X) || !ncol(X) || nrow(X) != length(y) || any(!is.finite(X)) ||
        is.null(colnames(X)) || anyNA(colnames(X)) || anyDuplicated(colnames(X))) {
      stop("Supply a finite model matrix with named, distinct coefficients for each order.")
    }
    if (!length(y) || any(!is.finite(y) | y < -1e-12 | y > 1 + 1e-12)) {
      stop("Each set must have a defined fractional response or zero working weight.")
    }
    if (!is.numeric(w) || length(w) != length(y) || any(!is.finite(w) | w < 0) || !any(w > 0)) {
      stop("Supply finite nonnegative raw proposal ratios with positive total weight.")
    }
    Z <- matrix(0, length(y), length(controls), dimnames = list(NULL, controls))
    if ("inclusion" %in% controls) Z[, "inclusion"] <- 1
    if ("response" %in% controls) {
      Z[, "response"] <- if (type == "relative") y else a$mean_richness - a$shared
    }
    list(idx = d$idx, X = X, y = pmax(0, pmin(1, y)), trial = trial, w = w,
         Z = Z, order = d$k)
  })
}

.global_raw <- function(a, theta) a$X *
  (a$trial * (a$y - stats::plogis(drop(a$X %*% theta))))

.global_cv_pilot <- function(a, theta, n) {
  m <- nrow(a$X)
  T <- .global_raw(a, theta) * a$w
  Z <- a$Z * a$w
  q <- ncol(Z)
  d <- ncol(T)
  zproducts <- as.matrix(expand.grid(seq_len(q), seq_len(q)))
  ztproducts <- as.matrix(expand.grid(seq_len(q), seq_len(d)))
  values <- cbind(Z, T, Z[, zproducts[, 1L], drop = FALSE] * Z[, zproducts[, 2L], drop = FALSE],
                  Z[, ztproducts[, 1L], drop = FALSE] * T[, ztproducts[, 2L], drop = FALSE])
  sums <- .inclusion_sums(a$idx, values, n) / m
  nu <- seq_len(q)
  nt <- q + seq_len(d)
  zz <- q + d + seq_len(q * q)
  zt <- q + d + q * q + seq_len(q * d)
  lapply(seq_len(n), function(i) {
    V <- matrix(sums[i, zz], q) - tcrossprod(sums[i, nu])
    K <- matrix(sums[i, zt], q) - tcrossprod(sums[i, nu], sums[i, nt])
    # Standardisation uses pilot variance only. Constant and aliased
    # controls receive zero adjustment in their unidentified directions.
    scale <- sqrt(pmax(diag(V), 0))
    good <- which(scale > 0)
    coef <- matrix(0, q, d)
    if (length(good)) {
      s <- scale[good]
      ev <- eigen(V[good, good, drop = FALSE] / tcrossprod(s), symmetric = TRUE)
      keep <- which(ev$values > max(ev$values) * 1e-10)
      if (length(keep)) {
        U <- ev$vectors[, keep, drop = FALSE]
        coef[good, ] <- (U %*% (crossprod(U, K[good, , drop = FALSE] / s) / ev$values[keep])) / s
      }
    }
    coef
  })
}

#' Plan exact control variates for globally sampled site scores
#'
#' Learns site-specific adjustments for raw Horvitz--Thompson score
#' contributions from an independent pilot. The inclusion control averages
#' to k/n. The response-weighted inclusion control averages to k/n times
#' the exact conditional relative-abundance dissimilarity or the exact
#' mean-richness-minus-shared-richness incidence dissimilarity. The
#' latter works as a control even for the constant-weight incidence fit
#' because its mean is known without approximating a ratio.
#' @param x,env Community and environment matrices in the same site order.
#' @param orders Distinct orders from two through one less than the site count.
#' @param pilot_paths Independent list from \code{\link{sample_sets_coupled}}.
#' @param theta Pilot parameter in order and design-column order. Position
#'   uses mean units for the default design.
#' @param type Relative abundance or incidence.
#' @param weighting Constant working weights, or mean richness for incidence.
#' @param controls Any subset of inclusion and response. Use an empty vector
#'   to evaluate raw Horvitz--Thompson scores without a pilot.
#' @param dist_matrix,ref As in \code{\link{ugdm_path_data}}.
#' @param design_matrix Optional function from one order element of
#'   \code{\link{ugdm_path_data}} to a numerical model matrix.
#' @return A plan for evaluating independent production paths.
#' @export
ugdm_global_score_plan <- function(x, env, orders, pilot_paths = NULL, theta = NULL,
    type = c("relative", "incidence"), weighting = c("constant", "richness"),
    controls = c("inclusion", "response"), dist_matrix = NULL, ref = NULL, design_matrix = NULL) {
  type <- match.arg(type)
  weighting <- match.arg(weighting)
  if (length(controls)) {
    controls <- unique(match.arg(controls, c("inclusion", "response"), several.ok = TRUE))
  }
  x <- .conditional_community(x)
  n <- nrow(x)
  orders <- .conditioned_orders(orders, n - 1L)
  if (nrow(env) != n || (type == "relative" && (weighting == "richness" || any(rowSums(x) == 0)))) {
    stop("Relative abundance requires positive site totals and constant working weights, with aligned environment.")
  }
  if (type == "incidence" && weighting == "constant" && sum(rowSums(x) == 0) >= min(orders)) {
    stop("All-empty sets have undefined incidence response under constant working weights.")
  }
  if (!is.null(design_matrix) && !is.function(design_matrix)) stop("design_matrix must be a function.")
  if (is.null(ref)) ref <- ugdm_design(env, cbind(seq_len(n - 1L), seq.int(2L, n)),
    dist_matrix, dist_range = if (!is.null(dist_matrix)) .pair_distance_range(dist_matrix))
  get_data <- function(paths) .global_score_data(x, env, paths, orders, type, weighting, controls,
                                                dist_matrix, ref, design_matrix)
  expected <- coefficients <- layout <- NULL
  pilot_evaluations <- 0L
  if (length(controls)) {
    pilot <- get_data(pilot_paths)
    loop_theta <- .conditioned_theta(theta, pilot)
    expected <- lapply(orders, function(k) {
      mu <- matrix(1, n, length(controls), dimnames = list(NULL, controls))
      if ("response" %in% controls) {
        mu[, "response"] <- conditional_zeta_dissimilarity(x, k, type)$site
      }
      k / n * mu
    })
    coefficients <- Map(.global_cv_pilot, pilot, loop_theta, MoreArgs = list(n = n))
    layout <- lapply(pilot, function(a) colnames(a$X))
    pilot_evaluations <- sum(vapply(pilot, function(a) nrow(a$X), integer(1)))
  }
  list(n = n, data = get_data, expected = expected, coefficients = coefficients,
       layout = layout, pilot_evaluations = pilot_evaluations)
}

#' Globally sampled conditional site scores and pseudovalues
#'
#' Every sampled set contributes to every site it contains. Production
#' averages use the raw target/proposal ratios and the exact k/n
#' inclusion probability. Control coefficients remain fixed while theta
#' varies. Each site's conditional score receives equal weight in the
#' returned unconditional score. Coupled paths preserve differences
#' between orders.
#' @param plan Plan from \code{\link{ugdm_global_score_plan}}.
#' @param paths Production paths from \code{\link{sample_sets_coupled}}
#'   independently sampled from the pilot.
#' @param theta Numeric vector of coefficients across orders, or a list of
#'   per-order vectors, matching the plan.
#' @return Stacked conditional scores, unconditional score, and
#'   sites-by-coefficients pseudovalues with pilot and production costs.
#' @export
ugdm_global_site_scores <- function(plan, paths, theta) {
  blocks <- plan$data(paths)
  if (!is.null(plan$layout) && !identical(lapply(blocks, function(a) colnames(a$X)), plan$layout)) {
    stop("Pilot and production coefficient layouts must match.")
  }
  theta <- .conditioned_theta(theta, blocks)
  n <- plan$n
  result <- lapply(seq_along(blocks), function(j) {
    a <- blocks[[j]]
    m <- nrow(a$X)
    weighted_residual <- a$w * a$trial * (a$y - stats::plogis(drop(a$X %*% theta[[j]])))
    H <- score_inclusion_cpp(a$idx, a$X, weighted_residual, a$w, n)[seq_len(n), -1L, drop = FALSE] / m
    if (!is.null(plan$coefficients)) {
      delta <- .inclusion_sums(a$idx, a$Z * a$w, n) / m - plan$expected[[j]]
      adjustment <- vapply(seq_len(n), function(i)
        drop(crossprod(delta[i, ], plan$coefficients[[j]][[i]])), numeric(ncol(a$X)))
      H <- H - t(matrix(adjustment, ncol(a$X), n))
    }
    cond <- H * n / a$order
    global <- colMeans(cond)
    colnames(cond) <- names(global) <- paste0("order", a$order, ".", colnames(a$X))
    V <- sweep(cond, 2L, global) * a$order * (n - 1) / (n - a$order)
    list(score = global, conditional = cond, pseudovalues = sweep(V, 2L, global, "+"))
  })
  list(score = unlist(lapply(result, "[[", "score")),
       conditional = do.call(cbind, lapply(result, "[[", "conditional")),
       pseudovalues = do.call(cbind, lapply(result, "[[", "pseudovalues")),
       evaluations = sum(vapply(blocks, function(a) nrow(a$X), integer(1))),
       pilot_evaluations = plan$pilot_evaluations)
}
