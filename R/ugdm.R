#' Distance matrix from coordinates
#'
#' @param coords Two-column matrix or data frame (x/longitude, y/latitude).
#' @param lonlat If `TRUE` (default `FALSE`) coordinates are decimal degrees and the haversine
#'   distance in km is returned; otherwise Euclidean distance.
#' @return A symmetric sites x sites matrix.
#' @export
#' @examples
#' coords_dist(cbind(runif(5), runif(5)))
coords_dist <- function(coords, lonlat = FALSE) {
  coords <- as.matrix(coords)
  if (ncol(coords) != 2) stop("`coords` must have two columns.")
  if (!lonlat) return(as.matrix(stats::dist(coords)))
  rad <- pi / 180
  la <- coords[, 2] * rad
  lo <- coords[, 1] * rad
  a <- sin(outer(la, la, "-") / 2)^2 + outer(cos(la), cos(la)) * sin(outer(lo, lo, "-") / 2)^2
  2 * 6371 * asin(pmin(sqrt(a), 1))
}

.pair_distance_range <- function(dm) {
  range(stats::as.dist(dm))
}

#' Set-level design for the multi-site GDUM
#'
#' Extends the pairwise generalised dissimilarity uniqueness model (Hernandez-Carrasco et al.
#' 2026) to sets of `k` sites. For each predictor `v` the *directional* (spread) columns are the
#' mean pairwise absolute difference, over the sites of the set, of monotone I-spline bases of the
#' predictor (three columns, `splines2::isp` of degree 2 as in the pairwise model); for `k = 2`
#' this is exactly `|f(x_i) - f(x_j)|`. The *non-directional* (position) column is the sum over
#' the sites of the set of the standardised predictor (for `k = 2`, `x_i + x_j`), so that its
#' coefficient is a site-level effect with the same meaning at every order. Geographic or other
#' distance can be added as a further directional predictor through the mean pairwise distance
#' of the set.
#'
#' @param env Data frame or matrix of sites x predictors (numeric). Predictors are standardised but not transformed: transform strongly skewed predictors (for example proportions with many zeros) beforehand; a warning is given when |skewness| > 2.
#' @param idx Integer matrix of sets from [sample_sets()].
#' @param dist_matrix Optional sites x sites distance matrix (see [coords_dist()]).
#' @param dist_range Optional fixed lower and upper bounds for mean set distance. A
#'   reference design takes precedence. By default the bounds come from the design's sets.
#' @param ref Optional reference design. Its fitted environmental transformation
#'   and distance scaling are reused for prediction or comparison. References
#'   from older versions still supply the distance basis.
#' @param basis_sites Optional indices of sites used to estimate predictor
#'   centring, scaling and spline bounds. Defaults to every row of `env`.
#'   An environmental basis in `ref` takes precedence.
#' @return An object of class `ugdm_design`.
#' @export
#' @examples
#' sim <- simulate_metacommunity(n_sites = 40, n_species = 30, seed = 1)
#' idx <- sample_sets(40, 3, M = 500, seed = 1)
#' d <- ugdm_design(sim$env, idx, dist_matrix = coords_dist(sim$coords))
#' dim(d$Xdir)
ugdm_design <- function(env, idx, dist_matrix = NULL, ref = NULL, dist_range = NULL,
                        basis_sites = NULL) {
  env <- as.matrix(as.data.frame(env))
  if (!is.numeric(env) || any(!is.finite(env))) stop("`env` must contain finite numeric predictors.")
  idx <- as.matrix(idx)
  if (!nrow(idx)) stop("No sets have a defined response.")
  if (!is.numeric(idx) || ncol(idx) < 2L || anyNA(idx) ||
      any(idx < 1L | idx > nrow(env) | idx != floor(idx))) stop("Invalid site indices.")
  if (is.null(basis_sites)) basis_sites <- seq_len(nrow(env))
  transform <- ref$environment_basis
  if (is.null(transform)) {
    if (!is.numeric(basis_sites) || length(basis_sites) < 4L ||
        anyNA(basis_sites) || anyDuplicated(basis_sites) ||
        any(basis_sites != floor(basis_sites) | basis_sites < 1L | basis_sites > nrow(env))) {
      stop("`basis_sites` must contain at least four distinct, valid site indices.")
    }
    training <- scale(env[basis_sites, , drop = FALSE])
    if (any(!is.finite(training))) stop("Basis sites contain constant predictors.")
    if (qr(training)$rank < ncol(training)) {
      warning("Predictors are linearly dependent (for example proportions that sum to one); importance and position effects are not identifiable. Drop one predictor.")
    } else if (ncol(training) > 1L) {
      vif <- tryCatch(diag(solve(stats::cor(training))), error = function(e) rep(Inf, ncol(training)))
      if (any(vif > 10)) warning("High collinearity among predictors (VIF > 10 for: ",
                                 paste(colnames(training)[vif > 10], collapse = ", "),
                                 "); importance shares may be unstable.")
    }
    skew <- colMeans(training^3)
    if (any(abs(skew) > 2)) warning("Strongly skewed predictor(s) (|skewness| > 2): ",
                                     paste(colnames(training)[abs(skew) > 2], collapse = ", "),
                                     ". Transform them before fitting.")
    splines <- lapply(seq_len(ncol(training)), function(v) {
      z <- unname(training[, v])
      if (length(unique(z)) < 4L) stop("A predictor needs at least four distinct values among basis sites.")
      splines2::isp(z, df = 3, degree = 2)
    })
    transform <- list(center = attr(training, "scaled:center"),
                      scale = attr(training, "scaled:scale"), splines = splines)
  }
  vars <- colnames(env)
  if (!identical(vars, names(transform$center))) stop("Reference and design predictor names must agree.")
  Z <- scale(env, center = transform$center, scale = transform$scale)
  M <- nrow(idx)
  k <- ncol(idx)
  V <- ncol(Z)
  Xdir <- matrix(0, M, 3 * V)
  Xpos <- matrix(0, M, V)
  for (v in seq_len(V)) {
    basis <- transform$splines[[v]]
    bounds <- attr(basis, "Boundary.knots")
    z <- pmin(pmax(Z[, v], bounds[1]), bounds[2])
    B <- unclass(stats::predict(basis, z))[, 1:3, drop = FALSE]
    Xdir[, 3 * (v - 1) + 1:3] <- mean_pair_isp_cpp(B, z, idx)
    Xpos[, v] <- rowSums(matrix(Z[idx, v], M, k))
  }
  colnames(Xdir) <- paste0(rep(vars, each = 3), "_f", 1:3)
  colnames(Xpos) <- vars
  Xdist <- NULL
  dist_basis <- NULL
  if (!is.null(dist_matrix)) {
    dd <- meanpair_dist_cpp(unname(as.matrix(dist_matrix)), idx)
    if (is.null(ref) || is.null(ref$dist_basis)) {
      if (is.null(dist_range)) dist_range <- range(dd)
      if (!is.numeric(dist_range) || length(dist_range) != 2L ||
          any(!is.finite(dist_range)) || dist_range[1] > dist_range[2]) {
        stop("`dist_range` must contain two finite bounds in ascending order.")
      }
      d01 <- pmin(pmax((dd - dist_range[1]) / max(diff(dist_range), 1e-12), 0), 1)
      dist_basis <- splines2::isp(d01, df = 3, degree = 2, Boundary.knots = c(0, 1))
      Xdist <- unclass(dist_basis)[, 1:3, drop = FALSE]
    } else {
      dist_range <- ref$dist_range
      dist_basis <- ref$dist_basis
      d01 <- pmin(pmax((dd - dist_range[1]) / max(diff(dist_range), 1e-12), 0), 1)
      Xdist <- unclass(stats::predict(dist_basis, d01))[, 1:3, drop = FALSE]
    }
    colnames(Xdist) <- paste0("distance_f", 1:3)
  }
  structure(list(Xdir = Xdir, Xpos = Xpos, Xdist = Xdist, vars = vars, idx = idx, n = nrow(env), k = k, dist_range = dist_range, dist_basis = dist_basis,
                 environment_basis = transform), class = "ugdm_design")
}

#' Response for the multi-site GDUM
#'
#' Dissimilarity-type response in `[0, 1]` and binomial trials for a set of sites.
#'
#' @param x Sites x species matrix.
#' @param idx Integer matrix of sets.
#' @param component One of `"turnover"` (Simpson-type: `1 - shared / min richness`), `"total"`
#'   (Sorensen-type: `1 - shared / mean richness`), `"nestedness"` (nestedness-resultant),
#'   or, for abundance data, `"relative"` (relative-abundance dissimilarity),
#'   `"balanced"`, `"bray_curtis"`, `"gradient"` (the last three use raw counts).
#' @param abundance_trials Constant working weight for raw-count abundance
#'   responses (default 1). It is not an observed binomial trial count. A
#'   constant multiplier does not change fitted coefficients or D2.
#' @param working_weights `"default"` uses the original richness weights for
#'   incidence and unit weights for relative abundance. `"constant"` uses unit
#'   weights for incidence and relative abundance. Raw-count weights are set
#'   by `abundance_trials` in either mode.
#' @return A list with `y` (response), `trials`, and `keep` (logical, sets with a defined response).
#' @export
#' @examples
#' sim <- simulate_metacommunity(n_sites = 40, n_species = 30, seed = 1)
#' idx <- sample_sets(40, 3, M = 200, seed = 1)
#' str(ugdm_response(sim$comm, idx, "turnover"))
ugdm_response <- function(x, idx, component = c("turnover", "total", "nestedness", "relative", "balanced", "bray_curtis", "gradient"), abundance_trials = 1, working_weights = c("default", "constant")) {
  component <- match.arg(component)
  working_weights <- match.arg(working_weights)
  if (component == "relative") {
    x <- .as_abundance(x)
    idx <- .species_indices(idx, nrow(x))
    totals <- rowSums(x)
    keep <- rowSums(matrix(totals[idx] > 0, nrow(idx))) == ncol(idx)
    q <- x / ifelse(totals > 0, totals, 1)
    y <- rep(NA_real_, nrow(idx))
    if (any(keep)) {
      y[keep] <- pmax(0, 1 - min_abundance_cpp(q, idx[keep, , drop = FALSE]))
    }
    return(list(y = y, trials = rep(1, length(y)), keep = keep))
  }
  z <- zeta_sets(x, idx, abundance = component %in% c("balanced", "bray_curtis", "gradient"))
  if (component == "turnover") { y <- 1 - z$turnover; tr <- z$richness_min }
  if (component == "total") { y <- 1 - z$total; tr <- z$richness_mean }
  if (component == "nestedness") { y <- z$nestedness; tr <- z$richness_min }
  if (component == "balanced") { y <- 1 - z$balanced; tr <- rep(abundance_trials, nrow(z)) }
  if (component == "bray_curtis") { y <- 1 - z$bray_curtis; tr <- rep(abundance_trials, nrow(z)) }
  if (component == "gradient") { y <- z$gradient; tr <- rep(abundance_trials, nrow(z)) }
  y <- pmin(pmax(y, 0), 1)
  keep <- is.finite(y) & tr > 0
  if (working_weights == "constant" && component %in% c("turnover", "total", "nestedness")) {
    tr <- rep(1, length(y))
  }
  list(y = y, trials = pmax(tr, 1e-8), keep = keep)
}

.model_weights <- function(n, trials = NULL, weights = NULL, sampling_weights = NULL) {
  if (n == 0L) stop("No sets have a defined response.")
  for (w in list(trials, weights, sampling_weights)) {
    if (!is.null(w) && (!is.numeric(w) || length(w) != n || any(!is.finite(w) | w < 0))) {
      stop("Trial, set and sampling weights must be finite, non-negative and match the sets.")
    }
  }
  effective <- rep(1, n)
  for (w in list(trials, weights, sampling_weights)) {
    if (!is.null(w)) effective <- effective * w
  }
  if (any(!is.finite(effective)) || !any(effective > 0)) {
    stop("Combined model weights must be finite with at least one positive weight.")
  }
  effective
}

.binom_dev <- function(z, mu, w) {
  mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
  t1 <- t2 <- numeric(length(z))
  pos <- z > 0
  neg <- z < 1
  t1[pos] <- z[pos] * log(z[pos] / mu[pos])
  t2[neg] <- (1 - z[neg]) * log((1 - z[neg]) / (1 - mu[neg]))
  2 * sum(w * (t1 + t2))
}

# weighted binomial (logit) IRLS with all slopes <= 0 (callers flip signs); intercept free and profiled out
.fit_logit_nnls <- function(X, z, w, maxit = 80, tol = 1e-9) {
  w <- w / max(w)
  w <- w / mean(w)
  p <- ncol(X)
  observed <- z[w > 0]
  if (all(observed == observed[1])) {
    mu <- min(max(observed[1], 1e-12), 1 - 1e-12)
    dev <- .binom_dev(z, rep(mu, length(z)), w)
    return(list(intercept = stats::qlogis(mu), coef = numeric(p), dev = dev,
                null_dev = dev, D2 = NA_real_, converged = TRUE))
  }
  wm <- min(max(sum(w * z) / sum(w), 1e-6), 1 - 1e-6)
  null_dev <- .binom_dev(z, rep(wm, length(z)), w)
  a <- stats::qlogis(wm)
  b <- rep(0, p)
  eta <- rep(a, length(z))
  mu <- stats::plogis(eta)
  dev <- .binom_dev(z, mu, w)
  conv <- FALSE
  if (p == 0) return(list(intercept = a, coef = numeric(0), dev = dev, null_dev = null_dev, D2 = 0, converged = TRUE))
  for (it in seq_len(maxit)) {
    v <- pmax(mu * (1 - mu), 1e-9)
    wt <- w * v
    zz <- eta + (z - mu) / v
    sw <- sum(wt)
    xm <- colSums(wt * X) / sw
    zm <- sum(wt * zz) / sw
    A <- sweep(X, 2, xm) * sqrt(wt)
    G <- crossprod(A)
    cv <- crossprod(A, (zz - zm) * sqrt(wt))
    R <- tryCatch(chol(G + diag(1e-10 * mean(diag(G)), p)), error = function(e) NULL)
    if (is.null(R)) break
    nb <- nnls::nnnpls(R, backsolve(R, cv, transpose = TRUE), con = rep(-1, p))$x
    na <- zm - sum(xm * nb)
    step <- 1
    repeat {
      b2 <- b + step * (nb - b)
      a2 <- a + step * (na - a)
      e2 <- a2 + drop(X %*% b2)
      m2 <- stats::plogis(e2)
      d2 <- .binom_dev(z, m2, w)
      if ((is.finite(d2) && d2 <= dev + 1e-9) || it == 1 || step < 1e-4) break
      step <- step / 2
    }
    ddev <- abs(dev - d2) / (abs(d2) + 0.1)
    a <- a2; b <- b2; eta <- e2; mu <- m2; dev <- d2
    if (ddev < tol) { conv <- TRUE; break }
  }
  list(intercept = a, coef = b, dev = dev, null_dev = null_dev, D2 = 1 - dev / null_dev, converged = conv)
}

.design_matrix <- function(d, vars = NULL, distance = TRUE) {
  V <- length(d$vars)
  sel <- if (is.null(vars)) seq_len(V) else match(vars, d$vars)
  width <- if (is.null(d$spread_width)) 3L else d$spread_width
  dcols <- as.vector(sapply(sel, function(v) width * (v - 1) + seq_len(width)))
  Xdir <- d$Xdir[, dcols, drop = FALSE]
  if (distance && !is.null(d$Xdist)) Xdir <- cbind(Xdir, d$Xdist)
  list(Xdir = Xdir, Xpos = d$Xpos[, sel, drop = FALSE], sel = sel)
}

#' Fit the multi-site GDUM
#'
#' Binomial (logit) model of the set-level dissimilarity response with monotone non-decreasing
#' directional (spread) effects of every predictor (and optionally distance) and free-sign
#' non-directional (position) effects. The intercept is free.
#'
#' @param design An object from [ugdm_design()].
#' @param y Response in `[0, 1]` (see [ugdm_response()]).
#' @param trials Response-specific working weights (same length as `y`).
#' @param weights Optional additional set weights (e.g. from [set_weights()]).
#' @param sampling_weights Optional base-to-proposal corrections from
#'   [sample_sets_species()]. The fitting criterion uses the product of
#'   `trials`, `weights` and `sampling_weights`.
#' @param vars Optional character vector: fit only these predictors.
#' @param distance Include the distance term if the design has one.
#' @param free_sign If `TRUE`, the directional effects may be non-decreasing or non-increasing
#'   (each spline basis function enters with both signs), as is appropriate for responses without an
#'   expected direction along environmental gradients, such as nestedness and abundance gradients.
#'   If `FALSE` (default) the directional effects are non-negative, as in the pairwise GDUM.
#' @return An object of class `ugdm_fit` with `height` (summed spline coefficients per predictor),
#'   `lambda` (position effect per standardised unit of the summed predictor), `D2`, `deviance`.
#'   A constant response gives an intercept-only fit with undefined `D2`.
#' @export
#' @examples
#' sim <- simulate_metacommunity(n_sites = 50, n_species = 40, seed = 1)
#' idx <- sample_sets(50, 2)
#' d <- ugdm_design(sim$env, idx)
#' r <- ugdm_response(sim$comm, idx, "turnover")
#' f <- ugdm_fit(d, r$y[r$keep], r$trials[r$keep])
#' f$height
ugdm_fit <- function(design, y, trials, weights = NULL, vars = NULL, distance = TRUE, free_sign = FALSE, sampling_weights = NULL) {
  if (nrow(design$Xpos) != length(y) || length(y) != length(trials)) stop("`y`, `trials` and the design must refer to the same sets.")
  dm <- .design_matrix(design, vars, distance)
  w <- .model_weights(length(y), trials, weights, sampling_weights)
  X <- if (free_sign) cbind(-dm$Xdir, dm$Xdir, dm$Xpos, -dm$Xpos) else cbind(-dm$Xdir, dm$Xpos, -dm$Xpos)
  f <- .fit_logit_nnls(X, y, w)
  p <- ncol(dm$Xdir)
  q <- ncol(dm$Xpos)
  cf <- f$coef
  if (free_sign) {
    beta <- -cf[seq_len(p)] + cf[p + seq_len(p)]
    lambda <- cf[2 * p + seq_len(q)] - cf[2 * p + q + seq_len(q)]
  } else {
    beta <- -cf[seq_len(p)]
    lambda <- cf[p + seq_len(q)] - cf[p + q + seq_len(q)]
  }
  vn <- design$vars[dm$sel]
  width <- if (is.null(design$spread_width)) 3L else design$spread_width
  height <- sapply(seq_along(vn), function(j) sum(beta[width * (j - 1) + seq_len(width)]))
  names(height) <- vn
  names(lambda) <- vn
  dist_height <- if (distance && !is.null(design$Xdist)) sum(beta[p - 2:0]) else NA_real_
  structure(list(beta = beta, height = height, lambda = lambda, dist_height = dist_height, intercept = f$intercept,
                 deviance = f$dev, null_deviance = f$null_dev, D2 = f$D2, converged = f$converged, vars = vn,
                 has_distance = distance && !is.null(design$Xdist)), class = "ugdm_fit")
}

#' @export
print.ugdm_fit <- function(x, ...) {
  cat("Multi-site GDUM fit: D2 =", round(x$D2, 3), "\n")
  print(round(data.frame(height = x$height, lambda = x$lambda), 3))
  if (x$has_distance) cat("distance height:", round(x$dist_height, 3), "\n")
  invisible(x)
}

.wcov <- function(a, b, w) {
  w <- w / sum(w)
  sum(w * (a - sum(w * a)) * (b - sum(w * b)))
}

#' Importance of predictors in a fitted multi-site GDUM
#'
#' Variance partitioning of the linear predictor as in joint species distribution models
#' (HMSC): the contribution of each predictor (directional plus position terms) is the
#' covariance of its term with the full linear predictor divided by the variance of the full
#' linear predictor, so shares sum to one and covariance between correlated predictors is shared.
#'
#' @param fit A `ugdm_fit`.
#' @param design The design used to fit it.
#' @param weights Optional set weights used in the fit.
#' @param trials Response-specific working weights used in the fit.
#' @param sampling_weights Optional proposal corrections used in the fit, as in [ugdm_fit()].
#' @return Named numeric vector of shares (predictors and, if present, `distance`).
#' @export
ugdm_importance <- function(fit, design, weights = NULL, trials = NULL, sampling_weights = NULL) {
  dm <- .design_matrix(design, fit$vars, fit$has_distance)
  V <- length(fit$vars)
  w <- .model_weights(nrow(dm$Xpos), trials, weights, sampling_weights)
  width <- if (is.null(design$spread_width)) 3L else design$spread_width
  eta <- vapply(seq_len(V), function(j) {
    cols <- width * (j - 1) + seq_len(width)
    drop(dm$Xdir[, cols, drop = FALSE] %*% fit$beta[cols]) + dm$Xpos[, j] * fit$lambda[j]
  }, numeric(length(w)))
  dim(eta) <- c(length(w), V)
  if (fit$has_distance) eta <- cbind(eta, drop(dm$Xdir[, width * V + 1:3, drop = FALSE] %*% fit$beta[width * V + 1:3]))
  tot <- rowSums(eta)
  vt <- .wcov(tot, tot, w)
  out <- if (vt > 0) vapply(seq_len(ncol(eta)), function(j) .wcov(eta[, j], tot, w) / vt, numeric(1)) else rep(NA_real_, ncol(eta))
  names(out) <- c(fit$vars, if (fit$has_distance) "distance")
  out
}

#' Shapley decomposition of deviance explained among groups of predictors
#'
#' Fits the model to every subset of predictor groups and averages each group's marginal gain in
#' deviance explained (D2) over all orderings. Contributions sum exactly to the D2 of the full
#' model and shared variance between correlated groups is divided among them.
#'
#' @param design A `ugdm_design`.
#' @param y,trials,weights,free_sign,sampling_weights As in [ugdm_fit()].
#' @param groups Named list; each element is a character vector of predictor names (use
#'   `"distance"` for the distance term). Every group can contain directional and position terms.
#' @return A named numeric vector of contributions to D2 with attribute `D2` (full model).
#' @export
#' @examples
#' sim <- simulate_metacommunity(n_sites = 40, n_species = 40, seed = 1)
#' idx <- sample_sets(40, 2)
#' d <- ugdm_design(sim$env, idx, dist_matrix = coords_dist(sim$coords))
#' r <- ugdm_response(sim$comm, idx, "turnover")
#' ugdm_shapley(d, r$y, r$trials, groups = list(env1 = "env1", env2 = "env2", distance = "distance"))
ugdm_shapley <- function(design, y, trials, groups, weights = NULL, free_sign = FALSE, sampling_weights = NULL) {
  G <- length(groups)
  if (G < 2 || G > 8) stop("Between 2 and 8 groups are supported.")
  w <- .model_weights(length(y), trials, weights, sampling_weights)
  if (all(y[w > 0] == y[which(w > 0)[1]])) {
    return(structure(stats::setNames(rep(NA_real_, G), names(groups)), D2 = NA_real_))
  }
  d2 <- numeric(2^G)
  for (m in seq_len(2^G - 1)) {
    sel <- which(bitwAnd(m, 2^(0:(G - 1))) > 0)
    vs <- unlist(groups[sel])
    env_vars <- setdiff(vs, "distance")
    use_dist <- "distance" %in% vs
    Xdir <- Xpos <- matrix(0, length(y), 0)
    if (length(env_vars)) {
      dm <- .design_matrix(design, env_vars, FALSE)
      Xdir <- dm$Xdir
      Xpos <- dm$Xpos
    }
    if (use_dist) Xdir <- cbind(Xdir, design$Xdist)
    X <- if (free_sign) cbind(-Xdir, Xdir, Xpos, -Xpos) else cbind(-Xdir, Xpos, -Xpos)
    d2[m + 1] <- .fit_logit_nnls(X, y, w)$D2
  }
  phi <- numeric(G)
  for (i in seq_len(G)) for (m in 0:(2^G - 1)) if (bitwAnd(m, 2^(i - 1)) == 0) {
    s <- sum(bitwAnd(m, 2^(0:(G - 1))) > 0)
    phi[i] <- phi[i] + factorial(s) * factorial(G - s - 1) / factorial(G) * (d2[bitwOr(m, 2^(i - 1)) + 1] - d2[m + 1])
  }
  names(phi) <- names(groups)
  attr(phi, "D2") <- d2[2^G]
  phi
}

#' Bayesian-bootstrap site weights
#'
#' Dirichlet(1, ..., 1) weights over sites with mean one, one column per draw. Using the same
#' matrix at every order gives *paired* draws, so that importance or ranks can be compared across
#' orders draw by draw.
#'
#' @param n Number of sites.
#' @param ndraw Number of draws.
#' @param seed Optional random seed.
#' @return An `n x ndraw` matrix.
#' @export
#' @examples
#' w <- site_weights(10, 5, seed = 1)
#' colMeans(w)
site_weights <- function(n, ndraw, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  g <- matrix(stats::rgamma(n * ndraw, 1), n, ndraw)
  n * sweep(g, 2, colSums(g), "/")
}

#' Set weights from site weights
#'
#' The weight of a set is the product of the weights of its sites.
#'
#' @param wsite Numeric vector of site weights.
#' @param idx Integer matrix of sets.
#' @return Numeric vector with one weight per set.
#' @export
set_weights <- function(wsite, idx) {
  apply(matrix(wsite[idx], nrow(idx)), 1, prod)
}
