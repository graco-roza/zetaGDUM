.context_penalty_basis <- function(basis) {
  bounds <- attr(basis, "Boundary.knots")
  cuts <- sort(unique(c(bounds, attr(basis, "knots"))))
  width <- diff(bounds)
  half <- diff(cuts) / 2
  middle <- utils::head(cuts, -1L) + half
  nodes <- rep(middle, 3L) + rep(c(-sqrt(3 / 5), 0, sqrt(3 / 5)), each = length(half)) * half
  weights <- rep(c(5 / 9, 8 / 9, 5 / 9), each = length(half)) * half
  deriv <- stats::predict(basis, nodes, derivs = 1L)
  values <- stats::predict(basis, nodes)
  slope <- 12 * drop(crossprod(values, weights * ((nodes - bounds[1L]) / width - .5))) / width
  list(derivative = crossprod(deriv * sqrt(weights)), null_slope = slope, width = width)
}

.context_geometry <- function(z, spline, design, mean_z, bounds, u) {
  width <- bounds[2L, ] - bounds[1L, ]
  per_predictor <- lapply(seq_len(ncol(z)), function(j) {
    at <- 3L * (j - 1L) + seq_len(3L)
    base <- design$Xdir[, at, drop = FALSE]
    set <- base * u[, j]
    if (design$k == 2L) {
      pair <- set
      context <- matrix(0, nrow(base), 3L, dimnames = dimnames(base))
    } else {
      midpoint <- outer(z[, j] / 2, z[, j] / 2, "+")
      bounded <- pmin(pmax((midpoint - bounds[1L, j]) / width[j], 0), 1)
      moment <- function(position) vapply(at, function(v) {
        spread <- abs(outer(spline[, v], spline[, v], "-"))
        meanpair_dist_cpp(position * spread, design$idx)
      }, numeric(nrow(base)))
      pair <- matrix(moment(bounded), nrow(base), 3L, dimnames = dimnames(base))
      pair <- pmin(pmax(pair, 0), base)
      context <- base * mean_z[, j] - matrix(moment(midpoint), nrow(base), 3L)
    }
    columns <- cbind(base - pair, pair, base - set, set)
    colnames(columns) <- paste0(rep(colnames(base), 4L),
                                rep(c("_pair_L", "_pair_U", "_set_L", "_set_U"), each = 3L))
    list(pair = pair, context = context, columns = columns, effective = set - pair)
  })
  assemble <- function(field) do.call(cbind, lapply(per_predictor, "[[", field))
  list(dir = assemble("columns"), pair = assemble("pair"),
       context = assemble("context"), effective = assemble("effective"))
}

#' Information in pair-midpoint and multisite contextual spread
#'
#' The diagnostic uses weighted singular-value decompositions on scaled
#' columns. Residualisation is diagnostic only. Fits keep the original
#' nonnegative endpoint constraints. For prediction outside the training
#' range, use the bounded context, which is the difference between the
#' bounded set and pair interpolators. The design also stores the exact
#' raw context in standardised predictor units.
#' @param design An extended design with contextual spread.
#' @param weights Finite nonnegative working information weights per set.
#' @param tolerance Relative singular-value threshold.
#' @return Ranks of nested weighted designs, canonical singular spectra,
#'   and fractions of pair and context column norms remaining after projection
#'   on their respective preceding model.
#' @export
ugdm_contextual_diagnostics <- function(design, weights = rep(1, nrow(design$idx)),
                                        tolerance = 1e-8) {
  if (!inherits(design, "ugdm_interaction_design") || design$spread != "contextual") {
    stop("Supply an interaction design with contextual spread.")
  }
  if (length(weights) != nrow(design$idx) || any(!is.finite(weights) | weights < 0) ||
      !any(weights > 0)) stop("Supply finite nonnegative weights with positive total.")
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1) {
    stop("Supply a singular-value threshold between zero and one.")
  }
  W <- sqrt(weights / sum(weights))
  main <- cbind(intercept = 1, design$base, design$geo, design$pos)
  pair <- design$pair_midpoint - design$base / 2
  context <- design$bounded_context
  decompose <- function(X) {
    a <- X * W
    scales <- sqrt(colSums(a^2))
    # Leave negligible raw columns small so standardising numerical
    # roundoff in an exactly symmetric position cannot create a direction.
    scales[scales < max(scales) * tolerance] <- max(scales)
    fit <- svd(sweep(a, 2L, scales, "/"), nv = 0L)
    use <- fit$d > max(fit$d) * tolerance
    list(Q = fit$u[, use, drop = FALSE], rank = sum(use), spectrum = fit$d)
  }
  residual <- function(X, Q) {
    x <- X * W
    denom <- sqrt(colSums(x^2))
    out <- sqrt(colSums((x - Q %*% crossprod(Q, x))^2))
    stats::setNames(ifelse(denom > tolerance, out / denom, NA_real_), colnames(X))
  }
  additive <- decompose(main)
  midpoint <- decompose(cbind(main, pair))
  set <- decompose(cbind(main, pair, context))
  list(rank = c(additive = additive$rank, midpoint = midpoint$rank, context = set$rank),
       rank_increment = c(midpoint = midpoint$rank - additive$rank,
                          context = set$rank - midpoint$rank),
       singular_values = list(additive = additive$spectrum, midpoint = midpoint$spectrum,
                              context = set$spectrum),
       pair_residual_fraction = residual(pair, additive$Q),
       context_residual_fraction = residual(context, midpoint$Q),
       context_rms = sqrt(colSums((design$contextual * W)^2)), tolerance = tolerance)
}

#' Repeated site partitions for adaptive MS-GDUM validation
#'
#' Without coordinates, partitions use independent site permutations.
#' With projected coordinates, each partition uses contiguous bands
#' perpendicular to a random spatial direction. Every site appears once
#' per partition and block sizes differ by at most one. Spatial bands
#' are a sensitivity design and do not imply independence between blocks.
#' @param n Number of observed sites.
#' @param blocks Number of validation blocks per partition.
#' @param repetitions Number of independent partition draws.
#' @param coords Optional numeric matrix of projected site coordinates,
#'   indexed by site.
#' @param seed Optional partition seed.
#' @return Integer matrix of site-by-partition assignments.
#' @export
ugdm_site_partitions <- function(n, blocks = 5L, repetitions = 10L,
                                coords = NULL, seed = NULL) {
  n <- .set_design_integer(n, "n", 2L)
  blocks <- .set_design_integer(blocks, "blocks", 2L, n)
  repetitions <- .set_design_integer(repetitions, "repetitions")
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    if (!is.numeric(coords) || nrow(coords) != n || ncol(coords) < 2L ||
        any(!is.finite(coords))) stop("Supply finite projected coordinates for every site.")
  }
  if (!is.null(seed)) set.seed(seed)
  groups <- rep(seq_len(blocks), length.out = n)
  groups <- sort(groups)
  out <- vapply(seq_len(repetitions), function(j) {
    index <- if (is.null(coords)) sample.int(n) else
      order(drop(coords %*% stats::rnorm(ncol(coords))))
    groups[order(index)]
  }, integer(n))
  colnames(out) <- paste0("partition", seq_len(repetitions))
  out
}
