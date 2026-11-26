.conditioned_orders <- function(orders, n) {
  if (!length(orders) || anyDuplicated(orders)) stop("Supply distinct zeta orders.")
  vapply(orders, .zeta_order_integer, integer(1), n = n, minimum = 2L, USE.NAMES = FALSE)
}

#' Site-conditioned uniform combination paths
#'
#' Each independent path begins with its focal site and a uniformly random
#' prefix of distinct other sites. Matching prefixes across orders retains
#' cross-order covariance. Every evaluation is counted, including the same
#' unordered combination appearing in more than one focal-site stratum.
#' @param n Number of sites.
#' @param orders Distinct zeta orders, at least two and at most n.
#' @param draws Positive paths per site, scalar or one per site.
#' @param coupled Share the remaining sites' permutation across orders.
#' @param seed Optional random seed.
#' @return A list in the layout accepted by [ugdm_path_data()], with focal
#'   site indices, per-site draw counts and total set evaluations as attributes.
#' @export
sample_sets_conditioned <- function(n, orders, draws = 100L, coupled = TRUE, seed = NULL) {
  n <- .set_design_integer(n, "n", 2L)
  orders <- .conditioned_orders(orders, n)
  if (length(draws) == 1L) draws <- rep(draws, n)
  if (length(draws) != n) stop("Supply a path count for each focal site or a scalar.")
  draws <- vapply(draws, .set_design_integer, integer(1), name = "draws", USE.NAMES = FALSE)
  if (!is.null(seed)) set.seed(seed)
  focal <- rep(seq_len(n), draws)
  paths <- function(k) {
    remaining <- if (k == 2L) matrix(sample.int(n - 1L, length(focal), replace = TRUE), ncol = 1L) else
      .path_permutations(n - 1L, k - 1L, length(focal), FALSE)
    cbind(focal, remaining + (remaining >= focal))
  }
  if (coupled) path <- paths(max(orders))
  out <- lapply(orders, function(k) {
    idx <- if (coupled) path[, seq_len(k), drop = FALSE] else paths(k)
    list(idx = idx, sampling_weights = rep(1, length(focal)), exact = FALSE)
  })
  structure(stats::setNames(out, orders), focal = focal, draws = draws,
            coupled = coupled, evaluations = length(focal) * length(orders))
}

.focal_predictor_moments <- function(site, k) {
  n <- nrow(site)
  pair <- (k - 1) / (n - 1)
  nonincident <- pair * (k - 2) / (n - 2)
  pos <- (site + pair * sweep(-site, 2L, colSums(site), "+")) / k
  dir <- vapply(seq_len(ncol(site)), function(j) {
    index <- order(site[, j])
    z <- site[index, j]
    sorted <- z * (2 * seq_len(n) - n) + sum(z) - 2 * cumsum(z)
    sums <- numeric(n)
    sums[index] <- sorted
    (pair * sums + nonincident * (sum(sums) / 2 - sums)) / choose(k, 2)
  }, numeric(n))
  list(position = pos, spread = dir)
}

.conditioned_data <- function(x, env, samples, distance, ref, type, weighting, controls, design_matrix) {
  focal <- attr(samples, "focal")
  n <- nrow(x)
  richness <- rowSums(x > 0)
  path <- ugdm_path_data(x, env, samples, distance, ref)
  feature <- path[[1L]]$design
  variables <- if ("predictors" %in% controls) .environment_for_moments(env, feature) else NULL
  orders <- lapply(path, function(a) {
    d <- a$design
    if (is.null(design_matrix)) {
      X <- cbind(intercept = 1, d$Xdir, d$Xdist, d$Xpos / d$k)
      colnames(X) <- c("intercept", colnames(d$Xdir), colnames(d$Xdist), paste0(d$vars, ".position"))
    } else X <- as.matrix(design_matrix(a))
    y <- a[[type]]
    w <- if (weighting == "richness") a$mean_richness else rep(1, length(y))
    # Empty sets under richness weighting have exactly zero score mass.
    y[!is.finite(y) & w == 0] <- 0
    if (!is.numeric(X) || nrow(X) != length(y) || !ncol(X) ||
        any(!is.finite(X)) || any(!is.finite(y) | y < -1e-12 | y > 1 + 1e-12)) {
      stop("Each sampled set needs a finite design and a defined fractional response.")
    }
    if (is.null(colnames(X))) colnames(X) <- paste0("coef", seq_len(ncol(X)))
    if (anyDuplicated(colnames(X))) stop("Use distinct design coefficient names per order.")
    Z <- matrix(numeric(), length(y), 0L)
    known <- matrix(numeric(), n, 0L)
    if (length(intersect(controls, c("response", "focal_incidence")))) {
      shared <- conditional_zeta_dissimilarity(x, d$k, type)
      if ("response" %in% controls) {
        Z <- cbind(Z, response = if (type == "relative") y else a$mean_richness - a$shared)
        known <- cbind(known, response = shared$site)
      }
      if ("focal_incidence" %in% controls) {
        r <- richness[focal]
        Z <- cbind(Z, focal_incidence = ifelse(r > 0, 1 - a$shared / pmax(r, 1), 0))
        known <- cbind(known, focal_incidence =
          ifelse(richness > 0, 1 - shared$shared$site / pmax(richness, 1), 0))
      }
    }
    if (!is.null(variables)) {
      moments <- .focal_predictor_moments(cbind(variables$spline, variables$z), d$k)
      spread <- seq_len(ncol(d$Xdir))
      position <- ncol(d$Xdir) + seq_along(d$vars)
      pos <- d$Xpos / d$k
      colnames(pos) <- paste0(d$vars, ".mean")
      Z <- cbind(Z, d$Xdir, pos)
      known <- cbind(known, moments$spread[, spread, drop = FALSE],
                     moments$position[, position, drop = FALSE])
      colnames(known) <- colnames(Z)
    }
    list(X = X, y = pmin(pmax(y, 0), 1), weights = w,
         control = Z, expectation = known, order = d$k)
  })
  list(orders = orders, focal = focal, counts = attr(samples, "draws"),
       n = n, ref = feature, evaluations = attr(samples, "evaluations"))
}

.conditioned_theta <- function(theta, parts) {
  sizes <- vapply(parts, function(p) ncol(p$X), integer(1))
  if (is.list(theta)) theta <- unlist(theta, use.names = FALSE)
  if (!is.numeric(theta) || length(theta) != sum(sizes) || any(!is.finite(theta))) {
    stop("Supply one finite coefficient vector per order or the concatenated vectors.")
  }
  unname(split(theta, rep(seq_along(parts), sizes)))
}

.conditioned_contrast <- function(contrast, size) {
  if (is.null(contrast)) return(NULL)
  if (is.null(dim(contrast))) contrast <- matrix(contrast, nrow = 1L)
  contrast <- as.matrix(contrast)
  if (!is.numeric(contrast) || ncol(contrast) != size || !nrow(contrast) ||
      any(!is.finite(contrast))) stop("Supply one contrast column per stacked coefficient.")
  contrast
}

.conditioned_raw <- function(part, theta) {
  part$X * (part$weights * (part$y - stats::plogis(drop(part$X %*% theta))))
}

.conditioned_pilot <- function(plan, theta) {
  theta <- .conditioned_theta(theta, plan$orders)
  strata <- split(seq_along(plan$focal), plan$focal)
  fits <- Map(function(part, b) {
    score <- .conditioned_raw(part, b)
    A <- lapply(strata, function(i) {
      if (!ncol(part$control)) return(matrix(0, 0, ncol(part$X)))
      Z <- part$control[i, , drop = FALSE]
      decomp <- qr(sweep(Z, 2L, colMeans(Z)))
      coef <- qr.coef(decomp, score[i, , drop = FALSE])
      coef[is.na(coef)] <- 0
      coef
    })
    corrected <- score
    if (ncol(part$control)) {
      for (i in seq_along(strata)) {
        rows <- strata[[i]]
        Z <- sweep(part$control[rows, , drop = FALSE], 2L, part$expectation[i, ])
        corrected[rows, ] <- score[rows, , drop = FALSE] - Z %*% A[[i]]
      }
    }
    list(A = A, corrected = corrected)
  }, plan$orders, theta)
  list(A = lapply(fits, "[[", "A"), corrected = lapply(fits, "[[", "corrected"))
}

.conditioned_counts <- function(total, minimum, priority) {
  n <- length(priority)
  excess <- total - n * minimum
  if (!any(priority > 0)) priority <- rep(1, n)
  target <- excess * priority / sum(priority)
  count <- floor(target)
  remainder <- excess - sum(count)
  if (remainder > 0) {
    i <- order(target - count, decreasing = TRUE)[seq_len(remainder)]
    count[i] <- count[i] + 1
  }
  as.integer(minimum + count)
}

#' Prepare independently piloted site-conditioned score integration
#'
#' Uses uniform sets conditional on each focal site and holds the predictor
#' reference fixed across orders and the independent pilot and production
#' samples. Automatic allocation divides the production path budget according
#' to pilot within-site SDs in the requested pseudovalue contrasts, with a
#' per-site floor. Sites keep equal weight in the unconditional score.
#'
#' The response control is exact for relative-abundance dissimilarity. For
#' incidence it uses mean richness times normalised dissimilarity, namely
#' mean richness minus shared richness, regardless of working weights.
#' A focal-richness incidence control is defined to be constant zero at empty
#' focal sites. Its target differs from the set-level normalised model
#' response. Predictor controls use exact mean-position and additive
#' mean-pairwise spread expectations from the reference basis.
#'
#' @param x,env Community and environment matrices with aligned site order.
#' @param orders Distinct orders from two through one less than the site count.
#' @param draws Production paths per site. In adaptive mode this is a scalar
#'   average per site, with total budget nrow(x) * draws.
#' @param theta Pilot parameter in the same convention as
#'   [ugdm_conditioned_scores()]. Required for controls or adaptive allocation.
#' @param type Relative-abundance or mean-richness-normalised total incidence.
#' @param weighting Unit working weights or mean richness for incidence.
#' @param controls Any combination of "response", "predictors" and, for
#'   incidence only, "focal_incidence".
#' @param pilot_draws Independent pilot paths per site, used when fitting
#'   controls or adapting the allocation.
#' @param adaptive Estimate an allocation for production draws from the pilot.
#' @param minimum_draws Minimum paths per site for adaptive allocation.
#' @param allocation_contrast Optional contrast matrix with one column per
#'   stacked coefficient. The default uses all coordinates in original units.
#' @param dist_matrix,ref As in [ugdm_path_data()].
#' @param design_matrix Optional function accepting one order element from
#'   [ugdm_path_data()] and returning its numerical model matrix.
#' @param coupled,seed As in [sample_sets_conditioned()].
#' @return A fixed score plan for [ugdm_conditioned_scores()], including
#'   separate pilot and production costs.
#' @export
ugdm_conditioned_score_plan <- function(x, env, orders, draws = 100L, theta = NULL,
    type = c("relative", "incidence"), weighting = c("constant", "richness"),
    controls = character(), pilot_draws = 40L, adaptive = FALSE, minimum_draws = 10L,
    allocation_contrast = NULL, dist_matrix = NULL, ref = NULL, design_matrix = NULL,
    coupled = TRUE, seed = NULL) {
  type <- match.arg(type)
  weighting <- match.arg(weighting)
  if (length(controls)) {
    controls <- match.arg(controls, c("response", "predictors", "focal_incidence"), several.ok = TRUE)
  }
  x <- .conditional_community(x)
  n <- nrow(x)
  orders <- .conditioned_orders(orders, n - 1L)
  if (type == "relative") {
    if (weighting == "richness" || "focal_incidence" %in% controls) {
      stop("Richness working weights and focal incidence controls require incidence data.")
    }
    if (any(rowSums(x) <= 0)) stop("Relative abundance requires positive site totals.")
  }
  if (type == "incidence" && weighting == "constant" && sum(rowSums(x) == 0) >= min(orders)) {
    stop("All-empty sets have undefined incidence response under constant working weights.")
  }
  if (!is.null(seed)) set.seed(seed)
  build <- function(count) .conditioned_data(x, env, sample_sets_conditioned(n, orders, count, coupled),
    dist_matrix, ref, type, weighting, controls, design_matrix)
  fit <- NULL
  pilot_cost <- 0
  if (adaptive || length(controls)) {
    pilot_draws <- .set_design_integer(pilot_draws, "pilot_draws", 2L)
    pilot <- build(pilot_draws)
    pilot_cost <- pilot$evaluations
    ref <- pilot$ref
    fit <- .conditioned_pilot(pilot, theta)
    if (adaptive) {
      draws <- .set_design_integer(draws, "draws")
      minimum_draws <- .set_design_integer(minimum_draws, "minimum_draws")
      if (draws < minimum_draws) stop("Average draws must be at least the allocation floor.")
      transformed <- Map(function(z, k) z * k * (n - 1) / (n - k), fit$corrected, orders)
      paired <- do.call(cbind, transformed)
      L <- .conditioned_contrast(allocation_contrast, ncol(paired))
      if (!is.null(L)) paired <- paired %*% t(L)
      groups <- split(seq_len(nrow(paired)), pilot$focal)
      s <- vapply(groups, function(i) sqrt(sum(apply(paired[i, , drop = FALSE], 2L, stats::var))),
                  numeric(1))
      draws <- .conditioned_counts(n * draws, minimum_draws, s)
    }
  }
  plan <- build(draws)
  plan$coefficients <- if (length(controls)) fit$A else NULL
  plan$ref <- NULL
  plan$strata <- split(seq_along(plan$focal), plan$focal)
  plan$pilot_evaluations <- pilot_cost
  plan
}

#' Evaluate matched site-conditioned scores and pseudovalues
#'
#' Estimates conditional scores by averaging over the independent paths for
#' each focal site. The equal average of conditional scores estimates the
#' global score. The exact complete-set deletion identity then gives
#' pseudovalues. For a contrast, each order's amplified path scores enter
#' the contrast before averaging over focal paths.
#'
#' @param plan A fixed integration plan from [ugdm_conditioned_score_plan()].
#' @param theta Concatenated or list of per-order coefficients in model-matrix
#'   column order. The default uses intercept, spread, optional distance and
#'   position in mean units.
#' @param contrast Optional contrast matrix with one column per stacked
#'   coefficient, or vector for a single contrast.
#' @return Global score, site conditional scores and pseudovalues, optional
#'   site contrast, and counts including every focal site's draw.
#' @export
ugdm_conditioned_scores <- function(plan, theta, contrast = NULL) {
  theta <- .conditioned_theta(theta, plan$orders)
  strata <- plan$strata
  rows <- lapply(seq_along(plan$orders), function(j) {
    block <- plan$orders[[j]]
    score <- .conditioned_raw(block, theta[[j]])
    if (!is.null(plan$coefficients)) for (i in seq_along(strata)) {
      ix <- strata[[i]]
      Z <- sweep(block$control[ix, , drop = FALSE], 2L, block$expectation[i, ])
      score[ix, ] <- score[ix, , drop = FALSE] - Z %*% plan$coefficients[[j]][[i]]
    }
    colnames(score) <- paste(block$order, colnames(block$X), sep = ":")
    score
  })
  average <- function(z) rowsum(z, plan$focal) / plan$counts
  conditional <- do.call(cbind, lapply(rows, average))
  G <- colMeans(conditional)
  multiplier <- rep(vapply(plan$orders, function(b) b$order * (plan$n - 1) /
    (plan$n - b$order), numeric(1)), vapply(rows, ncol, integer(1)))
  pseudovalues <- sweep(sweep(sweep(conditional, 2L, G), 2L, multiplier, "*"), 2L, G, "+")
  L <- .conditioned_contrast(contrast, length(G))
  paired <- NULL
  if (!is.null(L)) {
    direct <- sweep(do.call(cbind, rows), 2L, multiplier, "*") %*% t(L)
    focal <- average(direct)
    paired <- sweep(sweep(focal, 2L, colMeans(focal)), 2L, drop(L %*% G), "+")
  }
  list(score = G, conditional = conditional, pseudovalues = pseudovalues, contrast = paired,
       counts = plan$counts, evaluations = plan$evaluations, pilot_evaluations = plan$pilot_evaluations)
}
