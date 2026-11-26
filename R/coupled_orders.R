.path_permutations <- function(n, k, M, balanced) {
  if (!balanced) return(t(vapply(seq_len(M), function(i) sample.int(n, k), integer(k))))
  idx <- matrix(0L, M, k)
  for (start in seq.int(1L, M, by = n)) {
    m <- min(n, M - start + 1L)
    labels <- sample.int(n)
    offsets <- if (m == n) seq_len(n) - 1L else sample.int(n, m) - 1L
    idx[seq.int(start, length.out = m), ] <- matrix(
      labels[(outer(offsets, seq_len(k) - 1L, "+") %% n) + 1L], m, k)
  }
  idx
}

#' Coupled site combinations across zeta orders
#'
#' Uniform paths use common permutation prefixes. Balanced paths use all
#' cyclic starts of a random site permutation per full block. Species
#' mixtures share their branch and inverse-CDF uniforms and use a common
#' random ranking restricted to the selected species' occupied sites.
#' Species-mixture paths need not be nested. Exhaustive orders have no
#' random row correspondence with other orders.
#' @param x Numeric community matrix, sites by species.
#' @param orders Distinct orders from two through the number of sites.
#' @param M Number of random paths.
#' @param sampler Uniform, species-informed, or cyclic site-balanced.
#' @param alpha Targeted mixture fraction, scalar or one value per supplied order.
#' @param coupled Share randomness across orders. Independent mode uses the
#'   same order-specific marginal designs.
#' @param exact_limit Enumerate when there are at most this many sets and
#'   enumeration costs at most M set evaluations. Zero disables enumeration.
#' @param seed Optional random seed.
#' @return Named list of order-specific indices, raw proposal ratios and
#'   sampling metadata. Random draws at each order use row number as path ID.
#' @export
sample_sets_coupled <- function(x, orders, M = 10000L,
                                sampler = c("uniform", "species", "site"),
                                alpha = .5, coupled = TRUE, exact_limit = 50000L, seed = NULL) {
  x <- .conditional_community(x)
  n <- nrow(x)
  M <- .set_design_integer(M, "M")
  exact_limit <- .set_design_integer(exact_limit, "exact_limit", 0L)
  sampler <- match.arg(sampler)
  if (!length(orders) || anyDuplicated(orders)) stop("Supply distinct orders.")
  orders <- vapply(orders, .zeta_order_integer, integer(1), n = n, minimum = 2L)
  if (length(alpha) == 1L) alpha <- rep(alpha, length(orders))
  if (length(alpha) != length(orders)) stop("Supply one alpha per order or a scalar.")
  for (a in alpha) .check_alpha(a)
  if (!is.null(seed)) set.seed(seed)
  exact <- lchoose(n, orders) <= log(min(M, exact_limit)) + 1e-12
  P <- (x > 0) * 1L
  if (coupled && any(!exact)) {
    path <- .path_permutations(n, if (sampler == "species" && any(alpha[!exact] > 0)) n else max(orders[!exact]),
                               M, sampler == "site")
    if (sampler == "species") {
      branch <- stats::runif(M)
      select <- stats::runif(M)
    }
  }
  result <- lapply(seq_along(orders), function(j) {
    k <- orders[j]
    if (exact[j]) {
      idx <- t(utils::combn(n, k))
      return(list(idx = idx, sampling_weights = rep(1, nrow(idx)), exact = TRUE))
    }
    if (sampler == "species" && !coupled) {
      return(c(.sample_sets_mixture(P, k, M, alpha[j]), list(exact = FALSE)))
    }
    if (sampler != "species") {
      idx <- if (coupled) path[, seq_len(k), drop = FALSE] else
        .path_permutations(n, k, M, sampler == "site")
      return(list(idx = idx, sampling_weights = rep(1, M), exact = FALSE))
    }
    plan <- .sampling_plan(P, k, NULL)
    a <- if (is.finite(plan$log_expected_shared)) alpha[j] else 0
    targeted <- branch < a
    idx <- path[, seq_len(k), drop = FALSE]
    if (any(targeted)) {
      eligible <- which(plan$species_probability > 0)
      cm <- cumsum(plan$species_probability[eligible])
      cm[length(cm)] <- 1
      s <- eligible[findInterval(select[targeted], cm) + 1L]
      rows <- which(targeted)
      for (r in seq_along(s)) {
        row <- rows[r]
        idx[row, ] <- utils::head(path[row, ][P[path[row, ], s[r]] > 0], k)
      }
    }
    z <- shared_species_cpp(P, idx)
    list(idx = idx, sampling_weights = .mixture_weights(z, a, plan$log_expected_shared),
         shared = z, targeted = targeted, alpha = a, exact = FALSE)
  })
  structure(stats::setNames(result, orders), sampler = sampler, coupled = coupled)
}

#' Compute matched MS-GDUM designs and responses from paths
#'
#' Reuses running intersections, abundance minima, richness sums, position
#' sums and pairwise distance totals across nested prefixes. Orders without
#' a common prefix (including exact orders) are evaluated separately.
#' @param x,env Sites by species and sites by predictor matrices.
#' @param samples Output of [sample_sets_coupled()].
#' @param dist_matrix Optional site-by-site distance matrix.
#' @param ref Optional fixed reference design from [ugdm_design()].
#' @return Per-order design, matched normalised incidence and relative
#'   abundance responses, shared richness and raw sampling weights.
#' @export
ugdm_path_data <- function(x, env, samples, dist_matrix = NULL, ref = NULL) {
  x <- .conditional_community(x)
  n <- nrow(x)
  if (nrow(env) != n || !length(samples)) stop("Provide matching sites and a nonempty set sample.")
  if (!is.null(dist_matrix) && (!identical(dim(dist_matrix), c(n, n)) ||
      !is.numeric(dist_matrix) || any(!is.finite(dist_matrix)) || !isSymmetric(dist_matrix))) {
    stop("Supply a finite symmetric site-by-site distance matrix.")
  }
  if (is.null(ref)) ref <- ugdm_design(env, cbind(seq_len(n - 1L), seq.int(2L, n)),
    dist_matrix, dist_range = if (!is.null(dist_matrix)) .pair_distance_range(dist_matrix))
  if (!is.null(dist_matrix) && is.null(ref$dist_basis)) stop("Reference requires a distance basis.")
  e <- .environment_for_moments(env, ref)
  p <- ncol(e$z)
  totals <- rowSums(x)
  q <- x / ifelse(totals > 0, totals, 1)
  k <- vapply(samples, function(s) ncol(s$idx), integer(1))
  idx <- lapply(samples, function(s) .species_indices(s$idx, n))
  kernel <- function(i) path_kernels_cpp(i, x, q, cbind(e$spline, e$z), e$z,
    if (is.null(dist_matrix)) matrix(numeric(), 0, 0) else as.matrix(dist_matrix))
  values <- vector("list", length(samples))
  exact <- vapply(samples, function(s) isTRUE(s$exact), TRUE)
  random <- which(!exact)
  groups <- split(random, vapply(idx[random], nrow, integer(1)))
  for (j in groups) {
    j <- j[order(k[j])]
    values[j] <- kernel(idx[j])
  }
  for (i in which(exact)) values[i] <- kernel(idx[i])
  out <- lapply(seq_along(k), function(i) {
    a <- values[[i]]
    d <- ref
    d$idx <- idx[[i]]
    d$k <- unname(k[i])
    d$Xdir <- a[, 5L + seq_len(3L * p), drop = FALSE]
    colnames(d$Xdir) <- colnames(ref$Xdir)
    d$Xlinear <- a[, 5L + 3L * p + seq_len(p), drop = FALSE]
    colnames(d$Xlinear) <- ref$vars
    d$Xpos <- a[, 5L + 4L * p + seq_len(p), drop = FALSE]
    colnames(d$Xpos) <- ref$vars
    d$Xdist <- if (!is.null(dist_matrix)) {
      z <- pmin(pmax((a[, 5L] - ref$dist_range[1]) / max(diff(ref$dist_range), 1e-12), 0), 1)
      b <- unclass(stats::predict(ref$dist_basis, z))[, 1:3, drop = FALSE]
      colnames(b) <- colnames(ref$Xdist)
      b
    }
    list(design = d, incidence = a[, 3L], relative = a[, 4L], shared = a[, 1L],
         mean_richness = a[, 2L], sampling_weights = samples[[i]]$sampling_weights,
         exact = isTRUE(samples[[i]]$exact))
  })
  stats::setNames(out, k)
}
