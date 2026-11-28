.st3_coalitions <- function() {
  masks <- 0:15
  powers <- c(1L, 2L, 4L, 8L)
  size <- rowSums(sapply(powers, function(b) bitwAnd(masks, b) != 0L))
  contrast <- matrix(0, 16L, 4L, dimnames = list(NULL, st3_processes))
  for (j in seq_along(powers)) {
    missing <- masks[bitwAnd(masks, powers[j]) == 0L]
    w <- factorial(size[missing + 1L]) * factorial(3L - size[missing + 1L]) / factorial(4L)
    contrast[missing + 1L, j] <- -w
    contrast[missing + powers[j] + 1L, j] <- w
  }
  contrast
}

st3_apply <- function(plan, mask, coalition) {
  replace <- bitwAnd(mask, coalition) != 0L
  a <- plan$before
  a[replace] <- plan$after[replace]
  a
}

st3_shapley <- function(plan, sets, responses = .st3_responses, exact_limit = 32L,
                       draws = 32L, seed = 1L, retain_draws = FALSE, limit = 1e5) {
  if (!is.character(responses) || !length(responses) || anyDuplicated(responses) ||
      !all(responses %in% .st3_responses)) stop("Select distinct supported diversity responses.")
  n <- length(plan$sites)
  m <- nrow(plan$alignment)
  if (is.null(sets)) stop("Supply a matrix of spatial sets.")
  k <- ncol(sets)
  sets <- .st_sets(n, k, sets, limit)
  if (nrow(sets) > limit) stop("Spatial combination budget exceeded.")
  idx <- .st3_site_indices(sets, n, m)
  allocations <- st3_allocations(plan, exact_limit, draws, seed)
  before <- .st3_rectangle_state(plan$before, idx, responses)
  evaluate <- function(state) {
    z <- c(state$metrics, .st3_temporal_metrics(before, state, idx, responses))
    do.call(cbind, z[responses])
  }
  after <- .st3_rectangle_state(plan$after, idx, responses)
  F0 <- evaluate(before)
  F1 <- evaluate(after)
  shape <- c(nrow(sets), length(responses), 4L)
  phi_mean <- dispersion <- array(0, shape, dimnames = list(NULL, responses, st3_processes))
  total <- length(allocations$masks)
  valid <- is.finite(F0) & is.finite(F1)
  valid_count <- matrix(0L, nrow(F0), ncol(F0))
  retained <- if (retain_draws) vector("list", total) else NULL
  weights <- .st3_coalitions()
  flat0 <- as.vector(F0)
  flat1 <- as.vector(F1)
  err <- 0
  for (j in seq_len(total)) {
    mask <- allocations$masks[[j]]
    active <- sum(unique(as.vector(mask)))
    ids <- bitwAnd(0:15, active)
    values <- lapply(unique(ids), function(i) {
      if (i == 0L) return(flat0)
      if (i == active) return(flat1)
      as.vector(evaluate(.st3_rectangle_state(st3_apply(plan, mask, i), idx, responses)))
    })
    V <- do.call(cbind, values[match(ids, unique(ids))])
    good <- rowSums(!is.finite(V)) == 0L
    valid <- valid & matrix(good, nrow(valid))
    valid_count <- valid_count + matrix(good, nrow(valid))
    V[!good, ] <- 0
    attribution <- vapply(seq_len(4L), function(j) {
      upper <- which(weights[, j] > 0)
      lower <- upper - 2L^(j - 1L)
      drop((V[, upper, drop = FALSE] - V[, lower, drop = FALSE]) %*% weights[upper, j])
    }, numeric(nrow(V)))
    attribution <- matrix(attribution, ncol = 4L)
    err <- max(c(err, abs(rowSums(attribution)[good] - (flat1 - flat0)[good])))
    phi <- array(attribution, shape, dimnames = dimnames(phi_mean))
    deviation <- phi - phi_mean
    phi_mean <- phi_mean + deviation / j
    dispersion <- dispersion + deviation * (phi - phi_mean)
    if (retain_draws) {
      phi[!array(rep(good, 4L), shape)] <- NA_real_
      retained[[j]] <- phi
    }
  }
  phi <- phi_mean
  sd <- if (total > 1L) sqrt(pmax(dispersion, 0) / (total - 1L)) else dispersion
  se <- if (allocations$exact) sd * 0 else sd / sqrt(total)
  compatible <- array(rep(valid, 4L), shape)
  phi[!compatible] <- sd[!compatible] <- se[!compatible] <- NA_real_
  if (retain_draws) retained <- lapply(retained, function(a) {a[!compatible] <- NA_real_; a})
  list(contributions = phi, allocation_sd = sd, allocation_mc_se = se,
       draws = retained, before = F0, after = F1, delta = F1 - F0, valid = valid,
       allocation_support = valid_count / total, allocations = total, exact = allocations$exact,
       log_allocation_count = plan$log_allocation_count, max_efficiency_error = err,
       sites = lapply(seq_len(nrow(sets)), function(i) plan$sites[sets[i, ]]),
       alignment = plan$alignment)
}

st3_crossed <- function(x) {
  result <- x$contributions[, c("occupancy", "numerical"), , drop = FALSE]
  support <- x$valid[, "occupancy"] & x$valid[, "numerical"] & x$valid[, "abundance"]
  value <- aperm(result, c(1, 3, 2))
  value[!support, , ] <- NA_real_
  list(contributions = value, total_by_process = apply(value, c(1L, 2L), sum),
       delta = ifelse(support, x$delta[, "abundance"], NA_real_), valid = support)
}

st3_signed_summary <- function(phi, weights = NULL) {
  if (!is.matrix(phi) || ncol(phi) < 1L) stop("Supply a matrix of process contributions.")
  if (is.null(weights)) weights <- rep(1, nrow(phi))
  if (length(weights) != nrow(phi) || any(!is.finite(weights) | weights < 0)) stop("Invalid weights.")
  p <- rowSums(pmax(phi, 0))
  n <- rowSums(pmax(-phi, 0))
  gross <- p + n
  net <- p - n
  use <- complete.cases(phi) & is.finite(gross) & weights > 0
  means <- function(a) colSums(a[use, , drop = FALSE] * weights[use]) / sum(weights[use])
  list(by_combination = data.frame(differentiation = p, homogenisation = n, net, gross),
       importance = .st_ratio(abs(phi), gross),
       mean_positive = means(pmax(phi, 0)), mean_negative = means(pmax(-phi, 0)),
       mean_absolute = means(abs(phi)), absolute_mean = abs(means(phi)), support = sum(use))
}
