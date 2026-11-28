st2_windows <- function(times, temporal_order = 1L, starts = NULL, tolerance = 1e-8) {
  if (length(times) < 2L || any(!is.finite(times)) || any(diff(times) <= 0)) {
    stop("Use finite and strictly increasing survey times.")
  }
  m <- .st_one(temporal_order, length(times) %/% 2L, "temporal_order")
  if (is.null(starts)) {
    i <- seq_len(length(times) - 2L * m + 1L)
    starts <- cbind(i, i + m)
  }
  if (!is.matrix(starts) || ncol(starts) != 2L || !nrow(starts) || !is.numeric(starts) ||
      any(!is.finite(starts) | starts < 1 | starts != floor(starts) | starts + m - 1L > length(times)) ||
      any(starts[, 2L] - starts[, 1L] < m) || anyDuplicated(as.data.frame(starts))) {
    stop("Supply distinct pairs of non-overlapping chronological consecutive windows.")
  }
  a <- outer(starts[, 1L], seq_len(m) - 1L, "+")
  b <- outer(starts[, 2L], seq_len(m) - 1L, "+")
  at <- matrix(times[a], nrow(a))
  bt <- matrix(times[b], nrow(b))
  if (any(apply(abs((at - at[, 1L]) - (bt - bt[, 1L])), 1L, max) > tolerance)) {
    stop("Paired windows must have comparable within-window survey spacing.")
  }
  list(from = a, to = b, starts = starts, lag = rowMeans(bt) - rowMeans(at),
       temporal_order = m, times = times)
}

st2_basis <- function(x, env, training_times = NULL) {
  p <- .st_panel(x)
  env <- .st_environment(p, env)
  if (is.null(training_times)) training_times <- seq_along(p$times)
  training_times <- .st_orders(training_times, length(p$times), "training_times")
  flat <- .st_flat_env(env[, training_times, , drop = FALSE])
  zetaGDUM::ugdm_design(flat, matrix(c(1L, 2L), 1L))$environment_basis
}

# Equal-probability exhaustive or sampled sets use unit weights. A sampled set
# drawn with replacement using a known proposal q(I) receives inverse-q weight.
st2_transition_design <- function(x, env, spatial_order, temporal_order = 1L, sets = NULL,
                                  windows = NULL, basis = NULL, proposal = NULL,
                                  limit = 1e5) {
  p <- .st_panel(x)
  e <- .st_environment(p, env)
  n <- dim(e)[1L]
  k <- .st_one(spatial_order, n, "spatial_order")
  sets <- .st_sets(n, k, sets, limit)
  if (is.null(windows)) windows <- st2_windows(p$times, temporal_order)
  if (!identical(windows$times, p$times) || windows$temporal_order != temporal_order) {
    stop("Windows must match the panel times and temporal order.")
  }
  if (as.double(nrow(sets)) * nrow(windows$from) * k > limit) stop("Local design row budget exceeded.")
  if (!is.null(proposal) && (!is.numeric(proposal) || length(proposal) != nrow(sets) ||
                           any(!is.finite(proposal) | proposal <= 0 | proposal > 1))) {
    stop("Supply positive draw probabilities for each explicit combination.")
  }
  if (is.null(basis)) basis <- st2_basis(p, e, sort(unique(c(windows$from, windows$to))))
  vars <- dimnames(e)[[3L]]
  if (!identical(vars, names(basis$center))) stop("Predictors must match the training basis.")
  nw <- nrow(windows$from)
  v <- length(vars)
  # Both blocks contain site-level window means in original environmental units.
  window_environment <- function(w) {
    means <- lapply(seq_len(nrow(w)), function(j) {
      matrix(colMeans(aperm(e[, w[j, ], , drop = FALSE], c(2, 1, 3)), dims = 1L), n, v)
    })
    out <- do.call(rbind, means)
    colnames(out) <- vars
    out
  }
  flat <- rbind(window_environment(windows$from), window_environment(windows$to))
  z <- scale(flat, basis$center, basis$scale)
  rows <- expand.grid(set = seq_len(nrow(sets)), window = seq_len(nw))
  idx <- sets[rows$set, , drop = FALSE]
  initial <- idx + n * (rows$window - 1L)
  final <- initial + n * nw
  nr <- nrow(rows)
  row_id <- rep(seq_len(nr), each = k)
  gather <- function(a) as.vector(t(a))
  z0 <- z[gather(initial), , drop = FALSE]
  change <- z[gather(final), , drop = FALSE] - z0
  names_for <- function(a, label) {
    colnames(a) <- paste(label, colnames(a), sep = ".")
    a
  }
  if (k > 1L) {
    shared <- zetaGDUM::ugdm_design(flat, rbind(initial, final), ref = list(environment_basis = basis))
    b <- seq_len(nr)
    p0 <- shared$Xpos[b, , drop = FALSE]
    dp <- shared$Xpos[b + nr, , drop = FALSE] - p0
    s0 <- shared$Xdir[b, , drop = FALSE]
    ds <- shared$Xdir[b + nr, , drop = FALSE] - s0
    extra <- cbind(names_for(z0 - p0[row_id, , drop = FALSE] / k, "deviation"),
                   names_for(change - dp[row_id, , drop = FALSE] / k, "deviation_change"),
                   names_for(s0[row_id, , drop = FALSE], "spread"),
                   names_for(ds[row_id, , drop = FALSE], "spread_change"))
  } else {
    p0 <- z0
    dp <- change
    extra <- NULL
  }
  out <- cbind(names_for(p0[row_id, , drop = FALSE], "position"),
               names_for(dp[row_id, , drop = FALSE], "position_change"),
               extra, log_lag = log(windows$lag[rows$window[row_id]]),
               names_for((p0 * dp)[row_id, , drop = FALSE], "position_by_change"))
  weights <- if (is.null(proposal)) rep(1, nrow(sets)) else {
    w <- exp(min(log(proposal)) - log(proposal))
    w / mean(w)
  }
  rows$weight <- weights[rows$set]
  structure(list(X = out, rows = rows, row_id = row_id, site = gather(idx), sets = idx,
                 windows = windows, basis = basis, k = k, sites = dimnames(p$x)[[1L]],
                 times = p$times, weight = rows$weight[row_id]), class = "st2_design")
}

st2_states <- function(x, design) {
  p <- .st_panel(x)
  if (!identical(design$sites, dimnames(p$x)[[1L]]) || !identical(design$times, p$times)) {
    stop("Panel and design site and time axes must match.")
  }
  build <- function(w) {
    a <- lapply(seq_len(nrow(w)), function(j) .st2_minimum(p$x, w[j, ]))
    do.call(rbind, a)
  }
  windows <- design$rows$window[design$row_id]
  i <- design$site + dim(p$x)[1L] * (windows - 1L)
  from <- build(design$windows$from)[i, , drop = FALSE]
  to <- build(design$windows$to)[i, , drop = FALSE]
  dimnames(from) <- dimnames(to) <- list(NULL, dimnames(p$x)[[3L]])
  z <- from > 0
  y <- to > 0
  initial_count <- rowsum(z * 1L, design$row_id, reorder = FALSE)
  other <- if (design$k > 1L) initial_count[design$row_id, , drop = FALSE] - z else NULL
  list(z = z, y = y, other = other, before = from, after = to)
}

.st2_terms <- function(design, interaction) {
  terms <- colnames(design$X)
  if (!interaction) terms <- terms[!startsWith(terms, "position_by_change.")]
  if (length(unique(design$X[, "log_lag"])) == 1L) terms <- setdiff(terms, "log_lag")
  terms
}

# Collapse identical Bernoulli covariate rows to exact success and trial totals.
# Each physical row represents a site inside one combination and time contrast.
st2_risk_data <- function(design, states, event, terms, pooled = FALSE) {
  event <- match.arg(event, c("gain", "loss"))
  z <- states$z
  y <- states$y
  X <- cbind(intercept = 1, design$X[, terms, drop = FALSE])
  if (pooled) {
    z <- rowsum(z * 1L, design$row_id, reorder = FALSE) == design$k
    y <- rowsum(y * 1L, design$row_id, reorder = FALSE) == design$k
    first <- !duplicated(design$row_id)
    X <- X[first, , drop = FALSE]
    local <- which(first)
  } else {
    local <- seq_len(nrow(X))
  }
  opportunity <- if (event == "gain") !z else z
  success <- opportunity & if (event == "gain") y else !y
  if (design$k == 1L || pooled) {
    return(list(X = X, y = rowSums(success), n = rowSums(opportunity), local = local,
                weight = design$weight[local]))
  }
  collection <- lapply(0:(design$k - 1L), function(j) {
    category <- states$other == j
    list(X = cbind(X, other_occupancy = j / (design$k - 1L)),
         y = rowSums(success & category), n = rowSums(opportunity & category))
  })
  list(X = do.call(rbind, lapply(collection, `[[`, "X")),
       y = unlist(lapply(collection, `[[`, "y")), n = unlist(lapply(collection, `[[`, "n")),
       local = rep(local, design$k), weight = rep(design$weight, design$k))
}

.st2_fit_grouped <- function(group, ridge, multiplier = 1, diagnostics = TRUE) {
  w <- group$weight * multiplier
  keep <- group$n > 0 & w > 0
  A <- group$X[keep, , drop = FALSE]
  if (diagnostics) {
    singular <- if (nrow(A)) svd(A * sqrt(group$n[keep] * w[keep]), nu = 0, nv = 0)$d else numeric()
    q <- if (nrow(A)) qr(A) else NULL
    aliased <- if (!is.null(q) && q$rank < ncol(A)) {
      colnames(A)[q$pivot[seq.int(q$rank + 1L, ncol(A))]]
    } else character()
  }
  f <- .st_fit_risk(A, group$y[keep] * w[keep], group$n[keep] * w[keep], ridge)
  if (diagnostics) {
    f$singular_values <- singular
    f$aliased <- aliased
    f$opportunities <- sum(group$n[keep])
    f$events <- sum(group$y[keep])
    f$weighted_opportunities <- sum(group$n[keep] * w[keep])
  }
  f
}

st2_fit_transition <- function(x, design, interaction = FALSE, ridge = 0.1, pooled = FALSE) {
  if (!inherits(design, "st2_design") || !is.numeric(ridge) || length(ridge) != 1L ||
      !is.finite(ridge) || ridge <= 0) stop("Supply a Phase II design and a positive ridge.")
  states <- st2_states(x, design)
  terms <- .st2_terms(design, interaction)
  if (pooled) terms <- terms[!grepl("^deviation", terms)]
  events <- c("gain", "loss")
  groups <- setNames(lapply(events, function(e) st2_risk_data(design, states, e, terms, pooled)), events)
  fit <- lapply(groups, .st2_fit_grouped, ridge = ridge)
  structure(list(fit = fit, terms = terms, design = design, groups = groups, ridge = ridge,
                 pooled = pooled, species = colnames(states$z)), class = "st2_transition")
}

st2_predict_transition <- function(object, design = object$design, states) {
  if (!identical(design$basis, object$design$basis) || design$k != object$design$k ||
      design$windows$temporal_order != object$design$windows$temporal_order ||
      !identical(colnames(states$z), object$species) || nrow(states$z) != nrow(design$X)) {
    stop("Prediction requires the training basis, orders and species pool with aligned states.")
  }
  X <- cbind(intercept = 1, design$X[, object$terms, drop = FALSE])
  if (object$pooled) {
    H0 <- rowsum(states$z * 1L, design$row_id, reorder = FALSE) == design$k
    X <- X[!duplicated(design$row_id), , drop = FALSE]
    g <- .st_predict_risk(object$fit$gain, X)
    l <- .st_predict_risk(object$fit$loss, X)
    return(list(shared = ifelse(H0, 1 - l, g), initial_shared = H0))
  }
  score <- function(event) {
    f <- object$fit[[event]]
    if (design$k == 1L || f$status != "estimated") {
      return(matrix(.st_predict_risk(f, X), nrow(X), ncol(states$z)))
    }
    beta <- f$coefficients
    eta <- drop(X %*% beta[colnames(X)])
    plogis(eta + beta[["other_occupancy"]] * states$other / (design$k - 1L))
  }
  list(gain = score("gain"), loss = score("loss"))
}

st2_transition_summary <- function(prediction, design, states, occupancy = TRUE) {
  before <- rowsum(states$z * 1L, design$row_id, reorder = FALSE) == design$k
  if (is.null(prediction$shared)) {
    p <- ifelse(states$z, 1 - prediction$loss, prediction$gain)
    shared <- matrix(1, nrow(before), ncol(before))
    for (j in seq_len(design$k)) shared <- shared * p[seq.int(j, nrow(p), design$k), , drop = FALSE]
  } else {
    shared <- prediction$shared
  }
  events <- list(final_shared = shared, initial_shared = before,
                 gains = ifelse(before, 0, shared), losses = ifelse(before, 1 - shared, 0))
  colnames(events$final_shared) <- colnames(states$z)
  if (!occupancy || !is.null(prediction$shared)) return(events)
  pairs <- split(seq_len(nrow(design$X)), design$row_id)
  values <- vapply(pairs, function(i) {
    vapply(seq_len(ncol(states$z)), function(s) {
      z <- states$z[i, s]
      g <- prediction$gain[i[!z], s]
      l <- prediction$loss[i[z], s]
      if (anyNA(c(g, l))) return(rep(NA_real_, 4L))
      result <- st2_occupancy_probability(g, l)
      unlist(result[c("expansion", "contraction", "pure_redistribution", "redistribution")])
    }, numeric(4L))
  }, matrix(0, 4L, ncol(states$z)))
  for (j in 1:4) {
    event <- c("expansion", "contraction", "pure_redistribution", "redistribution")[j]
    events[[event]] <- t(matrix(values[j, , ], ncol(states$z), nrow(before)))
  }
  events
}
