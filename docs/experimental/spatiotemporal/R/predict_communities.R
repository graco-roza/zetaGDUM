.st2_draw_amount <- function(parameters, present, u) {
  a <- matrix(0, nrow(present), ncol(present))
  if (anyNA(present) || anyNA(parameters$location[present]) || anyNA(parameters$scale[present])) {
    stop("A selected prediction requires an identified occurrence and positive-abundance fit.")
  }
  location <- parameters$location[present]
  spread <- parameters$scale[present]
  p <- pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, u[present]))
  if (parameters$family == "lognormal") {
    a[present] <- exp(location + spread * stats::qnorm(p))
  } else {
    mu <- exp(location)
    zero <- stats::dnbinom(0, mu = mu, size = spread)
    a[present] <- pmax(1, stats::qnbinom(zero + (1 - zero) * p, mu = mu, size = spread))
  }
  a
}

st2_predict_communities <- function(transition, positive, design, states, draws = 20L,
                                     uniform_presence = NULL, uniform_abundance = NULL) {
  if (design$windows$temporal_order != 1L) stop("Use survey trajectories for higher temporal abundance orders.")
  if (is.null(transition$gain) || is.null(transition$loss)) stop("Use local gain and loss probabilities.")
  p <- ifelse(states$z, 1 - transition$loss, transition$gain)
  abundance <- st2_predict_positive_parameters(positive, design, states)
  shape <- c(dim(p), draws)
  if (is.null(uniform_presence)) uniform_presence <- array(stats::runif(prod(shape)), shape)
  if (is.null(uniform_abundance)) uniform_abundance <- array(stats::runif(prod(shape)), shape)
  if (!identical(as.numeric(dim(uniform_presence)), as.numeric(shape)) ||
      !identical(as.numeric(dim(uniform_abundance)), as.numeric(shape))) stop("Uniform draws must match prediction dimensions.")
  lapply(seq_len(draws), function(j) {
    result <- .st2_draw_amount(abundance, p > uniform_presence[, , j], uniform_abundance[, , j])
    colnames(result) <- colnames(states$z)
    result
  })
}

# The first observed assemblage anchors both enabled and disabled forecasts.
# Original survey histories are generated anew at every subsequent step.
st2_predict_trajectory <- function(object, positive, x, env, uniform_presence = NULL,
                                   uniform_abundance = NULL) {
  p <- .st_panel(x)
  d <- dim(p$x)
  if (!identical(dimnames(p$x)[[3L]], object$species)) stop("Species pool must match M2 training.")
  design <- st2_transition_design(p, env, 1L, basis = object$design$basis)
  transition <- st2_predict_dynamic(object, design, env)
  design <- st2_neighbourhood_design(design, env, object$neighbourhood)
  shape <- c(d[1L], d[2L] - 1L, d[3L])
  if (is.null(uniform_presence)) uniform_presence <- array(stats::runif(prod(shape)), shape)
  if (is.null(uniform_abundance)) uniform_abundance <- array(stats::runif(prod(shape)), shape)
  if (!identical(dim(uniform_presence), shape) || !identical(dim(uniform_abundance), shape)) stop("Incorrect trajectory innovation dimensions.")
  a <- p$x
  as_site_species <- function(z) matrix(z, d[1L], d[3L], dimnames = dimnames(a)[c(1L, 3L)])
  for (j in seq_len(d[2L] - 1L)) {
    rows <- which(design$rows$window == j)
    before <- as_site_species(a[, j, ])
    z <- before > 0
    step <- design
    step$X <- design$X[rows, , drop = FALSE]
    state <- list(z = z, before = before)
    ab <- st2_predict_positive_parameters(positive, step, state)
    risk <- ifelse(z, 1 - transition$loss[rows, , drop = FALSE], transition$gain[rows, , drop = FALSE])
    a[, j + 1L, ] <- .st2_draw_amount(ab, risk > as_site_species(uniform_presence[, j, ]),
                                      as_site_species(uniform_abundance[, j, ]))
  }
  st_panel(a, p$times)
}

# Incidence occupancy of windows in an already generated trajectory. No use
# of future observations to estimate a fitted dynamic transition probability.
st2_observed_projection <- function(x, design) {
  events <- st2_states(x, design)
  st2_transition_summary(list(shared = rowsum(events$y * 1L, design$row_id, reorder = FALSE) == design$k),
                        design, events, occupancy = FALSE)
}

st2_shared_coherence <- function(smaller, larger, tolerance = 1e-10) {
  a <- smaller$design
  b <- larger$design
  if (!identical(a$sites, b$sites) || !identical(a$times, b$times) ||
      !identical(dimnames(smaller$shared)[[2L]], dimnames(larger$shared)[[2L]])) stop("Compare matched axes and species.")
  tested <- violations <- 0
  maximum <- 0
  for (j in seq_len(nrow(b$sets))) {
    w <- b$rows$window[j]
    within <- vapply(seq_len(nrow(a$windows$from)), function(i) {
      all(a$windows$from[i, ] %in% b$windows$from[w, ]) &&
        all(a$windows$to[i, ] %in% b$windows$to[w, ])
    }, logical(1))
    cand <- which(a$rows$window %in% which(within))
    for (i in cand) {
      if (!all(a$sets[i, ] %in% b$sets[j, ])) next
      difference <- larger$shared[j, ] - smaller$shared[i, ]
      tested <- tested + sum(is.finite(difference))
      violations <- violations + sum(difference > tolerance, na.rm = TRUE)
      maximum <- max(c(maximum, difference), na.rm = TRUE)
    }
  }
  c(tested = tested, violations = violations, maximum_excess = maximum)
}

st2_overlap_coherence <- function(design, states, prediction, tolerance = 1e-10) {
  p <- ifelse(states$z, 1 - prediction$loss, prediction$gain)
  key <- paste(design$site, design$rows$window[design$row_id])
  groups <- split(seq_along(key), key)
  gaps <- unlist(lapply(groups[lengths(groups) > 1L], function(i) {
    apply(p[i, , drop = FALSE], 2L, function(z) if (all(is.finite(z))) diff(range(z)) else NA_real_)
  }))
  c(tested = sum(is.finite(gaps)), violations = sum(gaps > tolerance, na.rm = TRUE),
    maximum_gap = max(c(0, gaps), na.rm = TRUE))
}

st2_beta_accounting <- function(a, b, sets, components = c("total", "turnover", "nestedness",
                                                           "relative", "balanced", "bray_curtis", "gradient")) {
  if (!is.matrix(a) || !identical(dim(a), dim(b)) || any(!is.finite(a) | !is.finite(b) | a < 0 | b < 0)) {
    stop("Use matching finite nonnegative site-by-species matrices.")
  }
  responses <- lapply(components, function(component) {
    evaluate <- function(z) {
      v <- zetaGDUM::ugdm_response(z, sets, component)
      ifelse(v$keep, v$y, NA_real_)
    }
    z0 <- evaluate(a)
    z1 <- evaluate(b)
    zg <- evaluate(pmax(a, b))
    zl <- evaluate(pmin(a, b))
    gain <- ((zg - z0) + (z1 - zl)) / 2
    loss <- ((zl - z0) + (z1 - zg)) / 2
    data.frame(component, set = seq_len(nrow(sets)), before = z0, after = z1, change = z1 - z0,
               gain, loss, differentiating = pmax(gain, 0) + pmax(loss, 0),
               homogenising = pmin(gain, 0) + pmin(loss, 0))
  })
  do.call(rbind, responses)
}

st2_tatsumi <- function(a, b, index = c("sorensen", "bray-curtis"), components = c("four", "six")) {
  index <- match.arg(index)
  components <- match.arg(components)
  if (!is.matrix(a) || nrow(a) != 2L || !identical(dim(a), dim(b))) stop("Tatsumi pairwise accounting uses exactly two sites.")
  if (!requireNamespace("ecopart", quietly = TRUE)) stop("Install the authors' ecopart package for this comparator.")
  ecopart::ecopart.pair(a, b, index = index, components = components)
}
