# Phase II. Domains and conventions are inherited from the balanced panel.
st2_community <- function(x) {
  p <- .st_panel(x)
  d <- dim(p$x)
  a <- matrix(p$x, d[1] * d[2], d[3])
  total <- rowSums(a)
  list(panel = p, X = p$x > 0,
       Q = array(a / ifelse(total > 0, total, NA_real_), d, dimnames(p$x)),
       D = array(zetaGDUM::dominance_matrix(a, p = 0.8) > 0, d, dimnames(p$x)),
       empty = matrix(total == 0, d[1], d[2]))
}

.st2_time_sets <- function(times, m, mode, budget) {
  if (mode == "all") return(.st_indices(length(times), m, budget))
  count <- length(times) - m + 1L
  if (count > budget) stop("Temporal combination budget exceeded.")
  outer(seq_len(count), seq_len(m) - 1L, "+")
}

.st2_minimum <- function(a, times) {
  n <- dim(a)[1]
  s <- dim(a)[3]
  out <- matrix(a[, times[1L], ], n, s)
  for (j in times[-1L]) out[] <- pmin(out, a[, j, ])
  out
}

.st2_weights <- function(x) {
  total <- rowSums(x)
  x / ifelse(total == 0, 1, total)
}

# Returns exact order summaries and species-level contributions on the same grid.
st2_zeta <- function(x, spatial_max = dim(.st_panel(x)$x)[1],
                     temporal_max = NULL,
                     mode = c("all", "consecutive"), times = NULL,
                     empty = c("error", "undefined"), time_limit = 1e5) {
  mode <- match.arg(mode)
  empty <- match.arg(empty)
  z <- st2_community(x)
  d <- dim(z$X)
  if (is.null(times)) times <- seq_len(d[2])
  times <- .st_orders(times, d[2], "times")
  if (any(diff(times) <= 0) || (mode == "consecutive" && any(diff(times) != 1L))) {
    stop("Time indices must be chronological and consecutive for consecutive windows.")
  }
  ks <- seq_len(.st_one(spatial_max, d[1], "spatial_max"))
  if (is.null(temporal_max)) temporal_max <- length(times)
  ms <- seq_len(.st_one(temporal_max, length(times), "temporal_max"))
  count <- if (mode == "all") choose(length(times), ms) else length(times) - ms + 1L
  if (sum(count) > time_limit) stop("Temporal combination budget exceeded.")
  if (empty == "error" && any(z$empty[, times])) {
    stop("Relative-abundance zeta is undefined for site-time observations with zero total abundance.")
  }
  n <- d[1]
  weights <- outer(seq_len(n), ks, function(j, k) {
    exp(lchoose(n - j, k - 1L) - lchoose(n, k))
  })
  K <- length(ks)
  G <- length(ks) * length(ms)
  present <- abundant <- dominant <- matrix(0, G, d[3], dimnames = list(NULL, dimnames(z$X)[[3]]))
  surface <- expand.grid(spatial_order = ks, temporal_order = ms)
  surface$window <- mode
  for (m in ms) {
    idx <- .st2_time_sets(times, m, mode, time_limit)
    idx[] <- times[idx]
    block <- (m - 1L) * K + ks
    span <- z$panel$times[idx[, m]] - z$panel$times[idx[, 1L]]
    surface[block, c("time_sets", "minimum_span", "mean_span", "maximum_span")] <-
      list(nrow(idx), min(span), mean(span), max(span))
    for (j in seq_len(nrow(idx))) {
      selected <- idx[j, ]
      xp <- colSums(.st2_minimum(z$X, selected))
      dp <- colSums(.st2_minimum(z$D, selected))
      present[block, ] <- present[block, ] + outer(ks, xp, function(k, o) .st_choose_prob(o, k, n))
      dominant[block, ] <- dominant[block, ] + outer(ks, dp, function(k, o) .st_choose_prob(o, k, n))
      if (any(z$empty[, selected])) {
        abundant[block, ] <- NA_real_
      } else if (!is.na(abundant[block[1L], 1L])) {
        q <- .st2_minimum(z$Q, selected)
        ordered <- matrix(apply(q, 2L, sort), n, d[3])
        abundant[block, ] <- abundant[block, ] + crossprod(weights, ordered)
      }
    }
    present[block, ] <- present[block, ] / nrow(idx)
    dominant[block, ] <- dominant[block, ] / nrow(idx)
    abundant[block, ] <- abundant[block, ] / nrow(idx)
  }
  contrib <- list(presence = present, abundance = abundant, dominance = dominant)
  long <- lapply(names(contrib), function(type) {
    out <- surface
    out$response <- type
    out$zeta <- rowSums(contrib[[type]])
    out$effective_species <- .st_ratio(out$zeta^2, rowSums(contrib[[type]]^2))
    st_zeta_retention(out)
  })
  list(surface = do.call(rbind, long), contributions = contrib,
       conditional_minimum = .st_ratio(abundant, present),
       weights = lapply(contrib, .st2_weights))
}

st2_zeta_rectangle <- function(x, sites, times) {
  z <- st2_community(x)
  sites <- .st_orders(sites, dim(z$X)[1], "sites")
  times <- .st_orders(times, dim(z$X)[2], "times")
  if (any(diff(times) <= 0)) stop("Time indices must be chronological.")
  a <- z$X[sites, times, , drop = FALSE]
  b <- z$Q[sites, times, , drop = FALSE]
  c <- z$D[sites, times, , drop = FALSE]
  data.frame(species = dimnames(z$X)[[3]], presence = apply(a, 3L, all),
             abundance = apply(b, 3L, min), dominance = apply(c, 3L, all),
             row.names = NULL)
}
