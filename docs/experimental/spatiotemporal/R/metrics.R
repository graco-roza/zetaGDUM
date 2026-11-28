# Arrays are sites x ordered survey times x species. Keep the declared species
# pool, including columns that are zero throughout the panel.
st_panel <- function(x, times = seq_len(dim(x)[2]), missing = c("error", "complete_sites")) {
  missing <- match.arg(missing)
  d <- dim(x)
  if (length(d) != 3L || any(d == 0L) || !(is.numeric(x) || is.logical(x)) ||
      any(is.nan(x) | is.infinite(x) | x < 0, na.rm = TRUE)) {
    stop("Supply a non-negative sites x times x species array (NA for missing).")
  }
  if (!is.numeric(times) || length(times) != d[2] || any(!is.finite(times)) ||
      any(diff(times) <= 0)) stop("Survey times must be finite and strictly increasing.")
  labels <- dimnames(x)
  if (is.null(labels)) labels <- vector("list", 3)
  if (is.null(labels[[1]])) labels[[1]] <- paste0("site", seq_len(d[1]))
  if (is.null(labels[[2]])) labels[[2]] <- as.character(times)
  if (is.null(labels[[3]])) labels[[3]] <- paste0("species", seq_len(d[3]))
  if (any(vapply(labels, function(z) anyNA(z) || any(!nzchar(z)) || anyDuplicated(z) > 0,
                 logical(1)))) stop("Panel labels must be unique and non-empty within each axis.")
  dimnames(x) <- labels
  keep <- !apply(is.na(x), 1L, any)
  if (!all(keep) && missing == "error") stop("Balanced panel required. Use complete_sites explicitly.")
  if (!any(keep)) stop("No sites have complete observations across the panel.")
  structure(list(x = x[keep, , , drop = FALSE], times = times,
                 dropped_sites = labels[[1]][!keep]), class = "st_panel")
}

.st_panel <- function(x) if (inherits(x, "st_panel")) x else st_panel(x)

.st_orders <- function(x, max, arg) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x) | x < 1 | x > max | x != floor(x)) ||
      anyDuplicated(x)) stop(arg, " must contain distinct integer indices within the axis.")
  as.integer(x)
}

.st_one <- function(x, max, arg) {
  if (length(x) != 1L) stop(arg, " must be scalar.")
  .st_orders(x, max, arg)
}

.st_indices <- function(n, k, limit) {
  if (choose(n, k) > limit) stop("Combination budget exceeded. Reduce orders or supply explicit sets.")
  t(utils::combn(n, k))
}

.st_sets <- function(n, k, sets, limit) {
  if (is.null(sets)) return(.st_indices(n, k, limit))
  if (!is.matrix(sets) || !is.numeric(sets) || !nrow(sets) || ncol(sets) != k ||
      any(!is.finite(sets) | sets < 1 | sets > n | sets != floor(sets))) stop("Invalid sets matrix.")
  for (i in seq_len(nrow(sets))) .st_orders(sets[i, ], n, "sites in each set")
  sets
}

.st_ratio <- function(a, b) {
  result <- a / b
  result[!rep_len(!is.na(b) & b > 0, length(result))] <- NA_real_
  result
}
.st_choose_prob <- function(o, k, n) exp(lchoose(o, k) - lchoose(n, k))

.st_slice <- function(panel, time) {
  x <- panel$x
  matrix(x[, time, ], dim(x)[1], dim(x)[3], dimnames = dimnames(x)[c(1, 3)])
}

.st_history <- function(panel, sites, times) {
  x <- panel$x
  vapply(times, function(j) colSums(x[sites, j, , drop = FALSE] > 0,
                                  dims = 2L) == length(sites), logical(dim(x)[3]))
}

st_zeta <- function(x, spatial_order, temporal_order, sets = NULL, time_sets = NULL,
                    limit = 1e5) {
  p <- .st_panel(x)
  d <- dim(p$x)
  k <- .st_one(spatial_order, d[1], "spatial_order")
  m <- .st_one(temporal_order, d[2], "temporal_order")
  sets <- .st_sets(d[1], k, sets, limit)
  time_sets <- .st_sets(d[2], m, time_sets, limit)
  if (nrow(sets) * as.double(nrow(time_sets)) > limit) stop("Rectangle budget exceeded.")
  rows <- expand.grid(site_set = seq_len(nrow(sets)), time_set = seq_len(nrow(time_sets)))
  shared <- lapply(seq_len(nrow(rows)), function(r) {
    a <- p$x[sets[rows$site_set[r], ], time_sets[rows$time_set[r], ], , drop = FALSE]
    dimnames(p$x)[[3]][colSums(a > 0, dims = 2L) == k * m]
  })
  rows$spatial_order <- k
  rows$temporal_order <- m
  rows$zeta <- lengths(shared)
  rows$sites <- lapply(rows$site_set, function(i) dimnames(p$x)[[1]][sets[i, ]])
  rows$times <- lapply(rows$time_set, function(i) p$times[sort(time_sets[i, ])])
  rows$species <- shared
  rows
}

# Exact mean via the occupancy of each species after intersecting over times.
# Spatial combinations are never enumerated here.
st_zeta_decline <- function(x, spatial_max = dim(.st_panel(x)$x)[1],
                            temporal_max = dim(.st_panel(x)$x)[2], time_limit = 1e5) {
  p <- .st_panel(x)
  d <- dim(p$x)
  ks <- seq_len(.st_one(spatial_max, d[1], "spatial_max"))
  ms <- seq_len(.st_one(temporal_max, d[2], "temporal_max"))
  if (sum(choose(d[2], ms)) > time_limit) stop("Temporal combination budget exceeded.")
  estimates <- lapply(ms, function(m) {
    times <- .st_indices(d[2], m, time_limit)
    sums <- matrix(0, length(ks), d[3])
    for (i in seq_len(nrow(times))) {
      # Collapse the time axis within each site and species.
      shared_sites <- .st_slice(p, times[i, 1L]) > 0
      for (j in times[i, -1L]) shared_sites <- shared_sites & (.st_slice(p, j) > 0)
      occupancy <- colSums(shared_sites)
      for (k in ks) sums[k, ] <- sums[k, ] + .st_choose_prob(occupancy, k, d[1])
    }
    sums / nrow(times)
  })
  columns <- do.call(rbind, estimates)
  colnames(columns) <- dimnames(p$x)[[3]]
  surface <- expand.grid(spatial_order = ks, temporal_order = ms)
  surface$zeta <- rowSums(columns)
  surface$species_supported <- rowSums(columns > 0)
  surface$effective_species <- .st_ratio(surface$zeta^2, rowSums(columns^2))
  surface$rectangles <- choose(d[1], surface$spatial_order) * choose(d[2], surface$temporal_order)
  attr(surface, "species_contributions") <- columns
  surface
}

st_zeta_retention <- function(surface) {
  k <- surface$spatial_order
  m <- surface$temporal_order
  index <- paste(k, m)
  surface$spatial_retention <- .st_ratio(surface$zeta, surface$zeta[match(paste(k - 1L, m), index)])
  surface$temporal_retention <- .st_ratio(surface$zeta, surface$zeta[match(paste(k, m - 1L), index)])
  surface
}

st_decompose <- function(x, spatial_order = 1L, times = NULL, sets = NULL, limit = 1e5) {
  p <- .st_panel(x)
  d <- dim(p$x)
  if (is.null(times)) times <- seq_len(d[2])
  times <- .st_orders(times, d[2], "times")
  if (length(times) < 2L || any(diff(times) <= 0)) stop("Provide at least two increasing time indices.")
  k <- .st_one(spatial_order, d[1], "spatial_order")
  idx <- .st_sets(d[1], k, sets, limit)
  species <- dimnames(p$x)[[3]]
  first <- times[1]
  last <- tail(times, 1L)
  histories <- lapply(seq_len(nrow(idx)), function(i) {
    h <- matrix(.st_history(p, idx[i, ], times), d[3], length(times),
                dimnames = list(species, as.character(p$times[times])))
    h
  })
  events <- lapply(histories, function(h) {
    before <- h[, 1L]
    after <- h[, ncol(h)]
    full <- rowSums(h) == ncol(h)
    prior <- h[, -ncol(h), drop = FALSE]
    next_survey <- h[, -1L, drop = FALSE]
    list(persistent = species[before & after], gains = species[!before & after],
         losses = species[before & !after], continuous = species[full],
         interrupted = species[before & after & !full],
         intervening_only = species[!before & !after & rowSums(h) > 0],
         disappeared_during = species[rowSums(prior & !next_survey) > 0],
         gained_during = species[rowSums(!prior & next_survey) > 0],
         recolonised = species[apply(h, 1L, function(z) any(diff(which(z)) > 1L))])
  })
  out <- data.frame(site_set = seq_len(nrow(idx)), spatial_order = k,
                    start = p$times[first], end = p$times[last])
  for (label in names(events[[1]])) {
    members <- lapply(events, `[[`, label)
    out[[paste0(label, "_species")]] <- members
    out[[label]] <- lengths(members)
  }
  out$before <- out$persistent + out$losses
  out$after <- out$persistent + out$gains
  out$change <- out$gains - out$losses
  out$sites <- lapply(seq_len(nrow(idx)), function(i) dimnames(p$x)[[1]][idx[i, ]])
  out$history <- histories
  out
}

st_occupancy <- function(x, before = 1L, after = dim(.st_panel(x)$x)[2],
                         orders = seq_len(dim(.st_panel(x)$x)[1])) {
  p <- .st_panel(x)
  d <- dim(p$x)
  before <- .st_one(before, d[2], "before")
  after <- .st_one(after, d[2], "after")
  if (before >= after) stop("before must precede after.")
  orders <- .st_orders(orders, d[1], "orders")
  a <- .st_slice(p, before)
  b <- .st_slice(p, after)
  old <- a > 0
  new <- b > 0
  keep <- old & new
  o0 <- colSums(old)
  o1 <- colSums(new)
  both <- colSums(keep)
  t0 <- colSums(a)
  t1 <- colSums(b)
  m0 <- .st_ratio(t0, o0)
  m1 <- .st_ratio(t1, o1)
  series <- data.frame(species = colnames(a), occupancy_before = o0, occupancy_after = o1,
                       occupancy_change = o1 - o0, persistent = both, gained = o1 - both,
                       lost = o0 - both, expansion = o1 > o0, contraction = o1 < o0,
                       redistribution = o1 > both & o0 > both,
                       replacement_sites = pmin(o1 - both, o0 - both),
                       redistribution_only = o1 == o0 & o1 > both,
                       regional_appearance = o0 == 0 & o1 > 0,
                       regional_disappearance = o1 == 0 & o0 > 0,
                       total_before = t0, total_after = t1, total_change = t1 - t0,
                       conditional_before = m0, conditional_after = m1, conditional_change = m1 - m0,
                       dominance_before = .st_ratio(t0, sum(t0)), dominance_after = .st_ratio(t1, sum(t1)),
                       persistent_conditional_change = .st_ratio(colSums((b - a) * keep), both),
                       abundance_gained_sites = colSums(b * !old),
                       abundance_lost_sites = colSums(a * !new),
                       abundance_change_persistent_sites = colSums((b - a) * keep))
  series$dominance_change <- series$dominance_after - series$dominance_before
  sites <- rownames(a)
  members <- function(mask) lapply(seq_len(d[3]), function(j) sites[mask[, j]])
  series$new_sites <- members(!old & new)
  series$lost_sites <- members(old & !new)
  series$persistent_sites <- members(keep)
  contributions <- do.call(rbind, lapply(orders, function(k) {
    start <- .st_choose_prob(o0, k, d[1])
    end <- .st_choose_prob(o1, k, d[1])
    shared <- .st_choose_prob(both, k, d[1])
    data.frame(species = colnames(a), spatial_order = k, before = start, after = end,
               persistent = shared, gains = end - shared, losses = start - shared,
               change = end - start, row.names = NULL)
  }))
  list(species = series, contributions = contributions, times = p$times[c(before, after)])
}

# Delegates the spatial pairwise definitions and empty-pair conventions to MS-GDUM.
st_beta <- function(x) {
  p <- .st_panel(x)
  do.call(rbind, lapply(seq_along(p$times), function(t) {
    b <- zetaGDUM::pair_partition(.st_slice(p, t))
    data.frame(time = p$times[t], site_a = dimnames(p$x)[[1]][b$i],
               site_b = dimnames(p$x)[[1]][b$j], b)
  }))
}

# Established multisite Sørensen partition of Baselga (2010). This is a
# site-set functional based on pairwise exclusives and regional redundancy.
st_beta_multisite <- function(x, sites = NULL, time = 1L) {
  p <- .st_panel(x)
  d <- dim(p$x)
  if (is.null(sites)) sites <- seq_len(d[1])
  sites <- .st_orders(sites, d[1], "sites")
  time <- .st_one(time, d[2], "time")
  if (length(sites) < 2L) stop("Beta diversity requires at least two sites.")
  a <- .st_slice(p, time)[sites, , drop = FALSE] > 0
  pairs <- zetaGDUM::pair_partition(a * 1L)
  redundancy <- sum(a) - sum(colSums(a) > 0)
  lo <- sum(pmin(pairs$b, pairs$c))
  hi <- sum(pmax(pairs$b, pairs$c))
  total <- .st_ratio(lo + hi, 2 * redundancy + lo + hi)
  replacement <- .st_ratio(lo, redundancy + lo)
  c(sorensen = total, simpson = replacement, nestedness = total - replacement)
}

# Two-event Shapley accounting for the observed pairwise changes. The
# counterfactuals flip observed gains or losses alone and have no causal status.
st_beta_change <- function(x, before = 1L, after = dim(.st_panel(x)$x)[2]) {
  p <- .st_panel(x)
  d <- dim(p$x)
  before <- .st_one(before, d[2], "before")
  after <- .st_one(after, d[2], "after")
  if (before >= after) stop("before must precede after.")
  a <- .st_slice(p, before) > 0
  b <- .st_slice(p, after) > 0
  parts <- function(y) zetaGDUM::pair_partition(y * 1L)
  start <- parts(a)
  end <- parts(b)
  gain <- parts(a | b)
  loss <- parts(a & b)
  do.call(rbind, lapply(c("sorensen", "simpson", "nestedness"), function(term) {
    s <- start[[term]]
    e <- end[[term]]
    g <- gain[[term]]
    l <- loss[[term]]
    data.frame(site_a = rownames(a)[start$i], site_b = rownames(a)[start$j],
               component = term, before = s, after = e, change = e - s,
               gain_contribution = ((g - s) + (e - l)) / 2,
               loss_contribution = ((l - s) + (e - g)) / 2)
  }))
}
