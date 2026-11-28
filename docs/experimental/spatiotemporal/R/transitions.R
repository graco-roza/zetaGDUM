.st2_event <- function(before, after) {
  list(P = before & after, G = !before & after, L = before & !after)
}

st2_directional <- function(x, sites, before, after) {
  z <- st2_community(x)
  sites <- .st_orders(sites, dim(z$X)[1], "sites")
  before <- .st_one(before, dim(z$X)[2], "before")
  after <- .st_one(after, dim(z$X)[2], "after")
  if (before >= after) stop("before must precede after.")
  matrix_at <- function(a, t) matrix(a[sites, t, ], length(sites), dim(a)[3])
  x0 <- matrix_at(z$X, before)
  x1 <- matrix_at(z$X, after)
  h0 <- colSums(x0) == length(sites)
  h1 <- colSums(x1) == length(sites)
  incidence <- .st2_event(h0, h1)
  d0 <- colSums(matrix_at(z$D, before)) == length(sites)
  d1 <- colSums(matrix_at(z$D, after)) == length(sites)
  dominance <- .st2_event(d0, d1)
  q0 <- apply(matrix_at(z$Q, before), 2L, min)
  q1 <- apply(matrix_at(z$Q, after), 2L, min)
  gain <- pmax(q1 - q0, 0)
  loss <- pmax(q0 - q1, 0)
  raw <- data.frame(species = dimnames(z$X)[[3]], shared_before = h0, shared_after = h1,
                    persistent = incidence$P, gain = incidence$G, loss = incidence$L,
                    min_before = q0, min_after = q1, abundance_persistent = pmin(q0, q1),
                    abundance_gain = gain, abundance_loss = loss,
                    abundance_gain_occ = incidence$G * q1, abundance_loss_occ = incidence$L * q0,
                    abundance_gain_num = incidence$P * gain, abundance_loss_num = incidence$P * loss,
                    dominance_persistent = dominance$P, dominance_gain = dominance$G,
                    dominance_loss = dominance$L, dominance_gain_occ = dominance$G & !h0,
                    dominance_loss_occ = dominance$L & !h1,
                    dominance_gain_num = dominance$G & h0, dominance_loss_num = dominance$L & h1)
  summarise <- function(events) {
    P <- sum(events$P)
    G <- sum(events$G)
    L <- sum(events$L)
    denominator <- 2 * P + G + L
    c(before = P + L, after = P + G, P = P, G = G, L = L, change = G - L,
      D_gain = .st_ratio(G, denominator), D_loss = .st_ratio(L, denominator),
      D_temporal = .st_ratio(G + L, denominator))
  }
  amount <- c(min_before = sum(q0), min_after = sum(q1), change = sum(q1) - sum(q0),
              colSums(raw[, grep("^abundance_", names(raw)), drop = FALSE]))
  list(species = raw, sites = dimnames(z$X)[[1]][sites], times = z$panel$times[c(before, after)],
       presence = summarise(incidence), abundance = amount,
       dominance = c(summarise(dominance), colSums(raw[, grep("^dominance_.*_(occ|num)$", names(raw))])))
}

st2_occupancy <- function(x, sites = NULL, before = 1L, after = dim(.st_panel(x)$x)[2L],
                          orders = NULL) {
  p <- .st_panel(x)
  n <- dim(p$x)[1]
  if (is.null(sites)) sites <- seq_len(n)
  if (is.null(orders)) orders <- seq_len(n)
  sites <- .st_orders(sites, n, "sites")
  regional <- st_occupancy(p, before, after, orders)
  q <- regional$species
  a <- .st_slice(p, before)[sites, , drop = FALSE] > 0
  b <- .st_slice(p, after)[sites, , drop = FALSE] > 0
  g <- colSums(!a & b)
  l <- colSums(a & !b)
  list(combination = data.frame(species = colnames(a), r0 = colSums(a), r1 = colSums(b),
                                gains = g, losses = l, expansion = pmax(g - l, 0),
                                contraction = pmax(l - g, 0), redistribution = pmin(g, l),
                                pure_redistribution = g == l & g > 0,
                                regional_appearance = q$regional_appearance,
                                regional_disappearance = q$regional_disappearance),
       sites = rownames(a), regional = regional)
}

st2_poisson_binomial <- function(probability) {
  if (any(!is.finite(probability) | probability < 0 | probability > 1)) {
    stop("Bernoulli probabilities must be finite and in [0,1].")
  }
  mass <- 1
  for (p in probability) mass <- c(mass * (1 - p), 0) + c(0, mass * p)
  mass
}

st2_occupancy_probability <- function(gain, loss) {
  g <- st2_poisson_binomial(gain)
  l <- st2_poisson_binomial(loss)
  a <- seq_along(g) - 1L
  b <- seq_along(l) - 1L
  mass <- outer(g, l)
  delta <- outer(a, b, "-")
  list(gain_mass = g, loss_mass = l, joint_mass = mass,
       expansion = sum(mass[delta > 0]), contraction = sum(mass[delta < 0]),
       pure_redistribution = sum(mass[delta == 0]) - mass[1L, 1L],
       redistribution = sum(mass * outer(a, b, pmin)),
       mean_expansion = sum(mass * pmax(delta, 0)),
       mean_contraction = sum(mass * pmax(-delta, 0)))
}
