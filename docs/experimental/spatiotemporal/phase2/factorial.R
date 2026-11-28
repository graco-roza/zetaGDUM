# True group labels are used only by the generator and for scoring.
st2_factorial_groups <- function() expand.grid(occupancy = c("low", "high"),
                                               autocorrelation = c("low", "high"),
                                               abundance = c("low", "high"),
                                               stringsAsFactors = FALSE)

st2_factorial_world <- function(seed, scenario = "null", target = 1L, n = 60L, T = 12L,
                               group_sizes = rep(12L, 8L), region = FALSE, count = FALSE) {
  set.seed(seed)
  groups <- st2_factorial_groups()
  if (scenario == "unequal") group_sizes <- c(20, 16, 15, 13, 10, 9, 7, 6)
  if (scenario == "correlated") group_sizes <- c(36, 4, 4, 4, 4, 4, 4, 36)
  group <- rep(seq_len(8L), group_sizes)
  S <- length(group)
  species <- paste0("sp", seq_len(S))
  rho <- c(.15, .75)[match(groups$autocorrelation[group], c("low", "high"))]
  pi <- c(.25, .75)[match(groups$occupancy[group], c("low", "high"))]
  mean_abundance <- c(2, 16)[match(groups$abundance[group], c("low", "high"))]
  labels <- list(site = paste0("site", seq_len(n)), time = as.character(seq_len(T)), species = species)
  e <- array(rnorm(n * T * 2L), c(n, T, 2L), list(labels$site, labels$time, c("climate", "soil")))
  for (t in seq.int(2L,T)) e[, t, ] <- .7 * e[, t - 1L, ] + sqrt(1 - .7^2) * e[, t, ]
  blocks <- rep(seq_len(6L), length.out = n)
  if (region) e[, , 1L] <- e[, , 1L] + (blocks - 3.5) * .35
  if (scenario == "trend") e[, , 1L] <- sweep(e[, , 1L], 2L, seq(0, 1.8, length.out = T), "+")
  shared <- rep(0, T)
  regional <- matrix(0, n, T)
  if (scenario == "shock") {
    e[, , 1L] <- sweep(e[, , 1L], 2L, rnorm(T, 0, 1.4), "+")
    shared <- rnorm(T, 0, 1.6)
    regional <- matrix(rnorm(6L*T, 0, .8), 6L)[blocks, , drop = FALSE]
  }
  state <- array(runif(n * T * S), c(n, T, S))
  error <- array(rnorm(n * T * S), c(n, T, S))
  selected <- as.numeric(group %in% target)
  if (scenario %in% c("additive", "interaction", "nonlinear", "trend")) selected[] <- 1
  if (scenario %in% c("multiple", "unequal", "correlated")) selected <- as.numeric(group %in% c(1, 4, 6, 7))
  if (scenario == "opposing") selected <- ifelse(group %in% c(1, 4, 6, 7), 1, -1)
  terms <- switch(scenario, gain = c(1.2, 0, 0), loss = c(0, -1.2, 0), abundance = c(0, 0, .65),
                  additive = c(.8, -.7, .35), multiple = c(.8, -.7, .35),
                  interaction = c(1.2, -1.2, 0), nonlinear = c(.9, -.8, .4),
                  opposing = c(1.2, -1.2, .65), unequal = c(.8, -.7, .35),
                  correlated = c(.8, -.7, .35), trend = c(.8, -.7, .35), c(0,0,0))
  generator <- list(e = e, state = state, error = error, labels = labels, pi = pi, rho = rho,
                    mean_abundance = mean_abundance, shared = shared, regional = regional,
                    scenario = scenario, selected = selected, terms = terms, count = count)
  enabled <- st2_factorial_replay(generator, TRUE)
  disabled <- st2_factorial_replay(generator, FALSE)
  list(enabled = enabled, disabled = disabled, env = e, blocks = blocks, group = group,
       species = data.frame(species, group, pi, rho, mean_abundance), generator = generator)
}

st2_factorial_replay <- function(world, enabled, initial = NULL, origin = 1L, final = dim(world$e)[2L]) {
  n <- dim(world$e)[1]
  S <- length(world$pi)
  a <- array(0, c(n, final - origin + 1L, S),
             list(world$labels$site, world$labels$time[origin:final], world$labels$species))
  initial_gain <- qlogis(world$pi * (1-world$rho))
  initial_loss <- qlogis((1-world$pi) * (1-world$rho))
  effect <- function(t) {
    dx <- if (t == 1L) rep(0, n) else world$e[, t, 1L] - world$e[, t-1L, 1L]
    if (t > 1L && world$scenario == "interaction") dx <- world$e[, t-1L, 1L] * dx
    if (world$scenario == "nonlinear") dx <- dx^2 - .6
    if (t == 1L) dx[] <- 0
    outer(dx, world$selected) * enabled
  }
  amounts <- function(t, term) {
    eta <- matrix(log(world$mean_abundance), n, S, byrow = TRUE) + world$terms[3L] * term
    if (!world$count) return(exp(eta - .5^2/2 + .5*world$error[, t, ]))
    mu <- exp(eta)
    q0 <- dnbinom(0, size = 2, mu = mu)
    pmax(1, qnbinom(q0 + pnorm(world$error[,t, ]) * (1-q0), size = 2, mu = mu))
  }
  if (is.null(initial)) {
    if (origin != 1L) stop("Supply an observed starting assemblage after the first survey.")
    a[, 1L, ] <- (world$state[, 1L, ] < rep(world$pi, each = n)) * amounts(1L, effect(1L))
  } else {
    a[, 1L, ] <- initial
  }
  if (final > origin) for (t in seq.int(origin + 1L, final)) {
    term <- effect(t)
    shock <- world$shared[t] + world$regional[, t]
    g <- plogis(matrix(initial_gain, n, S, byrow = TRUE) + world$terms[1L] * term + shock)
    l <- plogis(matrix(initial_loss, n, S, byrow = TRUE) + world$terms[2L] * term - shock)
    previous <- a[, t - origin, ] > 0
    a[, t - origin + 1L, ] <- (world$state[, t, ] < ifelse(previous, 1 - l, g)) * amounts(t, term)
  }
  st_panel(a, origin:final)
}

st2_subset_panel <- function(p, sites, times = seq_along(p$times)) {
  st_panel(p$x[sites, times, , drop = FALSE], p$times[times])
}

st2_draw_sets <- function(n, k, budget) {
  if (choose(n, k) <= budget) return(.st_indices(n, k, budget))
  t(matrix(vapply(seq_len(budget), function(i) sort(sample.int(n, k)), integer(k)), k, budget))
}

st2_disable_change <- function(d, which = "climate") {
  disabled <- grepl(paste0("^(position_change|position_by_change|deviation_change|spread_change)\\.",
                            which, "(_f[0-9]+)?$"), colnames(d$X))
  d$X[, disabled] <- 0
  d
}
