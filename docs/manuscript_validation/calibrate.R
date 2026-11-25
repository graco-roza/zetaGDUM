# Rscript docs/manuscript_validation/calibrate.R OUT N SIGNAL [REPS=32] [B=79] [M=10000] [DEPENDENCE=iid]
source("docs/manuscript_validation/common.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Install Rcpp to run calibration.")
cache <- Sys.getenv("CALIBRATION_CPP_CACHE", tempdir())
Rcpp::sourceCpp("docs/manuscript_validation/calibration_helpers.cpp", cacheDir = cache)
args <- commandArgs(TRUE)
out <- args[1L]
n <- as.integer(args[2L])
signal <- args[3L]
reps <- if (length(args) > 3L) as.integer(args[4L]) else 32L
B <- if (length(args) > 4L) as.integer(args[5L]) else 79L
M <- if (length(args) > 5L) as.integer(args[6L]) else 10000L
dependence <- if (length(args) > 6L) match.arg(args[7L], c("iid", "spatial")) else "iid"
stopifnot(n %in% c(70L, 100L, 150L, 250L), signal %in% c("null", "weak", "strong"))
order_values <- pmax(2L, c(floor(.02 * n + .5), floor(.05 * n + .5), floor(.10 * n + .5), floor(.20 * n + .5)))
if (n >= 150L) order_values <- unique(c(2L, order_values))
P <- matrix(1L, n, 1L)
anchor <- data.frame(g = seq(-3, 3, length.out = 31), c = 2 * cos(seq(0, 7, length.out = 31)))
ref <- suppressWarnings(ugdm_design(anchor, sample_sets(nrow(anchor), 2, M = 100, seed = 109)))
map <- rbind(spread.gradient = c(0, 1, 1, 1, 0, 0, 0, 0, 0),
             spread.control = c(0, 0, 0, 0, 1, 1, 1, 0, 0),
             position.gradient = c(0, 0, 0, 0, 0, 0, 0, 1, 0),
             position.control = c(0, 0, 0, 0, 0, 0, 0, 0, 1))
level <- c(null = 0, weak = .30, strong = 1)[[signal]]
cpu <- function() sum(proc.time()[1:2])
known <- function(k) c(-.3, level * (1 - .5 * k / n) * c(.45, .4, .35),
                        c(0, 0, 0), level * .09, 0)
model_data <- function(idx, beta, env, u) {
  d <- ugdm_design(env, idx, ref = ref)
  X <- cbind(1, d$Xdir, d$Xpos)
  mu <- plogis(drop(X %*% beta))
  y <- mu + .85 * pmin(mu, 1 - mu) * rowMeans(matrix(u[idx], nrow(idx), ncol(idx)))
  list(X = X, y = y, idx = idx, design = d)
}
estimate <- function(a, start = numeric(ncol(a$X)), keep = rep(TRUE, length(a$y)),
                      w = rep(1, length(a$y))) {
  calibration_fit(a$X, a$y, unname(start), keep, w)
}
as_effect <- function(beta) setNames(drop(map %*% beta), rownames(map))
compare_helpers <- function() {
  set.seed(1231)
  omega <- c(.3, 1, .7, 4, 3, 1.2)
  for (k in c(2L, 4L)) {
    sets <- t(combn(6L, k))
    key <- function(x) apply(x, 1L, paste, collapse = "-")
    target <- apply(sets, 1L, function(i) prod(omega[i]))
    target <- target / sum(target)
    simulated <- calibration_cp(omega, k, 50000L)
    seen <- tabulate(match(key(simulated), key(sets)), nrow(sets)) / nrow(simulated)
    stopifnot(max(abs(seen - target)) < .008)
  }
  weights <- exp(seq(-1, 1, length.out = 30))
  pairs <- t(combn(30, 2))
  expected <- apply(pairs, 1L, function(i) prod(weights[i]))
  expected <- expected / sum(expected)
  set.seed(1232)
  observed <- calibration_cp(weights, 2, 80000)
  key <- function(x) x[, 1L] * 32L + x[, 2L]
  proportion <- tabulate(match(key(observed), key(pairs)), nrow(pairs)) / nrow(observed)
  stopifnot(max(abs(proportion - expected)) < .0015)
  n0 <- 28L
  env <- data.frame(g = rnorm(n0), c = rnorm(n0))
  a <- model_data(sample_sets(n0, 4, 1200L), known(4), env, tanh(rnorm(n0)))
  b <- estimate(a)
  fit <- ugdm_fit(a$design, a$y, rep(1, length(a$y)), free_sign = TRUE)
  stopifnot(max(abs(b - c(fit$intercept, fit$beta, fit$lambda))) < 2e-4)
  for (i in c(1L, 8L, 17L)) {
    keep <- rowSums(a$idx == i) == 0
    d <- subset_design(a$design, keep)
    refit <- ugdm_fit(d, a$y[keep], rep(1, sum(keep)), free_sign = TRUE)
    stopifnot(max(abs(estimate(a, b, keep) - c(refit$intercept, refit$beta, refit$lambda))) < 2e-4)
  }
}
if (n == 70L && signal == "null" && dependence == "iid") compare_helpers()
rows <- diagnostics <- mc_rows <- changes <- list()
for (r in seq_len(reps)) {
  message(n, " ", signal, " community ", r, "/", reps)
  stamp <- 1000000L + 1000000L * (r + 100L * (match(signal, c("null", "weak", "strong")) - 1L) +
    300L * (match(n, c(70L, 100L, 150L, 250L)) - 1L))
  if (dependence == "spatial") stamp <- stamp - 100000000 - 50000000
  set.seed(stamp)
  env <- data.frame(g = rnorm(n), c = rnorm(n))
  u <- rnorm(n)
  if (dependence == "spatial") {
    coords <- matrix(runif(2L * n), n, 2L)
    env$g <- .65 * sin(2 * pi * coords[, 1L]) * sqrt(2) + sqrt(1 - .65^2) * env$g
    u <- drop(t(chol(exp(-as.matrix(dist(coords)) / .24))) %*% u)
  }
  u <- tanh(u)
  W <- site_weights(n, B, seed = stamp + 200000L)
  set.seed(stamp + 300000L)
  G <- matrix(rnorm(399L * n), 399L, n) / sqrt(n * (n - 1L))
  fits <- list()
  for (k in order_values) {
    beta <- known(k)
    truth <- as_effect(beta)
    s <- draw(P, k, M, 0, seed = stamp + 1000L * k)
    a <- model_data(s$idx, beta, env, u)
    f <- estimate(a)
    effect <- as_effect(f)
    mu <- plogis(drop(a$X %*% f))
    h <- crossprod(a$X * sqrt(mu * (1 - mu))) / nrow(a$X)
    ev <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
    exc <- omission(a$idx, s$sampling_weights, n)
    t0 <- cpu()
    delete <- t(vapply(exc$keep, function(keep) as_effect(estimate(a, f, keep)), numeric(length(effect))))
    jack <- cpu() - t0
    t0 <- cpu()
    score <- a$X * as.numeric(a$y - mu) / length(mu)
    ins <- inclusion(a$idx, cbind(rep(1 / length(mu), length(mu)), score), n)
    mass <- ins[, 1L]
    pseudovalue <- (n - 1L) * (ins[, -1L, drop = FALSE] - outer(mass, colSums(score))) / (1 - mass)
    infl <- sweep(pseudovalue, 2L, colMeans(pseudovalue), "-") %*% solve(h, t(map))
    mul <- sweep(G %*% infl, 2L, effect, "+")
    mult <- cpu() - t0
    t0 <- cpu()
    bb <- t(vapply(seq_len(B), function(j) {
      if (s$exact) {
        w <- internal(".product_sampling_weights")(W[, j], a$idx, s$sampling_weights)
        as_effect(estimate(a, f, w = w))
      } else {
        set.seed(stamp + 500000L + k * 1200L + j)
        z <- model_data(calibration_cp(W[, j], k, M), beta, env, u)
        as_effect(estimate(z, f))
      }
    }, numeric(length(effect))))
    bayes <- cpu() - t0
    methods <- list(bayesian = bb, jackknife = delete, multiplier = mul)
    fits[[as.character(k)]] <- list(point = effect, truth = truth, methods = methods)
    for (method in names(methods)) {
      i <- interval(effect, methods[[method]], method)
      rows[[length(rows) + 1L]] <- data.frame(n, signal, r, order = k, ratio = k / n,
        method, i, truth = unname(truth))
    }
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(n, signal, r, order = k,
      budget = M, exact = s$exact, converged = all(is.finite(f)),
      rank = sum(ev > max(ev) * ncol(a$X) * .Machine$double.eps),
      parameters = ncol(a$X), lambda_min = min(ev), condition = max(ev) / min(ev),
      omission_ess = min(exc$ess), jack_cpu = jack, multiplier_cpu = mult, bb_cpu = bayes)
    # Independent numerical replicates on the same known community at manuscript budgets.
    if (r == 1L && signal == "weak" && dependence == "iid") for (m in c(10000L, 50000L)) {
      for (j in seq_len(8L)) {
        pool <- draw(P, k, m, 0, seed = stamp + 700000L + m + k * 200L + j)
        value <- as_effect(estimate(model_data(pool$idx, beta, env, u), f))
        mc_rows[[length(mc_rows) + 1L]] <- data.frame(n, signal, order = k, budget = m,
          replicate = j, term = names(value), estimate = unname(value), exact = pool$exact)
      }
    }
  }
  lo <- fits[[1L]]
  for (k in order_values[-1L]) {
    hi <- fits[[as.character(k)]]
    for (method in names(hi$methods)) {
      i <- interval(hi$point - lo$point, hi$methods[[method]] - lo$methods[[method]], method)
      changes[[length(changes) + 1L]] <- data.frame(n, signal, r, from = order_values[1L], order = k,
        method, i, truth = unname(hi$truth - lo$truth))
    }
  }
}
for (name in c("rows", "diagnostics", "changes", "mc_rows")) {
  v <- get(name)
  if (length(v)) write_study(do.call(rbind, v), out, paste(if (dependence == "iid") "cal" else "spatial", n, signal, name, sep = "_"))
}
