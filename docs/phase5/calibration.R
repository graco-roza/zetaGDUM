source("docs/phase5/common.R")
Rcpp::sourceCpp("docs/manuscript_validation/calibration_helpers.cpp",
  cacheDir = Sys.getenv("CALIBRATION_CPP_CACHE", tempdir()))
args <- commandArgs(TRUE)
n <- as.integer(args[1])
signal <- match.arg(args[2], c("null", "weak", "strong", "spatial", "prior"))
R <- if (length(args) > 2) as.integer(args[3]) else 30L
stopifnot(n %in% c(70, 100, 150, 250))
orders <- c(3L, round(n * .11))
M <- 900L
B <- 79L
proposals <- if (length(args) > 3) as.integer(args[4]) else 150L
qlev <- c(.50,.80,.90,.95)
refenv <- data.frame(g = seq(-3,3,length.out=81))
ref <- ugdm_design(refenv, cbind(1:80,2:81))
theta <- switch(signal, null = c(-.3,0,0,-.3,0,0), weak = c(-.3,.25,.12,-.3,.15,.25),
  strong = c(-.3,.7,.35,-.3,.45,.7), spatial = c(-.3,.25,.12,-.3,.15,.25),
  prior = NULL)
key <- paste(n, signal, sep = "_")
if (proposals != 150L) key <- paste0("integration_", key, "_", proposals)
output <- diagnostics <- list()

for (r in seq_len(R)) {
  set.seed(517L + 10000L * n + 330L * r)
  truth <- if (is.null(theta)) rnorm(6L, sd = .35) else theta
  env <- data.frame(g = rnorm(n))
  u <- rnorm(n)
  if (signal == "spatial") {
    coords <- cbind(runif(n), runif(n))
    u <- drop(t(chol(exp(-as.matrix(dist(coords)) / .3))) %*% u)
    env$g <- .7 * sin(6 * coords[,1]) + .65 * env$g
  }
  u <- tanh(u)
  z <- (env$g - ref$environment_basis$center) / ref$environment_basis$scale
  D <- abs(outer(z, z, "-"))
  dat <- function(idx, j) {
    k <- ncol(idx)
    X <- cbind(1, internal("meanpair_dist_cpp")(D, idx), rowMeans(matrix(z[idx], nrow(idx))))
    co <- truth[(j - 1) * 3 + 1:3]
    mu <- plogis(drop(X %*% co))
    y <- mu + .85 * pmin(mu, 1 - mu) * rowMeans(matrix(u[idx], nrow(idx)))
    list(X = X, y = y, idx = idx)
  }
  est <- function(d, keep = rep(TRUE, length(d$y)), start = numeric(3)) {
    calibration_fit(d$X, d$y, unname(start), keep, rep(1, length(d$y)))
  }
  long <- sample_sets_coupled(matrix(1, n, 1), orders, M, exact_limit = 0)
  a <- lapply(seq_along(orders), function(j) dat(long[[j]]$idx, j))
  point <- unlist(lapply(a, est))
  t0 <- sum(proc.time()[1:2])
  de <- do.call(cbind, lapply(seq_along(a), function(j) {
    z <- a[[j]]
    t(vapply(seq_len(n), function(i) est(z, rowSums(z$idx == i) == 0L,
      point[(j - 1) * 3 + 1:3]), numeric(3)))
  }))
  deletion_seconds <- sum(proc.time()[1:2]) - t0
  G <- site_weights(n, max(B, 700))
  t0 <- sum(proc.time()[1:2])
  pseudo <- ugdm_pseudovalue_bootstrap(point, de, weights = G)
  raw <- ugdm_pseudovalue_bootstrap(point, de, weights = G, centered = FALSE)
  pseudo_seconds <- sum(proc.time()[1:2]) - t0
  t0 <- sum(proc.time()[1:2])
  cp <- do.call(cbind, lapply(seq_along(a), function(j) {
    t(vapply(seq_len(B), function(l) est(dat(calibration_cp(G[,l], orders[j], M), j)), numeric(3)))
  }))
  cp_seconds <- sum(proc.time()[1:2]) - t0
  fit_data <- lapply(a, function(z) {
    d <- ugdm_model_design(ugdm_design(env, z$idx, ref = ref), "B", env)
    stopifnot(max(abs(d$Xdir[,1] - z$X[,2])) < 1e-10)
    list(design = d, y = z$y)
  })
  results <- list(cp = list(draws = cp, weights = rep(1/B, B), seconds = cp_seconds),
    pseudo = list(draws = pseudo$draws, weights = rep(1/ncol(G), ncol(G)),
                  seconds = pseudo_seconds + deletion_seconds),
    pseudo_raw = list(draws = raw$draws, weights = rep(1/ncol(G), ncol(G)),
                      seconds = pseudo_seconds + deletion_seconds))
  run_jel <- function(sd, type, b = proposals) tryCatch(
    ugdm_jel(fit_data, prior_sd = sd, draws = b, method = type, seed = 999 + r),
    error = function(e) list(reason = conditionMessage(e)))
  results$jel <- run_jel(if (signal == "prior") .35 else 2, "EL")
  if (n == 100L) {
    results$etel <- run_jel(if (signal == "prior") .35 else 2, "ETEL")
    if (signal %in% c("weak","strong")) results$jel_shrunk <- run_jel(.5, "EL")
  }
  C <- rbind(lo_spread = c(0,1,0,0,0,0), lo_position = c(0,0,1,0,0,0),
    hi_spread = c(0,0,0,0,1,0), hi_position = c(0,0,0,0,0,1),
    spread_change = c(0,-1,0,0,1,0), position_change = c(0,0,-1,0,0,1))
  u0 <- drop(C %*% truth)
  centers <- drop(C %*% point)
  ses <- internal(".jackknife_sd")(de %*% t(C))
  results$jackknife <- list(draws = NULL, seconds = deletion_seconds)
  for (method in names(results)) {
    b <- results[[method]]
    mock <- is.null(b$draws) && method != "jackknife"
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(n, signal, replicate = r, method,
      ess = if (is.null(b$ess)) NA_real_ else b$ess,
      rejected = if (is.null(b$rejected)) NA_real_ else b$rejected,
      seconds = if (is.null(b$seconds)) NA_real_ else b$seconds,
      reason = if (mock) b$reason else "")
    for (j in seq_len(nrow(C))) {
      if (method == "jackknife") {
        lower <- centers[j] - qt((1 + qlev) / 2, n - 1) * ses[j]
        upper <- centers[j] + qt((1 + qlev) / 2, n - 1) * ses[j]
        summary <- c(NA, NA, NA, centers[j])
      } else if (mock) {
        lower <- upper <- rep(NA_real_, length(qlev))
        summary <- rep(NA_real_, 4)
      } else {
        z <- drop(b$draws %*% C[j, ])
        lower <- qweight(z, b$weights, (1 - qlev)/2)
        upper <- qweight(z, b$weights, (1 + qlev)/2)
        summary <- c(sum(b$weights[z > 0]), sum(b$weights[z > .1]),
            sum(b$weights[z <= u0[j]]), sum(b$weights * z))
      }
      output[[length(output) + 1L]] <- data.frame(n, signal, replicate = r, method,
        contrast = rownames(C)[j], truth = unname(u0[j]), order_low = orders[1], order_high = orders[2],
        point = unname(centers[j]), posterior_mean = unname(summary[4]), p_positive = unname(summary[1]),
        p_practical = unname(summary[2]), posterior_cdf_truth = unname(summary[3]), level = qlev,
        lo = lower, hi = upper)
    }
  }
}
write5(do.call(rbind, output), paste0("calibration_", key))
write5(do.call(rbind, diagnostics), paste0("diagnostics_", key))
