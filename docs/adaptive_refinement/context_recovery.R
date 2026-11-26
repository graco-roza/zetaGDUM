source("docs/adaptive_refinement/ecology_common.R")
args <- commandArgs(TRUE)
scenario <- match.arg(args[1L], c("linear", "redundant", "positive", "negative", "collinear",
  "spatial", "sparse"))
B <- if (length(args) >= 2L) as.integer(args[2L]) else 16L
n <- switch(scenario, collinear = 150L, spatial = 250L, sparse = 150L, 70L)
orders <- c(2L, 6L, if (n <= 70L) 10L else 14L)
DGP <- function(r) {
  set.seed(715083L + 229L * r + n)
  s <- sort(runif(n))
  z <- rnorm(n)
  shift <- 0
  if (scenario == "spatial") {
    R <- exp(-abs(outer(s, s, "-")) / .08) + diag(1e-7, n)
    z <- rt(n, 11) * .30 + seq(-1.3, 1.3, length.out = n)
    shift <- .38 * drop(t(chol(R)) %*% rnorm(n))
  }
  env <- data.frame(climate = pmax(-2.7, pmin(2.7, z)))
  if (scenario == "collinear") env$surrogate <- env$climate + rnorm(n, sd = .06)
  if (scenario == "sparse") {
    opt <- seq(-2.8, 2.8, length.out = 80L)
    rate <- exp(-.5 * (outer(env$climate, opt, "-") / .23)^2)
    x <- matrix(rpois(length(rate), 4 * rate + .02), n)
    x[rowSums(x) == 0, 1L] <- 1
  } else x <- matrix(1, n, 1L)
  base <- study_reference(env, cbind(seq_len(n - 1L), 2:n))
  E <- internal(".environment_for_moments")(env, base)
  zz <- E$z[, 1L]
  shift <- shift + .24 * rnorm(n)
  truth <- function(idx) {
    k <- ncol(idx)
    S <- matrix(zz[idx], nrow(idx))
    R <- mean(range(zz))
    H <- diff(range(zz)) / 2
    m <- rowMeans(S)
    w <- (S - R) / H
    beta <- if (scenario == "linear") NULL else c(.25, .68, 1.05)
    f <- if (is.null(beta)) .70 * zz else drop(E$spline[, 1:3] %*% beta)
    F <- matrix(f[idx], nrow(idx))
    p <- .22 * m
    if (scenario != "linear") p <- p - .85 * rowMeans(w^2) + .52 * rowMeans(w^3)
    comb <- utils::combn(k, 2L)
    h <- b <- numeric(nrow(idx))
    strength <- switch(scenario, positive = .42, negative = -.42, collinear = .35, spatial = .40, 0)
    for (j in seq_len(ncol(comb))) {
      s1 <- comb[1L, j]
      s2 <- comb[2L, j]
      d <- abs(F[, s1] - F[, s2]) / ncol(comb)
      pair <- (S[, s1] + S[, s2]) / 2
      b <- b + d
      h <- h + strength * (m - pair) / H * d
    }
    eta <- -.72 + p + b + h
    list(y = if (scenario == "sparse") ugdm_response(x, idx, "relative")$y else
      plogis(eta + rowMeans(matrix(shift[idx], nrow(idx)))), eta = eta, interaction = h,
      sharing = if (scenario == "sparse")
        mean(internal("shared_species_cpp")(1L * (x > 0), idx) > 0) else NA_real_)
  }
  list(x = x, env = env, ref = base, truth = truth, coords = cbind(s, sin(pi * s)))
}
form <- data.frame(model = c("A", "B", "C", "D", "D_weak", "D_derivative", "D_strong",
    "D_coef", "D_all", "D_null"),
  position = c("linear", rep("nonlinear", 9L)),
  spread = c("constant", "constant", "varying", rep("contextual", 7L)),
  spread_smooth = c(rep(0, 8), .002, .002),
  position_smooth = c(rep(0, 8), .002, .002),
  interaction_shrink = c(rep(0, 4), 1e-6, .005, .15, .005, .005, .005),
  ridge = c(rep(0, 4), rep(1e-6, 6)),
  spread_null = c(rep(0, 9), .01),
  position_null = c(rep(0, 9), .01),
  interaction_penalty = c(rep("coefficients", 4), rep("derivative", 3), "coefficients",
    rep("derivative", 2)),
  penalty_scale = c(rep("native", 4), rep("range", 3), "native", rep("range", 2)))
design_metrics <- function(fit, d, y, calculate = TRUE) {
  if (!calculate) return(c(rank = NA_real_, parameters = NA_real_, condition = NA_real_))
  stats <- study_diag(fit, d, y)
  c(rank = stats$rank, parameters = stats$parameters, condition = stats$condition)
}
run <- function(r) {
  q <- DGP(r)
  ix <- sample_sets_coupled(q$x, orders, 1900L, sampler = "uniform", exact_limit = 0L)
  ev <- sample_sets_coupled(q$x, orders, 2800L, sampler = "uniform", exact_limit = 0L)
  result <- list()
  for (j in seq_along(orders)) {
    k <- orders[j]
    ixj <- ix[[j]]$idx
    evj <- ev[[j]]$idx
    ref <- q$ref
    y <- q$truth(ixj)$y
    actual <- q$truth(evj)
    des <- lapply(c("constant", "constant", "varying", "contextual"),
      function(s) study_d(q$env, ixj, ref, "nonlinear", s))
    des[[1L]] <- study_d(q$env, ixj, ref, "linear", "constant")
    assess <- lapply(seq_len(4L), function(i)
      study_d(q$env, evj, ref, form$position[i], form$spread[i]))
    score <- study_contrast(assess[[2L]], neighbours = 30L)
    true_eta <- drop(score$L %*% actual$eta)
    true_response <- drop(score$L %*% plogis(actual$eta))
    true_context <- drop(score$L[3L, , drop = FALSE] %*% actual$interaction)
    if (scenario == "sparse") true_eta[] <- true_response[] <- true_context <- NA_real_
    context_rank <- study_condition(des[[4L]])
    start <- NULL
    for (i in seq_len(nrow(form))) {
      shape <- min(i, 4L)
      opt <- if (i <= 4L) list() else as.list(form[i, 4:11])
      fit <- study_fit(des[[shape]], y, opt, if (i > 4L) start else NULL)
      if (shape == 4L && fit$converged) start <- fit$coef
      mu <- ugdm_predict(fit, assess[[shape]])
      c <- ugdm_interaction_contrast(fit, assess[[shape]], score$L)
      diag <- design_metrics(fit, des[[shape]], y, i <= 4L)
      observed <- study_rows(form$model[i], c)
      contextual <- NA_real_
      if (shape == 4L && context_rank["incremental_rank"] == ncol(des[[4L]]$bounded_context) &&
          context_rank["residual_context"] > .05) {
        width <- ncol(des[[4L]]$dir)
        widths <- seq_len(ncol(q$env))
        b <- matrix(fit$coef[1L + seq_len(width)], 12L, length(widths))
        delta <- as.numeric(b[10:12, , drop = FALSE] - b[7:9, , drop = FALSE])
        contextual <- drop(score$L[3L, , drop = FALSE] %*%
          assess[[4L]]$bounded_context %*% delta)
      }
      result[[length(result) + 1L]] <- data.frame(scenario, replicate = r, n, order = k,
        convergence = fit$converged, active = sum(fit$coef == 0), in_sample_D2 = fit$D2,
        test_RMSE = sqrt(mean((actual$y - mu)^2)),
        reference_RMSE = if (scenario == "sparse") NA_real_ else sqrt(mean((plogis(actual$eta) - mu)^2)),
        support_gap = score$gap, sharing = actual$sharing, observed,
        true_link_spread = true_eta[1L], true_link_position = true_eta[2L], true_link_joint = true_eta[3L],
        true_response_spread = true_response[1L], true_response_position = true_response[2L],
        true_response_joint = true_response[3L], true_context = true_context,
        fitted_context = contextual,
        context_incremental_rank = context_rank["incremental_rank"],
        context_orthogonal_norm = context_rank["residual_context"],
        rank = diag["rank"], parameters = diag["parameters"], condition = diag["condition"],
        row.names = NULL)
    }
  }
  do.call(rbind, result)
}
mode <- if (length(args) >= 3L) match.arg(args[3L], c("recovery", "selection")) else "recovery"
if (mode == "recovery") {
  all <- lapply(seq_len(B), function(r) {
    message(scenario, " ", r, "/", B)
    run(r)
  })
  write_study(do.call(rbind, all), "docs/adaptive_refinement/results", paste0("context_", scenario))
} else {
  source("docs/adaptive_refinement/ecology_selection_sim.R")
}
