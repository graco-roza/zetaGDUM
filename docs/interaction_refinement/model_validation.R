# Run from the repository root after installing the development package.
library(zetaGDUM)
args <- commandArgs(TRUE)
mode <- if (length(args)) match.arg(args[1L], c("all", "recovery", "allocation", "curves")) else "all"
B <- if (length(args) > 1L) as.integer(args[2L]) else 24L
stopifnot(is.finite(B), B > 0)
out <- "docs/interaction_refinement/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
orders <- c(2L, 6L, 12L)
scenarios <- c("linear", "nonlinear_position", "curved_turnover", "positive", "negative",
               "collinear", "spatial", "niches_narrow", "niches_wide")
models <- data.frame(name = c("original", "nonlinear_position", "varying_spread", "joint",
                              "spread_smooth", "position_smooth", "interaction_shrink", "ridge", "joint_penalty"),
  position = c("linear", "nonlinear", "linear", rep("nonlinear", 6L)),
  spread = c("constant", "constant", rep("varying", 7L)),
  spread_smooth = c(rep(0, 4), .01, 0, 0, 0, .01),
  position_smooth = c(rep(0, 5), .005, 0, 0, .005),
  interaction_shrink = c(rep(0, 6), .005, 0, .005),
  ridge = c(rep(0, 7), .001, .001))

layout <- function(r, scenario) {
  set.seed(90263L + 811L * r + 417L * match(scenario, scenarios))
  n <- if (scenario == "spatial") 100L else 70L
  s <- sort(runif(n))
  z <- if (scenario == "spatial") 2.4 * (s - .5) + rnorm(n, 0, .20) else rnorm(n)
  z <- pmax(-2.3, pmin(2.3, z))
  env <- data.frame(climate = z)
  if (scenario == "collinear") env$collinear <- z + rnorm(n, 0, .075)
  site_error <- if (scenario == "spatial") {
    drop(t(chol(exp(-abs(outer(s, s, "-")) / .14) + diag(1e-7, n))) %*% rnorm(n)) * .25
  } else rnorm(n, 0, .28)
  niche <- startsWith(scenario, "niches")
  x <- NULL
  if (niche) {
    width <- if (scenario == "niches_narrow") .36 else .85
    optimum <- seq(-2.5, 2.5, length.out = 54L)
    expected <- exp(-outer(z, optimum, "-")^2 / (2 * width^2)) *
      rep(exp(rnorm(length(optimum), 0, .8)), each = n)
    spatial_error <- drop(t(chol(exp(-abs(outer(s, s, "-")) / .12) + diag(1e-7, n))) %*% rnorm(n)) * .26
    expected <- expected * exp(spatial_error * rep(rep(c(-1, 1), 27), each = n))
    x <- matrix(rpois(length(expected), 5 * expected + .025), n)
    x[rowSums(x) == 0, 1L] <- 1
  }
  ref <- suppressWarnings(ugdm_design(env, cbind(seq_len(n - 1L), 2:n)))
  build <- function(idx, position = "nonlinear", spread = "varying") {
    ugdm_interaction_design(suppressWarnings(ugdm_design(env, idx, ref = ref)), env, position, spread)
  }
  zs <- scale(env$climate, center = ref$environment_basis$center[1L],
              scale = ref$environment_basis$scale[1L])
  directional_basis <- as.matrix(predict(ref$environment_basis$splines[[1L]], drop(zs)))
  linear_coefficients <- unname(stats::lm.fit(cbind(1, directional_basis), drop(zs))$coefficients[-1L])
  stopifnot(max(abs(stats::lm.fit(cbind(1, directional_basis), drop(zs))$residuals)) < 1e-10)
  position_coef <- if (scenario %in% c("nonlinear_position", "positive", "negative", "collinear", "spatial")) {
    c(.28, -.95, .70)
  } else c(.28, 0, 0)
  beta <- if (scenario %in% c("linear", "nonlinear_position")) {
    linear_coefficients * (.65 / mean(linear_coefficients))
  } else c(.3, .65, 1.65)
  delta <- switch(scenario, positive = 1.8, negative = -1.8, spatial = 1.25, collinear = 1.8, 0)
  response <- function(idx) {
    a <- build(idx)
    p <- drop(a$pos[, 1:3, drop = FALSE] %*% position_coef)
    f <- drop(a$base[, 1:3, drop = FALSE] %*% beta)
    h <- delta * (a$u[, 1L] - .5) * f
    eta <- -.78 + p + f + h
    if (niche) {
      y <- ugdm_response(x, idx, "relative")$y
      return(list(y = y, eta = rep(NA_real_, length(y)), p = NA_real_, f = NA_real_, h = NA_real_))
    }
    error <- rowMeans(matrix(site_error[idx], nrow(idx)))
    # The site term is shared across sets. Deterministic conditional responses
    # isolate recovery from independent binomial noise at the combination level.
    list(y = plogis(eta + error), eta = eta, p = p, f = f, h = h)
  }
  list(env = env, n = n, ref = ref, build = build, response = response, s = s,
       position_coef = position_coef, beta = beta, delta = delta)
}
config_contrast <- function(a) {
  z <- cbind(drop(a$mean_position[, 1L]), sqrt(drop(a$position_variance[, 1L])))
  rank <- scale(z)
  bounds <- apply(rank, 2L, quantile, probs = c(.30, .70))
  target <- expand.grid(position = bounds[, 1L], spread = bounds[, 2L])
  count <- min(30L, nrow(z) %/% 4L)
  W <- matrix(0, 4L, nrow(z))
  gap <- numeric(4L)
  used <- integer()
  for (j in seq_len(nrow(target))) {
    distance <- sqrt(rowSums(sweep(rank, 2L, unlist(target[j, ]))^2))
    available <- setdiff(order(distance), used)
    use <- available[seq_len(count)]
    W[j, use] <- 1 / length(use)
    gap[j] <- mean(distance[use])
    used <- c(used, use)
  }
  list(w = drop(c(1,-1,-1,1) %*% W), support_gap = max(gap))
}
components <- function(fit, a) {
  b <- fit$coef[seq_len(fit$nd) + 1L]
  p <- drop(a$pos %*% utils::tail(fit$coef, ncol(a$pos)))
  h <- rep(0, length(p))
  f <- if (fit$spread == "constant") drop(a$dir %*% b) else {
    base <- h
    for (j in seq_along(a$vars)) {
      jj <- 6L * (j - 1L) + seq_len(6L)
      lower <- b[jj[1:3]]
      upper <- b[jj[4:6]]
      d <- a$base[, 3L * (j - 1L) + seq_len(3L), drop = FALSE]
      base <- base + drop(d %*% ((upper + lower) / 2))
      h <- h + drop(d %*% (upper - lower)) * (a$u[, j] - .5)
    }
    base
  }
  list(f = f, p = p, h = h)
}
folds <- function(s, n) as.integer(cut(rank(s), breaks = seq(0, n, length.out = 6L)))
refit <- function(q, idx, y, evaluation, m, stats = TRUE) {
  a <- q$build(idx, m$position, m$spread)
  f <- do.call(ugdm_fit_interactions, c(list(a, y), as.list(m[4:7])))
  b <- q$build(evaluation, m$position, m$spread)
  if (!stats) return(list(fit = f, design = b))
  list(fit = f, design = b, diagnostics = ugdm_interaction_diagnostics(f, a, y))
}

if (mode %in% c("all", "recovery")) {
  rows <- list()
  for (scenario in scenarios) {
    message("Recovery: ", scenario)
    for (r in seq_len(B)) {
      q <- layout(r, scenario)
      split <- sample.int(q$n, ceiling(.28 * q$n))
      train <- setdiff(seq_len(q$n), split)
      for (k in orders) {
        pool <- matrix(train[sample_sets(length(train), k, 2400L)], ncol = k)
        evaluation <- matrix(split[sample_sets(length(split), k, 1500L)], ncol = k)
        target <- q$response(evaluation)
        reference <- q$build(evaluation)
        C <- config_contrast(reference)
        y <- q$response(pool)$y
        # Hold out the evaluation sites from model preprocessing as well.
        qfit <- q
        qfit$ref <- suppressWarnings(ugdm_design(q$env, pool, basis_sites = train))
        qfit$build <- function(idx, position = "nonlinear", spread = "varying") {
          ugdm_interaction_design(ugdm_design(q$env, idx, ref = qfit$ref), q$env, position, spread)
        }
        candidate <- c(1:9)
        if (r <= min(8L, B) && scenario != "collinear" && floor(length(train) / 5L) >= k) {
          selection <- suppressWarnings(ugdm_select_interactions(q$env[train, , drop = FALSE],
            folds(q$s[train], length(train)), k, function(idx) q$response(matrix(train[idx], ncol = k))$y,
            models[4:9, 4:7], M = 600L, seed = 189L + r))
          # At high orders, each blocked fold may have fewer sites than k.
          if (!is.na(selection$selected)) candidate <- c(candidate, 3L + selection$selected)
        } else selection <- NULL
        for (j in seq_along(candidate)) {
          m <- models[candidate[j], ]
          z <- refit(qfit, pool, y, evaluation, m)
          f <- z$fit
          d <- z$design
          ds <- z$diagnostics
          eta <- ugdm_predict(f, d, "link")
          v <- components(f, d)
          rows[[length(rows) + 1L]] <- data.frame(scenario, r, k, n = q$n,
            model = if (j <= 9L) m$name else "site_selected",
            penalty = m$name,
            selected_folds = if (j > 9 && !is.null(selection)) sum(selection$selected_without_fold == selection$selected, na.rm = TRUE) else NA,
            converged = f$converged, rank = ds$rank, parameters = ds$parameters,
            condition = ds$condition, active = length(ds$active), support_gap = C$support_gap,
            D2 = f$D2, interaction_D2 = unname(ds$deviance_allocation[3L]),
            spread_D2 = unname(ds$deviance_allocation[1L]), position_D2 = unname(ds$deviance_allocation[2L]),
            pred_RMSE = sqrt(mean((plogis(eta) - target$y)^2)),
            truth_RMSE = sqrt(mean((plogis(eta) - plogis(target$eta))^2)),
            spread_SD = sd(v$f), true_spread_SD = sd(target$f),
            position_SD = sd(v$p), true_position_SD = sd(target$p),
            link = sum(C$w * eta), true_link = sum(C$w * target$eta),
            response = sum(C$w * plogis(eta)), true_response = sum(C$w * plogis(target$eta)),
            explicit = sum(C$w * v$h), true_explicit = sum(C$w * target$h),
            position_residual = if (length(ds$nonlinear_position_residual)) min(ds$nonlinear_position_residual, na.rm = TRUE) else NA,
            interaction_residual = if (length(ds$interaction_residual)) min(ds$interaction_residual, na.rm = TRUE) else NA)
        }
      }
    }
  }
  result <- do.call(rbind, rows)
  con <- gzfile(file.path(out, "model_recovery.csv.gz"), "wt", compression = 9)
  write.csv(result, con, row.names = FALSE)
  close(con)
}

if (mode %in% c("all", "allocation")) {
  runs <- list()
  for (scenario in c("curved_turnover", "positive", "spatial")) {
    for (r in seq_len(min(B, 12L))) {
      q <- layout(r, scenario)
      comparison <- lapply(orders, function(k) {
        idx <- sample_sets(q$n, k, 4500L)
        list(idx = idx, C = config_contrast(q$build(idx)))
      })
      for (M in c(1000L, 6000L)) for (p in 1:3) {
        paths <- sample_sets_coupled(matrix(1, q$n, 1L), orders, M = M, exact_limit = 0L)
        for (j in seq_along(orders)) {
          idx <- paths[[j]]$idx
          test <- comparison[[j]]
          z <- refit(q, idx, q$response(idx)$y, test$idx, models[9L, ], FALSE)
          cont <- ugdm_interaction_contrast(z$fit, z$design, test$C$w)
          runs[[length(runs) + 1L]] <- data.frame(scenario, r, M, p, k = orders[j],
            link = cont$link, response = cont$response,
            explicit = sum(test$C$w * components(z$fit, z$design)$h))
        }
      }
    }
  }
  write.csv(do.call(rbind, runs), file.path(out, "model_pool_precision.csv"), row.names = FALSE)
}

if (mode %in% c("all", "curves")) {
  grid_curves <- function(basis, grid, beta, delta, position_coef) {
    z <- (grid - basis$center[1L]) / basis$scale[1L]
    b <- basis$splines[[1L]]
    bounds <- attr(b, "Boundary.knots")
    mid <- mean(bounds)
    half <- diff(bounds) / 2
    B <- predict(b, pmax(bounds[1L], pmin(bounds[2L], z)))
    t <- (z - mid) / half
    u <- pmin(1, pmax(-1, t))
    m <- cbind(z, u^2 + 2 * u * (t-u), u^3 + 3 * u^2 * (t-u))
    curves <- cbind(directional = drop(B %*% beta), position = drop(m %*% position_coef),
                    endpoint_difference = drop(B %*% delta))
    sweep(curves, 2L, colMeans(curves))
  }
  rows <- list()
  for (scenario in c("linear", "nonlinear_position", "curved_turnover", "positive", "negative", "spatial")) {
    for (r in seq_len(min(B, 12L))) {
      q <- layout(r, scenario)
      ref <- q$ref$environment_basis
      grid <- seq(quantile(q$env$climate, .10), quantile(q$env$climate, .90), length.out = 51)
      target <- grid_curves(ref, grid, q$beta, q$delta * q$beta, q$position_coef)
      for (k in orders) {
        idx <- sample_sets(q$n, k, 2400L)
        y <- q$response(idx)$y
        for (m in c(1:4, 9L)) {
          a <- q$build(idx, models$position[m], models$spread[m])
          f <- do.call(ugdm_fit_interactions, c(list(a, y), as.list(models[m, 4:7])))
          lower <- unname(f$coef[2:4])
          upper <- if (a$spread == "varying") unname(f$coef[5:7]) else lower
          position <- c(unname(tail(f$coef, ncol(a$pos))), rep(0, 3L - ncol(a$pos)))
          z <- grid_curves(ref, grid, (lower+upper)/2, upper-lower, position)
          RMSE <- sqrt(colMeans((z-target)^2))
          rows[[length(rows)+1L]] <- data.frame(scenario, r, k, model = models$name[m],
            directional_RMSE = RMSE["directional"], position_RMSE = RMSE["position"],
            endpoint_RMSE = RMSE["endpoint_difference"])
        }
      }
    }
  }
  write.csv(do.call(rbind, rows), file.path(out, "model_curve_recovery.csv"), row.names = FALSE)
}
