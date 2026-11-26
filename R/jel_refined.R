.jel_setup <- function(data, prior_sd, score, method, adjusted) {
  if (!length(data)) stop("Provide at least one order.")
  data <- lapply(data, function(d) if (is.null(d$design)) d else list(d))
  if (any(vapply(data, function(group) any(vapply(group, function(d)
    !is.null(d$fit) && !isTRUE(d$fit$free_sign), TRUE)), TRUE))) {
    stop("JEL requires regular unconstrained fits.")
  }
  prep <- lapply(data, .prepare_score_pools)
  first <- lapply(prep, function(p) p[[1L]])
  n <- first[[1L]]$n
  width <- vapply(first, function(p) ncol(p$X), integer(1))
  d <- sum(width)
  if (d > 12L || n <= d + 1L || any(vapply(first, function(p) p$n != n, TRUE))) {
    stop("Use a common site index and at most 12 joint coefficients.")
  }
  if (length(prior_sd) == 1L) prior_sd <- rep(prior_sd, d)
  if (length(prior_sd) != d || any(!is.finite(prior_sd) | prior_sd <= 0)) stop("Invalid prior SD.")
  if (adjusted && method != "EL") stop("The adjusted likelihood requires ordinary EL.")
  coef_names <- make.unique(unlist(lapply(first, function(p) paste0("order", p$k, ".", colnames(p$X)))))
  designs <- lapply(prep, function(pool) list(
    X = do.call(rbind, lapply(pool, function(p) p$X)),
    y = unlist(lapply(pool, function(p) p$y), use.names = FALSE),
    w = unlist(lapply(pool, function(p) p$effective), use.names = FALSE) / length(pool)))
  fits <- lapply(designs, function(p) {
    if (all(p$y[p$w > 0] == p$y[which(p$w > 0)[1L]])) stop("JEL requires nonconstant responses.")
    fit <- suppressWarnings(stats::glm.fit(p$X, p$y, weights = p$w,
      family = stats::binomial(), control = stats::glm.control(epsilon = 1e-10, maxit = 100)))
    if (!fit$converged || fit$rank != ncol(p$X) || any(!is.finite(fit$coefficients))) {
      stop("Unrestricted point fit is nonregular.")
    }
    fit
  })
  index <- split(seq_len(d), rep(seq_along(width), width))
  point <- stats::setNames(unlist(lapply(fits, function(f) f$coefficients), use.names = FALSE), coef_names)
  J <- function(theta) do.call(cbind, lapply(seq_along(prep), function(j)
    .pool_scores(prep[[j]], theta[index[[j]]], score)))
  moments <- J(point)
  influence <- do.call(cbind, lapply(seq_along(designs), function(j) {
    p <- designs[[j]]
    mu <- fits[[j]]$fitted.values
    H <- crossprod(p$X * sqrt(p$w * mu * (1 - mu)))
    z <- moments[, index[[j]], drop = FALSE]
    sweep(z, 2L, colMeans(z)) %*% solve(H)
  }))
  V <- crossprod(influence) / (n * (n - 1))
  ev <- eigen(V, symmetric = TRUE)
  if (max(ev$values) <= 0) stop("Zero joint site influence.")
  floor <- max(ev$values) * 1e-8
  W <- sweep(ev$vectors, 2L, sqrt(pmax(ev$values, floor)), "*")
  log_density <- function(z, certify = FALSE) {
    theta <- point + drop(W %*% z)
    el <- ugdm_empirical_likelihood(J(theta), method, adjusted = adjusted, certify = certify)
    list(value = el$loglik + sum(stats::dnorm(theta, sd = prior_sd, log = TRUE)),
      status = if (el$ok) "valid" else el$reason)
  }
  list(log_density = log_density, point = point, W = W, d = d,
       initial_regularized = sum(ev$values < floor))
}

.jel_mode <- function(target) {
  f <- function(z) {
    p <- target$log_density(z)$value
    if (is.finite(p)) -p else 1e15
  }
  origin <- numeric(target$d)
  if (!is.finite(target$log_density(origin)$value)) {
    why <- target$log_density(origin, certify = TRUE)$status
    stop("Initial empirical likelihood is unavailable: ", why)
  }
  opt <- stats::optim(origin, f, method = "BFGS",
    control = list(maxit = 100, reltol = 1e-9, ndeps = rep(1e-4, target$d)))
  if (!is.finite(target$log_density(opt$par)$value)) {
    opt$par <- origin
    opt$convergence <- 100L
  }
  H <- stats::optimHess(opt$par, f, control = list(ndeps = rep(.002, target$d)))
  ev <- eigen(H, symmetric = TRUE)
  # Clipping affects the proposal only. All draws still use the exact likelihood.
  vals <- pmin(pmax(ev$values, 1e-3), 1e6)
  scale <- tcrossprod(sweep(ev$vectors, 2L, 1 / sqrt(vals), "*"))
  list(center = opt$par, U = chol(scale), convergence = opt$convergence,
       curvature = ev$values, clipped = sum(vals != ev$values))
}

.jel_component <- function(center, U, weight, df) list(center = center, U = U, weight = weight, df = df)

.jel_importance <- function(N, parts, target, certify) {
  d <- target$d
  index <- sample.int(length(parts), N, replace = TRUE,
                      prob = vapply(parts, function(a) a$weight, numeric(1)))
  z <- matrix(0, N, d)
  for (j in seq_along(parts)) {
    rows <- which(index == j)
    m <- length(rows)
    if (!m) next
    a <- parts[[j]]
    t <- matrix(stats::rnorm(m * d), m, d) / sqrt(stats::rchisq(m, a$df) / a$df)
    z[rows, ] <- sweep(t %*% a$U, 2L, a$center, "+")
  }
  props <- vapply(parts, function(a) {
    u <- forwardsolve(t(a$U), t(sweep(z, 2L, a$center)))
    df <- a$df
    log(a$weight) + lgamma((df + d) / 2) - lgamma(df / 2) -
      d / 2 * log(df * pi) - sum(log(diag(a$U))) - (df + d) / 2 * log1p(colSums(u^2) / df)
  }, numeric(N))
  largest <- apply(props, 1L, max)
  logq <- largest + log(rowSums(exp(props - largest)))
  eval <- lapply(seq_len(N), function(i) target$log_density(z[i, ], certify))
  logw <- vapply(eval, function(a) a$value, numeric(1)) - logq
  if (!any(is.finite(logw))) stop("No importance proposals passed the empirical likelihood.")
  w <- exp(logw - max(logw))
  list(z = z, logw = logw, weights = w / sum(w),
       status = vapply(eval, function(a) a$status, character(1)))
}

.jel_chain_diagnostics <- function(chains) {
  n <- nrow(chains[[1L]])
  d <- ncol(chains[[1L]])
  half <- n %/% 2L
  split <- unlist(lapply(chains, function(a) list(a[seq_len(half), , drop = FALSE],
    a[n - rev(seq_len(half)) + 1L, , drop = FALSE])), recursive = FALSE)
  m <- length(split)
  means <- vapply(split, colMeans, numeric(d))
  variances <- vapply(split, function(x) apply(x, 2L, stats::var), numeric(d))
  W <- rowMeans(variances)
  B <- apply(means, 1L, stats::var) * half
  rhat <- ifelse(W > 0, sqrt(((half - 1) * W / half + B / half) / W), Inf)
  ess <- vapply(seq_len(d), function(j) {
    ac <- lapply(chains, function(x) {
      if (stats::var(x[, j]) <= 0) return(NULL)
      as.numeric(stats::acf(x[, j], lag.max = min(400L, n %/% 2L), plot = FALSE)$acf)[-1L]
    })
    if (any(lengths(ac) == 0L)) return(0)
    rho <- Reduce("+", ac) / length(ac)
    pairs <- rowSums(matrix(utils::head(rho, length(rho) %/% 2L * 2L), ncol = 2L, byrow = TRUE))
    stop_at <- match(TRUE, pairs <= 0, nomatch = length(pairs) + 1L) - 1L
    min(m * half, m * half / (1 + 2 * sum(pairs[seq_len(stop_at)])))
  }, numeric(1))
  list(rhat = rhat, ess = ess)
}

.jel_mcmc <- function(target, mode, draws, warmup, chains, adapt_covariance = TRUE) {
  d <- target$d
  acceptance <- numeric(chains)
  samples <- lapply(seq_len(chains), function(j) {
    z <- mode$center
    lp <- target$log_density(z)$value
    for (i in seq_len(40L)) {
      candidate <- z + drop(stats::rnorm(d) %*% mode$U)
      value <- target$log_density(candidate)$value
      if (is.finite(value)) {
        z <- candidate
        lp <- value
        break
      }
    }
    step <- 2.38 / sqrt(d)
    proposal <- mode$U
    modal <- crossprod(mode$U)
    history <- matrix(0, warmup, d)
    out <- matrix(0, draws, d)
    accepted <- 0L
    for (i in seq_len(warmup + draws)) {
      candidate <- z + step * drop(stats::rnorm(d) %*% proposal)
      p <- target$log_density(candidate)$value
      take <- log(stats::runif(1L)) < p - lp
      if (take) {
        z <- candidate
        lp <- p
      }
      if (i <= warmup) {
        step <- step * exp((take - .24) / sqrt(i + 15))
        history[i, ] <- z
        if (adapt_covariance && i %% 50L == 0L && i >= max(150L, 25L * d)) {
          recent <- history[seq.int(i %/% 4L + 1L, i), , drop = FALSE]
          proposal <- chol(.90 * stats::cov(recent) + .10 * modal)
        }
      } else {
        out[i - warmup, ] <- z
        accepted <- accepted + take
      }
    }
    acceptance[j] <<- accepted / draws
    out
  })
  list(samples = samples, acceptance = acceptance)
}

.jel_mixture_pilot <- function(target, mode, n) {
  run <- .jel_mcmc(target, mode, draws = n, warmup = max(150L, n %/% 2L), chains = 4L)
  z <- do.call(rbind, run$samples)
  center <- colMeans(z)
  S <- stats::cov(z)
  modal <- crossprod(mode$U)
  U <- chol(S + .03 * modal)
  scaled <- t(forwardsolve(t(U), t(sweep(z, 2L, center))))
  labels <- stats::kmeans(scaled, centers = 3, nstart = 3, iter.max = 40, algorithm = "Lloyd")$cluster
  components <- lapply(split(seq_len(nrow(z)), labels), function(idx) {
    u <- z[idx, , drop = FALSE]
    C <- if (nrow(u) > ncol(u)) stats::cov(u) else S
    .jel_component(colMeans(u), chol(C + .05 * S + .03 * modal),
                   .80 * nrow(u) / nrow(z), 8)
  })
  components <- c(components, list(
    .jel_component(mode$center, mode$U, .06, 8),
    .jel_component(center, U, .09, 5),
    .jel_component(center, 2.5 * U, .05, 3)))
  list(parts = components, diagnostic = .jel_chain_diagnostics(run$samples))
}

#' Curvature-whitened joint MS-GDUM JEL
#'
#' Uses the exact same site empirical likelihood in mode search, adaptively
#' proposed importance integration and independent random-walk MCMC.
#' Importance adaptation uses a pilot discarded from the final estimates.
#' @param data As in [ugdm_jel()]. An order may instead provide a list of
#'   independent set pools on its fixed set of sites.
#' @param prior_sd,method,seed As in [ugdm_jel()].
#' @param score Site score construction as in [ugdm_pool_scores()].
#' @param integration Adaptive importance sampling or independent chains.
#' @param draws Per-batch importance proposals or retained iterations per chain.
#' @param batches Number of independent final importance batches.
#' @param pilot Discarded importance draws or short-chain iterations per chain
#'   for adapting the proposal.
#' @param pilot_kernel Use an importance pilot or use short preliminary
#'   chains to adapt a multicomponent proposal. These chains are separate
#'   from any MCMC reference and their draws are discarded.
#' @param warmup,chains MCMC warmup per chain and number of chains.
#' @param adapt_covariance Adjust the random-walk covariance during warmup.
#'   The frozen proposal is used after warmup.
#' @param adjusted Use the adjusted ordinary EL as a separate target.
#'   Assess support before using this option.
#' @param certify Attempt a strict hull separation certificate on failures
#'   during importance integration.
#' @return Posterior draws and weights, batch/chain labels, mode and
#'   curvature diagnostics and Monte Carlo diagnostics. A successful
#'   computation does not establish ecological calibration.
#' @export
ugdm_jel_refined <- function(data, prior_sd = 2, method = c("EL", "ETEL"),
    score = c("centered", "conditional", "deletion", "inclusion"),
    integration = c("adaptive", "mcmc"), draws = 1000L, batches = 2L,
    pilot = 200L, warmup = 700L, chains = 4L,
    pilot_kernel = c("importance", "chains"), adapt_covariance = TRUE,
    adjusted = FALSE, certify = FALSE, seed = NULL) {
  method <- match.arg(method)
  score <- match.arg(score)
  integration <- match.arg(integration)
  pilot_kernel <- match.arg(pilot_kernel)
  start <- proc.time()[["elapsed"]]
  draws <- .set_design_integer(draws, "draws", 20L)
  if (!is.null(seed)) set.seed(seed)
  target <- .jel_setup(data, prior_sd, score, method, adjusted)
  mode <- .jel_mode(target)
  convert <- function(z) {
    x <- sweep(z %*% t(target$W), 2L, target$point, "+")
    colnames(x) <- names(target$point)
    x
  }
  if (integration == "adaptive") {
    batches <- .set_design_integer(batches, "batches", 2L)
    pilot <- .set_design_integer(pilot, "pilot", 20L)
    pilot_diag <- NULL
    pilot_ess <- NA_real_
    if (pilot_kernel == "chains") {
      tune <- .jel_mixture_pilot(target, mode, pilot)
      parts <- tune$parts
      pilot_diag <- tune$diagnostic
    } else {
      normal <- .jel_component(mode$center, mode$U, .85, 8)
      tails <- .jel_component(mode$center, 2.5 * mode$U, .15, 4)
      tune <- .jel_importance(pilot, list(normal, tails), target, FALSE)
      center <- colSums(tune$z * tune$weights)
      S <- crossprod(sweep(tune$z, 2L, center) * sqrt(tune$weights))
      U <- chol(.85 * S + .15 * crossprod(mode$U))
      parts <- list(.jel_component(center, U, .65, 8),
                    .jel_component(mode$center, mode$U, .27, 8),
                    .jel_component(mode$center, 2.5 * mode$U, .08, 4))
      pilot_ess <- 1 / sum(tune$weights^2)
    }
    run <- lapply(seq_len(batches), function(j) .jel_importance(draws, parts, target, certify))
    z <- do.call(rbind, lapply(run, function(r) r$z))
    logw <- unlist(lapply(run, function(r) r$logw))
    weight <- exp(logw - max(logw))
    weight <- weight / sum(weight)
    diag <- list(ess = 1 / sum(weight^2), batch_ess = vapply(run, function(r) 1 / sum(r$weights^2), 0.0),
      likelihood_status = table(unlist(lapply(run, function(r) r$status))), max_weight = max(weight),
      pilot_ess = pilot_ess, pilot_chain_diagnostics = pilot_diag)
    group <- rep(seq_len(batches), each = draws)
  } else {
    chains <- .set_design_integer(chains, "chains", 2L)
    warmup <- .set_design_integer(warmup, "warmup", 1L)
    run <- .jel_mcmc(target, mode, draws, warmup, chains, adapt_covariance)
    z <- do.call(rbind, run$samples)
    weight <- rep(1 / nrow(z), nrow(z))
    diag <- c(list(acceptance = run$acceptance),
      .jel_chain_diagnostics(lapply(run$samples, convert)))
    group <- rep(seq_len(chains), each = draws)
  }
  c(list(draws = convert(z), weights = weight, group = group,
    estimate = target$point, mode = drop(convert(matrix(mode$center, nrow = 1L))),
    mode_convergence = mode$convergence, curvature = mode$curvature, clipped = mode$clipped,
    initial_regularized = target$initial_regularized,
    method = method, integration = integration, adjusted = adjusted,
    seconds = proc.time()[["elapsed"]] - start), diag)
}
