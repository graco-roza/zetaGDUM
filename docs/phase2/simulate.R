# Run via Rscript docs/phase2/simulate.R OUTPUT STAGE [REPLICATES].
library(zetaGDUM)
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1L] else "docs/phase2/results"
stage <- if (length(args) > 1L) args[2L] else "all"
repetitions <- if (length(args) > 2L) as.integer(args[3L]) else 32L
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, name) write.csv(x, file.path(output, paste0(name, ".csv")), row.names = FALSE)
internal <- function(name) get(name, asNamespace("zetaGDUM"))
draw_sets <- internal(".inference_sets")
fit_sets <- internal(".inference_fit")
effects <- internal(".order_effects")
effect_info <- internal(".order_information")
product_weights <- internal(".product_sampling_weights")
kish <- internal(".kish")
groups <- list(gradient = "gradient", control = "control")
fit_quiet <- function(...) suppressWarnings(fit_sets(...))

generate_sites <- function(n, signal, structure, spatial, seed, breadth_override = NULL, occupancy_shape = NULL) {
  set.seed(seed)
  coords <- matrix(runif(2 * n), n, 2)
  gradient <- rnorm(n)
  if (spatial) gradient <- 0.65 * sin(2 * pi * coords[, 1]) * sqrt(2) + sqrt(1 - .65^2) * gradient
  env <- data.frame(gradient = gradient, control = rnorm(n))
  breadth <- if (is.null(breadth_override)) switch(structure, broad = 1.45, narrow = .48, patchy = .7) else breadth_override
  ns <- 60L
  opt <- runif(ns, -1.6, 1.6)
  if (is.null(occupancy_shape)) occupancy_shape <- switch(structure, broad = c(6, 2), narrow = c(2, 5), patchy = c(1.2, 3.2))
  occ <- rbeta(ns, occupancy_shape[1], occupancy_shape[2])
  contrast <- exp(-outer(gradient, opt, "-")^2 / (2 * breadth^2))
  if (structure == "patchy") {
    centers <- runif(ns)
    contrast <- contrast * (abs(outer(coords[, 1], centers, "-")) <= .35)
  }
  strength <- c(null = 0, weak = .4, strong = .9)[[signal]]
  prob <- sweep((1 - strength) + strength * contrast, 2, occ, "*")
  comm <- matrix(rbinom(n * ns, 1, as.vector(prob)) * (1 + rpois(n * ns, 3)), n, ns)
  # Prevent an undefined relative or incidence denominator without adding ubiquitous taxa.
  empty <- which(rowSums(comm) == 0)
  if (length(empty)) comm <- cbind(comm, diag(n)[, empty, drop = FALSE])
  u <- rnorm(n)
  if (spatial) u <- drop(t(chol(exp(-as.matrix(dist(coords)) / .24))) %*% u)
  list(comm = comm, env = env, P = (comm > 0) * 1L, u = tanh(u), coords = coords)
}

anchor <- data.frame(gradient = seq(-3, 3, length.out = 21), control = cos(seq(0, 7, length.out = 21)) * 2)
anchor_design <- suppressWarnings(ugdm_design(anchor, sample_sets(nrow(anchor), 2, M = 100, seed = 10)))

orders_for <- function(n) unique(c(2L, round(.3 * n), round(.6 * n), round(.85 * n), n - 1L))
scenarios <- data.frame(
  n = c(18, 18, 36, 36, 36, 72, 72, 36),
  signal = c("null", "weak", "null", "weak", "strong", "weak", "strong", "null"),
  structure = c("broad", "narrow", "narrow", "patchy", "broad", "narrow", "patchy", "patchy"),
  spatial = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))

known_basis <- function(k, n, signal) {
  level <- c(null = 0, weak = .30, strong = 1)[[signal]]
  beta <- level * c(.45, .40, .35, .13, .11, .09)
  beta[1:3] <- beta[1:3] * (1 - .5 * k / n)
  lambda <- level * c(.09, -.025)
  list(beta = beta, lambda = lambda, truth = c(spread.gradient = sum(beta[1:3]),
        spread.control = sum(beta[4:6]), position.gradient = lambda[1], position.control = lambda[2]))
}

known_sample <- function(sim, sets, b) {
  i <- sets$idx
  d <- ugdm_design(sim$env, i, ref = anchor_design)
  mu <- plogis(-.3 + drop(d$Xdir %*% b$beta + d$Xpos %*% b$lambda))
  # The multiplicative factor gives a bounded response when mu is near a boundary.
  # Use min(mu, 1-mu) for a strict bound throughout the supported design.
  y <- mu + .85 * pmin(mu, 1 - mu) * rowMeans(matrix(sim$u[i], nrow(i), ncol(i)))
  list(design = d, y = y, w = sets$sampling_weights, trials = rep(1, length(y)))
}

make_fit <- function(a, free_sign, w = a$w) {
  ugdm_fit(a$design, a$y, a$trials, sampling_weights = w, free_sign = free_sign)
}

safe_effect <- function(a, w = a$w) {
  f <- make_fit(a, TRUE, w)
  h <- effect_info(a$design, f, w)
  if (!f$converged || h$rank != h$parameters || !is.finite(f$D2)) {
    return(effects(f) * NA_real_)
  }
  effects(f)
}

intervals <- function(point, draws, method) {
  n <- nrow(draws)
  if (method == "jackknife") {
    se <- sqrt((n - 1) / n * colSums(sweep(draws, 2, colMeans(draws), "-")^2))
    span <- qt(.975, n - 1L) * se
    return(rbind(point - span, point + span))
  }
  apply(draws, 2, function(z) {
    if (!all(is.finite(z))) return(rep(NA_real_, 2))
    quantile(z, c(.025, .975), names = FALSE)
  })
}

calibrate <- function() {
  rows <- contrasts <- perm_rows <- identification <- finite <- list()
  B <- 39L
  requested <- Sys.getenv("ZETAGDUM_CALIBRATION_SCENARIOS")
  ids <- if (nzchar(requested)) as.integer(strsplit(requested, ",", fixed = TRUE)[[1L]]) else seq_len(nrow(scenarios))
  if (anyNA(ids) || anyDuplicated(ids) || any(!ids %in% seq_len(nrow(scenarios)))) {
    stop("Select distinct calibration scenarios from 1 to ", nrow(scenarios), ".")
  }
  for (scenario in ids) {
    config <- scenarios[scenario, ]
    rownames(config) <- NULL
    message("Calibration scenario ", scenario, " of ", nrow(scenarios))
    n <- config$n
    for (r in seq_len(repetitions)) {
      sim <- generate_sites(n, config$signal, config$structure, config$spatial, 10000L * scenario + r)
      W <- site_weights(n, B, seed = 81000L + r)
      fits <- list()
      for (k in orders_for(n)) {
        b <- known_basis(k, n, config$signal)
        sets <- draw_sets(sim$P, k, 1600L, 0, exact_limit = 4500L, seed = 91000L + 73L * r + k)
        a <- known_sample(sim, sets, b)
        center <- safe_effect(a)
        uncertainty <- ugdm_site_uncertainty(a$design, a$y, a$trials, a$w, free_sign = TRUE,
                                            ndraw = 399L, seed = 42000L + r)
        direct <- product <- matrix(NA_real_, B, length(center), dimnames = list(NULL, names(center)))
        for (j in seq_len(B)) {
          w <- product_weights(W[, j], a$design$idx, a$w)
          product[j, ] <- safe_effect(a, w)
          if (sets$exact) {
            direct[j, ] <- product[j, ]
          } else {
            set.seed(73200L + 100000L * r + 200L * k + j)
            weighted <- known_sample(sim, draw_sets(sim$P, k, 1600L, 0, omega = W[, j]), b)
            direct[j, ] <- safe_effect(weighted)
          }
        }
        null <- matrix(NA_real_, 399L, length(center))
        draws <- list(product_bb = product, direct_bb = direct,
                      jackknife = uncertainty$delete, multiplier = if (is.null(uncertainty$multiplier)) null else uncertainty$multiplier,
                      multiplier_raw = if (is.null(uncertainty$multiplier_raw)) null else uncertainty$multiplier_raw)
        point <- make_fit(a, TRUE)
        constrained <- make_fit(a, FALSE)
        identification[[length(identification) + 1L]] <- data.frame(
          scenario, config, replicate = r, order = k, exact = sets$exact,
          rank = uncertainty$diagnostics$rank, dimension = uncertainty$diagnostics$parameters,
          condition = uncertainty$diagnostics$condition, converged = point$converged,
          constrained_active = sum(constrained$beta <= 1e-7),
          min_delete_ess = min(uncertainty$delete_ess), D2 = point$D2)
        fits[[as.character(k)]] <- list(center = center, truth = b$truth, draws = draws)
        for (method in names(draws)) {
          ci <- intervals(center, draws[[method]], method)
          rows[[length(rows) + 1L]] <- data.frame(scenario, config, replicate = r, order = k, method,
            exact = sets$exact, effect = names(center), truth = unname(b$truth), estimate = unname(center),
            low = ci[1, ], high = ci[2, ], row.names = NULL)
        }
        if (sets$exact) {
          sampled <- vapply(seq_len(6L), function(j) {
            s <- draw_sets(sim$P, k, 500L, 0, seed = 620000L + 1000L * r + k * 10L + j)
            safe_effect(known_sample(sim, s, b))
          }, numeric(length(center)))
          finite[[length(finite) + 1L]] <- data.frame(scenario, config, replicate = r, order = k,
            effect = rep(names(center), 6L), reference = rep(unname(center), 6L),
            sampled = as.vector(sampled), draws = rep(seq_len(6L), each = length(center)))
        }
        # Global constrained permutation test, valid only under exchangeable iid null sites.
        if (!config$spatial && k %in% orders_for(n)[c(1L, 3L, 5L)]) {
          # A deterministic set pool is independent of the environmental site labels.
          set.seed(93000L + r + 100L * k)
          perm <- replicate(39L, {
            d <- ugdm_design(sim$env[sample.int(n), , drop = FALSE], sets$idx, ref = anchor_design)
            a$design <- d
            make_fit(a, FALSE)$D2
          })
          pval <- if (all(is.finite(c(constrained$D2, perm)))) (1 + sum(perm >= constrained$D2)) / 40 else NA_real_
          perm_rows[[length(perm_rows) + 1L]] <- data.frame(scenario, config, replicate = r, order = k, p = pval)
        }
      }
      for (j in seq.int(2L, length(fits))) {
        before <- fits[[j - 1L]]
        after <- fits[[j]]
        delta <- after$center - before$center
        truth <- after$truth - before$truth
        for (method in names(after$draws)) {
          ci <- intervals(delta, after$draws[[method]] - before$draws[[method]], method)
          contrasts[[length(contrasts) + 1L]] <- data.frame(scenario, config, replicate = r,
            from = as.integer(names(fits)[j - 1L]), to = as.integer(names(fits)[j]), method,
            effect = names(delta), truth = unname(truth), estimate = unname(delta), low = ci[1, ], high = ci[2, ],
            row.names = NULL)
        }
      }
    }
  }
  write_result(do.call(rbind, rows), "calibration")
  write_result(do.call(rbind, contrasts), "changes")
  if (length(perm_rows)) write_result(do.call(rbind, perm_rows), "permutation")
  write_result(do.call(rbind, identification), "identification")
  write_result(do.call(rbind, finite), "finite_reference")
}

precise <- function() {
  scenarios <- data.frame(n = c(24, 48, 48), signal = c("weak", "strong", "strong"),
    structure = c("narrow", "broad", "patchy"), spatial = c(FALSE, FALSE, TRUE))
  precision <- diagnostics <- pilot <- control <- thresholds <- list()
  for (s in seq_len(nrow(scenarios))) {
    config <- scenarios[s, ]
    rownames(config) <- NULL
    message("Precision scenario ", s)
    sim <- generate_sites(config$n, config$signal, config$structure, config$spatial, 1300L + s)
    ks <- round(c(.16, .48, .8) * config$n)
    for (k in ks) {
      # Evaluate the pilot against fresh fits at every candidate fraction.
      a <- suppressWarnings(ugdm_sampling_pilot(sim$comm, sim$env, k, M = 1000, nrep = 10, component = "total", seed = 2560L + k))
      pilot[[length(pilot) + 1L]] <- data.frame(scenario = s, config, order = k,
                                               selected = a$alpha, a$candidates)
      for (alpha in a$candidates$alpha) {
        p <- list(comm = sim$comm, env = sim$env, P = sim$P)
        z <- suppressWarnings(internal(".precision_run")(p, k, 4000, 12L, alpha, "total", "default", FALSE, NULL,
                                                        240000L + 500L * k + 1000L * alpha))
        mat <- z$effects[, grepl("^(spread|position)\\.", colnames(z$effects)), drop = FALSE]
        sec <- median(z$diagnostics$seconds)
        v <- sum(apply(mat, 2, var))
        usable <- with(z$diagnostics, all(converged & rank == parameters & is.finite(D2)))
        control[[length(control) + 1L]] <- data.frame(scenario = s, config, order = k, alpha, selected = a$alpha,
          variance = v, seconds = sec, cost = if (usable) v * sec else NA_real_,
          converged = mean(z$diagnostics$converged), shared = median(z$diagnostics$shared),
          ess = median(z$diagnostics$ess))
      }
    }
    z <- suppressWarnings(zeta_order_precision(sim$comm, sim$env, ks, nrep = 8L, alpha = "adaptive",
      pilot_nrep = 10L, reference_M = 50000L, exact_limit = 30000L, component = "total", groups = groups, seed = 6350L + s))
    precision[[s]] <- data.frame(scenario = s, config, z$precision)
    diagnostics[[s]] <- data.frame(scenario = s, config, z$diagnostics)
    thresholds[[s]] <- data.frame(scenario = s, config, retention = z$retention_break, confirmed = z$confirmed_budget)
  }
  write_result(do.call(rbind, precision), "precision")
  write_result(do.call(rbind, diagnostics), "precision_diagnostics")
  write_result(do.call(rbind, pilot), "pilot")
  write_result(do.call(rbind, control), "pilot_validation")
  write_result(do.call(rbind, thresholds), "thresholds")
}

focused <- function() {
  sim <- generate_sites(48L, "strong", "broad", FALSE, 18817L)
  assess <- c("spread.gradient", "position.gradient", "importance.gradient", "shapley.gradient")
  z <- suppressWarnings(zeta_order_precision(sim$comm, sim$env, c(8L, 23L),
    nrep = 12L, alpha = .25, reference_M = 50000L, component = "total",
    groups = groups, seed = 5277L, assess = assess, stop_after_confirmation = TRUE))
  rows <- z$precision[z$precision$effect %in% assess, ]
  write_result(data.frame(confirmed_budget = z$confirmed_budget, rows), "focused")
}

ecological <- function() {
  rows <- list()
  for (s in seq_len(nrow(scenarios))) {
    config <- scenarios[s, ]
    rownames(config) <- NULL
    message("Ecological scenario ", s)
    for (r in seq_len(repetitions)) {
      n <- config$n
      sim <- generate_sites(n, config$signal, config$structure, config$spatial, 800000L + s * 10000L + r)
      for (k in orders_for(n)) {
        seed <- 93000L + r + k
        s0 <- draw_sets(sim$P, k, 4000, .5, exact_limit = 6000, seed = seed)
        f <- fit_quiet(sim$comm, sim$env, s0, "total", "default", FALSE, groups, anchor_design)
        share <- internal(".species_zeta_presence")(sim$comm, k)
        cv <- if (k <= round(.4 * n) && r <= 12L) suppressWarnings(ugdm_validate(
          sim$comm, sim$env, order = k, component = "total", M = 1000L, M_test = 2000L,
          nrep = 2L, test_frac = .4, seed = seed, sampling = "species",
          test_sampling = "conventional"))$summary else NULL
        info <- f$information
        rows[[length(rows) + 1L]] <- data.frame(scenario = s, config, replicate = r, order = k, exact = s0$exact,
          retention_break = zeta_break(zeta_curve(sim$comm, kmax = n)),
          expected_shared = sum(share), observed_sharing = mean(s0$shared > 0),
          probability_sharing = mean(s0$sampling_weights * (s0$shared > 0)),
          ess = kish(f$weights), condition = info$condition, rank = info$rank, dimension = info$parameters,
          active = info$active, converged = f$fit$converged, D2 = f$fit$D2,
          D2_test = if (is.null(cv) || !any(is.finite(cv$D2_test))) NA_real_ else median(cv$D2_test, na.rm = TRUE),
          RMSE_test = if (is.null(cv)) NA_real_ else median(cv$rmse_test), as.list(f$effects))
      }
    }
  }
  write_result(do.call(rbind, rows), "ecological")
}

sensitivity <- function() {
  # Vary occupancy and niche breadth independently, with the same signal and sampling design.
  grid <- expand.grid(n = c(24L, 60L), breadth = c(.48, 1.45),
                      occupancy = c("sparse", "common"), spatial = c(FALSE, TRUE),
                      stringsAsFactors = FALSE)
  rows <- list()
  for (s in seq_len(nrow(grid))) {
    x <- grid[s, ]
    rownames(x) <- NULL
    for (r in seq_len(repetitions)) {
      shape <- if (x$occupancy == "sparse") c(2, 5) else c(6, 2)
      sim <- generate_sites(x$n, "weak", "broad", x$spatial, 523000L + 1200L * s + r,
                            breadth_override = x$breadth, occupancy_shape = shape)
      for (k in round(c(.25, .6, .9) * x$n)) {
        sets <- draw_sets(sim$P, k, 2000, .5, exact_limit = 2000, seed = 624000L + 100L * r + k)
        f <- fit_quiet(sim$comm, sim$env, sets, "total", "default", FALSE, groups, anchor_design)
        imp <- f$effects[c("importance.gradient", "importance.control")]
        rows[[length(rows) + 1L]] <- data.frame(x, replicate = r, order = k, ratio = k / x$n,
          sharing = mean(sets$sampling_weights * (sets$shared > 0)), ess = kish(f$weights),
          identified = f$information$rank == f$information$parameters && f$fit$converged && is.finite(f$fit$D2),
          gradient_first = unname(imp[1] > imp[2]), as.list(f$effects))
      }
    }
  }
  write_result(do.call(rbind, rows), "sensitivity")
}

additive <- function() {
  n <- 18L
  R <- 500L
  B <- 399L
  set.seed(8320L)
  w <- site_weights(n, R * B)
  x <- matrix(rnorm(n * R), n, R)
  rows <- list()
  for (r in seq_len(R)) {
    a <- x[, r]
    ww <- w[, (r - 1L) * B + seq_len(B), drop = FALSE]
    full <- colSums(ww)
    mass <- ww * rep(full, each = n) - ww^2
    inclusion <- sweep(mass, 2, colSums(mass), "/")
    pair <- colSums(a * inclusion)
    omit <- sweep(1 / ww, 2, colSums(1 / ww), "/")
    near <- (sum(a) - colSums(a * omit)) / (n - 1L)
    draws <- list(pair, near, rep(mean(a), B))
    for (j in seq_along(draws)) {
      interval <- quantile(draws[[j]], c(.025, .975))
      rows[[length(rows) + 1L]] <- data.frame(replicate = r, order = c(2L, n - 1L, n)[j],
        bb_se = sd(draws[[j]]), bb_covered = interval[1] <= 0 && 0 <= interval[2],
        jack_se = if (j < 3L) sd(a) / sqrt(n) else NA_real_,
        jack_covered = if (j < 3L) abs(mean(a)) <= qt(.975, n - 1) * sd(a) / sqrt(n) else NA,
        exact_sampling_se = 1 / sqrt(n))
    }
  }
  write_result(do.call(rbind, rows), "additive")
}

selected <- if (stage == "all") c("calibrate", "precise", "focused", "ecological", "sensitivity", "additive") else stage
for (one in selected) switch(one, calibrate = calibrate(), precise = precise(),
  focused = focused(), ecological = ecological(), sensitivity = sensitivity(), additive = additive(), stop("Unknown stage: ", one))
