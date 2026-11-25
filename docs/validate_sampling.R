#!/usr/bin/env Rscript
# Run from the repository root with the development version of zetaGDUM installed.
# Arguments are the output directory, then the optional paper_mee directory.
library(zetaGDUM)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else tempdir()
paper <- if (length(args) > 1L) args[2] else "../paper_mee"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
methods <- c("conventional", "species")
groups <- list(gradient1 = "env1", gradient2 = "env2")
median_finite <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
write_result <- function(x, name) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

community <- function(structure, n, seed) {
  sim <- switch(structure,
    broad = simulate_metacommunity(n, 65, breadth = 1.5, breadth_sd = 0.2,
                                  abundance_mean = 1.1, spatial = 0.5, seed = seed),
    narrow = simulate_metacommunity(n, 85, breadth = 0.53, breadth_sd = 0.35,
                                   abundance_mean = 1.2, spatial = 0.75, seed = seed),
    patchy = {
      pool <- rarity_species_pool(90, "Sr-Nn-La", p_dominant = 0.75, seed = seed)
      simulate_metacommunity(n, species = pool, spatial = 0.9, seed = seed + 10000L)
    })
  keep <- rowSums(sim$comm) > 0
  sim$comm <- sim$comm[keep, , drop = FALSE]
  sim$env <- sim$env[keep, , drop = FALSE]
  sim
}

coefficients <- function(fit) c(intercept = fit$intercept,
                                setNames(unname(fit$height), paste0("height_", names(fit$height))),
                                setNames(unname(fit$lambda), paste0("lambda_", names(fit$lambda))))
exact_fit <- function(sim, k, component) {
  idx <- t(combn(nrow(sim$comm), k))
  r <- ugdm_response(sim$comm, idx, component)
  d <- ugdm_design(sim$env, idx)
  f <- ugdm_fit(d, r$y, r$trials)
  list(fit = f, pshare = mean(zeta_sets(sim$comm, idx, FALSE)$shared > 0))
}

# Finite-community parameter target. Bootstrap coverage here is the frequency
# with which intervals contain the same community's exhaustive uniform fit.
recover <- function() {
  rows <- list()
  for (structure in c("broad", "narrow", "patchy")) {
    for (rep in seq_len(20)) {
      sim <- community(structure, 18, 100 + rep)
      for (k in c(3, 6, 10)) {
        if (nrow(sim$comm) < k) next
        for (component in c("total", "relative")) {
          exact <- exact_fit(sim, k, component)
          if (!is.finite(exact$fit$D2)) next
          truth <- coefficients(exact$fit)
          for (method in methods) {
            t0 <- proc.time()[["elapsed"]]
            r <- zeta_orders(sim$comm, sim$env, component = component, orders = k, M = 220,
                             ndraw = 49, seed = 200 + rep, sampling = method)
            sec <- proc.time()[["elapsed"]] - t0
            f <- r$point[[1]]
            # Heights and position effects are returned by the bootstrap pipeline.
            draw <- cbind(r$height[, 1, , drop = FALSE][, 1, ],
                          r$lambda[, 1, , drop = FALSE][, 1, ])
            interval <- apply(draw, 2, quantile, c(.025, .975), na.rm = TRUE)
            names(truth) <- names(coefficients(f))
            rows[[length(rows) + 1L]] <- data.frame(structure, rep, k, component, method,
              parameter = names(truth), truth = unname(truth), estimate = unname(coefficients(f)),
              covered = c(NA, truth[-1] >= interval[1, ] & truth[-1] <= interval[2, ]),
              width = c(NA, interval[2, ] - interval[1, ]), pshare = exact$pshare,
              d2 = f$D2, d2_exact = exact$fit$D2, seconds = sec,
              ess = r$point_diagnostics[1, "ess"], shared = r$point_diagnostics[1, "n_shared"],
              bootstrap_ess = median_finite(r$draw_diagnostics[, 1, "ess"]),
              bootstrap_shared = median_finite(r$draw_diagnostics[, 1, "n_shared"]))
          }
        }
      }
    }
  }
  write_result(do.call(rbind, rows), "recovery")
}

# Known-coefficient response with bounded, shared-site errors. The artificial
# response isolates parameter recovery and repeated-experiment coverage from
# misspecification of an ecological dissimilarity model.
known <- function() {
  rows <- list()
  beta <- c(0.35, 0.55, 0.25, 0.24, 0.35, 0.20)
  lambda <- c(0.07, -0.04)
  intercept <- 0.4
  true <- c(intercept, sum(beta[1:3]), sum(beta[4:6]), lambda)
  for (structure in c("broad", "narrow", "patchy")) {
    for (rep in seq_len(24)) {
      sim <- community(structure, 46, 1100 + rep)
      n <- nrow(sim$comm)
      set.seed(1200 + rep)
      u <- runif(n, -1, 1)
      W <- site_weights(n, 59)
      for (k in c(4, 14)) {
        for (method in methods) {
          set.seed(1300 + rep)
          one <- function(omega = NULL) {
            m <- 350
            if (method == "species") {
              s <- sample_sets_species(sim$comm, k, m, omega = omega)
              i <- s$idx
              w <- s$sampling_weights
            } else {
              i <- if (is.null(omega)) sample_sets(n, k, m) else sample_sets_weighted(omega, k, m)
              w <- NULL
            }
            d <- suppressWarnings(ugdm_design(sim$env, i))
            mu <- plogis(intercept + drop(d$Xdir %*% beta + d$Xpos %*% lambda))
            y <- mu + 0.8 * mu * (1 - mu) * rowMeans(matrix(u[i], nrow(i), k))
            coefficients(ugdm_fit(d, y, rep(1, m), sampling_weights = w))
          }
          point <- one()
          draws <- vapply(seq_len(ncol(W)), function(b) one(W[, b]), numeric(length(true)))
          interval <- apply(draws, 1, quantile, c(.025, .975))
          rows[[length(rows) + 1L]] <- data.frame(structure, rep, k, method,
              parameter = names(point), truth = true, estimate = unname(point),
              covered = true >= interval[1, ] & true <= interval[2, ],
              width = interval[2, ] - interval[1, ])
        }
      }
    }
  }
  write_result(do.call(rbind, rows), "known")
}

# High-order experiment. Each time trial executes the same end-to-end work
# (point and nine Bayesian-bootstrap fits with group Shapley) at its own M.
# Choose M by interpolation between separate timing pilots at M = 200 and 800.
high <- function() {
  rows <- list()
  for (structure in c("broad", "narrow", "patchy")) {
    message("Timing ", structure)
    for (k in c(3, 8, 16, 24)) {
      reference <- community(structure, 100, 501)
      fit_run <- function(sim, method, m, seed) {
        sec <- system.time(r <- suppressWarnings(zeta_orders(sim$comm, sim$env, component = "total",
          orders = k, M = m, ndraw = 9, seed = seed, groups = groups, sampling = method)))[["elapsed"]]
        list(r = r, sec = sec)
      }
      pilots <- lapply(methods, function(method) {
        vapply(c(200, 800), function(m) median(vapply(1:3, function(i) {
          fit_run(reference, method, m, 600 + i)$sec
        }, numeric(1))), numeric(1))
      })
      time_cap <- mean(vapply(pilots, function(t) mean(t), numeric(1)))
      for (rep in seq_len(6)) {
        sim <- community(structure, 100, 700 + rep)
        for (j in seq_along(methods)) {
          method <- methods[j]
          t <- pilots[[j]]
          budget <- round(200 + (time_cap - t[1]) * 600 / max(t[2] - t[1], 0.005))
          for (mode in c("equal_draws", "matched_time")) {
            m <- if (mode == "equal_draws") 500 else max(80, min(1800, budget))
            a <- fit_run(sim, method, m, 800 + rep)
            v <- suppressWarnings(ugdm_validate(sim$comm, sim$env, order = k, component = "total",
                M = m, M_test = 2000, nrep = 2, test_frac = 0.4, seed = 900 + rep, sampling = method,
                test_sampling = "conventional"))$summary
            p <- a$r$point_diagnostics[1, ]
            bd <- a$r$draw_diagnostics[, 1, ]
            rows[[length(rows) + 1L]] <- data.frame(structure, k, rep, method, mode, m, seconds = a$sec,
                time_cap, pshare = mean(zeta_sets(sim$comm, sample_sets(nrow(sim$comm), k, 4000, seed = 911),
                                                  abundance = FALSE)$shared > 0),
                occupancy_median = median(colMeans(sim$comm > 0)),
                occupancy_max = max(colMeans(sim$comm > 0)),
                d2 = a$r$point[[1]]$D2, ess = p[["ess"]], shared = p[["n_shared"]],
                alpha = p[["alpha"]], bootstrap_ess = median_finite(bd[, "ess"]),
                bootstrap_shared = median_finite(bd[, "n_shared"]),
                d2_test = median_finite(v$D2_test), rmse_test = median_finite(v$rmse_test))
          }
        }
      }
    }
  }
  write_result(do.call(rbind, rows), "high_orders")
}

# Independent throughput comparison with the same hard CPU budget for both
# samplers. Use the largest cumulative batch whose point fit completed in time.
deadline <- function() {
  rows <- list()
  budget <- 0.12
  cpu <- function() sum(proc.time()[1:2])
  for (structure in c("narrow", "patchy")) for (k in c(16, 24)) for (rep in seq_len(6)) {
    sim <- community(structure, 100, 1700 + rep)
    n <- nrow(sim$comm)
    set.seed(1800 + rep)
    held <- sort(sample.int(n, round(0.4 * n)))
    train <- setdiff(seq_len(n), held)
    x <- sim$comm[train, , drop = FALSE]
    test <- matrix(held[sample_sets(length(held), k, 3000, seed = 1900 + rep)], ncol = k)
    test_y <- ugdm_response(sim$comm, test, "total")
    test_d <- suppressWarnings(ugdm_design(sim$env, test))
    for (method in methods) {
      set.seed(2000 + rep)
      idx <- matrix(integer(), 0, k)
      weights <- numeric()
      best <- NULL
      t0 <- cpu()
      repeat {
        if (method == "species") {
          s <- sample_sets_species(x, k, 256)
        } else {
          s <- list(idx = sample_sets(nrow(x), k, 256), sampling_weights = rep(1, 256))
        }
        idx <- rbind(idx, matrix(train[s$idx], ncol = k))
        weights <- c(weights, s$sampling_weights)
        r <- ugdm_response(sim$comm, idx, "total")
        d <- suppressWarnings(ugdm_design(sim$env, idx))
        f <- ugdm_fit(d, r$y, r$trials, sampling_weights = weights)
        elapsed <- cpu() - t0
        if (elapsed > budget) break
        best <- list(idx = idx, weights = weights, fit = f, seconds = elapsed)
      }
      if (is.null(best)) stop("No fit completed within the common CPU budget.")
      z <- zeta_sets(sim$comm, best$idx, abundance = FALSE)$shared
      p <- ugdm_predict(best$fit, test_d)
      rows[[length(rows) + 1L]] <- data.frame(structure, k, rep, method, budget,
        m = nrow(best$idx), seconds = best$seconds, share = sum(z > 0),
        ess = sum(best$weights)^2 / sum(best$weights^2), d2 = best$fit$D2,
        rmse = sqrt(weighted.mean((test_y$y - p)^2, test_y$trials)))
    }
  }
  write_result(do.call(rbind, rows), "deadline")
}

birds <- function() {
  folder <- file.path(paper, "data_birds")
  read <- function(suffix) read.csv(file.path(folder, paste0("N63TBIRD-", suffix, ".csv")),
                                    row.names = 1, check.names = FALSE)
  x <- as.matrix(read("S"))
  e <- read("E")[rownames(x), c("AprMay", "DJF", "Urb", "Br", "Co", "Op")]
  d <- read("C")[, c("Long", "Lat")]
  keep <- complete.cases(e) & rowSums(x) > 0
  x <- x[keep, colSums(x) > 0, drop = FALSE]
  e <- e[keep, , drop = FALSE]
  d <- d[keep, , drop = FALSE]
  e[c("Urb", "Op")] <- lapply(e[c("Urb", "Op")], function(p) log1p(100 * p))
  list(comm = x, env = e, coords = d, groups = list(climate = c("AprMay", "DJF"),
                land_cover = c("Urb", "Br", "Co", "Op"), distance = "distance"))
}

plants <- function() {
  envir <- new.env()
  data("southwest", package = "gdm", envir = envir)
  a <- envir$southwest
  sites <- sort(unique(a$site))
  species <- sort(unique(as.character(a$species)))
  x <- matrix(0L, length(sites), length(species))
  x[cbind(match(a$site, sites), match(as.character(a$species), species))] <- 1L
  e <- a[match(sites, a$site), c("bio5", "bio6", "bio15", "bio19", "Long", "Lat")]
  keep <- complete.cases(e) & rowSums(x) > 0
  list(comm = x[keep, , drop = FALSE], env = e[keep, 1:4], coords = e[keep, 5:6],
       groups = list(temperature = c("bio5", "bio6"), precipitation = c("bio15", "bio19"), distance = "distance"))
}

butterflies <- function(path) {
  a <- readRDS(path)
  x <- as.matrix(a$comm)
  x[is.na(x) | x < 0] <- 0
  predictors <- a$predictors
  stopifnot(identical(as.character(rownames(x)), as.character(predictors$site)))
  vars <- c("Temp", "Prec", "hfp_1000", "modis_1000", "het_1000")
  env <- predictors[, vars]
  coords <- predictors[, c("x", "y")]
  keep <- complete.cases(env, coords) & rowSums(x) > 0
  x <- x[keep, , drop = FALSE]
  x <- x[, colSums(x) > 0, drop = FALSE]
  env <- env[keep, , drop = FALSE]
  for (v in vars) {
    y <- env[[v]]
    if (mean(((y - mean(y)) / sd(y))^3) > 1) env[[v]] <- log1p(y - min(y))
  }
  list(comm = x, env = env, coords = coords[keep, , drop = FALSE],
       groups = list(climate = c("Temp", "Prec"), land_cover = setdiff(vars, c("Temp", "Prec")),
                     distance = "distance"))
}

empirical <- function(only_butterflies = FALSE) {
  rows <- cv <- list()
  studies <- list()
  if (!only_butterflies) {
    if (file.exists(file.path(paper, "data_birds", "N63TBIRD-S.csv"))) studies$birds <- birds()
    if (requireNamespace("gdm", quietly = TRUE)) studies$plants <- plants()
  }
  butterfly_path <- Sys.getenv("ZETAGDUM_BUTTERFLIES")
  if (only_butterflies && !file.exists(butterfly_path)) stop("Set ZETAGDUM_BUTTERFLIES to the cleaned BDM RDS.")
  if (file.exists(butterfly_path)) studies$butterflies <- butterflies(butterfly_path)
  for (ds in names(studies)) {
    sim <- studies[[ds]]
    message("Fitting ", ds, ", sites ", nrow(sim$comm), ", species ", ncol(sim$comm))
    kb <- zeta_break(zeta_curve(sim$comm))
    orders <- if (ds == "birds") c(2, 8, 15, 28) else sort(unique(c(2, round((2 + kb) / 2), kb)))
    components <- if (ds == "birds") c("turnover", "balanced", "relative") else if (ds == "plants") c("total", "turnover") else c("turnover", "balanced")
    for (component in components) {
      for (method in methods) {
        r <- suppressWarnings(zeta_orders(sim$comm, sim$env, sim$coords, lonlat = TRUE,
              component = component, orders = orders, M = 2000, ndraw = 39, seed = 1,
              groups = sim$groups, sampling = method))
        for (i in seq_along(orders)) {
          sh <- apply(r$shapley[, i, , drop = FALSE], 3, median_finite)
          rows[[length(rows) + 1L]] <- data.frame(study = ds, component, method, order = orders[i],
            sites = nrow(sim$comm), species = ncol(sim$comm), break_order = kb,
            d2 = r$point[[i]]$D2, bootstrap_d2 = median_finite(r$D2[, i]),
            ess = r$point_diagnostics[i, "ess"], shared = r$point_diagnostics[i, "n_shared"],
            group = names(sh), shapley = unname(sh))
          if (orders[i] > round(nrow(sim$comm) * 0.3)) next
          v <- suppressWarnings(ugdm_validate(sim$comm, sim$env, sim$coords, lonlat = TRUE,
               component = component, order = orders[i], M = 2000, nrep = 10, seed = 100,
               sampling = method, test_sampling = "conventional"))$summary
          cv[[length(cv) + 1L]] <- data.frame(study = ds, component, method, order = orders[i], v)
        }
      }
    }
  }
  suffix <- if (only_butterflies) "_butterflies" else ""
  if (length(rows)) write_result(do.call(rbind, rows), paste0("empirical", suffix))
  if (length(cv)) write_result(do.call(rbind, cv), paste0("empirical_cv", suffix))
}

stages <- if (length(args) > 2L) args[3:length(args)] else c("recover", "known", "high", "deadline", "empirical")
for (stage in stages) {
  message("Stage ", stage)
  switch(stage, recover = recover(), known = known(), high = high(), deadline = deadline(), empirical = empirical(),
         butterflies = empirical(only_butterflies = TRUE),
         stop("Unknown validation stage: ", stage))
}
