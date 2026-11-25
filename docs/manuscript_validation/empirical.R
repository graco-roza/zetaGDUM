# Rscript docs/manuscript_validation/empirical.R OUT CASE [MC=8] [BB=100] [MODE=all]
source("docs/manuscript_validation/common.R")
args <- commandArgs(TRUE)
out <- args[1L]
case <- args[2L]
mc <- if (length(args) > 2L) as.integer(args[3L]) else 8L
B <- if (length(args) > 3L) as.integer(args[4L]) else 100L
mode <- if (length(args) > 4L) args[5L] else "all"
meta <- read_manifest()
meta <- meta[match(case, meta$id), , drop = FALSE]
d <- load_dd2(meta, Sys.getenv("DD2_DATA"))
n <- nrow(d$x)
p <- (d$x > 0) * 1L
range <- internal(".pair_distance_range")(d$dist)
orders <- as.integer(meta[c("low", "mid", "high")])
W <- site_weights(n, B, seed = 42006)
stages <- c("total", "relative")
responses <- function(s) cbind(total = ugdm_response(d$x, s$idx, "total", working_weights = "constant")$y,
                               relative = ugdm_response(d$x, s$idx, "relative", working_weights = "constant")$y)
context <- function(s) {
  design <- suppressWarnings(ugdm_design(d$env, s$idx, d$dist, dist_range = range))
  list(d = design, y = responses(s), w = s$sampling_weights,
       s = s, contrast = contrast_basis(d$env, design))
}
fit_all <- function(a, diagnose = FALSE) lapply(stages, function(response) {
  empirical_estimate(a$d, a$y[, response], a$w, a$contrast, diagnose)
})
as_draws <- function(l, j) do.call(rbind, lapply(l, function(x) x[[j]]$estimate))
target_sharing <- function(s) sum(s$sampling_weights * (s$shared > 0)) / sum(s$sampling_weights)
uncertainty <- mc_rows <- diagnostics <- timings <- list()
for (alpha in c(0, .5)) {
  label <- if (alpha == 0) "conventional" else "species"
  fits <- list()
  for (k in orders) {
    message(case, " ", label, " order ", k, " ", mode)
    for (M in c(10000L, 50000L)) {
      if (mode == "inference" && M != 10000L) next
      if (mode == "mc" && M != 50000L) next
      exact <- k == 2L && choose(n, 2L) <= M
      point_context <- context(draw(p, k, M, alpha, seed = 530000L + k))
      point <- fit_all(point_context, TRUE)
      nd <- if (mode == "inference" || exact) 1L else mc
      if (mode != "inference" || M == 10000L) for (r in seq_len(nd)) {
        sample <- if (r == 1L) point_context else context(draw(p, k, M, alpha,
          seed = 540000L + k + 1000L * r + M))
        rec <- if (r == 1L) point else fit_all(sample, TRUE)
        excluded <- matrix(sample$w, ncol = 1L)
        in_w <- inclusion(sample$s$idx, excluded, n)[, 1L]
        in_w2 <- inclusion(sample$s$idx, excluded^2, n)[, 1L]
        min_ess <- min((sum(sample$w) - in_w)^2 / (sum(sample$w^2) - in_w2))
        for (j in seq_along(stages)) {
          diagnostics[[length(diagnostics) + 1L]] <- data.frame(case, sampler = label, order = k,
            budget = M, replicate = r, response = stages[j], sets = nrow(sample$s$idx), exact,
            target_share = target_sharing(sample$s), proposal_share = mean(sample$s$shared > 0),
            expected_shared = sum(internal(".species_zeta_presence")(d$x, k)), ess = kish(sample$w),
            min_exclusion_ess = min_ess, rec[[j]]$diagnostics)
          mc_rows[[length(mc_rows) + 1L]] <- data.frame(case, sampler = label, order = k, budget = M,
            replicate = r, response = stages[j], term = names(rec[[j]]$estimate),
            estimate = unname(rec[[j]]$estimate))
        }
      }
      if (M != 10000L || mode == "mc") next
      a <- point_context
      exc <- omission(a$s$idx, a$w, n)
      s0 <- proc.time()[["elapsed"]]
      deleted <- lapply(exc$keep, function(keep) {
        o <- list(d = subset_design(a$d, keep), y = a$y[keep, , drop = FALSE], w = a$w[keep],
                  contrast = a$contrast)
        fit_all(o)
      })
      jk_time <- proc.time()[["elapsed"]] - s0
      s0 <- proc.time()[["elapsed"]]
      boot <- lapply(seq_len(B), function(b) {
        if (exact) {
          z <- a
          z$w <- internal(".product_sampling_weights")(W[, b], a$s$idx, a$w)
        } else {
          z <- context(draw(p, k, M, alpha, omega = W[, b], seed = 635000L + k + 101L * b))
        }
        fit_all(z)
      })
      bb_time <- proc.time()[["elapsed"]] - s0
      timings[[length(timings) + 1L]] <- data.frame(case, sampler = label, order = k,
        jk_time, bb_time, n_delete = n, n_boot = B)
      fits[[as.character(k)]] <- lapply(seq_along(stages), function(j) list(
        center = point[[j]]$estimate, bayesian = as_draws(boot, j), jackknife = as_draws(deleted, j)))
      # Long-form omission identifiers are joined across orders before forming changes.
      message("  min omission ESS ", round(min(exc$ess)), "; active coefficients ",
              paste(vapply(point, function(h) h$diagnostics$active, numeric(1L)), collapse = "/"))
    }
  }
  if (length(fits)) {
    for (j in seq_along(stages)) for (k in orders) {
      q <- fits[[as.character(k)]][[j]]
      if (k == orders[1]) from_values <- NA_integer_ else
        from_values <- unique(c(orders[1], orders[match(k, orders) - 1L]))
      for (from in unique(c(NA_integer_, from_values))) for (method in c("bayesian", "jackknife")) {
        if (is.na(from)) {
          point <- q$center
          draws <- q[[method]]
        } else {
          low <- fits[[as.character(from)]][[j]]
          point <- q$center - low$center
          draws <- q[[method]] - low[[method]]
        }
        i <- interval(point, draws, method)
        uncertainty[[length(uncertainty) + 1L]] <- data.frame(case, sampler = label,
          order = k, from, response = stages[j], method, i)
      }
    }
    # Relative-minus-incidence at every order and for the breakpoint-minus-pair change.
    for (k in orders) {
      q <- fits[[as.character(k)]]
      for (from in if (k == orders[1L]) NA_integer_ else c(NA_integer_, orders[1L])) {
        for (method in c("bayesian", "jackknife")) {
          point <- q[[2L]]$center - q[[1L]]$center
          draws <- q[[2L]][[method]] - q[[1L]][[method]]
          if (!is.na(from)) {
            a <- fits[[as.character(from)]]
            point <- point - (a[[2L]]$center - a[[1L]]$center)
            draws <- draws - (a[[2L]][[method]] - a[[1L]][[method]])
          }
          uncertainty[[length(uncertainty) + 1L]] <- data.frame(case, sampler = label,
            order = k, from, response = "relative_minus_total", method, interval(point, draws, method))
        }
      }
    }
  }
}
prefix <- paste(case, mode, sep = "_")
for (key in c("uncertainty", "mc_rows", "diagnostics", "timings")) {
  rows <- get(key)
  if (length(rows)) write_study(do.call(rbind, rows), out, paste(prefix, key, sep = "_"))
}
