source("docs/phase3/common.R")
set.seed(39284)

ensure_nonempty <- function(x) {
  empty <- which(rowSums(x) == 0)
  if (length(empty)) {
    singleton <- matrix(0L, nrow(x), length(empty))
    singleton[cbind(empty, seq_along(empty))] <- 1L
    x <- cbind(x, singleton)
  }
  x
}

run_cell <- function(s, k, M, R, ref_M, scenario) {
  n <- nrow(s$comm)
  ref_idx <- if (ref_M == 0L) t(combn(n, k)) else iid_sets(n, k, ref_M)
  basis <- ugdm_design(s$env, ref_idx)
  out <- list()
  audit_pairs <- iid_sets(n, 2, 160L)
  for (type in c("relative", "incidence")) {
    response_kind <- if (type == "relative") type else "total"
    rr <- ugdm_response(s$comm, ref_idx, response_kind)
    fr <- ugdm_fit(basis, rr$y, rr$trials)
    # All methods are evaluated on a common design to isolate fitting error.
    score_fit <- function(fit) {
      imp <- ugdm_importance(fit, basis, trials = rr$trials)
      local <- ugdm_local_distinctiveness(s$comm, fit, basis,
        if (type == "relative") "relative" else "normalised")
      c(fit$height, fit$lambda, imp,
        loss = logistic_loss(rr$y, ugdm_predict(fit, basis), rr$trials),
        prediction = ugdm_predict(fit, basis), local_spread = local$spread, local_position = local$position)
    }
    target <- score_fit(fr)
    # Enumerate 2 variable coalitions for comparable weighted deviance Shapley values.
    shapley <- function(full, d, r, idx, w, method, ms) {
      null <- logistic_loss(rr$y, rep(weighted.mean(rr$y, rr$trials), length(rr$y)), rr$trials)
      entropy <- logistic_loss(rr$y, rr$y, rr$trials)
      if (null - entropy < 1e-12) return(rep(NA_real_, 2))
      full_gain <- (null - logistic_loss(rr$y, ugdm_predict(full, basis), rr$trials)) / (null - entropy)
      one <- vapply(seq_along(s$env), function(j) {
        if (method == "conventional") {
          f <- ugdm_fit(d, r$y, r$trials, vars = names(s$env)[j], sampling_weights = w)
        } else {
          single <- ugdm_design(s$env[j], idx)
          nm <- c("intercept", colnames(single$Xdir), colnames(single$Xpos))
          single_ms <- ms
          single_ms$response <- ms$response[nm]
          single_ms$predictor <- ms$predictor[nm]
          f <- ugdm_fit_analytical(s$comm, s$env[j], single, type, w, method = method, moments = single_ms)
        }
        (null - logistic_loss(rr$y, ugdm_predict(f, basis), rr$trials)) / (null - entropy)
      }, numeric(1))
      .5 * (one + full_gain - rev(one))
    }
    ref_phi <- shapley(fr, basis, rr, ref_idx, NULL, "conventional", NULL)
    tm <- system.time(ms <- ugdm_crossproducts(s$comm, s$env, basis, type))[["elapsed"]]
    ref_predict <- grepl("^prediction", names(target))
    ref_local <- grepl("^local_", names(target))
    for (repeat_id in seq_len(R)) {
      set.seed(938 + n * 100000L + k * 1000L + repeat_id)
      methods <- c("iid", "species", "site", "pair", "environment", "hybrid", "control_variate")
      if (k > 2L && n > 3 && any(n %% (2:floor(sqrt(n))) == 0)) methods <- setdiff(methods, "pair")
      for (method in methods) {
        analytical <- method %in% c("hybrid", "control_variate")
        set.seed(938 + n * 100000L + k * 1000L + repeat_id)
        elapsed <- system.time({
          w <- rep(1, M)
          if (method == "species") {
            sets <- sample_sets_species(s$comm, k, M)
            idx <- sets$idx
            w <- sets$sampling_weights
          } else if (method %in% c("site", "pair", "environment")) {
            idx <- sample_sets_balanced(n, k, M, method, ceiling(rank(s$env[[1]]) * 4 / n))
          } else idx <- iid_sets(n, k, M)
          design <- ugdm_design(s$env, idx, ref = basis)
          r <- ugdm_response(s$comm, idx, response_kind)
          fit <- if (analytical) ugdm_fit_analytical(s$comm, s$env, design, type, w,
            method = method, moments = ms) else ugdm_fit(design, r$y, r$trials, sampling_weights = w)
        })[["elapsed"]]
        value <- score_fit(fit)
        phi <- shapley(fit, design, r, idx, w, if (analytical) method else "conventional", ms)
        diag <- tabulate(idx, n)
        drop_w <- zetaGDUM:::.inclusion_sums(idx, cbind(w, w^2), n)
        w0 <- sum(w) - drop_w[, 1]
        w2 <- sum(w^2) - drop_w[, 2]
        info <- zetaGDUM:::.order_information(design, fit, w * r$trials)
        shared <- zeta_sets(s$comm, idx, abundance = FALSE)$shared > 0
        present <- matrix(FALSE, nrow(idx), n)
        present[cbind(rep(seq_len(nrow(idx)), ncol(idx)), as.vector(idx))] <- TRUE
        pairs <- colSums(present[, audit_pairs[, 1L]] & present[, audit_pairs[, 2L]])
        out[[length(out) + 1L]] <- data.frame(scenario, n, k, M, repeat_id, type, method,
          used = if (analytical) fit$analytical$used else method,
          converged = fit$converged, rank = info$rank, parameters = info$parameters,
          active = info$active, sharing = weighted.mean(shared, w), nonzero = sum(shared),
          ess = sum(w)^2 / sum(w^2), min_delete_ess = min(ifelse(w2 > 0, w0^2 / w2, 0)),
          inclusion_cv = sd(diag) / mean(diag), pair_cv = sd(pairs) / mean(pairs),
          elapsed, moments_seconds = if (analytical) tm else 0,
          design_kib = as.numeric(object.size(design)) / 1024,
          moments_kib = if (analytical) as.numeric(object.size(ms)) / 1024 else 0,
          spread_error = unname(value[1] - target[1]), position_error = unname(value[3] - target[3]),
          importance_error = unname(value[5] - target[5]), phi_error = phi[1] - ref_phi[1],
          same_importance_rank = (value[5] > value[6]) == (target[5] > target[6]),
          prediction_rmse = sqrt(mean((value[ref_predict] - target[ref_predict])^2)),
          local_shapley_rmse = sqrt(mean((value[ref_local] - target[ref_local])^2)),
          loss_difference = unname(value["loss"] - target["loss"]), sampled_D2 = fit$D2)
      }
    }
  }
  do.call(rbind, out)
}

args <- commandArgs(TRUE)
mode <- if (length(args)) args[[1]] else "all"
stopifnot(mode %in% c("all", "exact", "scale"))
if (mode != "scale") {
  sim <- simulate_metacommunity(13, 35, breadth = 0.85, spatial = 0.25, seed = 49)
  sim$comm <- ensure_nonempty(sim$comm)
  result <- do.call(rbind, lapply(c(3, 5), function(k) {
    do.call(rbind, lapply(c(120, 600), function(M) run_cell(sim, k, M, 24, 0, "small_exhaustive")))
  }))
  write_result(result, "fitting_exact.csv")
}
if (mode != "exact") {
  specs <- read.csv("docs/manuscript_validation/communities.csv")
  out <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    n <- specs$sites[i]
    sim <- simulate_metacommunity(n, specs$species[i], breadth = c(.35,.35,.45,.35,.55)[i],
                                  abundance_mean = c(1.4,1,1,1,1)[i],
                                  spatial = c(.75,.2,.6,.35,.85)[i], seed = 741 + i)
    sim$comm <- ensure_nonempty(sim$comm)
    out[[i]] <- run_cell(sim, specs$high[i], 1500, 10, 30000, specs$community[i])
  }
  write_result(do.call(rbind, out), "fitting_scale.csv")
}
