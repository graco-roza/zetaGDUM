.st2_group_sums <- function(matrix, group) {
  cbind(overall = rowSums(matrix), t(rowsum(t(matrix), group, reorder = TRUE)))
}

st2_group_surface <- function(p, group, kmax, mmax) {
  z <- st2_zeta(p, kmax, mmax, mode = "consecutive", empty = "undefined")
  do.call(rbind, lapply(names(z$contributions), function(response) {
    values <- .st2_group_sums(z$contributions[[response]], group)
    a <- expand.grid(spatial_order = seq_len(kmax), temporal_order = seq_len(mmax),
                     group = 0:8, KEEP.OUT.ATTRS = FALSE)
    data.frame(a, response, value = as.vector(values))
  }))
}

st2_factorial_effect <- function(sim, replicate, kmax = 3L, mmax = 3L, budget = 30L, draws = 4L,
                                 interaction = sim$generator$scenario == "interaction") {
  p <- sim$enabled
  n <- dim(p$x)[1]
  T <- dim(p$x)[2]
  test_sites <- sort(sample.int(n, n %/% 3L))
  train <- st2_subset_panel(p, setdiff(seq_len(n), test_sites))
  test <- st2_subset_panel(p, test_sites)
  baseline <- st2_subset_panel(sim$disabled, test_sites)
  ref <- st2_basis(train, sim$env)
  dynamic_design <- st2_transition_design(train, sim$env, 1L, basis = ref)
  m2 <- st2_fit_dynamic(train, dynamic_design, interaction = interaction)
  a2 <- st2_fit_positive(train, dynamic_design, hierarchical = TRUE, interaction = interaction)
  a <- st2_group_surface(test, sim$group, kmax, mmax)
  b <- st2_group_surface(baseline, sim$group, kmax, mmax)
  truth <- a$value - b$value
  shape <- c(length(test_sites), T - 1L, length(sim$group))
  draws_m2 <- vapply(seq_len(draws), function(j) {
    u <- array(runif(prod(shape)), shape)
    ua <- array(runif(prod(shape)), shape)
    prediction <- st2_predict_trajectory(m2, a2, baseline, sim$env, u, ua)
    reference <- st2_predict_trajectory(m2, a2, baseline, sim$env, u, ua, transform = st2_disable_change)
    st2_group_surface(prediction, sim$group, kmax, mmax)$value -
      st2_group_surface(reference, sim$group, kmax, mmax)$value
  }, numeric(nrow(a)))
  a$value <- NULL
  result <- list(cbind(a, architecture = "M2", estimand = "trajectory_surface", truth,
                       estimate = rowMeans(draws_m2), prediction_mc_sd = apply(draws_m2, 1L, sd),
                       n_windows = NA_integer_, stringsAsFactors = FALSE))
  for (m in seq_len(mmax)) {
    windows <- st2_windows(test$times, m)
    # Comparable conditional estimand for set models: both simulated outcome
    # windows start from the same last observed assemblage in the initial window.
    replays <- lapply(seq_len(nrow(windows$from)), function(w) {
      origin <- tail(windows$from[w, ], 1L)
      final <- tail(windows$to[w, ], 1L)
      before <- sim$disabled$x[, origin, ]
      list(enabled = st2_factorial_replay(sim$generator, TRUE, before, origin, final)$x[, -1L, , drop = FALSE],
           disabled = st2_factorial_replay(sim$generator, FALSE, before, origin, final)$x[, -1L, , drop = FALSE])
    })
    for (k in seq_len(kmax)) {
      d0 <- st2_transition_design(train, sim$env, k, m,
                                   sets = st2_draw_sets(dim(train$x)[1], k, budget), basis = ref)
      d1 <- st2_transition_design(baseline, sim$env, k, m, windows = windows,
                                   sets = st2_draw_sets(length(test_sites), k, budget), basis = ref)
      z <- st2_states(baseline, d1)
      w <- d1$rows$window[d1$row_id]
      observed <- lapply(c("enabled", "disabled"), function(state) {
        ys <- do.call(rbind, lapply(replays, function(r) .st2_minimum(r[[state]], seq_len(m))))
        ys[test_sites[d1$site] + n * (w - 1L), , drop = FALSE]
      })
      rectangle <- function(a) {
        if (m == 1L) return(st2_rectangle_quantities(a, d1, beta = FALSE)$species)
        list(presence = (rowsum((a > 0) * 1L, d1$row_id, reorder = FALSE) == k) * 1.0)
      }
      observed <- lapply(observed, rectangle)
      dd <- st2_disable_change(d1)
      fit <- list(M0 = st2_fit_transition(train, d0, interaction = interaction, pooled = TRUE),
                  M1 = st2_fit_transition(train, d0, interaction = interaction))
      prediction <- lapply(fit, function(model) lapply(list(d1, dd), function(d) {
        st2_predict_transition(model, d, z)
      }))
      prediction$M2 <- lapply(list(identity, st2_disable_change), function(transform) {
        st2_m2_window(m2, baseline, sim$env, d1, transform)
      })
      for (name in names(prediction)) {
        estimates <- prediction[[name]]
        shared <- lapply(estimates, function(v) st2_transition_summary(v, d1, z, FALSE)$final_shared)
        effects <- list(presence = .st2_group_sums(shared[[1L]] - shared[[2L]], sim$group))
        draw_sd <- list()
        if (m == 1L && name == "M1") {
          positive <- st2_fit_positive(train, d0, interaction = interaction)
          dshape <- c(dim(z$z), draws)
          u <- array(runif(prod(dshape)), dshape)
          ua <- array(runif(prod(dshape)), dshape)
          v1 <- st2_predict_communities(estimates[[1L]], positive, d1, z, draws, u, ua)
          v0 <- st2_predict_communities(estimates[[2L]], positive, dd, z, draws, u, ua)
          for (response in c("abundance", "dominance")) {
            est <- lapply(seq_len(draws), function(b) {
              delta <- rectangle(v1[[b]])[[response]] - rectangle(v0[[b]])[[response]]
              .st2_group_sums(delta, sim$group)
            })
            effects[[response]] <- Reduce("+", est) / draws
            draw_sd[[response]] <- apply(vapply(est, colMeans, numeric(9L)), 1L, sd)
          }
        }
        for (response in names(effects)) {
          t0 <- .st2_group_sums(observed[[1L]][[response]] - observed[[2L]][[response]], sim$group)
          estimate <- effects[[response]]
          observed_support <- is.finite(t0[, 1L])
          support <- observed_support & is.finite(estimate[, 1L])
          mean_supported <- function(y) if (any(support)) colMeans(y[support, , drop = FALSE]) else rep(NA_real_, 9L)
          result[[length(result) + 1L]] <- data.frame(spatial_order = k, temporal_order = m,
            group = 0:8, response, architecture = name, estimand = "conditional_window",
            truth = unname(mean_supported(t0)), estimate = unname(mean_supported(estimate)),
            prediction_mc_sd = if (is.null(draw_sd[[response]])) NA_real_ else draw_sd[[response]],
            n_windows = sum(support), row.names = NULL)
        }
      }
    }
  }
  cbind(scenario = sim$generator$scenario, replicate, do.call(rbind, result))
}
