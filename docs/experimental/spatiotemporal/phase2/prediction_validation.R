st2_m2_window <- function(model, panel, env, design, transform = identity) {
  # Predict only transitions actually required by test window pairs.
  p <- .st_panel(panel)
  starts <- sort(unique(unlist(lapply(seq_len(nrow(design$windows$from)), function(j) {
    seq.int(tail(design$windows$from[j, ], 1L), tail(design$windows$to[j, ], 1L) - 1L)
  }))))
  windows <- st2_windows(p$times, 1L, cbind(starts, starts + 1L))
  local <- transform(st2_transition_design(p, env, 1L, windows = windows, basis = model$design$basis))
  fitted <- st2_predict_dynamic(model, local, env)
  n <- dim(p$x)[1L]
  S <- dim(p$x)[3L]
  g <- array(fitted$gain, c(n, length(starts), S))
  l <- array(fitted$loss, c(n, length(starts), S))
  profiles <- lapply(seq_len(nrow(design$windows$from)), function(w) {
    J0 <- design$windows$from[w, ]
    J1 <- design$windows$to[w, ]
    at <- tail(J0, 1L)
    prob <- matrix(p$x[, at, ] > 0, n, S)
    for (t in seq.int(at, J1[1L] - 1L)) {
      i <- match(t, starts)
      prob <- prob * (1 - l[, i, ]) + (1 - prob) * g[, i, ]
    }
    for (t in J1[-1L]) prob <- prob * (1 - l[, match(t - 1L, starts), ])
    prob
  })
  if (design$windows$temporal_order == 1L) {
    mapping <- match(design$windows$from[design$rows$window[design$row_id], 1L], starts)
    row <- design$site + n * (mapping - 1L)
    return(lapply(fitted, function(a) a[row, , drop = FALSE]))
  }
  idx <- design$site + n * (design$rows$window[design$row_id] - 1L)
  prob <- do.call(rbind, profiles)[idx, , drop = FALSE]
  list(gain = prob, loss = 1 - prob)
}

st2_rectangle_quantities <- function(a, design, beta = TRUE) {
  if (nrow(a) != nrow(design$X)) stop("Prediction rows must match the design.")
  k <- design$k
  n <- nrow(design$sets)
  idx <- matrix(seq_len(nrow(a)), n, k, byrow = TRUE)
  X <- a > 0
  D <- zetaGDUM::dominance_matrix(a) > 0
  total <- rowSums(a)
  Q <- a / ifelse(total > 0, total, NA_real_)
  reduce <- function(low, min_fun) {
    z <- low[idx[, 1L], , drop = FALSE]
    for (j in seq_len(k)[-1L]) z <- min_fun(z, low[idx[, j], , drop = FALSE])
    z
  }
  species <- list(presence = reduce(X, `&`) * 1.0, abundance = reduce(Q, pmin),
                  dominance = reduce(D, `&`) * 1.0)
  out <- as.data.frame(lapply(species, rowSums))
  if (beta && k > 1L) {
    for (component in c("total", "turnover", "nestedness", "relative", "bray_curtis", "balanced", "gradient")) {
      q <- zetaGDUM::ugdm_response(a, idx, component)
      out[[component]] <- ifelse(q$keep, q$y, NA_real_)
    }
  }
  list(species = species, rectangle = out)
}

st2_scores <- function(prob, observed, selected) {
  a <- prob[selected]
  y <- observed[selected]
  keep <- is.finite(a)
  a <- a[keep]
  y <- y[keep]
  p <- pmax(1e-10, pmin(1 - 1e-10, a))
  c(n = length(a), observed_rate = mean(y), predicted_rate = mean(a),
    log_loss = -mean(ifelse(y, log(p), log1p(-p))), brier = mean((a - y)^2))
}

st2_score_predictions <- function(q, design, states, name) {
  k <- design$k
  H0 <- q$initial_shared
  H1 <- rowsum(states$y * 1L, design$row_id, reorder = FALSE) == k
  out <- list()
  add <- function(metric, values) data.frame(model = name, response = metric, as.list(values), row.names = NULL)
  out[[1]] <- add("shared_gain", st2_scores(q$final_shared, H1, !H0))
  out[[2]] <- add("shared_loss", st2_scores(1 - q$final_shared, !H1, H0))
  if (!is.null(q$local_gain)) {
    out[[3]] <- add("local_gain", st2_scores(q$local_gain, states$y, !states$z))
    out[[4]] <- add("local_loss", st2_scores(q$local_loss, !states$y, states$z))
  }
  if (!is.null(q$expansion)) {
    g <- rowsum((!states$z & states$y) * 1L, design$row_id, reorder = FALSE)
    l <- rowsum((states$z & !states$y) * 1L, design$row_id, reorder = FALSE)
    all <- matrix(TRUE, nrow(H0), ncol(H0))
    out[[5]] <- add("expansion", st2_scores(q$expansion, g > l, all))
    out[[6]] <- add("contraction", st2_scores(q$contraction, l > g, all))
    out[[7]] <- add("pure_redistribution", st2_scores(q$pure_redistribution, g == l & g > 0, all))
    out[[8]] <- add("balanced_redistribution",
                    c(n = sum(is.finite(q$redistribution)), observed_rate = mean(pmin(g, l)),
                      predicted_rate = mean(q$redistribution, na.rm = TRUE), log_loss = NA_real_,
                      brier = mean((q$redistribution - pmin(g, l))^2, na.rm = TRUE)))
  }
  do.call(rbind, out)
}

st2_validate_split <- function(sim, split, order_pairs, set_budget = 36L, draws = 4L, interaction = FALSE) {
  p <- sim$enabled
  n <- dim(p$x)[1]
  T <- dim(p$x)[2]
  sites <- seq_len(n)
  times <- seq_len(T)
  fold <- if (split %in% c("sites", "joint")) sort(sample.int(n, n %/% 3L)) else which(sim$blocks %in% c(5, 6))
  test_sites <- if (split == "forward") sites else fold
  train_sites <- if (split == "forward") sites else setdiff(sites, fold)
  cut <- ceiling(T * .65)
  train_times <- if (split %in% c("forward", "joint")) seq_len(cut) else times
  test_times <- if (split %in% c("forward", "joint")) seq.int(cut, T) else times
  train <- st2_subset_panel(p, train_sites, train_times)
  test <- st2_subset_panel(p, test_sites, test_times)
  train_env <- sim$env[, train_times, , drop = FALSE]
  test_env <- sim$env[, test_times, , drop = FALSE]
  ref <- st2_basis(train, train_env)
  dl <- st2_transition_design(train, train_env, 1L, basis = ref)
  m2 <- st2_fit_dynamic(train, dl, interaction = interaction)
  a2 <- st2_fit_positive(train, dl, hierarchical = TRUE, interaction = interaction)
  test_local <- st2_transition_design(test, test_env, 1L, basis = ref)
  states_local <- st2_states(test, test_local)
  m2_local <- st2_predict_dynamic(m2, test_local, test_env)
  m2_draws <- st2_predict_communities(m2_local, a2, test_local, states_local, draws)
  sizes <- list()
  scores <- quantities <- coherence <- calibration <- list()
  for (order_row in seq_len(nrow(order_pairs))) {
    k <- order_pairs[order_row, 1L]
    m <- order_pairs[order_row, 2L]
    if (m > min(length(train_times), length(test_times)) %/% 2L) next
    dt <- st2_transition_design(train, train_env, k, m, sets = st2_draw_sets(length(train_sites), k, set_budget), basis = ref)
    ds <- st2_transition_design(test, test_env, k, m, sets = st2_draw_sets(length(test_sites), k, set_budget), basis = ref)
    z <- st2_states(test, ds)
    model <- list(M0 = st2_fit_transition(train, dt, pooled = TRUE),
                  M1 = st2_fit_transition(train, dt),
                  M1_interaction = st2_fit_transition(train, dt, interaction = TRUE))
    p <- lapply(model, st2_predict_transition, design = ds, states = z)
    p$M2 <- st2_m2_window(m2, test, test_env, ds)
    alpha <- if (m == 1L) st2_fit_positive(train, dt) else NULL
    observed_shared <- rowsum(z$y * 1L, ds$row_id, reorder = FALSE) == k
    observed_amount <- if (m == 1L) st2_rectangle_quantities(z$after, ds)$rectangle else NULL
    for (name in names(p)) {
      q <- st2_transition_summary(p[[name]], ds, z, occupancy = name != "M0" && k > 1L)
      if (!is.null(p[[name]]$gain)) {
        q$local_gain <- p[[name]]$gain
        q$local_loss <- p[[name]]$loss
      }
      result <- st2_score_predictions(q, ds, z, name)
      H0 <- q$initial_shared
      for (event in c("gain", "loss")) {
        use <- if (event == "gain") !H0 else H0
        rate <- if (event == "gain") q$final_shared[use] else 1 - q$final_shared[use]
        observed <- if (event == "gain") observed_shared[use] else !observed_shared[use]
        good <- is.finite(rate)
        if (any(good)) {
          bin <- pmin(10L, floor(10 * rate[good]) + 1L)
          a <- aggregate(cbind(predicted = unname(rate[good]), observed = unname(observed[good]), n = rep(1L, sum(good))) ~ bin,
                         FUN = sum)
          calibration[[length(calibration)+1L]] <- cbind(split, spatial_order = k, temporal_order = m,
                                                        model = name, event, a)
        }
      }
      scores[[length(scores) + 1L]] <- cbind(split, spatial_order = k, temporal_order = m, result)
      presence_error <- rowSums(q$final_shared) - rowSums(observed_shared)
      values <- data.frame(response = "presence", mae = mean(abs(presence_error), na.rm = TRUE),
                           rmse = sqrt(mean(presence_error^2, na.rm = TRUE)), support = sum(is.finite(presence_error)))
      if (m == 1L && name %in% c("M1", "M2")) {
        forecasts <- if (name == "M1") {
          st2_predict_communities(p[[name]], alpha, ds, z, draws)
        } else {
          idx <- ds$site + length(test_sites) * (ds$rows$window[ds$row_id] - 1L)
          lapply(m2_draws, function(v) v[idx, , drop = FALSE])
        }
        rectangle <- lapply(forecasts, function(v) st2_rectangle_quantities(v, ds)$rectangle)
        for (metric in setdiff(names(observed_amount), "presence")) {
          v <- vapply(rectangle, `[[`, numeric(nrow(ds$sets)), metric)
          prediction <- rowMeans(v) # Common defined support in all predictive draws.
          error <- prediction - observed_amount[[metric]]
          values <- rbind(values, data.frame(response = metric, mae = mean(abs(error), na.rm = TRUE),
                                              rmse = sqrt(mean(error^2, na.rm = TRUE)), support = sum(is.finite(error))))
        }
      }
      quantities[[length(quantities) + 1L]] <- cbind(split, spatial_order = k, temporal_order = m, model = name, values)
      if (name %in% c("M0", "M1", "M2")) {
        key <- paste(name, k, m, sep = ":")
        sizes[[key]] <- list(design = ds, shared = q$final_shared)
        if (name != "M0" && k > 1L) {
          cc <- st2_overlap_coherence(ds, z, p[[name]])
          coherence[[length(coherence) + 1L]] <- data.frame(split, model = name, spatial_order = k,
            temporal_order = m, comparison = "overlap", tested = unname(cc[1L]), violations = unname(cc[2L]), maximum_excess = unname(cc[3L]),
            row.names = NULL)
        }
      }
    }
  }
  # A common panel-origin trajectory supplies the higher-order nesting benchmark.
  forecast <- st2_predict_trajectory(m2, a2, test, test_env)
  for (k in 1:3) for (m in 1:3) for (name in c("M0", "M1", "M2")) {
    small <- sizes[[paste(name, k, m, sep = ":")]]
    if (is.null(small)) next
    for (axis in c("space", "time")) {
      large <- sizes[[paste(name, k + as.integer(axis == "space"), m + as.integer(axis == "time"), sep = ":")]]
      if (is.null(large)) next
      if (name == "M2") {
        small$shared <- st2_observed_projection(forecast, small$design)$final_shared
        large$shared <- st2_observed_projection(forecast, large$design)$final_shared
      }
      ch <- st2_shared_coherence(small, large)
      coherence[[length(coherence) + 1L]] <- data.frame(split, model = name, spatial_order = k,
        temporal_order = m, comparison = axis, tested = unname(ch[1L]), violations = unname(ch[2L]), maximum_excess = unname(ch[3L]),
        row.names = NULL)
    }
  }
  list(scores = do.call(rbind, scores), quantities = do.call(rbind, quantities), coherence = do.call(rbind, coherence),
       calibration = do.call(rbind, calibration), dynamic_status = c(gain = m2$fit$gain$status, loss = m2$fit$loss$status))
}
