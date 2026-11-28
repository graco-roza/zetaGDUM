.st2_local_blocks <- function(design, blocks) {
  if (length(blocks) != length(design$sites) || anyNA(blocks)) stop("Supply one block per site.")
  matrix(as.character(blocks[design$sets]), nrow(design$sets), ncol(design$sets))
}

.st2_moving_counts <- function(n, length) {
  start <- sample.int(n - length + 1L, ceiling(n / length), replace = TRUE)
  path <- as.vector(outer(seq_len(length) - 1L, start, "+"))
  tabulate(head(path, n), n)
}

st2_transition_uncertainty <- function(model, event, blocks, draws = 40L,
                                       temporal_length = 2L, exact_delete = FALSE,
                                       terms = names(model$fit[[event]]$coefficients)) {
  design <- model$design
  f <- model$fit[[event]]
  if (f$status != "estimated") return(NULL)
  g <- model$groups[[event]]
  all_terms <- names(f$coefficients)
  if (!all(terms %in% all_terms)) stop("Choose fitted coefficient names.")
  keep <- g$n > 0
  X <- g$X[keep, , drop = FALSE]
  local <- g$local[keep]
  indexing <- design$row_id[local]
  w <- g$weight[keep]
  n <- g$n[keep]
  y <- g$y[keep]
  p <- plogis(drop(X %*% f$coefficients))
  score <- X * (w * (y - n * p))
  curvature_weight <- w * n * p * (1 - p)
  h <- crossprod(X * sqrt(curvature_weight)) + diag(model$ridge, ncol(X))
  inverse <- solve(h)
  block_members <- .st2_local_blocks(design, blocks)[indexing, , drop = FALSE]
  site_members <- .st2_local_blocks(design, seq_along(design$sites))[indexing, , drop = FALSE]
  estimate <- function(multiplier, exact = FALSE) {
    if (exact) {
      full <- numeric(length(keep))
      full[keep] <- multiplier
      return(.st2_fit_grouped(g, model$ridge, full, diagnostics = FALSE)$coefficients[terms])
    }
    update <- inverse %*% colSums(score * (multiplier - 1))
    f$coefficients[terms] + drop(update)[match(terms, all_terms)]
  }
  delete_variance <- function(membership) {
    ids <- unique(as.vector(membership))
    if (length(ids) < 3L) return(rep(NA_real_, length(terms)))
    values <- vapply(ids, function(id) {
      retained <- rowSums(membership == id) == 0
      if (exact_delete) return(estimate(retained, TRUE))
      # Newton update at the full fit, with curvature on retained rows.
      dropped <- X[!retained, , drop = FALSE]
      h_retained <- h - crossprod(dropped * sqrt(curvature_weight[!retained]))
      shift <- solve(h_retained, -colSums(score[!retained, , drop = FALSE]))
      f$coefficients[terms] + shift[match(terms, all_terms)]
    }, numeric(length(terms)))
    values <- matrix(values, nrow = length(terms))
    sqrt((length(ids) - 1) * rowMeans((values - rowMeans(values))^2))
  }
  horizon <- nrow(design$windows$from)
  if (temporal_length > horizon || temporal_length < 1L || draws < 2L) stop("Invalid bootstrap length or draw count.")
  w_index <- design$rows$window[indexing]
  chrono <- order(design$windows$starts[, 2L], design$windows$starts[, 1L])
  ids <- unique(as.vector(block_members))
  # Every record in a combination receives the product over its distinct
  # spatial blocks, preserving all site/species histories within a draw.
  unique_member <- lapply(seq_len(ncol(block_members)), function(j) {
    if (j == 1L) return(rep(TRUE, nrow(block_members)))
    rowSums(block_members[, seq_len(j - 1L), drop = FALSE] == block_members[, j]) == 0
  })
  member_index <- matrix(match(block_members, ids), nrow(block_members))
  temporal <- crossed <- matrix(NA_real_, length(terms), draws)
  for (b in seq_len(draws)) {
    t <- integer(horizon)
    t[chrono] <- .st2_moving_counts(horizon, temporal_length)
    time <- t[w_index]
    temporal[, b] <- estimate(time)
    block_count <- tabulate(sample.int(length(ids), length(ids), replace = TRUE), length(ids))
    space <- rep(1, length(time))
    for (j in seq_len(ncol(block_members))) {
      take <- unique_member[[j]]
      space[take] <- space[take] * block_count[member_index[take, j]]
    }
    crossed[, b] <- estimate(time * space)
  }
  se <- cbind(naive = unname(f$working_independent_se[terms]), delete_site = delete_variance(site_members),
              delete_spatial_block = delete_variance(block_members),
              moving_temporal = apply(temporal, 1L, stats::sd),
              crossed = apply(crossed, 1L, stats::sd))
  rownames(se) <- terms
  list(estimate = f$coefficients[terms], se = se, exact_delete = exact_delete, bootstrap_draws = draws)
}

# Profile one group of design columns against all nuisance columns.
.st2_profile_information <- function(model, event, prefix) {
  f <- model$fit[[event]]
  if (f$status != "estimated") return(NA_real_)
  g <- model$groups[[event]]
  p <- plogis(drop(g$X %*% f$coefficients))
  weight <- g$n * g$weight * p * (1 - p)
  X <- g$X * sqrt(weight)
  tested <- startsWith(colnames(X), prefix)
  if (!any(tested)) return(NA_real_)
  conditional <- qr.resid(qr(X[, !tested, drop = FALSE]), X[, tested, drop = FALSE])
  sum(conditional^2) / (sum(g$n * g$weight) * sum(tested))
}

# Uniform sampled-rectangle uncertainty is computational error conditional on
# the panel. Species weights and means are exact, independently of sampling.
st2_order_diagnostic <- function(x, env, spatial_max, temporal_max, budget = 60L,
                                 mc_replicates = 3L, basis = NULL) {
  p <- .st_panel(x)
  d <- dim(p$x)
  if (temporal_max > d[2] %/% 2L || budget < 2L || mc_replicates < 2L) stop("Increase budget or reduce orders.")
  if (is.null(basis)) basis <- st2_basis(p, env)
  exact <- st2_zeta(p, spatial_max, temporal_max, mode = "consecutive", empty = "undefined")
  community <- st2_community(p)
  out <- species <- list()
  for (k in seq_len(spatial_max)) for (m in seq_len(temporal_max)) {
    started <- proc.time()[["elapsed"]]
    n <- d[1]
    exhaustive <- choose(n, k) <= budget
    subsets <- function() {
      if (exhaustive) .st_indices(n, k, budget) else
        t(matrix(vapply(seq_len(budget), function(i) sort(sample.int(n, k)), integer(k)), k, budget))
    }
    order_rows <- exact$surface$spatial_order == k & exact$surface$temporal_order == m
    surface <- exact$surface[order_rows, , drop = FALSE]
    index <- (m - 1L) * spatial_max + k
    ts <- .st2_time_sets(p$times, m, "consecutive", 1e5)
    rectangle <- function(sets) {
      z <- vapply(seq_len(nrow(sets) * nrow(ts)), function(r) {
        i <- ((r - 1L) %% nrow(sets)) + 1L
        j <- ((r - 1L) %/% nrow(sets)) + 1L
        a <- sets[i, ]; b <- ts[j, ]
        c(presence = sum(colSums(community$X[a,b,,drop=FALSE], dims = 2L) == k*m),
          abundance = sum(apply(community$Q[a,b,,drop=FALSE], 3L, min)),
          dominance = sum(colSums(community$D[a,b,,drop=FALSE], dims = 2L) == k*m))
      }, numeric(3L))
      c(colMeans(t(z)), rowMeans(z == 0))
    }
    sets <- subsets()
    sampled <- vapply(seq_len(mc_replicates), function(rep) rectangle(if (rep == 1L) sets else subsets()), numeric(6L))
    design <- st2_transition_design(p, env, k, m, sets = sets, basis = basis)
    dominant <- st_panel(community$D, p$times)
    models <- list(presence = st2_fit_transition(p, design),
                   dominance = st2_fit_transition(dominant, design))
    mem <- sum(vapply(models, function(x) as.numeric(utils::object.size(x)), numeric(1)))
    keep_row <- !duplicated(design$row_id)
    terms <- setdiff(.st2_terms(design, FALSE), grep("^deviation", colnames(design$X), value = TRUE))
    X <- cbind(intercept = 1, design$X[keep_row, terms, drop = FALSE])
    relative <- st_panel(ifelse(is.na(community$Q), 0, community$Q), p$times)
    site_min <- st2_states(relative, design)
    collapse <- function(a) {
      out <- a[keep_row, , drop = FALSE]
      for (j in seq_len(k)[-1L]) out <- pmin(out, a[which(keep_row) + j - 1L, , drop = FALSE])
      out
    }
    a <- collapse(site_min$before)
    b <- collapse(site_min$after)
    empty <- function(windows) {
      out <- vapply(seq_len(nrow(windows)), function(j) rowSums(community$empty[, windows[j, ], drop = FALSE]) > 0, logical(n))
      matrix(out, n)
    }
    w <- design$rows$window[design$row_id]
    undefined <- empty(design$windows$from)[cbind(design$site, w)] |
      empty(design$windows$to)[cbind(design$site, w)]
    valid_row <- rowsum(undefined * 1L, design$row_id, reorder = FALSE)[, 1L] == 0
    for (response in surface$response) for (event in c("gain","loss")) {
      i <- match(response, surface$response)
      if (response == "abundance") {
        # Diagnostic opportunity: a minimum can rise if < 1 and fall if > 0.
        # Profile a weighted Gaussian working model of the eligible change amounts.
        opportunities <- if (event == "gain") a < 1 else a > 0
        success <- if (event == "gain") b > a else b < a
        amount <- if (event == "gain") pmax(b - a, 0) else pmax(a - b, 0)
        weight <- rowSums(opportunities)
        adequate <- valid_row & weight > 0
        A <- X[adequate, , drop = FALSE] * sqrt(weight[adequate])
        y <- rowSums(amount * opportunities)[adequate] / sqrt(weight[adequate])
        rank <- if (nrow(A)) qr(A)$rank else 0L
        status <- if (!nrow(A)) "no_opportunities" else if (rank < ncol(X)) "aliased" else "working_gaussian"
        profile <- function(prefix) {
          tested <- startsWith(colnames(A), prefix)
          if (!any(tested) || status != "working_gaussian") return(NA_real_)
          errors <- qr.resid(qr(A), y)
          variance <- sum(errors^2) / (nrow(A) - rank)
          if (!is.finite(variance) || variance <= 0) return(NA_real_)
          Q <- qr.resid(qr(A[, !tested, drop = FALSE]), A[, tested, drop = FALSE])
          sum(Q^2) / (sum(weight[adequate]) * sum(tested) * variance)
        }
        count <- sum(weight[adequate])
        events <- sum(success[adequate, , drop = FALSE] & opportunities[adequate, , drop = FALSE])
        profile_pos <- profile("position.")
        profile_spread <- profile("spread.")
        position <- grep("^position\\.", colnames(A))[1L]
        if (status == "working_gaussian" && !is.na(position) && n > k && n >= 3L) {
          retained_sets <- design$sets[adequate, , drop = FALSE]
          deleted <- vapply(seq_len(n), function(j) {
            use <- rowSums(retained_sets == j) == 0
            if (!any(use)) return(NA_real_)
            fit <- stats::lm.fit(A[use, , drop = FALSE], y[use])
            if (fit$rank < ncol(A)) return(NA_real_)
            fit$coefficients[position]
          }, numeric(1))
          se <- sqrt((n - 1L) * mean((deleted - mean(deleted))^2))
        } else se <- NA_real_
      } else {
        fit <- models[[response]]
        f <- fit$fit[[event]]
        influence <- if (nrow(design$windows$from) > 1L) {
          st2_transition_uncertainty(fit, event, rep(seq_len(5L), length.out = n),
                                      draws = 12L, temporal_length = 2L,
                                      terms = grep("^position\\.", names(f$coefficients), value = TRUE)[1L])
        } else NULL
        status <- f$status
        rank <- f$rank
        count <- f$opportunities
        events <- f$events
        profile_pos <- .st2_profile_information(fit, event, "position.")
        profile_spread <- .st2_profile_information(fit, event, "spread.")
        se <- if (is.null(influence)) NA_real_ else influence$se[1L, "delete_site"]
      }
      out[[length(out) + 1L]] <- data.frame(surface[i, , drop = FALSE], event, opportunities = count,
        events, status, rank, profiled_position = profile_pos, profiled_spread = profile_spread,
        position_delete_site_se = se, zero_fraction = mean(sampled[i + 3L, ]),
        combination_mc_sd = stats::sd(sampled[i, ]), set_count = nrow(sets), exhaustive,
        elapsed_seconds = proc.time()[["elapsed"]] - started, fit_bytes = mem, row.names = NULL)
    }
    species[[length(species) + 1L]] <- do.call(rbind, lapply(names(exact$weights), function(response) {
      data.frame(spatial_order = k, temporal_order = m, response,
                 species = colnames(exact$weights[[response]]),
                 weight = exact$weights[[response]][index, ], row.names = NULL)
    }))
  }
  list(orders = do.call(rbind, out), weights = do.call(rbind, species))
}
