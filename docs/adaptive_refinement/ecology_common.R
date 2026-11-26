source("docs/manuscript_validation/common.R")

# Ordered cuts of projections of the geographic distance matrix. Each held-out
# group is a contiguous slab along a geographic direction; sizes differ by at
# most one. All candidate penalties use the same training and held-out paths.
study_blocks <- function(distance, partitions = 3L, groups = 4L) {
  coords <- stats::cmdscale(distance, k = 2L, add = TRUE)$points
  angle <- seq(0, pi, length.out = partitions + 1L)[seq_len(partitions)]
  lapply(angle, function(a) {
    z <- drop(coords %*% c(cos(a), sin(a)))
    cut(rank(z, ties.method = "first"), breaks = seq(0, length(z), length.out = groups + 1L),
        include.lowest = TRUE, labels = FALSE)
  })
}

study_idx <- function(sites, k, M) {
  matrix(sites[sample_sets(length(sites), k, M)], ncol = k)
}

study_contrast <- function(d, climate = 1L, neighbours = 30L) {
  environmental <- cbind(d$mean_position[, climate], sqrt(d$position_variance[, climate]))
  centre <- colMeans(environmental)
  scale <- apply(environmental, 2L, stats::sd)
  n <- nrow(environmental)
  if (n < 4L || any(!is.finite(scale) | scale < 1e-12)) {
    return(list(L = matrix(NA_real_, 3L, n), gap = NA_real_))
  }
  standardised <- sweep(sweep(environmental, 2L, centre), 2L, scale, "/")
  percentiles <- apply(standardised, 2L, stats::quantile, probs = c(.30, .70), names = FALSE)
  goal <- expand.grid(position = percentiles[, 1L], spread = percentiles[, 2L])
  n_each <- min(neighbours, n %/% 4L)
  W <- matrix(0, 4L, n)
  used <- integer()
  gap <- numeric(4L)
  for (j in seq_len(nrow(goal))) {
    distance <- sqrt(rowSums(sweep(standardised, 2L, unlist(goal[j, ]))^2))
    near <- utils::head(setdiff(order(distance), used), n_each)
    W[j, near] <- 1 / n_each
    used <- c(used, near)
    gap[j] <- mean(distance[near])
  }
  A <- rbind(spread = c(-.5, -.5, .5, .5), position = c(-.5, .5, -.5, .5),
             joint = c(1, -1, -1, 1))
  list(L = A %*% W, gap = max(gap))
}

study_fit <- function(d, y, candidate, start = NULL) {
  args <- c(list(design = d, y = y), as.list(candidate))
  if (!is.null(start) && "start" %in% names(formals(ugdm_fit_interactions))) args$start <- start
  do.call(ugdm_fit_interactions, args)
}

study_d <- function(env, idx, base, position, spread, geo = NULL) {
  ugdm_interaction_design(ugdm_design(env, idx, geo, ref = base),
    env, position = position, spread = spread)
}

study_diag <- function(fit, design, y) {
  ugdm_interaction_diagnostics(fit, design, y)
}

study_reference <- function(env, idx, geo = NULL, train = seq_len(nrow(env))) {
  dist_range <- if (is.null(geo)) NULL else range(as.dist(geo[train, train, drop = FALSE]))
  ugdm_design(env, idx, geo, basis_sites = train, dist_range = dist_range)
}

study_rows <- function(form, result) {
  data.frame(form, link_spread = result$link[1L], link_position = result$link[2L],
    link_joint = result$link[3L], response_spread = result$response[1L],
    response_position = result$response[2L], response_joint = result$response[3L])
}

study_policy <- function() {
  settings <- rbind(
    c(0, 0, 0, 0, 0, 0), c(1e-6, 1e-6, 1e-6, 1e-6, 0, 0),
    c(.01, 1e-6, 1e-6, 1e-6, 0, 0), c(1e-6, .01, 1e-6, 1e-6, 0, 0),
    c(1e-6, 1e-6, .01, 1e-6, 0, 0), c(1e-6, 1e-6, .25, 1e-6, 0, 0),
    c(.01, .01, .01, 1e-6, 0, 0), c(.01, .01, .01, 1e-6, .01, .01),
    c(.01, .01, .01, .01, 0, 0))
  colnames(settings) <- c("spread_smooth", "position_smooth", "interaction_shrink", "ridge",
    "spread_null", "position_null")
  as.data.frame(settings)
}

study_new <- function(settings) {
  c(as.list(settings), list(interaction_penalty = "derivative", penalty_scale = "range"))
}

study_link_parts <- function(f, d) {
  layout <- internal(".interaction_layout")(d, f$has_distance)
  eta <- drop(layout$X %*% f$coef)
  intercept <- unname(f$coef[1L])
  position <- drop(layout$X[, layout$pos, drop = FALSE] %*% f$coef[layout$pos])
  list(eta = eta, intercept = intercept, position = position,
       spread = eta - intercept - position)
}

study_shapley <- function(additive, full, baseline_design, full_design) {
  p <- study_link_parts(additive, baseline_design)
  base <- plogis(p$intercept)
  combined <- plogis(p$eta)
  spread <- .5 * (plogis(p$intercept + p$spread) - base +
    combined - plogis(p$intercept + p$position))
  prediction <- ugdm_predict(full, full_design)
  cbind(predicted = prediction, base = base, spread = spread,
    position = combined - base - spread, modification = prediction - combined)
}

study_condition <- function(d, climate = 1L) {
  base <- cbind(1, d$base, d$pos, d$geo)
  scale <- function(X) {
    nrm <- sqrt(colSums(X^2))
    nrm[nrm <= 1e-14] <- 1
    sweep(X, 2L, nrm, "/")
  }
  v <- if (!is.null(d$contextual)) d$contextual else {
    at <- 6L * (climate - 1L) + (1:3)
    d$dir[, at + 3L, drop = FALSE] - d$dir[, at, drop = FALSE]
  }
  other <- cbind(base, d$pair_midpoint)
  residual <- qr.resid(qr(scale(other)), scale(v))
  projected <- sqrt(colSums(residual^2))
  rank <- function(X) {
    val <- svd(scale(X), nu = 0L, nv = 0L)$d
    sum(val > max(val) * 1e-8)
  }
  c(base_rank = rank(base), pair_rank = rank(other), incremental_rank = rank(cbind(other, v)) -
    rank(other), raw_context_SD = median(apply(v, 2L, stats::sd)),
    residual_context = median(projected))
}

study_choose <- function(fold_results, field = "validation_loss") {
  agg <- lapply(split(fold_results, fold_results$candidate), function(a) {
    L <- mean(a[[field]])
    stability <- sqrt(sum(vapply(a[c("spread", "position", "joint")], var, 0.0)))
    c(loss = L, stability = stability)
  })
  agg <- do.call(rbind, agg)
  if (!any(is.finite(agg[, "loss"]))) stop("All regularisation candidates failed.")
  threshold <- min(agg[, "loss"]) * 1.05 + .001
  eligible <- which(agg[, "loss"] <= threshold & is.finite(agg[, "stability"]))
  selected <- eligible[which.min(agg[eligible, "stability"])]
  list(selected = selected, statistics = agg)
}
