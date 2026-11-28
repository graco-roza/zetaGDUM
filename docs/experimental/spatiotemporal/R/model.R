# Align environment by names when supplied, otherwise by retained panel order.
.st_environment <- function(panel, env) {
  if (!is.array(env) || length(dim(env)) != 3L || !is.numeric(env) ||
      !dim(env)[3] || dim(env)[2] != length(panel$times)) {
    stop("Environment must be a numeric sites x times x predictors array.")
  }
  dn <- dimnames(env)
  target <- dimnames(panel$x)
  if (!is.null(dn[[1]])) {
    if (anyDuplicated(dn[[1]]) || !all(target[[1]] %in% dn[[1]])) stop("Environment site labels do not match.")
    env <- env[target[[1]], , , drop = FALSE]
  }
  if (!is.null(dn[[2]])) {
    if (anyDuplicated(dn[[2]]) || !setequal(target[[2]], dn[[2]])) stop("Environment time labels do not match.")
    env <- env[, target[[2]], , drop = FALSE]
  }
  if (dim(env)[1] != dim(panel$x)[1] || any(!is.finite(env))) {
    stop("Environment must be finite and aligned to the retained panel.")
  }
  vars <- dimnames(env)[[3]]
  if (is.null(vars)) vars <- paste0("env", seq_len(dim(env)[3]))
  if (anyNA(vars) || any(!nzchar(vars)) || anyDuplicated(vars)) stop("Invalid predictor names.")
  dimnames(env) <- list(target[[1]], target[[2]], vars)
  env
}

.st_flat_env <- function(env) {
  matrix(env, ncol = dim(env)[3], dimnames = list(NULL, dimnames(env)[[3]]))
}

st_design <- function(x, env, spatial_order = 1L, transitions = NULL, sets = NULL,
                      dist_matrix = NULL, ref = NULL, limit = 1e5) {
  p <- .st_panel(x)
  e <- .st_environment(p, env)
  n <- dim(e)[1]
  nt <- dim(e)[2]
  k <- .st_one(spatial_order, n, "spatial_order")
  idx <- .st_sets(n, k, sets, limit)
  if (nt < 2L) stop("A transition design needs at least two survey times.")
  if (is.null(transitions)) transitions <- cbind(seq_len(nt - 1L), seq.int(2L, nt))
  if (!is.matrix(transitions) || !is.numeric(transitions) || ncol(transitions) != 2L ||
      !nrow(transitions) || any(!is.finite(transitions) | transitions != floor(transitions) |
                                  transitions < 1 | transitions > nt) ||
      any(transitions[, 1] >= transitions[, 2]) || anyDuplicated(as.data.frame(transitions))) {
    stop("Transitions must be distinct pairs of increasing time indices.")
  }
  if (as.double(nrow(idx)) * nrow(transitions) > limit) stop("Transition budget exceeded.")
  flat <- .st_flat_env(e)
  rows <- expand.grid(set = seq_len(nrow(idx)), transition = seq_len(nrow(transitions)))
  sets_by_row <- idx[rows$set, , drop = FALSE]
  from <- transitions[rows$transition, 1L]
  to <- transitions[rows$transition, 2L]
  initial <- sets_by_row + n * (from - 1L)
  final <- sets_by_row + n * (to - 1L)
  if (!is.null(dist_matrix)) {
    if (k == 1L || !is.matrix(dist_matrix) || !is.numeric(dist_matrix) ||
        !identical(dim(dist_matrix), c(n, n)) ||
        any(!is.finite(dist_matrix) | dist_matrix < 0) ||
        !isTRUE(isSymmetric(unname(dist_matrix))) || any(diag(dist_matrix) != 0)) {
      stop("Distance must be a finite symmetric site distance matrix for k >= 2.")
    }
    if (!is.null(rownames(dist_matrix)) || !is.null(colnames(dist_matrix))) {
      sites <- dimnames(p$x)[[1]]
      if (!identical(rownames(dist_matrix), sites) || !identical(colnames(dist_matrix), sites)) {
        stop("Distance matrix labels must match panel site order.")
      }
    }
  }
  center <- ref$center
  spread <- NULL
  if (k > 1L) {
    transformed <- zetaGDUM::ugdm_design(flat, rbind(initial, final),
                                        ref = ref$native)
    pos <- transformed$Xpos
    spread <- transformed$Xdir
    reference <- list(native = list(environment_basis = transformed$environment_basis),
                      center = transformed$environment_basis$center,
                      scale = transformed$environment_basis$scale)
  } else {
    z <- if (is.null(center)) scale(flat) else scale(flat, center, ref$scale)
    if (!all(is.finite(z))) stop("Predictors require nonzero training variance.")
    pos <- z[c(initial, final), , drop = FALSE]
    reference <- list(center = attr(z, "scaled:center"), scale = attr(z, "scaled:scale"))
  }
  if (!identical(names(reference$center), colnames(flat))) stop("Reference predictors must match.")
  i <- seq_len(nrow(rows))
  j <- i + nrow(rows)
  label <- function(x, prefix) {
    colnames(x) <- paste0(prefix, ".", colnames(x))
    x
  }
  position0 <- pos[i, , drop = FALSE]
  position1 <- pos[j, , drop = FALSE]
  change <- position1 - position0
  terms <- cbind(label(position0, "position_initial"), label(position1, "position_final"),
                 label(change, "position_change"),
                 label(position0 * change, "position_by_change"))
  magnitude <- vapply(seq_len(ncol(flat)), function(v) {
    v0 <- matrix(flat[initial, v], nrow(rows), k)
    v1 <- matrix(flat[final, v], nrow(rows), k)
    rowMeans(abs(v1 - v0)) / reference$scale[v]
  }, numeric(nrow(rows)))
  magnitude <- matrix(magnitude, nrow(rows), ncol(flat), dimnames = list(NULL, colnames(flat)))
  terms <- cbind(terms, label(magnitude, "change_magnitude"))
  if (k > 1L) {
    initial_spread <- spread[i, , drop = FALSE]
    final_spread <- spread[j, , drop = FALSE]
    terms <- cbind(terms, label(initial_spread, "spread_initial"), label(final_spread, "spread_final"),
                   label(final_spread - initial_spread, "spread_change"))
  }
  if (!is.null(dist_matrix)) {
    pairs <- utils::combn(k, 2L)
    distances <- vapply(seq_len(ncol(pairs)), function(j) {
      dist_matrix[cbind(sets_by_row[, pairs[1, j]], sets_by_row[, pairs[2, j]])]
    }, numeric(nrow(rows)))
    terms <- cbind(terms, geographic_distance = rowMeans(matrix(distances, nrow(rows))))
  }
  terms <- cbind(terms, lag = p$times[to] - p$times[from])
  structure(list(X = terms, sets = sets_by_row, from = from, to = to, ref = reference,
                 sites = dimnames(p$x)[[1]], times = p$times), class = "st_design")
}

# Grouped Bernoulli working likelihood by risk set. Small ridge regularizes
# separation in mixed data. Rank deficiency and structural boundaries are explicit.
.st_fit_risk <- function(X, events, opportunities, ridge) {
  keep <- opportunities > 0
  Z <- X[keep, , drop = FALSE]
  y <- events[keep]
  w <- opportunities[keep]
  beta <- setNames(rep(NA_real_, ncol(X)), colnames(X))
  if (!length(y)) return(list(status = "no_opportunities", coefficients = beta, rank = 0L))
  rank <- qr(Z)$rank
  rate <- sum(y) / sum(w)
  if (rate %in% c(0, 1)) {
    return(list(status = "boundary", rate = rate, coefficients = beta, rank = rank))
  }
  if (rank < ncol(X)) return(list(status = "aliased", coefficients = beta, rank = rank))
  penalty <- rep(ridge, ncol(Z))
  objective <- function(b) {
    eta <- drop(Z %*% b)
    sum(w * (pmax(eta, 0) + log1p(exp(-abs(eta)))) - y * eta) +
      sum(penalty * b^2) / 2
  }
  gradient <- function(b) {
    drop(crossprod(Z, w * stats::plogis(Z %*% b) - y)) + penalty * b
  }
  fit <- stats::optim(c(stats::qlogis(rate), rep(0, ncol(X) - 1L)), objective, gradient,
                      method = "BFGS", control = list(maxit = 1000L, reltol = 1e-11))
  beta[] <- fit$par
  mu <- stats::plogis(drop(Z %*% beta))
  information <- crossprod(Z * sqrt(w * mu * (1 - mu)))
  # The independent-trial standard error is diagnostic only.
  working_se <- tryCatch(setNames(sqrt(diag(chol2inv(chol(information)))), colnames(Z)),
                         error = function(e) beta * NA_real_)
  list(status = if (fit$convergence == 0L) "estimated" else "nonconverged",
       coefficients = beta, rank = rank, rate = rate,
       working_independent_se = working_se,
       curvature = eigen(information / sum(w), symmetric = TRUE, only.values = TRUE)$values)
}

.st_predict_risk <- function(fit, X) {
  if (fit$status == "boundary") return(rep(fit$rate, nrow(X)))
  if (fit$status == "estimated") return(stats::plogis(drop(X %*% fit$coefficients)))
  rep(NA_real_, nrow(X))
}

.st_model_matrix <- function(design, terms) {
  if (!is.character(terms) || anyDuplicated(terms) || !all(terms %in% colnames(design$X))) {
    stop("Select unique terms from colnames(design$X).")
  }
  cbind(intercept = 1, design$X[, terms, drop = FALSE])
}

st_gdum <- function(x, design, gain_terms = character(), loss_terms = gain_terms,
                    pool = dimnames(.st_panel(x)$x)[[3]], ridge = 0.1, blocks = NULL) {
  p <- .st_panel(x)
  species <- dimnames(p$x)[[3]]
  if (!inherits(design, "st_design") || !identical(design$sites, dimnames(p$x)[[1]]) ||
      !identical(design$times, p$times)) stop("Design and panel axes must match.")
  if (!is.character(pool) || !length(pool) || anyDuplicated(pool) || !all(pool %in% species)) {
    stop("Use a fixed pool of distinct species labels from the panel.")
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge <= 0) {
    stop("ridge must be a positive scalar.")
  }
  counts <- vapply(seq_len(nrow(design$X)), function(i) {
    h <- matrix(.st_history(p, design$sets[i, ], c(design$from[i], design$to[i])),
                nrow = length(species))[match(pool, species), , drop = FALSE]
    a <- h[, 1]
    b <- h[, 2]
    c(before = sum(a), persistent = sum(a & b), gain = sum(!a & b), loss = sum(a & !b))
  }, numeric(4L))
  counts <- as.data.frame(t(counts))
  counts$gain_opportunities <- length(pool) - counts$before
  counts$loss_opportunities <- counts$before
  X <- list(gain = .st_model_matrix(design, gain_terms), loss = .st_model_matrix(design, loss_terms))
  fit <- lapply(names(X), function(event) {
    .st_fit_risk(X[[event]], counts[[event]], counts[[paste0(event, "_opportunities")]], ridge)
  })
  names(fit) <- names(X)
  if (!is.null(blocks)) {
    if (length(blocks) != length(design$sites) || anyNA(blocks) ||
        length(unique(blocks)) < 3L) stop("Supply three or more site blocks in panel order.")
    memberships <- matrix(blocks[design$sets], nrow(design$sets))
    block_ids <- unique(blocks)
    for (event in names(fit)) {
      deleted <- t(vapply(block_ids, function(block) {
        keep <- rowSums(memberships == block) == 0
        .st_fit_risk(X[[event]][keep, , drop = FALSE], counts[[event]][keep],
                     counts[[paste0(event, "_opportunities")]][keep], ridge)$coefficients
      }, fit[[event]]$coefficients))
      rownames(deleted) <- block_ids
      fit[[event]]$delete_block <- deleted
      fit[[event]]$jackknife_se <- sqrt((nrow(deleted) - 1) / nrow(deleted) *
                                      colSums(sweep(deleted, 2, colMeans(deleted))^2))
    }
  }
  structure(list(fit = fit, counts = counts, design = design,
                 terms = list(gain = gain_terms, loss = loss_terms), pool = pool,
                 ridge = ridge), class = "st_gdum")
}

st_predict <- function(object, design = object$design, before = NULL) {
  if (is.null(before)) {
    if (!identical(design, object$design)) stop("Supply observed baseline counts for new design rows.")
    before <- object$counts$before
  }
  if (!inherits(design, "st_design") || !identical(design$ref, object$design$ref) ||
      ncol(design$sets) != ncol(object$design$sets)) stop("Reuse the fitted design basis and spatial order.")
  if (!is.numeric(before) || length(before) != nrow(design$X) ||
      any(!is.finite(before) | before != floor(before) | before < 0 | before > length(object$pool))) {
    stop("Supply an initial shared-species count per design row within the fitted pool.")
  }
  g <- .st_predict_risk(object$fit$gain, .st_model_matrix(design, object$terms$gain))
  l <- .st_predict_risk(object$fit$loss, .st_model_matrix(design, object$terms$loss))
  gains <- ifelse(before < length(object$pool), (length(object$pool) - before) * g, 0)
  losses <- ifelse(before > 0, before * l, 0)
  data.frame(gain_probability = g, loss_probability = l, persistent = before - losses,
             gains = gains, losses = losses, before = before, after = before + gains - losses,
             change = gains - losses)
}

# Spatial comparison path, using a common basis over the panel. No new
# definition of turnover, nestedness, relative abundance or uniqueness.
st_spatial_gdum <- function(x, env, spatial_order = 2L, component = "turnover",
                            sets = NULL, free_sign = FALSE, limit = 1e5) {
  p <- .st_panel(x)
  e <- .st_environment(p, env)
  n <- dim(e)[1]
  k <- .st_one(spatial_order, n, "spatial_order")
  if (k < 2L) stop("Spatial uniqueness requires at least two sites.")
  idx <- .st_sets(n, k, sets, limit)
  flat <- .st_flat_env(e)
  lapply(seq_along(p$times), function(j) {
    y <- zetaGDUM::ugdm_response(.st_slice(p, j), idx, component)
    design <- zetaGDUM::ugdm_design(flat, idx[y$keep, , drop = FALSE] + n * (j - 1L))
    list(time = p$times[j], response = y, design = design,
         fit = zetaGDUM::ugdm_fit(design, y$y[y$keep], y$trials[y$keep], free_sign = free_sign))
  })
}

# Conditional positive abundance on all occupied rows, or on the same
# persistently occupied site-species pairs at both endpoints.
st_abundance_gdum <- function(x, design, terms = character(), persistent_only = FALSE,
                              response = c("final", "log_ratio"), blocks = NULL) {
  response <- match.arg(response)
  if (response == "log_ratio") persistent_only <- TRUE
  p <- .st_panel(x)
  if (!inherits(design, "st_design") || ncol(design$sets) != 1L ||
      !identical(design$sites, dimnames(p$x)[[1]]) ||
      !identical(design$times, p$times)) stop("Use a matching order-one design.")
  X <- .st_model_matrix(design, terms)
  records <- lapply(seq_len(nrow(X)), function(i) {
    a <- p$x[design$sets[i, 1], design$from[i], ]
    b <- p$x[design$sets[i, 1], design$to[i], ]
    keep <- b > 0 & (!persistent_only | a > 0)
    data.frame(row = rep(i, sum(keep)), initial = a[keep], final = b[keep],
               species = dimnames(p$x)[[3]][keep])
  })
  rows <- do.call(rbind, records)
  if (!nrow(rows)) stop("No eligible positive abundance observations.")
  Z <- X[rows$row, , drop = FALSE]
  if (qr(Z)$rank < ncol(Z)) stop("Abundance terms are aliased.")
  y <- if (response == "final") rows$final else log(rows$final) - log(rows$initial)
  family <- if (response == "final") stats::Gamma(link = "log") else stats::gaussian()
  fit <- stats::glm.fit(Z, y, family = family)
  result <- list(fit = fit, data = rows, design = design, terms = terms,
                 persistent_only = persistent_only, response = response)
  if (!is.null(blocks)) {
    if (length(blocks) != length(design$sites) || anyNA(blocks) ||
        length(unique(blocks)) < 3L) stop("Supply three or more site blocks in panel order.")
    membership <- blocks[design$sets[rows$row, 1]]
    deleted <- t(vapply(unique(blocks), function(block) {
      keep <- membership != block
      A <- Z[keep, , drop = FALSE]
      if (!any(keep) || qr(A)$rank < ncol(Z)) return(rep(NA_real_, ncol(Z)))
      stats::glm.fit(A, y[keep], family = family)$coefficients
    }, numeric(ncol(Z))))
    colnames(deleted) <- colnames(Z)
    rownames(deleted) <- unique(blocks)
    result$delete_block <- deleted
    result$jackknife_se <- sqrt((nrow(deleted) - 1) / nrow(deleted) *
                                 colSums(sweep(deleted, 2, colMeans(deleted))^2))
  }
  result
}
