# Run from repo root. Use the same validated DD2 inputs as the prior analyses.
source("docs/adaptive_refinement/ecology_common.R")
arg <- commandArgs(TRUE)
if (length(arg) < 2L) stop("Supply the DD2 directory and a row of the community manifest.")
meta <- read_manifest()[as.integer(arg[2L]), ]
dat <- load_dd2(meta, arg[1L])
x <- dat$x
env <- dat$env
geo <- dat$dist
n <- nrow(x)
order_pair <- unlist(meta[c("mid", "high")], use.names = FALSE)
M <- if (length(arg) >= 3L) as.integer(arg[3L]) else 9000L
M_fold <- if (length(arg) >= 4L) as.integer(arg[4L]) else 500L
folds <- study_blocks(geo)
candidate <- study_policy()
methods <- c("occurrence", "relative")
get_y <- function(idx, kind) ugdm_response(x, idx, if (kind == "occurrence") "total" else "relative",
                                          working_weights = "constant")$y
set.seed(203112L + n)
paths <- function(m) lapply(sample_sets_coupled(x, order_pair, m,
  sampler = "uniform", exact_limit = 0L), "[[", "idx")
reference_idx <- paths(max(14000L, M))
fit_idx <- paths(M)
ref <- study_reference(env, fit_idx[[1L]], geo)
fold_data <- vector("list", length(folds) * 4L)
z <- 0L
for (partition in folds) for (hold in sort(unique(partition))) {
  train <- which(partition != hold)
  test <- which(partition == hold)
  k <- order_pair[1L]
  tr <- study_idx(train, k, M_fold)
  te <- study_idx(test, k, M_fold)
  train_ref <- study_reference(env, tr, geo, train)
  z <- z + 1L
  fold_data[[z]] <- list(train_idx = tr, test_idx = te,
    train = study_d(env, tr, train_ref, "nonlinear", "varying", geo),
    test = study_d(env, te, train_ref, "nonlinear", "varying", geo),
    common = study_d(env, reference_idx[[1L]][seq_len(min(2000L, nrow(reference_idx[[1L]]))), ],
      train_ref, "nonlinear", "varying", geo))
}
centred <- study_d(env, reference_idx[[1L]][seq_len(min(2000L, nrow(reference_idx[[1L]]))), ],
  ref, "nonlinear", "varying", geo)
contrast <- study_contrast(centred, neighbours = 20L)
site_output <- model_output <- cv_output <- list()
for (kind in methods) {
  message(meta$community, " ", kind, " smoothing")
  rows <- vector("list", length(fold_data) * nrow(candidate))
  for (i in seq_along(fold_data)) {
    a <- fold_data[[i]]
    y_train <- get_y(a$train_idx, kind)
    y_test <- get_y(a$test_idx, kind)
    start <- NULL
    for (j in seq_len(nrow(candidate))) {
      f <- study_fit(a$train, y_train, study_new(candidate[j, ]), start)
      if (f$converged) start <- f$coef
      mu <- ugdm_predict(f, a$test)
      l <- if (f$converged) loss(y_test, mu, rep(1 / length(mu), length(mu))) else Inf
      c <- ugdm_interaction_contrast(f, a$common, contrast$L)
      rows[[(i - 1L) * nrow(candidate) + j]] <- data.frame(partition = ceiling(i / 4L),
        fold = (i - 1L) %% 4L + 1L, candidate = j, validation_loss = l,
        converged = f$converged, spread = c$response[1L], position = c$response[2L],
        joint = c$response[3L])
    }
  }
  cv <- do.call(rbind, rows)
  # Repeated partitions overlap. This loss allowance is descriptive and has
  # no independent-block standard-error interpretation.
  selected <- study_choose(cv)
  chosen <- selected$selected
  cv_output[[kind]] <- cbind(meta[c("id", "community")][rep(1L, nrow(cv)), ],
    response = kind, selected = chosen, cv, candidate[cv$candidate, ])
  retention <- zeta_local_retention(x, if (kind == "occurrence") "occurrence" else "relative")
  path <- list()
  sites <- list()
  for (r in seq_along(order_pair)) {
    k <- order_pair[r]
    message(meta$community, " ", kind, " fits at order ", k)
    ix <- fit_idx[[r]]
    ev <- reference_idx[[r]]
    y <- get_y(ix, kind)
    y_ref <- get_y(ev, kind)
    sharing <- mean(internal("shared_species_cpp")(1L * (x > 0), ev) > 0)
    form <- data.frame(model = c("A", "B", "C", "C_adaptive", "D"),
      position = c("linear", rep("nonlinear", 4L)),
      spread = c("constant", "constant", "varying", "varying", "contextual"))
    evaluations <- function(z) study_d(env, ev, ref, form$position[z], form$spread[z], geo)
    additive_ref <- evaluations(2L)
    scores <- study_contrast(additive_ref)
    support <- retention$summary
    observed <- if (kind == "relative") 1 - retention$profile[, k] else
      internal(".conditional_sample_mean")(ev, cbind(y_ref), rep(1, length(y_ref)), n)[, 1L]
    results <- vector("list", nrow(form))
    for (j in seq_len(nrow(form))) {
      f <- form[j, ]
      d <- study_d(env, ix, ref, f$position, f$spread, geo)
      setup <- if (j < 4L) list() else study_new(candidate[chosen, ])
      fit <- study_fit(d, y, setup)
      diagnostic <- study_diag(fit, d, y)
      e <- evaluations(j)
      result <- ugdm_interaction_contrast(fit, e, scores$L)
      study <- study_rows(f$model, result)
      stats <- if (j >= 3L) study_condition(d) else rep(NA_real_, 5L)
      out <- data.frame(id = meta$id, community = meta$community, response = kind, order = k,
        set_count = nrow(ix), assessment_count = nrow(ev), selected = chosen,
        converged = fit$converged, gradient = fit$projected_gradient, D2 = fit$D2,
        response_SD = sd(y_ref), information_rank = diagnostic$rank,
        information_parameters = diagnostic$parameters, condition = diagnostic$condition,
        active = length(diagnostic$active),
        spread_allocation = unname(diagnostic$deviance_allocation["spread"]),
        position_allocation = unname(diagnostic$deviance_allocation["position"]),
        modification_allocation = unname(diagnostic$deviance_allocation["interaction"]),
        support_gap = scores$gap,
        sharing = sharing,
        study, base_rank = stats[1L], pair_rank = stats[2L], contextual_rank = stats[3L],
        contextual_SD = stats[4L], residual_context = stats[5L], row.names = NULL)
      results[[j]] <- out
      if (j == 2L) additive <- fit
      if (j %in% c(2L, 4L, 5L)) {
        values <- study_shapley(additive, fit, additive_ref, e)
        local <- internal(".conditional_sample_mean")(ev, values, rep(1, nrow(ev)), n)
        colnames(local) <- colnames(values)
        sites[[length(sites) + 1L]] <- data.frame(id = meta$id, community = meta$community,
          response = kind, order = k, model = f$model, site = seq_len(n),
          characteristic_retention = support$T, conditional_retention = support$T_plus,
          retention = retention$profile[, k], observed = observed, local,
          residual = observed - local[, "predicted"], row.names = NULL)
      }
    }
    path[[r]] <- do.call(rbind, results)
  }
  site_output[[kind]] <- do.call(rbind, sites)
  model_output[[kind]] <- do.call(rbind, path)
}
root <- "docs/adaptive_refinement/results"
write_study(do.call(rbind, cv_output), root, paste0("ecology_cv_", meta$id))
write_study(do.call(rbind, model_output), root, paste0("ecology_models_", meta$id))
write_study(do.call(rbind, site_output), root, paste0("ecology_sites_", meta$id))
