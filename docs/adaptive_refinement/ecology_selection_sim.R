# Called by context_recovery.R in selection mode.
candidate <- study_policy()
study_selection <- function(r) {
  q <- DGP(r)
  k <- 6L
  set.seed(3314L + 790L * r)
  split <- if (scenario == "spatial") q$coords else NULL
  partitions <- ugdm_site_partitions(n, 4L, 3L, split)
  eval_idx <- study_idx(seq_len(n), k, 1500L)
  eval_truth <- q$truth(eval_idx)
  reference <- study_d(q$env, eval_idx, q$ref, "nonlinear", "varying")
  C <- study_contrast(reference, neighbours = 20L)$L
  target <- drop(C %*% plogis(eval_truth$eta))
  summaries <- list()
  for (partition in seq_len(ncol(partitions))) for (hold in 1:4) {
    train <- which(partitions[, partition] != hold)
    test <- which(partitions[, partition] == hold)
    train_idx <- study_idx(train, k, 340L)
    test_idx <- study_idx(test, k, 340L)
    ref <- study_reference(q$env, train_idx, train = train)
    d <- study_d(q$env, train_idx, ref, "nonlinear", "varying")
    validation <- study_d(q$env, test_idx, ref, "nonlinear", "varying")
    common <- study_d(q$env, eval_idx, ref, "nonlinear", "varying")
    fit_response <- q$truth(train_idx)$y
    trial_response <- q$truth(test_idx)$y
    start <- NULL
    for (i in seq_len(nrow(candidate))) {
      f <- study_fit(d, fit_response, study_new(candidate[i, ]), start)
      if (f$converged) start <- f$coef
      mu <- ugdm_predict(f, validation)
      estimate <- ugdm_interaction_contrast(f, common, C)$response
      summaries[[length(summaries) + 1L]] <- data.frame(partition, hold, candidate = i,
        loss = if (f$converged) loss(trial_response, mu, rep(1 / length(mu), length(mu))) else Inf,
        spread = estimate[1L], position = estimate[2L], joint = estimate[3L])
    }
  }
  fold_runs <- do.call(rbind, summaries)
  selection <- study_choose(fold_runs, field = "loss")
  validation_paths <- function(partition) {
    if (!any(is.finite(fold_runs$loss[fold_runs$partition == partition]))) return(NA_integer_)
    study_choose(fold_runs[fold_runs$partition == partition, ], field = "loss")$selected
  }
  partition_choice <- vapply(seq_len(ncol(partitions)), validation_paths, integer(1L))
  evaluate <- study_idx(seq_len(n), k, 2100L)
  y <- q$truth(evaluate)$y
  d <- study_d(q$env, evaluate, q$ref, "nonlinear", "varying")
  per_fit <- lapply(seq_len(nrow(candidate)), function(i) {
    f <- study_fit(d, y, study_new(candidate[i, ]))
    e <- ugdm_interaction_contrast(f, reference, C)$response
    data.frame(candidate = i, spread_error = e[1L] - target[1L],
      position_error = e[2L] - target[2L], joint_error = e[3L] - target[3L],
      converged = f$converged)
  })
  fits <- do.call(rbind, per_fit)
  metric <- rowSums(fits[c("spread_error", "position_error", "joint_error")]^2)
  metric[!fits$converged] <- Inf
  cbind(data.frame(scenario, r, rule = selection$selected,
    loss_only = which.min(selection$statistics[, "loss"]), oracle = which.min(metric),
    number_agreeing = sum(partition_choice == selection$selected, na.rm = TRUE)),
    fits, selection$statistics)
}
rows <- lapply(seq_len(B), function(i) {
  message("selection ", scenario, " ", i, "/", B)
  study_selection(i)
})
write_study(do.call(rbind, rows), "docs/adaptive_refinement/results",
  paste0("selection_", scenario))
