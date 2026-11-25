# Rscript docs/manuscript_validation/holdout.R OUT [SPLITS=30]
source("docs/manuscript_validation/common.R")
args <- commandArgs(TRUE)
out <- args[1L]
splits <- if (length(args) > 1L) as.integer(args[2L]) else 30L
meta <- read_manifest()
rows <- list()
for (i in seq_len(nrow(meta))) {
  a <- load_dd2(meta[i, ], Sys.getenv("DD2_DATA"))
  message(meta$id[i], " held-out")
  n <- nrow(a$x)
  P <- (a$x > 0) * 1L
  for (r in seq_len(splits)) {
    set.seed(34000L + r)
    train <- sample.int(n, floor(.7 * n))
    test <- setdiff(seq_len(n), train)
    range <- internal(".pair_distance_range")(a$dist[train, train, drop = FALSE])
    for (k in as.integer(meta[i, c("low", "mid", "high")])) {
      s_test <- draw(P[test, , drop = FALSE], k, 2000L, 0, seed = 800000L + r * 100L + k)
      idx_test <- matrix(test[s_test$idx], ncol = k)
      for (alpha in c(0, .5)) {
        s_train <- draw(P[train, , drop = FALSE], k, 10000L, alpha, seed = 900000L + r * 100L + k)
        idx_train <- matrix(train[s_train$idx], ncol = k)
        # Estimate all predictor transformations from the training sites.
        d <- suppressWarnings(ugdm_design(a$env, idx_train, a$dist, dist_range = range, basis_sites = train))
        dt <- suppressWarnings(ugdm_design(a$env, idx_test, a$dist, ref = d))
        for (response in c("total", "relative")) {
          y <- ugdm_response(a$x, idx_train, response, working_weights = "constant")$y
          yt <- ugdm_response(a$x, idx_test, response, working_weights = "constant")$y
          w <- s_train$sampling_weights
          full <- ugdm_fit(d, y, rep(1, length(y)), sampling_weights = w)
          spread <- internal(".fit_logit_nnls")(-cbind(d$Xdir, d$Xdist), y, w)
          pred <- ugdm_predict(full, dt)
          sp <- plogis(spread$intercept - drop(cbind(dt$Xdir, dt$Xdist) %*% spread$coef))
          rmse <- sqrt(mean((pred - yt)^2))
          rmse_sp <- sqrt(mean((sp - yt)^2))
          rows[[length(rows) + 1L]] <- data.frame(case = meta$id[i], order = k, response,
            sampler = if (alpha == 0) "conventional" else "species", split = r,
            rmse, rmse_spread = rmse_sp, gain = rmse_sp - rmse,
            converged = full$converged && spread$converged)
        }
      }
    }
  }
}
write_study(do.call(rbind, rows), out, "holdout")
