source("docs/within_order/simulation_design.R")
args <- commandArgs(TRUE)
scenario <- args[1L]
k <- as.integer(args[2L])
set.seed(18041L + 40L * k + match(scenario, c("baseline","equal","correlated","position","sparse","spatial","broad")))
a <- lapply(seq_len(3L), function(j) {
  sites <- within_sites(3600L, scenario)
  idx <- sample_sets(nrow(sites$env), k, 85000L)
  lapply(c("incidence","relative"), function(response) within_data(sites, k, 85000L, response, idx = idx))
})
result <- do.call(rbind, lapply(seq_len(3L), function(j) {
  do.call(rbind, lapply(seq_len(2L), function(i) {
    x <- a[[j]][[i]]
    b <- a[[j %% 3L + 1L]][[i]]
    f <- within_fit(x)
    out <- ugdm_within_order_importance(f, x$design, x$y, x$trials,
      evaluation = b$design, evaluation_y = b$y, evaluation_trials = b$trials)
    data.frame(scenario, order = k, reference = j, response = c("incidence","relative")[i],
      predictor = names(out$predictors$absolute), importance = out$predictors$absolute,
      full = out$predictors$total, converged = f$converged)
  }))
}))
within_write(result, paste0("reference_", scenario, "_", k, ".csv"))
