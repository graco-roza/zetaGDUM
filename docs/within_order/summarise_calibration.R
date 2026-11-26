source("docs/within_order/simulation_design.R")
paths <- list.files(within_dir, "^calibration_.*[.]csv[.]gz$", full.names = TRUE)
if (!length(paths)) stop("Run calibration.R before summarising.")
retry_paths <- list.files(within_dir, "^score_setup_retry_.*[.]csv[.]gz$", full.names = TRUE)
data <- do.call(rbind, lapply(c(paths, retry_paths), function(path) {
  a <- read.csv(gzfile(path))
  # Keep the replicate and draw budgets visible.
  name <- basename(path)
  a$optimizer_retry <- startsWith(name, "score_setup_retry_")
  name <- sub("^score_setup_retry_", "calibration_", name)
  key <- strsplit(sub("[.]csv[.]gz$", "", name), "_", fixed = TRUE)[[1L]]
  for (field in c("score_pool_seconds", "covariance_condition", "information_condition"))
    if (!field %in% names(a)) a[[field]] <- NA_real_
  a$B <- as.integer(key[6L])
  a$draw_budget <- as.integer(key[7L])
  a$set_budget <- as.integer(key[8L])
  a
}))
keys <- c("scenario", "n", "order", "response", "B", "draw_budget", "set_budget", "method", "estimand")
data <- data[!duplicated(data[c(keys, "replicate")], fromLast = TRUE), ]
summary <- lapply(split(seq_len(nrow(data)), interaction(data[keys], drop = TRUE)), function(i) {
  q <- data[i, ]
  usable <- is.finite(q$lo95)
  a <- q[usable, ]
  val <- function(z) if (any(is.finite(z))) mean(z[is.finite(z)]) else NA_real_
  cbind(q[1L, keys], total = nrow(q), usable = nrow(a),
    optimizer_retries = sum(q$optimizer_retry), truth = val(q$truth), mean_draws = val(q$draws_ok),
    convergence = val(q$point_fit), response_sd = val(q$response_sd),
    time_seconds = val(q$seconds), point_seconds = val(q$point_elapsed),
    score_pool_seconds = val(q$score_pool_seconds),
    covariance_condition = val(log10(q$covariance_condition)),
    information_condition = val(log10(q$information_condition)),
    coverage80 = val(a$lo80 <= a$truth & a$truth <= a$hi80),
    coverage95 = val(a$lo95 <= a$truth & a$truth <= a$hi95),
    width80 = val(a$hi80 - a$lo80), width95 = val(a$hi95 - a$lo95),
    point_bias = val(q$point - q$truth),
    point_rmse = sqrt(val((q$point - q$truth)^2)),
    bias = val(a$median - a$truth), rmse = sqrt(val((a$median - a$truth)^2)),
    practical_positive = val(a$p_practical >= .95), positive = val(a$p_greater >= .95),
    practical_negative = val(a$p_reverse_practical >= .95), negative = val(a$p_greater <= .05),
    mean_leading = val(a$p_leading))
})
within_write(do.call(rbind, summary), "calibration_summary.csv")
leaders <- data[!is.na(data$p_leading), ]
if (nrow(leaders)) {
  split_runs <- split(leaders, interaction(leaders[c(keys[1:8], "replicate")], drop = TRUE))
  ranks <- do.call(rbind, lapply(split_runs, function(a) {
    tol <- 1e-10
    truth_leads <- a$truth >= max(a$truth) - tol
    observed <- a$p_leading >= max(a$p_leading) - tol
    data.frame(a[1, keys[1:8]], replicate = a$replicate[1], truth_tied = sum(truth_leads) > 1L,
      success = mean(truth_leads[observed]), declared_probability = max(a$p_leading),
      assigned_mass = sum(a$p_leading[truth_leads]))
  }))
  within_write(ranks, "calibration_rank.csv")
  ranks$probability_bin <- cut(ranks$declared_probability,
    breaks = c(0, .4, .5, .6, .7, .8, .9, 1), include.lowest = TRUE)
  rank_keys <- c(keys[1:8], "probability_bin")
  reliability <- lapply(split(ranks, ranks[rank_keys], drop = TRUE), function(a) {
    data.frame(a[1L, rank_keys], communities = nrow(a),
      mean_probability = mean(a$declared_probability),
      correct_leader = mean(a$success), mean_mass_on_reference_leader = mean(a$assigned_mass))
  })
  within_write(do.call(rbind, reliability), "calibration_rank_reliability.csv")
}
