source("docs/within_order/simulation_design.R")
args <- commandArgs(TRUE)
n <- as.integer(args[1L])
k <- as.integer(args[2L])
response <- match.arg(args[3L], c("incidence", "relative"))
scenario <- args[4L]
generating <- if (scenario == "spatial-naive") "spatial" else scenario
B <- as.integer(args[5L])
draws <- as.integer(args[6L])
M <- as.integer(args[7L])
methods <- if (length(args) > 7L) strsplit(args[8L], ",", fixed = TRUE)[[1L]] else
  c("bayesian_bootstrap", "gaussian_gmm", "multiplier")
reference_path <- file.path(within_dir, paste0("reference_", generating, "_", k, ".csv"))
stopifnot(file.exists(reference_path))
reference <- read.csv(reference_path)
truth <- tapply(reference$importance[reference$response == response],
                reference$predictor[reference$response == response], mean)
if (scenario == "equal") truth[c("strong", "weak")] <- mean(truth[c("strong", "weak")])
fixed <- within_reference()
set.seed(81444L + k)
validation <- within_data(within_sites(2600L, generating), k, 4500L, response, fixed)
threshold <- .01
run_replicate <- function(i) {
  set.seed(14018L + 1871L * i + 45L * n + k)
  S <- within_sites(n, generating)
  a <- within_data(S, k, M, response, fixed)
  f <- within_fit(a)
  score_started <- proc.time()[["elapsed"]]
  score_pools <- lapply(seq_len(2L), function(j) {
    x <- within_data(S, k, M, response, fixed)
    x[c("design", "y", "trials")]
  })
  score_seconds <- proc.time()[["elapsed"]] - score_started
  build <- function(idx, sites) {
    x <- within_data(list(env = S$env[sites, , drop = FALSE],
                         counts = S$counts[sites, , drop = FALSE]),
                     k, M, response, fixed, idx)
    x[c("design", "y", "trials")]
  }
  point_start <- proc.time()[["elapsed"]]
  point <- ugdm_within_order_importance(f, a$design, a$y, a$trials,
    evaluation = validation$design, evaluation_y = validation$y,
    evaluation_trials = validation$trials, failure = "warn")$predictors
  point_elapsed <- proc.time()[["elapsed"]] - point_start
  result <- lapply(methods, function(method) {
    elapsed <- proc.time()[["elapsed"]]
    res <- tryCatch(ugdm_within_order_resample(f, a$design, a$y, a$trials,
      evaluation = validation$design, evaluation_y = validation$y,
      evaluation_trials = validation$trials, method = method, M = M,
      draws = draws, build = build, threshold = threshold,
      score_pools = if (method %in% c("gaussian_gmm", "multiplier")) score_pools else NULL,
      blocks = if (scenario == "spatial") S$blocks else NULL,
      seed = 8271 + i), error = identity)
    sec <- proc.time()[["elapsed"]] - elapsed
    E <- diag(3L)
    index <- utils::combn(3L, 2L)
    # Rebuild row contrasts in the canonical predictor order.
    L <- rbind(E, E[index[1L, ], , drop = FALSE] - E[index[2L, ], , drop = FALSE])
    vars <- c("strong","weak","irrelevant")
    labels <- c(vars, paste(vars[index[1L, ]], vars[index[2L, ]], sep = " - "))
    expected <- drop(L %*% truth[vars])
    a_hat <- drop(L %*% point$absolute[vars])
    base <- data.frame(n, order = k, response, scenario, replicate = i, method,
      estimand = labels, truth = expected, point = a_hat, point_elapsed, seconds = sec, score_pool_seconds = if (method %in% c("gaussian_gmm", "multiplier")) score_seconds else 0,
      response_sd = sd(a$y), point_fit = f$converged, full = point$total)
    base$covariance_condition <- base$information_condition <- NA_real_
    if (inherits(res, "error")) {
      D <- matrix(NA_real_, 0L, 3L)
      error <- conditionMessage(res)
    } else {
      D <- as.matrix(res$predictor_draws)[, vars, drop = FALSE]
      if (!is.null(res$score_diagnostics)) {
        base$covariance_condition <- max(res$score_diagnostics$covariance_condition)
        base$information_condition <- max(res$score_diagnostics$information_condition)
      }
      error <- ""
    }
    valid <- rowSums(is.finite(D)) == 3L
    D <- D[valid, , drop = FALSE]
    base$draws_ok <- nrow(D)
    base$error <- error
    base$median <- base$lo80 <- base$hi80 <- base$lo95 <- base$hi95 <-
      base$p_greater <- base$p_practical <- base$p_reverse_practical <- base$p_leading <- NA_real_
    if (nrow(D) >= ceiling(.80 * draws)) {
      Q <- D %*% t(L)
      q <- t(apply(Q, 2L, quantile, c(.5,.1,.9,.025,.975), names = FALSE))
      base[, c("median","lo80","hi80","lo95","hi95")] <- q
      base$p_greater <- colMeans(Q > 0)
      base$p_practical <- colMeans(Q > threshold)
      base$p_reverse_practical <- colMeans(Q < -threshold)
      leaders <- D == apply(D, 1L, max)
      base$p_leading[seq_len(3L)] <- colMeans(leaders / rowSums(leaders))
    }
    base
  })
  if (i %% 10L == 0L) message("Replicate ", i, "/", B)
  do.call(rbind, result)
}
key <- paste("calibration", scenario, n, k, response, B, draws, M, sep = "_")
replicates <- seq_len(B)
if (identical(args[9L], "--retry-failed")) {
  original <- read.csv(gzfile(file.path(within_dir, paste0(key, ".csv.gz"))))
  failed <- grepl("^Coalition .* failed", original$error)
  replicates <- sort(unique(original$replicate[failed & original$method %in% methods]))
  if (!length(replicates)) stop("There are no recorded setup failures to retry.")
  key <- sub("^calibration", "score_setup_retry", key)
}
workers <- as.integer(Sys.getenv("MS_GDUM_CALIBRATION_WORKERS", "1"))
stopifnot(!is.na(workers), workers >= 1L)
out <- if (workers > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(replicates, run_replicate, mc.cores = workers, mc.set.seed = FALSE)
} else {
  lapply(replicates, run_replicate)
}
if (!all(vapply(out, is.data.frame, logical(1L)))) {
  stop("An ecological replicate failed outside the recorded fitting stages.")
}
within_write(do.call(rbind, out), paste0(key,".csv.gz"))
