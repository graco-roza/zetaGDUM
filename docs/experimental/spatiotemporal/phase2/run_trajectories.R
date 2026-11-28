root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase2_load.R"))
st2_load(root)
source(file.path(root, "phase2/factorial.R"))
out <- Sys.getenv("ST2_OUTPUT", file.path(root, "phase2/results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cases <- expand.grid(replicate = 1:4, split = c("sites", "regions", "forward", "joint"))
work <- function(idx) {
  case <- cases[idx, ]
  sim <- st2_factorial_world(60500L + case$replicate, "additive", region = TRUE)
  n <- dim(sim$enabled$x)[1L]
  test_sites <- switch(as.character(case$split), sites = sort(sample.int(n, n %/% 3L)),
                       regions = which(sim$blocks >= 5L), forward = seq_len(n),
                       joint = which(sim$blocks >= 5L))
  train_sites <- if (case$split == "forward") seq_len(n) else setdiff(seq_len(n), test_sites)
  train_times <- if (case$split %in% c("forward", "joint")) 1:8 else 1:12
  test_times <- if (case$split %in% c("forward", "joint")) 8:12 else 1:12
  a <- st2_subset_panel(sim$enabled, train_sites, train_times)
  b <- st2_subset_panel(sim$enabled, test_sites, test_times)
  e0 <- sim$env[, train_times, , drop = FALSE]
  e1 <- sim$env[, test_times, , drop = FALSE]
  train <- st2_transition_design(a, e0, 1)
  m2 <- st2_fit_dynamic(a, train)
  positive <- st2_fit_positive(a, train, hierarchical = TRUE)
  measure <- function(p) st2_zeta(p, 3, 3, mode = "consecutive", empty = "undefined")$surface
  true <- measure(b)
  draws <- vapply(1:8, function(draw) measure(st2_predict_trajectory(m2, positive, b, e1))$zeta,
                  numeric(nrow(true)))
  data.frame(replicate = case$replicate, split = case$split, true[c("spatial_order", "temporal_order", "response")],
             observed = true$zeta, predicted = rowMeans(draws),
             prediction_mc_sd = apply(draws, 1L, stats::sd), row.names = NULL)
}
jobs <- seq_len(nrow(cases))
workers <- as.integer(Sys.getenv("ST2_WORKERS", "3"))
ans <- if (.Platform$OS.type == "windows" || workers == 1L) lapply(jobs, work) else
  parallel::mclapply(jobs, work, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
if (any(vapply(ans, inherits, logical(1), "try-error"))) stop("A trajectory worker failed.")
z <- do.call(rbind, ans)
z$error <- z$predicted - z$observed
utils::write.csv(z, file.path(out, "trajectory_predictions.csv"), row.names = FALSE, na = "")
