root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase2_load.R"))
st2_load(root)
for (f in c("factorial", "prediction_validation")) source(file.path(root, "phase2", paste0(f, ".R")))
out <- Sys.getenv("ST2_OUTPUT", file.path(root, "phase2/results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
reps <- as.integer(Sys.getenv("ST2_PRED_REPS", "4"))
draws <- as.integer(Sys.getenv("ST2_DRAWS", "4"))
workers <- as.integer(Sys.getenv("ST2_WORKERS", "3"))
tasks <- expand.grid(replicate = seq_len(reps), scenario = c("null", "interaction", "additive", "opposing"),
                     split = c("sites", "regions", "forward", "joint"), stringsAsFactors = FALSE)
if (Sys.getenv("ST2_PILOT") == "1") tasks <- tasks[tasks$scenario %in% c("null", "interaction") & tasks$split %in% c("joint"), ]
grid <- as.matrix(expand.grid(spatial_order = 1:3, temporal_order = 1:3))
job <- function(i) {
  x <- tasks[i, ]
  sim <- st2_factorial_world(24000L + x$replicate, x$scenario, region = TRUE)
  ans <- st2_validate_split(sim, x$split, grid, draws = draws,
                            interaction = x$scenario == "interaction")
  fields <- c("scores", "quantities", "coherence", "calibration")
  result <- lapply(ans[fields], function(tab) if (nrow(tab)) cbind(replicate = x$replicate, scenario = x$scenario, tab) else NULL)
  result$status <- data.frame(x, model = "M2", event = names(ans$dynamic_status), status = unname(ans$dynamic_status))
  result
}
indices <- seq_len(nrow(tasks))
ans <- if (.Platform$OS.type == "windows" || workers == 1L) lapply(indices, job) else
  parallel::mclapply(indices, job, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
errors <- vapply(ans, inherits, logical(1), "try-error")
if (any(errors)) stop(paste(unlist(ans[errors]), collapse = "\n"))
for (name in names(ans[[1L]])) {
  utils::write.csv(do.call(rbind, lapply(ans, `[[`, name)), file.path(out, paste0("prediction_", name, ".csv")),
                   row.names = FALSE, na = "")
}
