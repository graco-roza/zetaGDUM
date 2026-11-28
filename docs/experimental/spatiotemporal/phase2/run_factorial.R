root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase2_load.R"))
st2_load(root)
for (f in c("factorial", "prediction_validation", "effect_validation")) {
  source(file.path(root, "phase2", paste0(f, ".R")))
}
out <- Sys.getenv("ST2_OUTPUT", file.path(root, "phase2/results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
reps <- as.integer(Sys.getenv("ST2_FACTOR_REPS", "4"))
draws <- as.integer(Sys.getenv("ST2_DRAWS", "4"))
workers <- as.integer(Sys.getenv("ST2_WORKERS", "3"))
primary <- expand.grid(scenario = c("gain", "loss", "abundance"), target = 1:8, stringsAsFactors = FALSE)
extras <- data.frame(scenario = c("null", "additive", "multiple", "interaction", "nonlinear", "opposing",
                                 "unequal", "correlated", "trend", "shock"), target = 0L)
scenarios <- rbind(primary, extras)
if (Sys.getenv("ST2_PILOT") == "1") scenarios <- scenarios[c(1L, 17L, 25L, 28L), ]
args <- merge(scenarios, data.frame(replicate = seq_len(reps)))
job <- function(i) {
  row <- args[i, ]
  s <- st2_factorial_world(15000 + row$replicate, row$scenario, row$target, region = row$scenario == "shock")
  result <- st2_factorial_effect(s, row$replicate, draws = draws)
  cbind(target = row$target, result)
}
indices <- seq_len(nrow(args))
result <- if (.Platform$OS.type == "windows" || workers == 1L) lapply(indices, job) else
  parallel::mclapply(indices, job, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
errors <- vapply(result, inherits, logical(1), "try-error")
if (any(errors)) stop(paste(unlist(result[errors]), collapse = "\n"))
z <- do.call(rbind, result)
z$error <- z$estimate - z$truth
utils::write.csv(z, file.path(out, "factorial_contrasts.csv"), row.names = FALSE, na = "")
