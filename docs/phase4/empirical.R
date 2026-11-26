source("docs/phase4/common.R")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3L)
meta <- read_manifest()[as.integer(args[2L]), ]
component <- match.arg(args[3L], c("total", "relative"))
data <- load_dd2(meta, args[1L])
# The first metadata predictor is locked before any response/model inspection.
env <- data$env[1]
grid <- sort(unique(c(meta$low, meta$mid, meta$high, ceiling(meta$sites * .20), meta$sites)))
message(meta$community, " ", component, " predictor ", names(env), " orders ", paste(grid, collapse = ","))
fit <- phase4_profile(data$x, env, grid, component, reference_M = 12000, seed = 413 + meta$sites)
key <- paste(meta$id, component, sep = "_")
stamp <- function(x) cbind(community = meta$community, n = meta$sites,
                            predictor = names(env), component = component, x)
phase4_write(stamp(cbind(fit$profile, retention_break = fit$retention_break,
                          highest_supported_order = fit$highest_supported_order)), "profile", key)
phase4_write(stamp(fit$precision), "precision", key)
phase4_write(stamp(fit$pilots), "pilot", key)
boot <- paired_bootstrap(data$x, env, grid[grid < meta$sites], component, 4000, 39, 195 + meta$sites)
phase4_write(stamp(boot), "bootstrap", key)
