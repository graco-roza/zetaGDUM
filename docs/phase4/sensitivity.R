source("docs/phase4/common.R")
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
a <- read_manifest()[as.integer(args[2L]), ]
d <- load_dd2(a, args[1L])
env <- d$env[1L]
grid <- sort(unique(c(3L, unlist(a[c("low", "mid", "high")], use.names = FALSE))))
message(a$community, " reference sensitivity: ", paste(grid, collapse = ","))
ans <- phase4_profile(d$x, env, grid, "relative", reference_M = 50000,
                      seed = 403 + a$sites)
stamp <- function(x) cbind(community = a$community, precursor = a$low, n = a$sites,
                            predictor = names(env), component = "relative", reference_M = 50000, x)
phase4_write(stamp(cbind(ans$profile, highest_supported_order = ans$highest_supported_order)),
             "sensitivity", a$id)
phase4_write(stamp(ans$precision), "sensitivity_precision", a$id)
