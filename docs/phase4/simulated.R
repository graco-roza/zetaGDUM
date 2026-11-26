source("docs/phase4/common.R")
args <- commandArgs(TRUE)
case <- if (length(args)) as.integer(args[1L]) else 1L
spec <- data.frame(n = c(70,100,150,250,72,18), species = c(40,100,250,350,19,35),
  breadth = c(.32,.6,.9,1.3,NA,.85), abundance = c(1,1.4,2,1.8,NA,2),
  gradient = c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE), spatial = c(.2,.5,.8,0,NA,.2))
stopifnot(case %in% seq_len(nrow(spec)))
a <- spec[case, ]
s <- if (case == 5L) support_case(a$n) else simulate_metacommunity(a$n, a$species,
  breadth = a$breadth, abundance_mean = a$abundance, spatial = a$spatial,
  use_env = c(a$gradient, FALSE), seed = 8810 + case)
s$comm <- safe_community(s$comm)
richness <- rowSums(s$comm > 0)
occupancy <- colMeans(s$comm > 0)
phase4_write(cbind(case, a, species_retained = sum(occupancy > 0),
  mean_richness = mean(richness), richness_sd = sd(richness),
  mean_occupancy = mean(occupancy[occupancy > 0]),
  rare_species_fraction = mean(occupancy[occupancy > 0] <= .1),
  median_positive_abundance = median(s$comm[s$comm > 0])), "composition", paste0("sim", case))
if (length(args) > 1L && args[2L] == "summary") quit(save = "no")
env <- s$env["env1"]
mid <- max(3L, round(a$n * .05))
high <- max(5L, round(a$n * .12))
orders <- sort(unique(c(2, mid, high, if (case == 5L) 5L, a$n)))
for (component in c("total", "relative")) {
  fit <- phase4_profile(s$comm, env, orders, component, reference_M = 9000,
                       seed = 489 + case)
  key <- paste0("sim", case, "_", component)
  stamp <- function(x) cbind(case, component, a[rep(1L, nrow(x)), ], x, row.names = NULL)
  phase4_write(stamp(cbind(fit$profile, retention_break = fit$retention_break,
                            highest_supported_order = fit$highest_supported_order)), "profile", key)
  phase4_write(stamp(fit$precision), "precision", key)
  phase4_write(stamp(fit$pilots), "pilot", key)
}
