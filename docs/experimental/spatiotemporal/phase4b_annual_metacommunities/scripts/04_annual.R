source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")
source(file.path(directory,"scripts","trajectory.R"))
load_kernel()
x <- readRDS(file.path(processed,"biotime-604.rds"))
r <- analyse_trajectory(x$panel,x$stratum,draws=64L,keep_sets=TRUE)
for (name in names(r)[vapply(r,is.data.frame,logical(1))]) write_result(r[[name]],name)
write_result(data.table(test=c("annual telescope","annual Shapley efficiency"),
                        error=c(r$telescope_error,r$efficiency_error)), "empirical_validation")
# Structured block deletions of whole years, evaluated on fixed regional classifications.
u <- r$annual[response %in% responses[1:4]]
u[, timeblock := (initial-min(initial)) %/% 3L]
u <- u[, .(net=mean(net)),by=.(k,response,process,timeblock)]
u <- u[, .(temporal_block_se=jackknife(vapply(timeblock,
  function(t) mean_valid(net[timeblock != t]), numeric(1))), blocks=.N),
  by=.(k,response,process)]
write_result(u,"temporal_block_sensitivity")
