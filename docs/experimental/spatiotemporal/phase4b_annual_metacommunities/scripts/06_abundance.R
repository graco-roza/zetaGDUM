source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")
source(file.path(directory,"scripts","trajectory.R"))
Rcpp::sourceCpp(file.path(phase,"phase4_swedish_lakes","scripts","rectangle_kernel.cpp"),
                cacheDir=file.path(processed,"cpp-cache"))
x <- readRDS(file.path(processed,"biotime-604.rds"))
p <- st_panel(x$counts,x$panel$times)
n <- dim(p$x)[1]
years <- p$times
metrics <- c(1L,2L,4L,7L,8L,9L,12L,13L)
metric_names <- c("zeta_presence","zeta_abundance","incidence","abundance","occupancy",
                  "numerical_heterogeneity","local_sorensen","local_relative_bray_curtis")
out <- list()
for (t in seq_len(length(years)-1L)) {
  plan <- st3_plan(p,t,t+1L)
  alloc <- st3_allocations(plan,draws=32,seed=5109+t)
  for (k in 1:6) {
    sets <- t(combn(n,k))
    z <- st4_evaluate_cpp(plan$before,plan$after,lapply(seq_len(nrow(sets)),function(j) sets[j,]),
                          alloc$masks,.st3_coalitions())
    use <- apply(is.finite(z$phi[,7:9,,,drop=FALSE]),1,all)
    good <- function(v) max(c(0,abs(v)),na.rm=TRUE)
    stopifnot(z$efficiency_error < 1e-11,
              good(z$phi[use,7,,]-z$phi[use,8,,]-z$phi[use,9,,])<1e-11,
              good(z$phi[,4,4,])<1e-11)
    for (r in seq_along(metrics)) {
      j <- metrics[r]
      valid <- if (j %in% 7:9) use else apply(is.finite(z$phi[,j,,,drop=FALSE]),1,all)
      for (i in 1:4) {
        v <- apply(z$phi[valid,j,i,,drop=FALSE],1,mean)
        draw <- apply(z$phi[valid,j,i,,drop=FALSE],4,mean_valid)
        out[[length(out)+1L]] <- data.table(k=k,initial=years[t],final=years[t+1L],
          response=metric_names[r],process=st3_processes[i],before=mean_valid(z$before[valid,j]),
          after=mean_valid(z$after[valid,j]),net=mean_valid(v),
          negative_magnitude=mean_valid(pmax(-v,0)),positive_magnitude=mean_valid(pmax(v,0)),
          allocation_se=if (alloc$exact) 0 else sd(draw)/sqrt(length(draw)),
          support=sum(valid),all_sets=nrow(sets),shapley_error=z$efficiency_error)
      }
    }
  }
}
write_result(rbindlist(out),"count_sensitivity")
