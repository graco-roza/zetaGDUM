source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")
source(file.path(directory,"scripts","trajectory.R"))
load_kernel()
x <- readRDS(file.path(processed,"biotime-604.rds"))
p <- x$panel
n <- dim(p$x)[1]
y <- p$times
contrasts <- list(first_annual=c(1,2),mid_annual=c(5,6),endpoint=c(1,length(y)))
out <- list()
for (name in names(contrasts)) {
  j <- contrasts[[name]]
  plan <- st3_plan(p,j[1],j[2])
  allocations <- st3_allocations(plan,draws=512L,seed=730L)
  for (k in 2:6) {
    sets <- t(combn(n,k))
    z <- evaluate(plan,sets,allocations)
    for (budget in c(32L,128L,512L)) {
      draw <- z$draws[,,,seq_len(min(budget,dim(z$draws)[4L])),drop=FALSE]
      q <- z
      q$draws <- draw
      q$phi <- apply(draw,c(1,2,3),mean)
      out[[length(out)+1L]] <- summarize_phi(q,sets,x$stratum,exhaustive=TRUE)[,
        `:=`(k=k,contrast=name,budget=budget,actual_draws=dim(draw)[4L],
            initial=y[j[1]],final=y[j[2]])]
    }
  }
}
write_result(rbindlist(out),"allocation_sensitivity")
