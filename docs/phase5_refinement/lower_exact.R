source("docs/phase5_refinement/common.R")
Rcpp::sourceCpp("docs/manuscript_validation/calibration_helpers.cpp",
  cacheDir = Sys.getenv("CALIBRATION_CPP_CACHE",tempdir()))
out <- list()
for (scenario in list(c(70,"null"),c(70,"weak"),c(150,"strong"))) {
  n <- as.integer(scenario[1])
  signal <- scenario[2]
  for (r in 1:24) {
    z <- j_layout(n,r,signal)
    a <- z$exact_at(1L)
    ii <- a$design$idx
    m <- length(a$y)
    solve <- function(keep,start=numeric(3)) calibration_fit(a$X,a$y,unname(start),keep,rep(1,m))
    point <- solve(rep(TRUE,m))
    deletions <- t(vapply(seq_len(n),function(j)solve(rowSums(ii==j)==0,point),numeric(3)))
    V <- internal(".jackknife_sd")(deletions)
    bb <- ugdm_pseudovalue_bootstrap(point,deletions,draws=599,seed=984+r)$draws
    for (j in 2:3) for (level in c(.5,.8,.9,.95)) {
      q <- qt((1+level)/2,n-1)*V[j]
      b <- quantile(bb[,j],c((1-level)/2,(1+level)/2),names=FALSE)
      for (method in c("jackknife","parameter_BB")) {
        ci <- if (method=="jackknife") point[j]+c(-q,q) else b
        out[[length(out)+1L]] <- data.frame(n,signal,r,effect=c("lower_spread","lower_position")[j-1L],
          method,level,order1=2L,sets=m,estimate=point[j],truth=z$truth[j],
          width=diff(ci),covered=z$truth[j]>=ci[1] && z$truth[j]<=ci[2])
      }
    }
  }
}
write_refinement(do.call(rbind,out),"exact_lower")
