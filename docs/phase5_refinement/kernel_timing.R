source("docs/phase5_refinement/common.R")
z <- j_layout(70L,2L,"weak")
out <- list()
for (M in c(900L,5000L,20000L,50000L)) {
  a <- z$pool(M)
  p <- lapply(a,function(d) internal(".prepare_score_pools")(list(d))[[1L]])
  theta <- replicate(50L,rnorm(6L,sd=.1)+z$truth)
  calc <- function(f) lapply(seq_len(ncol(theta)),function(t) lapply(seq_along(p),function(j)
    f(p[[j]],theta[(j-1L)*3L+1:3,t])))
  old <- function(p,theta) internal("score_inclusion_cpp")(p$idx,p$X,
    p$effective*(p$y-plogis(drop(p$X%*%theta))),p$w,p$n)[,-1L,drop=FALSE]
  fused <- function(p,theta) internal("site_score_kernel_cpp")(p$sites_by_row,p$design_by_row,
    p$y,p$effective,theta,p$n)
  time <- function() sum(proc.time()[1:2])
  tick <- time(); prev <- calc(old); slow <- (time()-tick)*20
  tick <- time(); next_score <- calc(fused); fast <- (time()-tick)*20
  error <- max(abs(unlist(prev)-unlist(next_score)))
  stopifnot(error<1e-10)
  out[[length(out)+1L]] <- data.frame(M,old_ms=slow,fused_ms=fast,error)
}
write_refinement(do.call(rbind,out),"kernel_timing")
