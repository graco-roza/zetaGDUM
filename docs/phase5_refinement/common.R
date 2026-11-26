source("docs/phase5/common.R")
write_refinement <- function(x, key) write_study(x,"docs/phase5_refinement/runs",key)

j_layout <- function(n, r, signal) {
  set.seed(3947L + as.integer(n) * 10000L + r * 140L)
  ref_env <- data.frame(g = seq(-3, 3, length.out = 81))
  ref <- ugdm_design(ref_env, cbind(1:80, 2:81))
  env <- data.frame(g = rnorm(n))
  u <- tanh(rnorm(n))
  z <- (env$g-ref$environment_basis$center)/ref$environment_basis$scale
  distance <- abs(outer(z,z,"-"))
  truth <- switch(signal, null=c(-.3,0,0,-.3,0,0),
    weak=c(-.3,.25,.12,-.3,.15,.25), strong=c(-.3,.7,.35,-.3,.45,.70))
  orders <- c(2L, if (n < 40L) 4L else round(.11*n))
  build <- function(idx,j) {
    k <- ncol(idx)
    spread <- internal("meanpair_dist_cpp")(distance,idx)
    pos <- rowSums(matrix(z[idx],nrow(idx)))
    X <- cbind(intercept=1,g_linear=spread,g.position=pos/k)
    eta <- drop(X %*% truth[(j-1L)*3L+1:3])
    mu <- plogis(eta)
    y <- mu + .85 * pmin(mu,1-mu) * rowMeans(matrix(u[idx],nrow(idx)))
    d <- list(idx=idx,n=n,k=k,Xdir=matrix(spread,ncol=1L,dimnames=list(NULL,"g_linear")),
      Xpos=matrix(pos,ncol=1L,dimnames=list(NULL,"g")),vars="g",spread_width=1L,
      environment_basis=ref$environment_basis)
    list(X=X,design=d,y=y)
  }
  gen <- function(M, sampler="uniform", exact=FALSE) {
    if (exact) return(lapply(seq_along(orders),function(j)
      build(t(combn(n,orders[j])),j)))
    s <- sample_sets_coupled(matrix(1,n,1),orders,M,sampler=sampler,exact_limit=0)
    lapply(seq_along(s),function(j) build(s[[j]]$idx,j))
  }
  list(pool=gen, exact_at=function(j) build(t(combn(n,orders[j])),j),
       truth=truth, orders=orders)
}

map6 <- rbind(lower_spread=c(0,1,0,0,0,0), lower_position=c(0,0,1,0,0,0),
  higher_spread=c(0,0,0,0,1,0), higher_position=c(0,0,0,0,0,1),
  paired_spread=c(0,-1,0,0,1,0), paired_position=c(0,0,-1,0,0,1))

centers <- function(a) sweep(a,2L,colMeans(a))
