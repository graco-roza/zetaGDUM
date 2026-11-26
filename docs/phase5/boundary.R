source("docs/phase5/common.R")
n <- 100L
R <- 36L
orders <- c(4L, 12L)
out <- list()
for (r in seq_len(R)) {
  set.seed(24953 + r)
  env <- data.frame(g = rnorm(n))
  u <- tanh(rnorm(n))
  x <- matrix(1,n,1)
  A <- ugdm_path_data(x, env, sample_sets_coupled(x,orders,1100,exact_limit=0))
  weights <- site_weights(n,799)
  theta <- c(0, .30)
  estimate <- function(X, Y, w, start=NULL) {
    internal(".variant_solution")(X,Y,w,matrix(0,3,3),2L,start)$coef
  }
  f <- vector("list",2)
  for (j in 1:2) {
    d <- ugdm_model_design(A[[j]]$design,"B")
    X <- cbind(1,d$Xdir,d$Xpos/d$k)
    mu <- plogis(-.3 + theta[j] * X[,2] + .25 * X[,3])
    y <- mu + .8*pmin(mu,1-mu)*rowMeans(matrix(u[d$idx],nrow(X)))
    w <- rep(1,length(y))
    a <- estimate(X,y,w)
    v <- vapply(seq_len(n),function(i) {
      use <- rowSums(d$idx == i) == 0L
      estimate(X[use,],y[use],w[use],a)[2]
    },0.0)
    f[[j]] <- list(point=a[2],delete=v)
  }
  all <- c(f,list(list(point=f[[2]]$point-f[[1]]$point,delete=f[[2]]$delete-f[[1]]$delete)))
  truth <- c(theta,.30)
  for (j in seq_along(all)) {
    a <- all[[j]]
    jk <- internal(".jackknife_sd")(matrix(a$delete,ncol=1))
    bb <- c(ugdm_pseudovalue_bootstrap(a$point,a$delete,weights=weights)$draws)
    out[[length(out)+1L]] <- data.frame(r, contrast=c("lower","upper","change")[j],truth=truth[j],
      point=a$point, active=a$point < 1e-8,
      jk_cover=abs(a$point-truth[j])<=qt(.975,n-1)*jk,
      pseudo_cover=quantile(bb,.025)<=truth[j] && quantile(bb,.975)>=truth[j],
      negative_draws=mean(bb<0), practical_positive=mean(bb>.1))
  }
}
write5(do.call(rbind,out),"boundary")
