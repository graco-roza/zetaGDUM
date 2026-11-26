source("docs/phase5_refinement/common.R")
n <- as.integer(commandArgs(TRUE)[1L])
R <- 24L
orders <- c(3L,round(.11*n))
out <- list()
mnames <- c("A","B","adjacent","curvature","ridge","combined")
for (kind in c("linear","nonlinear","collinear")) for (r in seq_len(R)) {
  set.seed(5012L+n*1000L+r)
  x <- rnorm(n)
  z <- if (kind=="collinear") x+.12*rnorm(n) else rnorm(n)
  env <- data.frame(g=x,h=z)
  pa <- matrix(1,n,1)
  dat <- suppressWarnings(ugdm_path_data(pa,env,sample_sets_coupled(pa,orders,2400,exact_limit=0)))
  ref <- suppressWarnings(ugdm_path_data(pa,env,sample_sets_coupled(pa,orders,14000,exact_limit=0),
                                     ref=dat[[1]]$design))
  u <- tanh(rnorm(n))
  for (j in seq_along(dat)) {
    d <- dat[[j]]$design
    basis <- d$environment_basis
    q <- (quantile(x,c(.25,.75),names=FALSE)-basis$center[1])/basis$scale[1]
    scale_spread <- .1*colSums(diff(predict(basis$splines[[1]],q)))
    truth_b <- c(.1,1.7,.12) * if (j==1L) 1 else .8
    slope <- if (j==1L) .42 else .72
    g <- function(a) if (kind=="linear") a$Xlinear[,1]*.75 else drop(a$Xdir[,1:3]%*%truth_b)
    mu <- function(a) plogis(-.35 + g(a) + slope*a$Xpos[,1]/a$k)
    b_ref <- if (kind=="linear") "B" else "A"
    dr <- ugdm_model_design(ref[[j]]$design,b_ref)
    ry <- mu(ref[[j]]$design)
    fit_r <- ugdm_fit_variant(dr,ry,model=b_ref)
    allocated <- ugdm_variant_diagnostics(fit_r,dr,ry)$relative_shapley[["position"]]
    truth <- c(spread=if (kind=="linear") .075*diff(q) else sum(scale_spread*truth_b),
      position=slope*diff(q)/d$k,relative_position=allocated)
    eta <- mu(d)
    y <- eta + .85*pmin(eta,1-eta)*rowMeans(matrix(u[d$idx],nrow(d$idx)))
    for (m in mnames) {
      model <- if (m %in% c("A","B")) m else "C"
      d0 <- ugdm_model_design(d,model)
      smooth <- if (m=="ridge") 0 else if (m=="adjacent") .001 else .008
      ridge <- if (m %in% c("ridge","combined")) .0001 else 0
      regularization <- if (m %in% c("curvature","combined")) "curvature" else "adjacent"
      f <- ugdm_fit_variant(d0,y,model=model,smooth=smooth,ridge=ridge,regularization=regularization)
      dg <- ugdm_variant_diagnostics(f,d0,y)
      b <- f$beta[1:3]
      b0 <- if (m=="B") .1*f$beta[1]*diff(q) else sum(b*scale_spread)
      measured <- c(b0,unname(f$lambda[1])*diff(q),dg$relative_shapley[["position"]])
      curvature <- if (m=="B") 0 else drop(crossprod(b,ugdm_curvature_penalty(d)[[1]]%*%b))
      out[[length(out)+1L]] <- data.frame(n,kind,r,order=d$k,model=m,
        effect=names(truth),estimate=measured,truth=unname(truth),curvature,
        condition=dg$condition,converged=f$converged)
    }
  }
}
write_refinement(do.call(rbind,out),paste0("curvature_",n))
