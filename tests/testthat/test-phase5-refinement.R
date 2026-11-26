test_that("integrated curvature measures the actual spline derivative", {
  env <- data.frame(a = seq(-2,2,length.out=20), b = seq(-1,1,length.out=20)^3)
  d <- ugdm_design(env, t(combn(20, 2)))
  out <- ugdm_curvature_penalty(d)
  expect_identical(names(out),names(env))
  for (j in seq_along(out)) {
    C <- out[[j]]
    stopifnot(inherits(d$environment_basis$splines[[j]],"ISpline"))
    b <- attr(d$environment_basis$splines[[j]],"Boundary.knots")
    beta <- c(.12,1.9,.30)
    f <- function(z) drop(splines2::iSpline(z,degree=2,Boundary.knots=b,derivs=2)%*%beta)^2
    expect_equal(drop(crossprod(beta,C%*%beta)),integrate(f,b[1],b[2])$value,tolerance=1e-10)
    expect_true(all(eigen(C,symmetric=TRUE)$values >= -1e-12))
    # Verify the null function, without assuming its coefficient pattern.
    ee <- eigen(C,symmetric=TRUE)
    z <- seq(b[1],b[2],length.out=25)
    f0 <- drop(splines2::iSpline(z,degree=2,Boundary.knots=b)%*%ee$vectors[,3])
    expect_lt(max(abs(resid(lm(f0~z)))),1e-10)
  }
  y <- plogis(.5+d$Xdir[,2]*.4-d$Xpos[,1]*.2)
  c <- ugdm_fit_variant(d,y,model="C",regularization="curvature")
  expect_true(c$converged)
  expect_true(all(c$beta>=0))
  expect_equal(c$penalty[2:4,2:4],.001*out[[1]]+diag(.0001,3),ignore_attr=TRUE)
})

test_that("pseudovalue diagnostics centre each independent pool", {
  x <- matrix(c(4,3,8,6,-1,10),3,2)
  y <- matrix(c(2,3,7,5,4,18),3,2)
  q <- ugdm_pseudovalue_noise(x,y)
  s <- function(a) sweep(a,2,colMeans(a))
  expect_equal(q$numerical_covariance,crossprod(s(x)-s(y))/12)
  expect_equal(q$pooled_jackknife,crossprod((s(x)+s(y))/2)/6)
  expect_equal(q$average_jackknife-q$pooled_jackknife,q$numerical_covariance/2)
  expect_equal(ugdm_pseudovalue_noise(sweep(x,2,c(-12,47),"+"),y),q)
})

test_that("complete conditional inclusion and matched deletions agree", {
  set.seed(83)
  n <- 9
  idx <- t(combn(n,3))
  env <- data.frame(g=rnorm(n))
  d <- ugdm_model_design(ugdm_design(env,idx),"B",env)
  y <- runif(nrow(idx))
  a <- list(design=d,y=y)
  b <- c(.1,.5,-.2)
  p <- ugdm_pool_scores(list(a,a),b,"conditional")
  o <- ugdm_score_pseudovalues(d,y,b,method="deletion")
  expect_equal(p,o$pseudovalues,tolerance=1e-11)
  w <- runif(nrow(idx),.5,3)
  X <- cbind(1,d$Xdir,d$Xpos/d$k)
  R <- X*(w*(y-plogis(drop(X%*%b))))
  mean <- colMeans(R)
  cond <- t(vapply(seq_len(n),function(i) colSums(R[rowSums(idx==i)>0,,drop=FALSE])*n/3/nrow(X),numeric(3)))
  v <- sweep(sweep(cond,2,mean)*3*(n-1)/(n-3),2,mean,"+")
  inc <- ugdm_score_pseudovalues(d,y,b,sampling_weights=w,method="inclusion")
  expect_equal(unname(inc$pseudovalues),unname(v),tolerance=1e-10)
})

test_that("adjusted JEL resolves a certified separating moment", {
  z <- matrix(2:22,ncol=1L)
  a <- ugdm_empirical_likelihood(z,certify=TRUE)
  expect_identical(a$reason,"zero strictly outside the moment hull")
  b <- ugdm_empirical_likelihood(z,adjusted=TRUE)
  x <- c(z,-log(length(z))/2*mean(z))
  expect_true(b$ok)
  expect_equal(sum(b$weights*x),0,tolerance=1e-7)
  expect_equal(sum(b$weights),1)
})

test_that("refined importance and independent chains target the same moments", {
  set.seed(94)
  n <- 24
  env <- data.frame(g=rnorm(n))
  x <- matrix(1,n,1)
  e <- .11*tanh(rnorm(n))
  a <- ugdm_path_data(x,env,sample_sets_coupled(x,c(2,4),400,exact_limit=0,seed=18))
  d <- lapply(a,function(q) {
    u <- ugdm_model_design(q$design,"B")
    y <- plogis(-.4+u$Xdir[,1]*.55+u$Xpos[,1]/u$k*.38) +
      rowMeans(matrix(e[u$idx],nrow(u$idx)))
    list(design=u,y=y)
  })
  i <- ugdm_jel_refined(d,draws=120,pilot=75,seed=8)
  m <- ugdm_jel_refined(d,integration="mcmc",draws=350,warmup=250,chains=3,seed=5)
  expect_equal(sum(i$weights),1)
  expect_equal(length(i$batch_ess),2)
  expect_equal(sum(i$likelihood_status),240)
  expect_gt(i$ess,5)
  expect_equal(length(m$rhat),6)
  expect_true(all(is.finite(m$rhat)))
  expect_lt(abs(sum(i$weights*i$draws[,3])-mean(m$draws[,3])),.15)
})
