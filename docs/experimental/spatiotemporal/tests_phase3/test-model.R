test_that("native environmental columns and common-penalty fitted totals agree", {
  sim <- st3_simulate("expansion_numerical",seed=66)
  p <- sim$panel
  w <- st3_pair_windows(p$times)
  sets <- st3_sample_sets(length(sim$blocks),2,40)
  d <- st3_environment_design(p,sim$env,sets,w)
  native <- zetaGDUM::ugdm_design(matrix(sim$env,ncol=2,dimnames=list(NULL,c("climate","soil"))),
                                 sets,ref=list(environment_basis=d$basis$native))
  expect_equal(unname(d$X[seq_len(nrow(sets)),"position.climate"]),unname(native$Xpos[,1L]))
  expect_equal(unname(d$X[seq_len(nrow(sets)),paste0("spread.climate_f",1:3)]),unname(native$Xdir[,1:3]))
  y <- st3_collect(p,sets,w,crossed=TRUE)
  fit <- st3_fit_process_gdum(d,y$y)
  total <- st3_fit_process_gdum(d,cbind(total=rowSums(y$y)))
  curves <- st3_predict_process_gdum(fit,d,y$before)
  expect_equal(curves$total,drop(st3_predict_process_gdum(total,d)$contributions),ignore_attr=TRUE)
  expect_equal(rowSums(y$y),y$after-y$before)
  expect_true("log_lag" %in% fit$dropped)
  expect_true(all(curves$extrapolation==FALSE))
  expect_true(fit$effective_df<=ncol(fit$X))
  contrast <- colMeans(fit$X[1:15,,drop=FALSE])-colMeans(fit$X[16:30,,drop=FALSE])
  u <- st3_contrast_uncertainty(fit,contrast,sim$blocks,draws=5)
  expect_equal(dim(u$se),c(8,3))
  expect_true(all(is.finite(u$se)))
})

test_that("designed opposing changes cancel net abundance dissimilarity", {
  for (seed in 1:3) {
    z <- st3_simulate("opposing",seed=seed)
    for (k in 2:6) {
      sets <- st3_sample_sets(length(z$blocks),k,45)
      a <- st3_collect(z$panel,sets,st3_pair_windows(z$panel$times)[1:3])
      expect_lt(max(abs(rowSums(a$y))),1e-12)
      expect_gt(mean(a$y[,"numerical"]),.05)
      expect_lt(mean(a$y[,"expansion"]),-.05)
      expect_equal(unname(a$y[,c("contraction","redistribution")]),matrix(0,nrow(a$y),2))
    }
  }
})

test_that("pure simulations agree with analytical zeta and all supported windows conserve diversity", {
  for (case in c("expansion","contraction","redistribution","numerical")) {
    s <- st3_simulate(case,n=8,surveys=8)
    p <- s$panel
    sets <- t(combn(8,3))
    plan <- st3_plan(p,1,2)
    x <- st3_shapley(plan,sets)
    expected <- sum(.st_choose_prob(colSums(plan$after>0),3,8)-.st_choose_prob(colSums(plan$before>0),3,8))
    expect_equal(mean(x$delta[,"zeta_presence"]),expected)
    expect_equal(mean(x$contributions[,"zeta_presence",case]),expected)
    expect_true(all(x$contributions[,,setdiff(st3_processes,case)]==0,na.rm=TRUE))
    for (m in 1:4) {
      w <- st3_pair_windows(p$times,m)[[1]]
      z <- st3_shapley(st3_plan(p,w$before,w$after),sets)
      expect_lt(z$max_efficiency_error,1e-12)
    }
  }
})

test_that("weighted fits use common support and flag invalid extrapolated totals", {
  sim <- st3_simulate("filtering",seed=900)
  p <- sim$panel
  s <- st3_sample_sets(20,2,14)
  w <- st3_pair_windows(p$times)[1:4]
  design <- st3_environment_design(p,sim$env,s,w)
  y <- cbind(expansion=rep(1.1,nrow(design$X)),numerical=rep(.4,nrow(design$X)))
  y[1,1] <- NA_real_
  fit <- st3_fit_process_gdum(design,y,weights=seq_len(nrow(y)))
  expect_false(fit$keep[1])
  scores <- st3_predict_process_gdum(fit,design,rep(.4,nrow(y)))
  expect_true(all(scores$incompatible))
  z <- design
  z$X[,"position.climate"] <- 100
  expect_true(all(st3_predict_process_gdum(fit,z)$extrapolation))
  contrast <- colMeans(fit$X[1:5,,drop=FALSE])-colMeans(fit$X[6:10,,drop=FALSE])
  total <- st3_fit_process_gdum(design,cbind(total=rowSums(y)),weights=seq_len(nrow(y)))
  expect_equal(dim(st3_contrast_uncertainty(total,contrast,sim$blocks,draws=5)$se),c(1,3))
})
