component_fixture <- function(k = 2L) {
  set.seed(42)
  n <- 31L
  env <- data.frame(phosphorus = runif(n, 10, 75), water = rnorm(n))
  x <- cbind(1, matrix(rpois(n * 9, 1), n))
  idx <- sample_sets(n, k, 200)
  list(x = x, env = env, idx = idx, y = ugdm_component_response(x, idx))
}

test_that("component partitions and coherent predictions close exactly", {
  a <- component_fixture(4L)
  y <- a$y
  expect_equal(y[, "total"], ugdm_response(a$x, a$idx, "relative")$y)
  expect_equal(y[, "occupancy"] + y[, "heterogeneity"], y[, "total"])
  d <- ugdm_component_design(a$env, a$idx)
  f <- ugdm_fit_components(d, y)
  p <- ugdm_predict_components(f, d)
  expect_true(f$converged)
  expect_true(all(as.matrix(p[c("shared","occupancy","heterogeneity")]) >= 0))
  expect_equal(rowSums(p[c("shared","occupancy","heterogeneity")]), rep(1, nrow(y)))
  g <- ugdm_fit_components(d, y, "independent")
  q <- ugdm_predict_components(g, d)
  expect_equal(p$total, q$total)
  expect_true(mean(abs(q$closure_error)) > 1e-6)
  idx <- a$idx[nrow(a$idx):1L, ]
  expect_error(ugdm_fit_components(ugdm_component_design(a$env, idx), y), "set order")
})

test_that("zero dissimilarity and share boundaries have exact predictions", {
  a <- component_fixture()
  d <- ugdm_component_design(a$env, a$idx)
  x <- matrix(1, nrow(a$x), 2)
  y <- ugdm_component_response(x, a$idx)
  f <- ugdm_fit_components(d,y)
  p <- ugdm_predict_components(f,d)
  expect_false(f$share_identified)
  expect_true(all(p$total == 0 & p$shared == 1 & p$heterogeneity == 0))
  expect_true(all(is.na(ugdm_assess_components(f,d,y,y)$metrics$D2)))
  x <- diag(nrow(a$x))
  y <- ugdm_component_response(x,a$idx)
  f <- ugdm_fit_components(d,y)
  p <- ugdm_predict_components(f,d)
  expect_true(f$share_identified)
  expect_true(all(p$occupancy == 1 & p$heterogeneity == 0))
  y <- cbind(total=c(0,.6,1),occupancy=c(0,.5,.7),heterogeneity=c(0,.1,.3))
  f <- zetaGDUM:::.component_fits(matrix(1,3,1),matrix(0,1,1),y,c(3,1,1),"coherent",integer())
  expect_equal(plogis(f$share$coef), .4/1.6, tolerance=1e-6)
})

test_that("interaction restrictions and height of the component hierarchy are coherent", {
  a <- component_fixture()
  main <- ugdm_component_design(a$env,a$idx)
  full <- ugdm_component_design(a$env,a$idx,"interaction",reference=main)
  p <- zetaGDUM:::.component_problem(main,.001,.01)
  q <- zetaGDUM:::.component_problem(full,.001,.01)
  keep <- setdiff(seq_len(ncol(q$X)), unlist(lapply(q$groups,"[[","interaction")))
  A <- zetaGDUM:::.component_transform(q,keep)
  expect_equal(unname(q$X %*% A), unname(p$X), tolerance=1e-12)
  expect_equal(crossprod(A,q$P %*% A), p$P, tolerance=1e-12)
  f <- ugdm_fit_components(full,a$y)
  h <- ugdm_component_importance(f,full,a$y)
  expect_length(h$diagnostics$failed_coalitions,0)
  expect_true("interaction" %in% h$hierarchy$term)
  expect_equal(tapply(h$hierarchy$absolute,h$hierarchy$response,sum),
    tapply(h$hierarchy$D2,h$hierarchy$response,unique), tolerance=1e-9)
  expect_equal(tapply(h$predictors$absolute,h$predictors$response,sum),
    tapply(h$predictors$D2,h$predictors$response,unique), tolerance=1e-9)
  g <- ugdm_component_importance(f,full,a$y,exact_max=0,permutations=10,seed=85)
  expect_equal(tapply(g$hierarchy$absolute,g$hierarchy$response,sum),
    tapply(g$hierarchy$D2,g$hierarchy$response,unique), tolerance=1e-9)
  expect_true(all(is.finite(g$hierarchy$mc_se)))
  q <- ugdm_component_diagnostics(f,full,a$y)
  expect_equal(nrow(q$interaction_residual),6L)
  expect_true(all(q$stages$information_rank <= q$stages$parameters))
})

test_that("surface support counts distinct training sets and masks extrapolation", {
  a <- component_fixture()
  d <- ugdm_component_design(a$env,a$idx)
  f <- ugdm_fit_components(d,a$y)
  shifted <- a$env
  shifted$phosphorus <- shifted$phosphorus + 1e3
  g <- ugdm_component_design(shifted,a$idx,reference=f)
  expect_false(any(ugdm_component_support(f,g)$supported))
  expect_true(all(is.na(ugdm_predict_components(f,g,mask=TRUE)$total)))
  s <- ugdm_component_surface(f,d,"phosphorus",bins=4,min_sets=1,min_sites=2)
  expect_equal(nrow(s),16L)
  expect_true(any(s$supported))
  expect_true(all(is.na(s$total[!s$supported])))
  idx <- matrix(rep(c(1L,2L),times=15L),ncol=2L,byrow=TRUE)
  d <- ugdm_component_design(a$env,idx)
  f <- ugdm_fit_components(d,ugdm_component_response(a$x,idx))
  expect_false(any(ugdm_component_support(f,d,min_sets=2,min_sites=2)$supported))
})

test_that("spatial validation and existing confirmation accept component fits", {
  a <- component_fixture()
  q <- ugdm_validate_components(a$x,a$env,fold=rep(1:3,length.out=nrow(a$x)),
    combinations=80,seed=12)
  expect_equal(nrow(q),3L*3L*2L*3L)
  expect_true(all(q$closure_RMSE[q$architecture=="coherent"] < 1e-14))
  d <- ugdm_component_design(a$env,a$idx)
  f <- ugdm_fit_components(d,a$y)
  v <- ugdm_within_order_diagnostic(f,d,a$y[,"occupancy"],NULL,NULL,
    heldout_D2=.15,min_response_sd=0,min_heldout_D2=0,max_mc_se=.03,
    max_width=.2,max_condition=1e7,component="occupancy")
  expect_identical(v$status,"prediction_only")
  weights <- seq_len(nrow(a$y))
  f <- ugdm_fit_components(d,a$y,sampling_weights=weights)
  v <- ugdm_within_order_diagnostic(f,d,a$y[,"occupancy"],NULL,NULL,
    heldout_D2=.15,min_response_sd=0,min_heldout_D2=0,max_mc_se=.03,
    max_width=.2,max_condition=1e7,component="occupancy")
  expect_equal(v$sampling_ESS,sum(weights)^2/sum(weights^2))
  bad <- a$x
  bad[1L,1L] <- NA
  expect_error(ugdm_validate_components(bad,a$env,fold=rep(1:3,length.out=nrow(a$x))), "complete")
})
