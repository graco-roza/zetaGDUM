test_that("coupled samplers have exact uniform and species-mixture marginals", {
  x <- 1L * cbind(c(1, 1, 1, 0, 0, 0), c(1, 1, 1, 1, 0, 1), c(0, 1, 1, 1, 1, 0))
  key <- function(s) apply(s, 1, function(r) paste(sort(r), collapse = "-"))
  for (method in c("uniform", "site", "species")) {
    s <- sample_sets_coupled(x, c(4, 2, 3), 24000L, method, alpha = c(.3,.2,.65),
                             exact_limit = 0, seed = 1729)
    for (j in seq_along(s)) {
      k <- ncol(s[[j]]$idx)
      all <- t(combn(6, k))
      share <- zeta_sets(x, all)$shared
      p <- rep(1 / nrow(all), nrow(all))
      h <- if (method == "species") p * (1 - s[[j]]$alpha + s[[j]]$alpha * share / mean(share)) else p
      id <- match(key(s[[j]]$idx), key(all))
      expect_lt(max(abs(tabulate(id, length(p)) / length(id) - h)), .017)
      expect_equal(s[[j]]$sampling_weights, (p / h)[id], tolerance = 1e-10)
      expect_lt(abs(mean(s[[j]]$sampling_weights * share[id]) - mean(share)), .04)
    }
    if (method != "species") expect_equal(s[["3"]]$idx, s[["4"]]$idx[, 1:3])
    if (method == "site") expect_true(all(tabulate(s[["4"]]$idx[1:6, ], 6) == 4L))
  }
  expect_true(sample_sets_coupled(x, 2, 20, exact_limit = 100)$`2`$exact)
})

test_that("incremental designs and responses agree with conventional evaluation", {
  sim <- simulate_metacommunity(14, 20, abundance_mean = 2, seed = 842)
  sim$comm <- cbind(2, sim$comm)
  dm <- coords_dist(sim$coords)
  s <- sample_sets_coupled(sim$comm, c(4, 2, 6), 150, exact_limit = 100, seed = 17)
  a <- ugdm_path_data(sim$comm, sim$env, s, dm)
  for (j in seq_along(a)) {
    idx <- s[[j]]$idx
    d <- ugdm_design(sim$env, idx, dm, ref = a[[j]]$design)
    expect_equal(d$Xdir, a[[j]]$design$Xdir, tolerance = 1e-10)
    expect_equal(d$Xpos, a[[j]]$design$Xpos, tolerance = 1e-10)
    expect_equal(d$Xdist, a[[j]]$design$Xdist, tolerance = 1e-10)
    expect_equal(ugdm_response(sim$comm, idx, "total")$y, a[[j]]$incidence, tolerance = 1e-10)
    expect_equal(ugdm_response(sim$comm, idx, "relative")$y, a[[j]]$relative, tolerance = 1e-10)
    expect_equal(ugdm_model_design(d, "B", sim$env)$Xdir,
                 ugdm_model_design(a[[j]]$design, "B")$Xdir, tolerance = 1e-10)
  }
  d <- a[[1]]$design
  r <- a[[1]]$relative
  for (model in c("A", "B", "C")) {
    dd <- ugdm_model_design(d, model)
    ff <- ugdm_fit_variant(dd, r, model = model, position = "sum")
    gg <- ugdm_fit_variant(dd, r, model = model, position = "mean")
    expect_equal(gg$position_coefficients, d$k * ff$position_coefficients)
    expect_equal(ugdm_predict(ff, dd), ugdm_predict(gg, dd))
    alloc <- ugdm_variant_diagnostics(ff, dd, r)
    expect_equal(sum(alloc$shapley), ff$D2)
    cf <- ugdm_configuration_effects(ff, dd, c(1, 2), c(3, 4))$effects
    expect_equal(cf$spread_response + cf$position_response, cf$dissimilarity_to - cf$dissimilarity_from)
  }
  expect_lt(max(abs(ugdm_predict(ugdm_fit(d, r, rep(1, length(r))), d) -
                  ugdm_predict(ugdm_fit_variant(d, r), d))), 1e-4)
})

test_that("score pseudovalues satisfy complete-set identities", {
  set.seed(92)
  n <- 10
  k <- 3
  env <- data.frame(g = rnorm(n))
  idx <- t(combn(n, k))
  d <- ugdm_model_design(ugdm_design(env, idx), "B", env)
  theta <- c(.2, .6, -.12)
  y <- runif(nrow(idx))
  X <- cbind(intercept = 1, d$Xdir, g.position = d$Xpos[, 1] / k)
  cur <- X * (y - plogis(drop(X %*% theta)))
  goods <- vapply(seq_len(n), function(i) colMeans(cur[rowSums(idx == i) == 0, ]), rep(0, 3))
  for (scheme in c("deletion", "centered", "inclusion")) {
    a <- ugdm_score_pseudovalues(d, y, theta, method = scheme)
    expect_equal(unname(a$deleted), unname(t(goods)), tolerance = 1e-11)
    expect_equal(colMeans(a$pseudovalues), a$full, tolerance = 1e-11)
    expect_equal(unname(a$pseudovalues),
                 sweep(-(n - 1) * unname(t(goods)), 2L, n * colMeans(cur), "+"), tolerance = 1e-11)
  }
})

test_that("Dirichlet covariance and empirical-likelihood constraints hold", {
  set.seed(744)
  a <- cbind(rnorm(12), rnorm(12))
  p <- ugdm_pseudovalue_bootstrap(c(0, 4), a, draws = 55000)
  expect_equal(p$bootstrap_covariance, p$jackknife_covariance * 11 / 13)
  expect_lt(max(abs(cov(p$draws) - p$bootstrap_covariance)), .09)
  m <- sweep(a, 2, colMeans(a) - c(.15, .1))
  for (method in c("EL", "ETEL")) {
    g <- ugdm_empirical_likelihood(m, method)
    expect_true(g$ok)
    expect_equal(sum(g$weights), 1)
    expect_equal(unname(colSums(m * g$weights)), c(0, 0), tolerance = 1e-7)
    expect_equal(g$loglik, sum(log(nrow(m) * g$weights)), tolerance = 1e-9)
    expect_false(ugdm_empirical_likelihood(matrix(1:12), method)$ok)
  }
})

test_that("numerical budget confirmation and ecological status are separate", {
  set.seed(9)
  f <- function(M, s) c(spread = mean(stats::rnorm(M, sd = if (s == "a") .06 else .12)))
  a <- ugdm_contrast_budget(f, c("a", "b"), c(100, 400, 1600),
    c(spread = .02), c(spread = 1e-8), c(spread = TRUE), nrep = 12, seed = 93)
  expect_equal(a$contrasts$tolerance, .02)
  expect_true(a$contrasts$passed)
  report <- ugdm_order_status(2:5, rep("spread", 4), c(TRUE, TRUE, FALSE, TRUE),
    c(TRUE, TRUE, TRUE, FALSE), c(TRUE, FALSE, TRUE, TRUE))
  expect_identical(report$profile$status, c("supported", "computationally unresolved",
    "ecologically uninformative", "statistical precision unresolved"))
})

test_that("joint JEL uses one weighted posterior on stacked site scores", {
  set.seed(181)
  n <- 24L
  env <- data.frame(g = rnorm(n))
  x <- matrix(1,n,1)
  u <- .15 * tanh(rnorm(n))
  p <- ugdm_path_data(x,env,sample_sets_coupled(x,c(2,4),450,exact_limit=0,seed=411))
  data <- lapply(p, function(a) {
    d <- ugdm_model_design(a$design,"B")
    eta <- -.3 + .47*d$Xdir[,1] + .25*d$Xpos[,1]/d$k
    list(design=d,y=plogis(eta) + rowMeans(matrix(u[d$idx],nrow(d$idx))))
  })
  out <- ugdm_jel(data,draws=220,seed=10)
  expect_equal(dim(out$draws),c(220L,6L))
  expect_equal(sum(out$weights),1)
  expect_gt(out$ess, 1)
  expect_true(all(is.finite(out$estimate)))
  data[[1]]$fit <- ugdm_fit_variant(data[[1]]$design,data[[1]]$y,model="B")
  expect_error(ugdm_jel(data, draws=20),"unconstrained")
})

test_that("constant responses have undefined deviance explained in all variants", {
  env <- data.frame(x = seq_len(8))
  idx <- t(combn(8, 2))
  d <- suppressWarnings(ugdm_design(env, idx))
  for (m in c("A", "B", "C")) {
    a <- ugdm_model_design(d, m, env)
    f <- ugdm_fit_variant(a, rep(1, nrow(idx)), model = m)
    expect_true(is.na(f$D2))
  }
  expect_error(ugdm_jel(list(list(design=d, y=rep(1,nrow(idx)))), draws=20),
               "nonconstant")
})

test_that("divergent species paths recompute the right set kernels", {
  sim <- simulate_metacommunity(14, 35, abundance_mean = 2, seed = 182)
  x <- cbind(sim$comm, 2)
  s <- sample_sets_coupled(x,c(2,4,6),85L,"species",alpha=c(.2,.8,.6),exact_limit=0,seed=821)
  a <- ugdm_path_data(x,sim$env,s)
  expect_true(any(a[[1]]$design$idx != a[[2]]$design$idx[,1:2]))
  for (j in seq_along(s)) {
    expect_equal(a[[j]]$shared, zeta_sets(x,s[[j]]$idx)$shared)
    expect_equal(a[[j]]$relative, ugdm_response(x,s[[j]]$idx,"relative")$y)
    expect_equal(a[[j]]$design$Xdir,
      ugdm_design(sim$env,s[[j]]$idx,ref=a[[j]]$design)$Xdir,tolerance=1e-10)
  }
})

test_that("linear spread supports exact cross-products and the site effect map", {
  s <- simulate_metacommunity(12, 26, abundance_mean = 2, seed = 111)
  x <- cbind(s$comm,2)
  idx <- t(combn(12,3))
  d <- ugdm_model_design(ugdm_design(s$env,idx,coords_dist(s$coords)),"B",s$env)
  r <- ugdm_response(x,idx,"relative")
  f <- ugdm_fit(d,r$y,r$trials,free_sign=TRUE)
  X <- zetaGDUM:::.order_terms(d,f)
  exact <- ugdm_crossproducts(x,s$env,d,"relative")$response
  use <- is.finite(exact)
  expect_equal(exact[use],colMeans(X*r$y)[use],tolerance=1e-10)
  map <- zetaGDUM:::.order_effect_map(d,f)
  expect_equal(drop(map %*% c(f$intercept,f$beta,f$lambda)),zetaGDUM:::.order_effects(f))
})
