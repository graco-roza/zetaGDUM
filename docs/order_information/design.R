library(zetaGDUM)

output <- "docs/order_information/results"
scenarios <- data.frame(name = c("narrow", "broad", "dominant", "correlated", "spatial", "null"),
  species = c(32, 64, 48, 80, 48, 40), breadth = c(.26, .85, .32, .50, .62, .65),
  occupancy = c(-.25, 1.8, -.8, .45, .75, .7), effect = c(1, .65, .80, .42, .85, 0),
  coupling = c(.6, -.30, .5, 1.1, .2, 0), noise = c(.28, .50, .04, .75, .45, .5),
  correlation = c(0, 0, 0, .83, .25, 0), spatial = c(0, 0, 0, 0, .6, 0))
scenario_index <- function(name) match(name, scenarios$name)
order_grid <- function(name) if (name %in% c("broad", "dominant")) 2:40 else 2:30

community <- function(n, name, seed) {
  set.seed(seed)
  d <- scenarios[scenario_index(name), ]
  Z <- matrix(rnorm(n * 2L), n)
  if (d$spatial > 0) for (i in 2:n) {
    Z[i, ] <- d$spatial * Z[i-1L, ] + sqrt(1-d$spatial^2) * Z[i, ]
  }
  Z[, 2] <- d$correlation * Z[, 1] + sqrt(1-d$correlation^2) * Z[, 2]
  env <- setNames(as.data.frame(2 * pnorm(Z) - 1), c("gradient", "secondary"))
  j <- seq_len(d$species)
  prevalence <- d$occupancy + 1.15 * qnorm((j-.5)/length(j))
  opt1 <- seq(-1,1,length.out=length(j))
  opt2 <- opt1[(seq_along(j) * 13L) %% length(j) + 1L]
  niche <- -d$effect * (outer(env$gradient, opt1, "-")^2 / d$breadth^2 +
        .36 * outer(env$secondary, opt2, "-")^2)
  occ <- plogis(sweep(niche, 2, prevalence, "+"))
  p <- matrix(runif(n * length(j)), n)
  noise <- matrix(rnorm(n * length(j)), n)
  q <- (p < occ) * exp(rep(d$coupling * prevalence, each=n) + d$noise * noise + .30*niche)
  background <- if (name == "dominant") {
    cbind(500 + rep(c(1, 1.01), length.out = n), 100)
  } else matrix(.5, n, 1)
  list(x = cbind(background, q), env = env,
       blocks = cut(seq_len(n), breaks = max(5L, n %/% 7L), labels=FALSE))
}

ref_env <- data.frame(gradient=seq(-1,1,length.out=101L), secondary=seq(1,-1,length.out=101L)[
  (seq_len(101L)*19L) %% 101L + 1L])
reference_basis <- ugdm_design(ref_env, cbind(1:100, 2:101))
layout <- getFromNamespace(".interaction_layout", "zetaGDUM")
solve_coalition <- getFromNamespace(".variant_solution", "zetaGDUM")
score_inputs <- getFromNamespace(".within_score_inputs", "zetaGDUM")
score_state <- getFromNamespace(".within_score_state", "zetaGDUM")
build <- function(a, k, M, response, seed = NULL, idx = NULL) {
  if (is.null(idx)) idx <- sample_sets(nrow(a$x), k, M, seed = seed)
  base <- ugdm_design(a$env, idx, ref = reference_basis)
  d <- ugdm_interaction_design(base,a$env, "nonlinear", "constant")
  list(design = d, y = ugdm_response(a$x, idx, response, working_weights = "constant")$y,
       X = layout(d, FALSE)$X)
}
fit_game <- function(a) {
  p <- ncol(a$X)
  # The package additive layout puts the intercept first, then directional and position terms.
  group1 <- c(2:4,8:10)
  group2 <- c(5:7,11:13)
  stopifnot(p == 13L)
  keep <- list(1L, c(1L,group1), c(1L,group2), seq_len(p))
  P <- diag(c(0,rep(.001,p-1L)))
  w <- rep(1,nrow(a$X))
  fit <- function(i) {
    J <- keep[[i]]
    solve_coalition(a$X[,J,drop=FALSE], a$y, w, P[J,J,drop=FALSE],
                    which(J %in% 2:7))
  }
  fits <- lapply(seq_along(keep), fit)
  list(fits=fits, keep=keep, w=w, P=P)
}

evaluate_game <- function(game, b, draws = NULL) {
  softplus <- function(z) pmax(z,0)+log1p(exp(-abs(z)))
  loss <- function(i) {
    coeff <- if (is.null(draws)) game$fits[[i]]$coef else draws[[i]]
    eta <- b$X[,game$keep[[i]],drop=FALSE] %*% coeff
    colMeans(softplus(eta) - b$y * eta)
  }
  D <- lapply(seq_along(game$keep), loss)
  saturated <- mean(ifelse(b$y > 0 & b$y < 1, -b$y*log(b$y)-(1-b$y)*log1p(-b$y),0))
  den <- D[[1L]]-saturated
  scale <- function(x) ifelse(den > .Machine$double.eps, x / den, NA_real_)
  a <- scale((D[[1L]]-D[[2L]]+D[[3L]]-D[[4L]])/2)
  c <- scale((D[[1L]]-D[[3L]]+D[[2L]]-D[[4L]])/2)
  cbind(gradient=a, secondary=c, difference=a-c, D2=scale(D[[1L]]-D[[4L]]))
}

score_draws <- function(a, pool, game, groups, B) {
  p <- list(X=a$X,y=a$y,weights=game$w,penalty=game$P)
  n <- length(groups)
  s <- score_inputs(p, a$design, n)
  pools <- list(s, score_inputs(list(X=pool$X,y=pool$y,weights=rep(1,length(pool$y))), pool$design,n))
  G <- max(groups)
  z <- matrix(rnorm(G*B), G, B)
  lapply(seq_along(game$keep), function(i) {
    J <- game$keep[[i]]
    state <- score_state(p,s,J,game$fits[[i]]$coef,which(J %in% 2:7),pools,n,groups,NULL)
    proposed <- state$beta + state$inv %*% (state$g + crossprod(state$scores,z))
    if (i == 1L) return(proposed)
    vapply(seq_len(B), function(j) state$project(proposed[,j]), numeric(length(J)))
  })
}
write_result <- function(x,path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  temporary <- tempfile(tmpdir=dirname(path))
  con <- if (endsWith(path,".gz")) gzfile(temporary,"wt") else file(temporary,"wt")
  write.csv(x, con, row.names=FALSE)
  close(con)
  stopifnot(file.rename(temporary,path))
}
