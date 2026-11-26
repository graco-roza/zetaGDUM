source("docs/phase5/common.R")
args <- commandArgs(TRUE)
n <- as.integer(args[1])
R <- if (length(args) > 1) as.integer(args[2]) else 24L
out <- list()
for (case in c("null", "linear_weak", "linear_strong", "nonlinear", "collinear")) {
  for (r in seq_len(R)) {
    set.seed(187 + n * 1000 + r)
    z <- rnorm(n)
    env <- data.frame(g = z, h = if (case == "collinear") z + .12 * rnorm(n) else rnorm(n))
    orders <- c(3L, round(.11 * n))
    x <- matrix(1, n, 1)
    s <- ugdm_path_data(x, env, sample_sets_coupled(x, orders, 2100, exact_limit = 0))
    u <- tanh(rnorm(n))
    for (j in seq_along(orders)) {
      d <- s[[j]]$design
      linear <- case %in% c("null","linear_weak","linear_strong")
      true_b <- switch(case, null = c(0,0,0), linear_weak = c(.18,0,0),
        linear_strong = c(.85,0,0), nonlinear = c(.1,1.7,.12), collinear = c(.1,1.7,.12))
      true_b <- true_b * (if(j == 1L) 1 else .8)
      lambda <- switch(case, null = 0, linear_weak = .14, linear_strong = .42, .42) *
        if (j == 1L) 1 else 1.8
      A <- if (linear) d$Xlinear[,1L] * true_b[1] else drop(d$Xdir[,1:3] %*% true_b)
      mu <- plogis(-.3 + A + lambda * d$Xpos[,1] / d$k)
      y <- mu + .85 * pmin(mu, 1 - mu) * rowMeans(matrix(u[d$idx], nrow(d$idx)))
      q <- (quantile(env$g,c(.25,.75),names=FALSE) - d$environment_basis$center[1]) /
        d$environment_basis$scale[1]
      truth <- c(spread = .1 * if (linear) diff(q) * true_b[1] else
                   sum(true_b * colSums(diff(predict(d$environment_basis$splines[[1]],q)))),
                 position = unname(diff(q) * lambda / d$k))
      for (m in c("A","B","C","C_smooth")) {
        actual <- if (m == "C_smooth") "C" else m
        dm <- ugdm_model_design(d, actual)
        t0 <- sum(proc.time()[1:2])
        f <- ugdm_fit_variant(dm, y, model = actual, smooth = .001, ridge = if (m=="C_smooth") 0 else .0001)
        fx <- effect5(list(design = dm, y = y, trials = rep(1,length(y)),sampling_weights = rep(1,length(y))),
                      f, env)
        out[[length(out) + 1L]] <- data.frame(n, case, r, order = d$k, model = m,
          effect = names(truth), truth = unname(truth), point = unname(fx[1:2]),
          pred_RMSE = sqrt(mean((ugdm_predict(f,dm) - mu)^2)), converged = f$converged,
          condition = unname(fx[8]), active = unname(fx[9]), cpu = sum(proc.time()[1:2]) - t0)
      }
    }
  }
}
write5(do.call(rbind, out), paste0("generating_", n))
