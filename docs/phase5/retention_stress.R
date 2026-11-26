source("docs/phase5/common.R")
n <- 250L
group <- ceiling(seq_len(n) / 4)
x <- cbind(core = 1, 1e7 * outer(group, seq_len(max(group)), "=="))
x[229:n,1] <- 0
x[229:n,-1] <- x[229:n,-1] * (1 + 1e-7)
env <- data.frame(gradient = seq_len(n) / n)
orders <- c(4L, 13L, 30L, 50L, n)
exact <- zeta_order_support(x, orders)
base <- ugdm_path_data(x, env, sample_sets_coupled(x, orders, 3500, seed = 482))
rows <- lapply(seq_along(orders), function(j) {
  a <- model_data(base[[j]], "relative", "B")
  f <- fit5(a, "B")
  d <- ugdm_variant_diagnostics(f, a$design, a$y, sampling_weights=a$sampling_weights)
  cbind(exact[j,], response_sd=sd(a$y), converged=f$converged, D2=f$D2,
        scaled_condition=d$condition, rank=d$rank, parameters=d$parameters)
})
write5(do.call(rbind,rows),"retention_stress")
