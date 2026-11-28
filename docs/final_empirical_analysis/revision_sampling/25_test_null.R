# 25_test_null.R: variance of the held-out response (test sites/combinations are generated from seeds, so they are reproduced exactly) -> skill = 1 - MSE / Var(test response), a null that does not depend on the training sampler.
# (Deviance explained against the TRAINING mean is not comparable between samplers: balanced training sets have a different mean total dissimilarity than the random test combinations.)
source("revision_sampling/11_known_rep.R")   # also defines run_meta helpers via 00_sim_common.R
ks <- unique(read.csv("revision_sampling/tables/sim_meta_sup.csv")[, c("n_sites", "env", "rep", "order")])
one <- function(i, sim) { n <- ks$n_sites[i]; dist <- ks$env[i]; rep <- ks$rep[i]; k <- ks$order[i]
  if (sim == "meta") { sd0 <- 1e6 * rep + 1e3 * match(dist, ENVS) + as.integer(log10(n) * 100); set.seed(sd0); gen_env(n, dist); set.seed(sd0 + 7); te <- gen_env(400L, dist); cte <- sim_counts(te$E$TP); idx_te <- rand_sets(400L, k, 1500L); y <- unclass(ugdm_component_response(cte, idx_te)) }
  else { sd0 <- 1e6 * rep + 1e3 * match(dist, ENVS) + as.integer(log10(n) * 100) + 3e8; set.seed(sd0); tr <- gen_env(n, dist); uT <- rnorm(n, 0, U_SD); uS <- rnorm(n, 0, U_SD); set.seed(sd0 + 7); te <- gen_env(400L, dist); uTe <- rnorm(400, 0, U_SD); uSe <- rnorm(400, 0, U_SD); idx_te <- rand_sets(400L, k, 1500L)
    y <- unclass(gen_resp(te$E, te$xy, idx_te, k, DDG[[as.character(k)]], KP[[as.character(k)]], uTe, uSe)$y) }
  data.frame(ks[i, ], test_mean_total = mean(y[, "total"]), test_var_total = mean((y[, "total"] - mean(y[, "total"]))^2), test_var_occupancy = mean((y[, "occupancy"] - mean(y[, "occupancy"]))^2), test_var_heterogeneity = mean((y[, "heterogeneity"] - mean(y[, "heterogeneity"]))^2)) }
for (sim in c("meta", "known")) write.csv(do.call(rbind, mclapply(seq_len(nrow(ks)), one, sim = sim, mc.cores = 6)), sprintf("revision_sampling/tables/sim_%s_testnull.csv", sim), row.names = FALSE)
