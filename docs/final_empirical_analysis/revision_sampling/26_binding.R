# 26_binding.R: which support criterion limits availability? For replicate 1 of series S (all sites x orders x environments x samplers), the maximum number of reference configurations that satisfy all criteria at every position at one spread value (max_common) is recomputed with
# each criterion REMOVED in turn. A criterion is "binding" when its removal alone lifts max_common to >= 10. Diagnostic only; the main results never use a relaxed criterion.
source("revision_sampling/00_sim_common.R"); source("revision_sampling/01b_support_diag.R"); NS <- c(100L, 250L, 500L, 1000L, 2500L, 5000L)
U <- expand.grid(n = NS, k = KS, env = ENVS, stringsAsFactors = FALSE)
one <- function(i) { n <- U$n[i]; dist <- U$env[i]; k <- U$k[i]; rep <- 1L; sd0 <- 1e6 * rep + 1e3 * match(dist, ENVS) + as.integer(log10(n) * 100); set.seed(sd0); tr <- gen_env(n, dist); out <- list()
  for (samp in SAMPS) { set.seed(sd0 + 100 * k + match(samp, SAMPS)); sm <- if (samp == "random") list(idx = rand_sets(n, k, 4000L)) else bal_sets(tr$E$TP, k, 4000L); if (nrow(sm$idx) < 200L) next
    d <- tryCatch(prepare_diag(tr$E, tr$xy, sm$idx, k), error = function(e) NULL); if (is.null(d)) next
    cnt <- function(use) { m <- 0L; for (j in seq_along(GRID)) { z <- rep(TRUE, d$nref); for (p in seq_along(POS)) { zz <- rep(TRUE, d$nref); if ("other" %in% use) zz <- zz & d$other_ok; if ("inside" %in% use) zz <- zz & d$CR$inside[, j, p]; if ("cell" %in% use) zz <- zz & d$CR$cell[, j, p]; if ("nn" %in% use) zz <- zz & d$CR$nn[, j, p]; z <- z & zz }; m <- max(m, sum(z)) }; m }
    A <- c("other", "inside", "cell", "nn"); base <- cnt(A); r <- data.frame(n_sites = n, order = k, env = dist, sampler = samp, max_common_all = base, without_nuisance_cells = cnt(setdiff(A, "other")), without_range = cnt(setdiff(A, "inside")), without_TP_cell = cnt(setdiff(A, "cell")), without_nn = cnt(setdiff(A, "nn")),
      pass_nuisance_cells = mean(d$other_ok), pass_range_mid = mean(d$CR$inside[, 11, 1] & d$CR$inside[, 11, 3]), pass_TP_cell_mid = mean(d$CR$cell[, 11, 1] & d$CR$cell[, 11, 3]), pass_nn_mid = mean(d$CR$nn[, 11, 1] & d$CR$nn[, 11, 3])); out[[samp]] <- r }
  do.call(rbind, out) }
R <- do.call(rbind, mclapply(seq_len(nrow(U)), one, mc.cores = 6)); write.csv(R, "revision_sampling/tables/sim_support_binding_criteria.csv", row.names = FALSE); print(dim(R))
