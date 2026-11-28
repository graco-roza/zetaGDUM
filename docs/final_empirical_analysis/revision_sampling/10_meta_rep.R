# 10_meta_rep.R: one metacommunity replicate unit = one site draw (n sites, environmental distribution, replicate) at one zeta order, evaluated for BOTH combination samplers on the SAME sites and communities (paired design).
# Truth = independent Monte Carlo of the package response functions on freshly simulated communities for the matched configurations of the replicate (not the fitted model). Unsupported contrasts are recorded, never imputed.
source("revision_sampling/00_sim_common.R")
run_meta <- function(n, dist, rep, k, C, series, nreps_truth = 25L) {
  sd0 <- 1e6 * rep + 1e3 * match(dist, ENVS) + as.integer(log10(n) * 100); set.seed(sd0); tr <- gen_env(n, dist); cnt <- sim_counts(tr$E$TP)
  set.seed(sd0 + 7); te <- gen_env(400L, dist); cte <- sim_counts(te$E$TP); idx_te <- rand_sets(400L, k, 1500L); yte <- ugdm_component_response(cte, idx_te)
  tag <- data.frame(series = series, n_sites = n, env = dist, rep = rep, order = k, n_combos_target = C)
  FIT <- list(); SUP <- list(); CON <- list()
  for (samp in SAMPS) { set.seed(sd0 + 100 * k + match(samp, SAMPS)); t0 <- Sys.time()
    sm <- if (samp == "random") list(idx = rand_sets(n, k, C), tries = NA_integer_) else bal_sets(tr$E$TP, k, C)
    idx <- sm$idx; row <- cbind(tag, sampler = samp, n_combos = nrow(idx), n_sites_used = length(unique(as.vector(idx))), balanced_tries = sm$tries)
    if (nrow(idx) < 200L) { SUP[[samp]] <- cbind(row, data.frame(supported = FALSE, reason = "fewer than 200 feasible balanced combinations")); next }
    res <- tryCatch({
      y <- ugdm_component_response(cnt, idx); p <- prepare_sim(tr$E, tr$xy, idx, k); y <- ugdm_component_response(cnt, p$dd$idx)
      dd <- list(additive = p$dd, spread_only = ugdm_component_design(tr$E, p$dd$idx, "spread_only", coordinates = tr$xy, lonlat = TRUE)); fits <- lapply(dd, function(z) ugdm_fit_components(z, y, "coherent")); ybar <- colSums(unclass(y) * (fits$additive$weights / sum(fits$additive$weights)))
      fr <- do.call(rbind, lapply(names(fits), function(m) { dte <- ugdm_component_design(te$E, idx_te, m, coordinates = te$xy, lonlat = TRUE, reference = dd[[m]]); cbind(row, model = m, converged = fits[[m]]$converged, as.data.frame(t(heldout(fits[[m]], dte, yte, ybar)))) }))
      sp <- cbind(row, tab_support(p)); sp$reason <- if (sp$supported) "supported" else if (sp$max_common_profiles < MINREF) "fewer than 10 common profiles at any spread" else "no run of >= 3 grid points and >= 0.20 log10 width with >= 10 common profiles"
      cr <- NULL
      if (sp$supported) { ids <- select_ids(p$valid[[sp$interval]], MAXPROF); sv <- GRID[parse_iv(sp$interval)[1]:parse_iv(sp$interval)[2]]; m0 <- ceiling(length(sv) / 2); sv <- sv[c(m0 - 2L, m0, m0 + 2L)]; tt <- truth_contrasts(p, ids, sv, nreps_truth)
        est <- do.call(rbind, lapply(names(fits), function(m) contrasts_for(predict_cells(p, dd[[m]], fits[[m]], m, ids, sv), m, sv))); cr <- merge(cbind(row, est), setNames(tt, c("contrast", paste0("true_", c("dT", "dO", "dH")))), by = "contrast"); cr$n_ids <- length(ids); cr$s_lo <- sv[1]; cr$s_hi <- sv[3]; cr$s_mid <- sv[2] }
      list(fr = fr, sp = sp, cr = cr) }, error = function(e) { message("ERR ", series, " n", n, " ", dist, " r", rep, " k", k, " ", samp, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) { SUP[[samp]] <- cbind(row, data.frame(supported = FALSE, reason = "fit or support error")); next }
    FIT[[samp]] <- res$fr; SUP[[samp]] <- res$sp; if (!is.null(res$cr)) CON[[samp]] <- res$cr; SUP[[samp]]$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs")) }
  list(fit = bind_rows(FIT), sup = bind_rows(SUP), con = bind_rows(CON)) }
