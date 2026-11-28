# 11_known_rep.R: known-parameter replicate unit. Responses are generated from an additive MS-GDUM with KNOWN link-scale coefficients (defined on a fixed population reference basis, so the truth is the same function of TP position/spread in every replicate),
# with binomial sampling noise (300 trials) and a site-level random effect (sd 0.25 on both the total and the heterogeneity-share logit; a combination's effect = mean of its sites' effects), so that sites, not only combinations, carry independent information.
# Truth contrasts = the same matched-cell prediction machinery applied to the TRUE coefficients on the replicate's own supported reference configurations. The fitted models never see the truth.
source("revision_sampling/00_sim_common.R")
set.seed(31); POP <- gen_env(4000L, "uniform"); DDG <- list(); KP <- list(); TRIALS <- 300L; U_SD <- .25
TARGET_TOTAL <- c(`2` = .45, `3` = .60, `5` = .75, `10` = .88, `20` = .95)
for (k in KS) { set.seed(k); i0 <- rand_sets(4000L, k, 2500L); DDG[[as.character(k)]] <- ugdm_component_design(POP$E, i0, "additive", coordinates = POP$xy, lonlat = TRUE); X <- layoutfn(DDG[[as.character(k)]], TRUE)$X; nm <- colnames(X)
  bT <- setNames(numeric(length(nm)), nm); bS <- bT; bT[c("TP_f1", "TP_f2", "TP_f3", "TN_f1", "TN_f2", "TN_f3", "Chla_f1", "Chla_f2", "Chla_f3")] <- c(1.2, .8, .4, .3, .2, .1, .2, .1, .05); bT[c("TP", "TP_quadratic", "TN", "Chla")] <- c(.45, -.2, .15, .10)
  bS[c("TP_f1", "TP_f2", "TP_f3")] <- c(-1, -.6, -.3); bS[c("TP", "TP_quadratic")] <- c(-.35, .1)
  oth <- drop(X[, -1, drop = FALSE] %*% bT[-1]); bT[1] <- uniroot(function(b) mean(plogis(b + oth)) - TARGET_TOTAL[[as.character(k)]], c(-12, 12))$root
  oS <- drop(X[, -1, drop = FALSE] %*% bS[-1]); bS[1] <- uniroot(function(b) mean(plogis(b + oS)) - .35, c(-12, 12))$root; KP[[as.character(k)]] <- list(bT = bT, bS = bS) }
cells_from_coef <- function(p, dd_ref, model, bT, bS, ids, svals) { G <- list(); for (pn in names(PCOL)) for (sv in svals) G[[paste(model, pn, sv)]] <- cell_parts(p, dd_ref, model, ids, log10(PCOL[[pn]]), sv, bT, bS); G }
gen_resp <- function(E, xy, idx, k, ddg, kp, ubar_sites_T, ubar_sites_S) { d <- ugdm_component_design(E, idx, "additive", coordinates = xy, lonlat = TRUE, reference = ddg); X <- layoutfn(d, TRUE)$X; eT <- drop(X %*% kp$bT); eS <- drop(X %*% kp$bS)
  mT <- rowMeans(matrix(ubar_sites_T[idx], nrow(idx))); mS <- rowMeans(matrix(ubar_sites_S[idx], nrow(idx))); yT <- rbinom(nrow(idx), TRIALS, plogis(eT + mT)) / TRIALS; yS <- rbinom(nrow(idx), TRIALS, plogis(eS + mS)) / TRIALS; H <- yT * yS
  list(y = structure(cbind(total = yT, occupancy = yT - H, heterogeneity = H), sets = d$idx, class = c("ugdm_component_response", "matrix", "array")), mu_T = plogis(eT), mu_S = plogis(eS)) }
run_known <- function(n, dist, rep, k, C, series) {
  sd0 <- 1e6 * rep + 1e3 * match(dist, ENVS) + as.integer(log10(n) * 100) + 3e8; set.seed(sd0); tr <- gen_env(n, dist); uT <- rnorm(n, 0, U_SD); uS <- rnorm(n, 0, U_SD)
  set.seed(sd0 + 7); te <- gen_env(400L, dist); uTe <- rnorm(400, 0, U_SD); uSe <- rnorm(400, 0, U_SD); idx_te <- rand_sets(400L, k, 1500L); ddg <- DDG[[as.character(k)]]; kp <- KP[[as.character(k)]]
  yte <- gen_resp(te$E, te$xy, idx_te, k, ddg, kp, uTe, uSe); tag <- data.frame(series = series, n_sites = n, env = dist, rep = rep, order = k, n_combos_target = C)
  FIT <- list(); SUP <- list(); CON <- list()
  for (samp in SAMPS) { set.seed(sd0 + 100 * k + match(samp, SAMPS)); t0 <- Sys.time()
    sm <- if (samp == "random") list(idx = rand_sets(n, k, C), tries = NA_integer_) else bal_sets(tr$E$TP, k, C); idx <- sm$idx
    row <- cbind(tag, sampler = samp, n_combos = nrow(idx), n_sites_used = length(unique(as.vector(idx))), balanced_tries = sm$tries)
    if (nrow(idx) < 200L) { SUP[[samp]] <- cbind(row, data.frame(supported = FALSE, reason = "fewer than 200 feasible balanced combinations")); next }
    res <- tryCatch({ p <- prepare_sim(tr$E, tr$xy, idx, k); g <- gen_resp(tr$E, tr$xy, p$dd$idx, k, ddg, kp, uT, uS); y <- g$y
      dd <- list(additive = p$dd, spread_only = ugdm_component_design(tr$E, p$dd$idx, "spread_only", coordinates = tr$xy, lonlat = TRUE)); fits <- lapply(dd, function(z) ugdm_fit_components(z, y, "coherent")); ybar <- colSums(unclass(y) * (fits$additive$weights / sum(fits$additive$weights)))
      orc <- dev_frac(unclass(yte$y)[, "total"], yte$mu_T, ybar[["total"]])
      fr <- do.call(rbind, lapply(names(fits), function(m) { dte <- ugdm_component_design(te$E, idx_te, m, coordinates = te$xy, lonlat = TRUE, reference = dd[[m]]); cbind(row, model = m, converged = fits[[m]]$converged, as.data.frame(t(heldout(fits[[m]], dte, yte$y, ybar))), D2_total_oracle = orc) }))
      sp <- cbind(row, tab_support(p)); sp$reason <- if (sp$supported) "supported" else if (sp$max_common_profiles < MINREF) "fewer than 10 common profiles at any spread" else "no run of >= 3 grid points and >= 0.20 log10 width with >= 10 common profiles"; cr <- NULL
      if (sp$supported) { ids <- select_ids(p$valid[[sp$interval]], MAXPROF); sv <- GRID[parse_iv(sp$interval)[1]:parse_iv(sp$interval)[2]]; m0 <- ceiling(length(sv) / 2); sv <- sv[c(m0 - 2L, m0, m0 + 2L)]
        tt <- contrasts_for(cells_from_coef(p, ddg, "additive", kp$bT, kp$bS, ids, sv), "additive", sv); tt <- tt[, c("contrast", "dT", "dO", "dH", "dEtaT", "dEtaS")]; names(tt)[-1] <- paste0("true_", names(tt)[-1])
        est <- do.call(rbind, lapply(names(fits), function(m) contrasts_for(predict_cells(p, dd[[m]], fits[[m]], m, ids, sv), m, sv))); cr <- merge(cbind(row, est), tt, by = "contrast"); cr$n_ids <- length(ids); cr$s_lo <- sv[1]; cr$s_hi <- sv[3]; cr$s_mid <- sv[2] }
      list(fr = fr, sp = sp, cr = cr) }, error = function(e) { message("ERR ", series, " n", n, " ", dist, " r", rep, " k", k, " ", samp, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) { SUP[[samp]] <- cbind(row, data.frame(supported = FALSE, reason = "fit or support error")); next }
    FIT[[samp]] <- res$fr; SUP[[samp]] <- res$sp; if (!is.null(res$cr)) CON[[samp]] <- res$cr; SUP[[samp]]$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs")) }
  list(fit = bind_rows(FIT), sup = bind_rows(SUP), con = bind_rows(CON)) }
