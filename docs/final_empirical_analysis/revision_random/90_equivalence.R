# 90_equivalence.R: numerical implementation validation of zetaGDUM against existing implementations (NOT a predictive validation; no cross-validation; no random site splits).
# EPA macroinvertebrates (BM) and zooplankton (ZP), shared pool, log10 TP, TN, Chl-a; identical sampled site combinations within each comparison.
#   E1  Additive MS-GDUM at k = 2 (ugdm_fit: spread I-splines + linear position + geographic-distance I-splines, logit link, free intercept) vs original GDUM (gdmmTMB::gdmm, binomial/logit, mono = TRUE).
#       Same likelihood and link: coefficient agreement is tested directly. Reference response: numerator = y * trials, denominator = trials.
#   E2  Spread-only MS-GDUM at k = 2 vs GDM (gdm::gdm; negative-exponential link mu = 1 - exp(-eta), weighted deviance). Environmental predictors only (see notes).
#   E3  Spread-only MS-GDUM at k > 2 (k = 3, 5) vs I-spline MS-GDM (zetadiv::Zeta.msgdm; binomial(log) of Simpson-normalised zeta, sign-constrained, set-wise mean pairwise I-spline differences).
# Two separate questions are reported for E2 and E3: (i) design agreement (the set-level I-spline spread columns built by zetaGDUM's kernel from the reference model's own basis vs the reference's design) and
# (ii) coefficient agreement under the reference's link and estimation criterion (temporary validation-specific fit: reference likelihood on the zetaGDUM-kernel design). zetaGDUM's native fits (logit link, own basis) are
# also compared with the reference but are labelled UNMATCHED: different link/likelihood/basis, so no coefficient equivalence is claimed for them. The default package behaviour is not modified.
# Usage (final/, ZG_LIB6 set): Rscript revision_random/90_equivalence.R
source("common.R"); source("lib_ugdm.R"); suppressMessages({library(gdmmTMB); library(gdm); library(zetadiv)})
# splines2 is deliberately NOT attached: it would mask gdmmTMB::isp (degree 3) inside the model formulas
OUT <- "revision_random/tables"; POOL <- "shared"; M <- 3000L; SEED <- 4242L
Ii <- zetadiv:::.Ii; kern <- function(B, z, idx) zetaGDUM:::mean_pair_isp_cpp(B, z, idx)
ramsay <- function(x, kv) { xs <- (x - kv[1]) / (kv[3] - kv[1]); ms <- (kv[2] - kv[1]) / (kv[3] - kv[1]); ts <- c(0, 0, ms, 1, 1); sapply(1:3, function(i) vapply(xs, function(u) if (u >= 1) 1 else if (u <= 0) 0 else Ii(i, 2, u, ts), 0)) }
cor0 <- function(a, b) if (length(a) < 3 || sd(a) == 0 || sd(b) == 0) NA_real_ else cor(a, b)
CO <- list(); CU <- list(); DS <- list(); SM <- list()
add_coef <- function(cmp, tx, k, cls, pred, idx, term, ref, zet, matched) CO[[length(CO) + 1]] <<- data.frame(comparison = cmp, taxon = tx, order = k, term_class = cls, predictor = pred, index = idx, term = term, reference = ref, zetaGDUM = zet, difference = zet - ref, matched_link_and_criterion = matched)
add_curve <- function(cmp, tx, k, pred, x, fref, fzet, matched) CU[[length(CU) + 1]] <<- data.frame(comparison = cmp, taxon = tx, order = k, predictor = pred, x_log10 = x, f_reference = fref, f_zetaGDUM = fzet, matched_link_and_criterion = matched)
summ <- function(cmp, tx, k, matched, claim, design_diff = NA, ...) { D <- do.call(rbind, CO[vapply(CO, function(z) z$comparison[1] == cmp && z$taxon[1] == tx && z$order[1] == k, logical(1))]); sp <- D$term_class == "spline"
  Hs <- D[sp, ] |> dplyr::group_by(predictor) |> dplyr::summarise(r = sum(reference), z = sum(zetaGDUM), .groups = "drop"); Cc <- do.call(rbind, CU[vapply(CU, function(z) z$comparison[1] == cmp && z$taxon[1] == tx && z$order[1] == k, logical(1))])
  SM[[length(SM) + 1]] <<- data.frame(comparison = cmp, taxon = tx, order = k, matched_link_and_criterion = matched, n_terms = nrow(D), max_abs_diff_all = max(abs(D$difference)), max_abs_diff_intercept = max(abs(D$difference[D$term_class == "intercept"])), max_abs_diff_spline = max(abs(D$difference[sp])),
    max_abs_diff_spline_sum = max(abs(Hs$z - Hs$r)), max_abs_diff_position = if (any(D$term_class == "position")) max(abs(D$difference[D$term_class == "position"])) else NA, max_abs_diff_distance = if (any(D$term_class == "distance")) max(abs(D$difference[D$term_class == "distance"])) else NA,
    cor_all_coefficients = cor0(D$reference, D$zetaGDUM), cor_spline_coefficients = cor0(D$reference[sp], D$zetaGDUM[sp]), cor_spline_sums = cor0(Hs$r, Hs$z), design_max_abs_diff = design_diff, curve_max_abs_diff = max(abs(Cc$f_reference - Cc$f_zetaGDUM)),
    curve_shape_cor_min = min(tapply(seq_len(nrow(Cc)), Cc$predictor, function(i) cor0(Cc$f_reference[i], Cc$f_zetaGDUM[i]))), ..., claim = claim, stringsAsFactors = FALSE) }
for (cc in c("BM", "ZP")) { d <- POOLS[[POOL]][[cc]]; e <- get_env(d$env); xy <- as.matrix(d$env[, c("lon", "lat")]); DM <- hav(xy); n <- nrow(e); vars <- names(e); tx <- cc; grid01 <- seq(0, 1, length.out = 101)
  # ============ E1: additive MS-GDUM (k = 2) vs gdmmTMB ============
  il <- sample_sets(n, 2, M, seed = SEED); r <- ugdm_response(d$x, il, "total"); k <- r$keep; y <- r$y[k]; den <- r$trials[k]
  dsg <- ugdm_design(e, il, dist_matrix = DM); f <- ugdm_fit(dsg, y, den); Z <- as.data.frame(scale(e, center = dsg$environment_basis$center, scale = dsg$environment_basis$scale))
  Xp <- data.frame(dist01 = pmin(pmax((DM[cbind(il[, 1], il[, 2])] - dsg$dist_range[1]) / diff(dsg$dist_range), 0), 1))
  m <- suppressWarnings(gdmm(Y_diss = y * den, Y_den = den, D = il[k, ], X = Z, X_pair = Xp[k, , drop = FALSE], diss_formula = ~ gdmmTMB::isp(TP) + gdmmTMB::isp(TN) + gdmmTMB::isp(Chla), uniq_formula = ~ TP + TN + Chla, pair_formula = ~ gdmmTMB::isp(dist01), mono = TRUE, mono_pair = TRUE, family = "binomial", link = "logit", control = list(eval.max = 20000, iter.max = 20000, rel.tol = 1e-14)))
  rp <- m$obj$report(); cmp <- "E1 additive MS-GDUM k=2 vs gdmmTMB (matched)"
  add_coef(cmp, tx, 2, "intercept", NA, 0, "intercept", rp$intercept, f$intercept, TRUE)
  for (vi in seq_along(vars)) for (j in 1:3) add_coef(cmp, tx, 2, "spline", vars[vi], j, paste0(vars[vi], "_I", j), rp$e_beta[3 * (vi - 1) + j], f$beta[3 * (vi - 1) + j], TRUE)
  for (j in 1:3) add_coef(cmp, tx, 2, "distance", "geographic distance", j, paste0("distance_I", j), rp$e_beta_p[j], f$beta[9 + j], TRUE)
  for (vi in seq_along(vars)) add_coef(cmp, tx, 2, "position", vars[vi], 1, paste0(vars[vi], "_position"), rp$lambda[vi], unname(f$lambda[vi]), TRUE)
  pl <- as.numeric(predict(m, CI = FALSE, type = "link")); eta <- f$intercept + drop(cbind(dsg$Xdir, dsg$Xdist)[k, ] %*% f$beta) + drop(dsg$Xpos[k, ] %*% f$lambda)
  DS[[length(DS) + 1]] <- data.frame(comparison = cmp, taxon = tx, order = 2, quantity = "linear predictor, all pairs", max_abs_diff = max(abs(pl - eta)), note = "gdmmTMB predict(type = 'link') vs zetaGDUM eta")
  # design agreement: gdmmTMB isp(x) basis (degree 2, no interior knot) vs zetaGDUM basis, pair differences
  des <- 0; for (vi in seq_along(vars)) { B1 <- unclass(gdmmTMB::isp(Z[[vars[vi]]])); B2 <- unclass(dsg$environment_basis$splines[[vi]])[, 1:3]; des <- max(des, max(abs(B1[, 1:3] - B2))) }
  DS[[length(DS) + 1]] <- data.frame(comparison = cmp, taxon = tx, order = 2, quantity = "I-spline basis (site level), gdmmTMB::isp vs zetaGDUM", max_abs_diff = des, note = "")
  for (vi in seq_along(vars)) { v <- vars[vi]; zg <- seq(min(Z[[v]]), max(Z[[v]]), length.out = 101); sp1 <- gdmmTMB::isp(Z[[v]]); sp2 <- dsg$environment_basis$splines[[vi]]; Bg <- unclass(predict(sp1, zg))[, 1:3]; Bz <- unclass(predict(sp2, zg))[, 1:3]
    add_curve(cmp, tx, 2, v, zg * dsg$environment_basis$scale[[v]] + dsg$environment_basis$center[[v]], drop(Bg %*% rp$e_beta[3 * (vi - 1) + 1:3]), drop(Bz %*% f$beta[3 * (vi - 1) + 1:3]), TRUE) }
  summ(cmp, tx, 2, TRUE, "Coefficient equivalence supported (same logit likelihood, same basis, same pairs)", des, converged_reference = m$opt$message, D2_zetaGDUM = f$D2)
  # ============ E2: spread-only MS-GDUM (k = 2) vs gdm ============
  il2 <- il[k, ]; tab <- data.frame(distance = y, weights = den, s1.xCoord = 0, s1.yCoord = 0, s2.xCoord = 0, s2.yCoord = 0); for (v in vars) tab[[paste0("s1.", v)]] <- e[il2[, 1], v]; for (v in vars) tab[[paste0("s2.", v)]] <- e[il2[, 2], v]; class(tab) <- c("gdmData", "data.frame")
  g <- suppressWarnings(gdm(tab, geo = FALSE)); kn <- g$knots; cg <- g$coefficients
  Xg <- do.call(cbind, lapply(seq_along(vars), function(vi) { B <- ramsay(e[[vars[vi]]], kn[3 * (vi - 1) + 1:3]); kern(B, e[[vars[vi]]], il2) })); colnames(Xg) <- paste0(rep(vars, each = 3), "_f", 1:3)
  wdev <- function(mu) { mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12); 2 * sum(den * (ifelse(y > 0, y * log(y / mu), 0) + ifelse(y < 1, (1 - y) * log((1 - y) / (1 - mu)), 0))) }
  obj <- function(p) wdev(1 - exp(-(p[1] + drop(Xg %*% p[-1])))) / sum(den); best <- NULL
  for (st in list(c(.5, rep(.05, 9)), c(1, rep(.2, 9)), c(.2, rep(.5, 9)))) {  # independent starts (gdm's own solution is NOT used as a start)
    o <- optim(st, obj, method = "L-BFGS-B", lower = c(1e-8, rep(0, 9)), control = list(maxit = 20000, factr = 10, pgtol = 0)); o <- optim(o$par, obj, method = "L-BFGS-B", lower = c(1e-8, rep(0, 9)), control = list(maxit = 20000, factr = 10, pgtol = 0)); if (is.null(best) || o$value < best$value) best <- o }
  cmp2 <- "E2 spread-only MS-GDUM k=2 vs gdm (matched: negative-exponential link, gdm deviance, gdm knots/basis)"; add_coef(cmp2, tx, 2, "intercept", NA, 0, "intercept", g$intercept, best$par[1], TRUE)
  for (vi in seq_along(vars)) for (j in 1:3) add_coef(cmp2, tx, 2, "spline", vars[vi], j, paste0(vars[vi], "_I", j), cg[3 * (vi - 1) + j], best$par[1 + 3 * (vi - 1) + j], TRUE)
  dev_gdm_params <- wdev(1 - exp(-(g$intercept + drop(Xg %*% cg)))); DS[[length(DS) + 1]] <- data.frame(comparison = cmp2, taxon = tx, order = 2, quantity = "weighted deviance: gdm reported / gdm coefficients on zetaGDUM-kernel design / validation-routine optimum", max_abs_diff = abs(dev_gdm_params - g$gdmdeviance), note = sprintf("gdm %.4f; gdm coefficients on kernel design %.4f; validation optimum %.4f", g$gdmdeviance, dev_gdm_params, best$value * sum(den)))
  # basis agreement: gdm fitted curve vs Ramsay reconstruction with gdm coefficients
  ie <- gdm::isplineExtract(g); bd <- 0
  for (vi in seq_along(vars)) { v <- vars[vi]; xg <- ie$x[, v]; fg <- ie$y[, v]; fr <- drop(ramsay(xg, kn[3 * (vi - 1) + 1:3]) %*% cg[3 * (vi - 1) + 1:3]); bd <- max(bd, max(abs(fg - fr)))
    fz <- drop(ramsay(xg, kn[3 * (vi - 1) + 1:3]) %*% best$par[1 + 3 * (vi - 1) + 1:3]); add_curve(cmp2, tx, 2, v, xg, fg, fz, TRUE) }
  DS[[length(DS) + 1]] <- data.frame(comparison = cmp2, taxon = tx, order = 2, quantity = "gdm fitted I-spline curve vs standard (Ramsay) I-spline reconstruction with gdm knots and coefficients", max_abs_diff = bd, note = "isplineExtract() plotting curve (evaluated by gdm on its own grid) vs standard reconstruction; the fitted design itself is reproduced exactly (see linear-predictor row)")
  des2 <- max(abs(g$intercept + drop(Xg %*% cg) - g$ecological)); DS[[length(DS) + 1]] <- data.frame(comparison = cmp2, taxon = tx, order = 2, quantity = "gdm linear predictor (ecological) vs intercept + zetaGDUM-kernel design on gdm's own basis (standard Ramsay I-splines, gdm knots) x gdm coefficients", max_abs_diff = des2, note = "design identity")
  summ(cmp2, tx, 2, TRUE, "Design identity and coefficient agreement under matched link/criterion (negative-exponential, gdm deviance); environmental predictors only", des2, deviance_gdm = g$gdmdeviance, deviance_validation_optimum = best$value * sum(den))
  # native (unmatched): zetaGDUM spread-only NNLS (logit, own basis, no interior knot) vs gdm
  dsg0 <- ugdm_design(e, il, dist_matrix = NULL); nat <- zetaGDUM:::.fit_logit_nnls(-dsg0$Xdir[k, ], y, den); bn <- -nat$coef; cmp2n <- "E2n spread-only MS-GDUM k=2 native fit vs gdm (UNMATCHED: logit link, own basis)"
  add_coef(cmp2n, tx, 2, "intercept", NA, 0, "intercept", g$intercept, nat$intercept, FALSE); for (vi in seq_along(vars)) for (j in 1:3) add_coef(cmp2n, tx, 2, "spline", vars[vi], j, paste0(vars[vi], "_I", j), cg[3 * (vi - 1) + j], bn[3 * (vi - 1) + j], FALSE)
  for (vi in seq_along(vars)) { v <- vars[vi]; xg <- ie$x[, v]; zg <- (xg - dsg0$environment_basis$center[[v]]) / dsg0$environment_basis$scale[[v]]; bk <- attr(dsg0$environment_basis$splines[[vi]], "Boundary.knots"); Bz <- unclass(predict(dsg0$environment_basis$splines[[vi]], pmin(pmax(zg, bk[1]), bk[2])))[, 1:3]; add_curve(cmp2n, tx, 2, v, xg, ie$y[, v], drop(Bz %*% bn[3 * (vi - 1) + 1:3]), FALSE) }
  summ(cmp2n, tx, 2, FALSE, "No coefficient equivalence claimed (different link, likelihood and basis)", NA)
  # ============ E3: spread-only MS-GDUM (k > 2) vs zetadiv ============
  pa <- as.data.frame((d$x > 0) * 1L); env_df <- as.data.frame(e)
  for (kk in c(3L, 5L)) { set.seed(SEED + kk); zm <- suppressWarnings(Zeta.msgdm(pa, env_df, xy = NULL, order = kk, sam = M, reg.type = "ispline", normalize = "Simpson", rescale.pred = FALSE, silent = TRUE, empty.row = 0)); set.seed(SEED + kk); idx <- t(replicate(M, sample(1:n, kk)))
    rt <- ugdm_response(d$x, idx, "turnover"); sim <- 1 - rt$y; resp_diff <- max(abs(sim[rt$keep] - zm$val[rt$keep]))
    Xz <- do.call(cbind, lapply(vars, function(v) { B <- ramsay(e[[v]], c(min(e[[v]]), median(e[[v]]), max(e[[v]]))); kern(B, e[[v]], idx) })); Pz <- as.matrix(zm$predictors[, 1:9]); dd <- max(abs(Xz - Pz))
    cmp3 <- sprintf("E3 spread-only MS-GDUM k=%d vs zetadiv (matched: binomial(log), Simpson zeta, sign constraints)", kk); cz <- coef(zm$model); dfz <- data.frame(zeta.val = zm$val, Xz); names(dfz)[-1] <- names(zm$predictors)[1:9]
    mod <- zetadiv::glm.cons(zeta.val ~ ., data = dfz, family = binomial("log"), method = "glm.fit.cons", cons = -1, cons.inter = -1, control = list()); cm <- coef(mod)
    add_coef(cmp3, tx, kk, "intercept", NA, 0, "intercept", cz[[1]], cm[[1]], TRUE); for (vi in seq_along(vars)) for (j in 1:3) add_coef(cmp3, tx, kk, "spline", vars[vi], j, paste0(vars[vi], "_I", j), cz[[1 + 3 * (vi - 1) + j]], cm[[1 + 3 * (vi - 1) + j]], TRUE)
    DS[[length(DS) + 1]] <- data.frame(comparison = cmp3, taxon = tx, order = kk, quantity = "zetadiv response (Simpson zeta) vs 1 - zetaGDUM turnover response, replicated sets", max_abs_diff = resp_diff, note = "sets reproduced by replaying zetadiv's sample() sequence")
    DS[[length(DS) + 1]] <- data.frame(comparison = cmp3, taxon = tx, order = kk, quantity = "set-level spread design: zetaGDUM kernel on zetadiv basis vs zetadiv predictors", max_abs_diff = dd, note = "rescale.pred = FALSE")
    for (vi in seq_along(vars)) { v <- vars[vi]; xg <- seq(min(e[[v]]), max(e[[v]]), length.out = 101); Bg <- ramsay(xg, c(min(e[[v]]), median(e[[v]]), max(e[[v]]))); add_curve(cmp3, tx, kk, v, xg, drop(Bg %*% (-cz[1 + 3 * (vi - 1) + 1:3])), drop(Bg %*% (-cm[1 + 3 * (vi - 1) + 1:3])), TRUE) }
    summ(cmp3, tx, kk, TRUE, "Design identical to zetadiv up to tolerance; coefficients reproduced under zetadiv's own likelihood", dd, response_max_abs_diff = resp_diff)
    # native (unmatched): ugdm design (own basis), logit, turnover (Simpson-type) response
    dk <- ugdm_design(e, idx); nat <- zetaGDUM:::.fit_logit_nnls(-dk$Xdir[rt$keep, ], rt$y[rt$keep], rt$trials[rt$keep]); bn <- -nat$coef; cmp3n <- sprintf("E3n spread-only MS-GDUM k=%d native fit vs zetadiv (UNMATCHED: logit/dissimilarity vs log/similarity, own basis)", kk)
    add_coef(cmp3n, tx, kk, "intercept", NA, 0, "intercept", -cz[[1]], nat$intercept, FALSE); for (vi in seq_along(vars)) for (j in 1:3) add_coef(cmp3n, tx, kk, "spline", vars[vi], j, paste0(vars[vi], "_I", j), -cz[[1 + 3 * (vi - 1) + j]], bn[3 * (vi - 1) + j], FALSE)
    for (vi in seq_along(vars)) { v <- vars[vi]; xg <- seq(min(e[[v]]), max(e[[v]]), length.out = 101); Bg <- ramsay(xg, c(min(e[[v]]), median(e[[v]]), max(e[[v]]))); zg <- (xg - dk$environment_basis$center[[v]]) / dk$environment_basis$scale[[v]]; Bn <- unclass(predict(dk$environment_basis$splines[[vi]], zg))[, 1:3]
      add_curve(cmp3n, tx, kk, v, xg, drop(Bg %*% (-cz[1 + 3 * (vi - 1) + 1:3])), drop(Bn %*% bn[3 * (vi - 1) + 1:3]), FALSE) }
    summ(cmp3n, tx, kk, FALSE, "No coefficient equivalence claimed (different link, response orientation and basis)", NA) } }
Cf <- do.call(rbind, CO); Cu <- do.call(rbind, CU); Dsg <- do.call(rbind, DS); Sm <- dplyr::bind_rows(SM)
write.csv(Cf, file.path(OUT, "eq_coefficients.csv"), row.names = FALSE); write.csv(Cu, file.path(OUT, "eq_fitted_curves.csv"), row.names = FALSE); write.csv(Dsg, file.path(OUT, "eq_design_agreement.csv"), row.names = FALSE); write.csv(Sm, file.path(OUT, "eq_validation_table.csv"), row.names = FALSE)
print(as.data.frame(Sm[, c("comparison", "taxon", "order", "max_abs_diff_all", "max_abs_diff_intercept", "max_abs_diff_spline", "cor_all_coefficients", "design_max_abs_diff", "curve_max_abs_diff")]), digits = 3); print(as.data.frame(Dsg[, c("comparison", "taxon", "order", "max_abs_diff")]), digits = 3)
