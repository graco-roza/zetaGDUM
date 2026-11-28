# 00_fig6_common.R: shared helpers for the redesigned Figure 6 (principal package: additive framework only). EPA shared lake pool (principal), full pools (sensitivity).
# Spread-only MS-GDM and additive MS-GDUM are both fitted with the coherent component architecture (ugdm_fit_components) on IDENTICAL site combinations, response definitions, I-spline bases, link functions, nonnegativity constraints and ridge penalty; they differ only by the TP/TN/Chl-a position terms.
source("common.R"); suppressMessages(library(RANN)); source("revision_fig6/01_support_machinery.R"); library(parallel)
TXL <- c(BM = "EPA macroinvertebrates", ZP = "EPA zooplankton"); ORD <- list(BM = c(2L, 5L), ZP = c(2L, 4L)); POOLSET <- c("shared", "full"); MODELS6 <- c(spread_only = "spread_only", additive = "additive"); MODEL_LAB <- c(spread_only = "Spread-only MS-GDM", additive = "Additive MS-GDUM")
PCOL <- c(Low = 10, Intermediate = 35, High = 100); MAXPROF <- 250L
get_d <- function(pool, tx) { d <- POOLS[[pool]][[tx]]; list(x = d$x, env = d$env) }
bvec <- function(fit, stage) drop(fit$transform %*% fit$fits[[stage]]$coef)
cell_parts <- function(p, dd_ref, model, ids, m_log, s, bT, bS, bounds = NULL) {
  ref <- p$ref[ids, , drop = FALSE]; sh <- p$shape[ids, , drop = FALSE]; k <- p$k; out <- list()
  for (b in split(seq_along(ids), ceiling(seq_along(ids) / 100L))) { ii <- as.vector(t(ref[b, , drop = FALSE])); en <- p$E[ii, , drop = FALSE]; xy <- p$xy[ii, , drop = FALSE]; en$TP <- as.vector(t(m_log + s * sh[b, , drop = FALSE])); idx <- matrix(seq_len(length(ii)), ncol = k, byrow = TRUE)
    zz <- ugdm_component_design(en, idx, model, coordinates = xy, lonlat = TRUE, reference = dd_ref); X <- layoutfn(zz, !is.null(zz$geo))$X; cn <- colnames(X); stopifnot(identical(cn, names(bT)))
    sp <- grep("^TP_f", cn); pos <- match(c("TP", "TP_quadratic", "TP_cubic"), cn); oth <- setdiff(seq_along(cn), c(sp, pos)); part <- function(bb, cols) if (length(cols)) drop(X[, cols, drop = FALSE] %*% bb[cols]) else rep(0, nrow(X))
    out[[length(out) + 1]] <- data.frame(profile = ids[b], oT = part(bT, oth), sT = part(bT, sp), pT = part(bT, pos), oS = part(bS, oth), sS = part(bS, sp), pS = part(bS, pos)) }
  do.call(rbind, out) }
resp_from <- function(eT, eS) { T <- plogis(eT); H <- T * plogis(eS); data.frame(total = T, heterogeneity = H, occupancy = T - H) }
select_ids <- function(ids, n = MAXPROF) if (length(ids) <= n) ids else ids[unique(round(seq(1, length(ids), length.out = n)))]
parse_iv <- function(s) as.integer(strsplit(s, ":", fixed = TRUE)[[1]])
best_interval <- function(keys) { if (!length(keys)) return(NA_character_); z <- t(vapply(keys, parse_iv, integer(2))); keys[order(-(z[, 2] - z[, 1]), z[, 1], z[, 2])[1]] }
# contrasts from profile-level parts: A) spread s1 -> s2 at fixed position; B) position m1 -> m2 at fixed spread (additive: Shapley decomposition into spline-slope (spread terms) and explicit position terms)
state_resp <- function(P, a_part_from, b_part_from) { # P: list of per-cell data frames; returns resp for chosen sources
  eT <- a_part_from$oT + a_part_from$sT + b_part_from$pT; eS <- a_part_from$oS + a_part_from$sS + b_part_from$pS; resp_from(eT, eS) }
contrast_A <- function(G, model, pos, s1, s2) { a <- G[[paste(model, pos, s1)]]; b <- G[[paste(model, pos, s2)]]; ra <- resp_from(a$oT + a$sT + a$pT, a$oS + a$sS + a$pS); rb <- resp_from(b$oT + b$sT + b$pT, b$oS + b$sS + b$pS); d <- 100 * (rb - ra)
  data.frame(model = model, position = pos, s_from = s1, s_to = s2, n_profiles = nrow(a), dT_pp = mean(d$total), dO_pp = mean(d$occupancy), dH_pp = mean(d$heterogeneity), sd_dT = sd(d$total), sd_dO = sd(d$occupancy), sd_dH = sd(d$heterogeneity), dEtaT = mean((b$oT + b$sT + b$pT) - (a$oT + a$sT + a$pT)), dEtaS = mean((b$oS + b$sS + b$pS) - (a$oS + a$sS + a$pS)),
    dEtaT_spread_terms = mean(b$sT - a$sT), dEtaT_position_terms = mean(b$pT - a$pT), dEtaS_spread_terms = mean(b$sS - a$sS), dEtaS_position_terms = mean(b$pS - a$pS)) }
contrast_B <- function(G, model, p1, p2, s) { a <- G[[paste(model, p1, s)]]; b <- G[[paste(model, p2, s)]]; r11 <- resp_from(a$oT + a$sT + a$pT, a$oS + a$sS + a$pS); r22 <- resp_from(b$oT + b$sT + b$pT, b$oS + b$sS + b$pS); d <- 100 * (r22 - r11)
  # Shapley over two factors: (a) spline-slope terms (spread terms evaluated at the second position) and (b) explicit position terms (evaluated at the second position)
  mix <- function(sp_from, pos_from) resp_from(a$oT + sp_from$sT + pos_from$pT, a$oS + sp_from$sS + pos_from$pS); r00 <- r11; r10 <- mix(b, a); r01 <- mix(a, b); phi_a <- 100 * ((r10 - r00) + (r22 - r01)) / 2; phi_b <- 100 * ((r01 - r00) + (r22 - r10)) / 2
  has_pos <- model == "additive"; data.frame(model = model, position_from = p1, position_to = p2, spread = s, n_profiles = nrow(a), dT_pp = mean(d$total), dO_pp = mean(d$occupancy), dH_pp = mean(d$heterogeneity), sd_dT = sd(d$total), sd_dO = sd(d$occupancy), sd_dH = sd(d$heterogeneity),
    spline_slope_dT = mean(phi_a$total), spline_slope_dO = mean(phi_a$occupancy), spline_slope_dH = mean(phi_a$heterogeneity), explicit_position_dT = if (has_pos) mean(phi_b$total) else 0, explicit_position_dO = if (has_pos) mean(phi_b$occupancy) else 0, explicit_position_dH = if (has_pos) mean(phi_b$heterogeneity) else 0,
    dEtaT = mean((b$oT + b$sT + b$pT) - (a$oT + a$sT + a$pT)), dEtaT_spline_slope = mean(b$sT - a$sT), dEtaT_explicit_position = mean(b$pT - a$pT), dEtaS = mean((b$oS + b$sS + b$pS) - (a$oS + a$sS + a$pS)), dEtaS_spline_slope = mean(b$sS - a$sS), dEtaS_explicit_position = mean(b$pS - a$pS)) }
comp_measures <- function(dT, dO, dH, thr = .05) { gross <- abs(dO) + abs(dH); comp <- (gross - abs(dT)) / 2; data.frame(gross_pp = gross, net_total_pp = abs(dT), compensated_pp = comp, opposing = (dO * dH < 0) & pmin(abs(dO), abs(dH)) >= thr, fraction_compensated = ifelse(gross >= 1, comp / gross, NA_real_)) }
