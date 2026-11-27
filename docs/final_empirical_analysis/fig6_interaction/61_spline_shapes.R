# 61_spline_shapes.R: fitted TP-spread spline shapes at different TP positions, for the total-dissimilarity stage and the heterogeneity-share stage,
# for every informative zeta order and every model of 60_cv_interactions.R (additive; explicit interaction, all predictors or TP only; penalties 0.01, 1e-4, 1e-6).
# All-site fits (6000 sampled combinations, seed 7000 + k, as in the existing Figure 6 surface fits). Package code unmodified (zetaGDUM f7e115f).
# Spread shape. For a pair of sites the spread term is sum_j beta_j |I_j(z1) - I_j(z2)| and the I-splines are monotone, so for k = 2 it equals |f(z1) - f(z2)| with
#   f(z) = sum_j beta_j I_j(z).  Additive model: beta_j = TP_fj (position independent).  Explicit interaction: the design interpolates two endpoint curves by the clipped mean set position u,
#   beta_j(u) = (1 - u) beta_{L,j} + u beta_{U,j}.  f is plotted (centred at the lowest TP) on the MODEL (logit) scale. The position terms (TP, TP_quadratic, TP_cubic) shift the level, not the shape.
# Link decomposition: model-scale spread increment delta = f(z80) - f(z20) (20th and 80th percentile of log10 TP) is applied at the fitted baseline linear predictors of the sets nearest each target position
#   (median of the 200 nearest training sets, both stages) and converted to the response scale: in the ADDITIVE model delta is identical at every position, so position differences on the response scale are
#   link-scale (nonlinear logistic) effects; in the INTERACTION models delta itself changes with position (explicit interaction).
# Usage (run from final/, ZG_LIB6 set):  POOL=shared|full [ONLY_JOBS=1] Rscript fig6_interaction/61_spline_shapes.R
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr); library(parallel)})
POOL <- Sys.getenv("POOL", "shared"); NC <- as.integer(Sys.getenv("NCORES", "9")); OUT <- "fig6_interaction/tables"; ONLY <- Sys.getenv("ONLY_JOBS", "")
D <- readRDS(sprintf("fig6_%s/out/prepared.rds", POOL))
get_env <- function(d) data.frame(TP = log10(d$env$TP), TN = log10(d$env$TN), Chla = log10(d$env$Chla)); get_xy <- function(d) as.matrix(d$env[, c("lon", "lat")])
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10, 12), ZP = c(2, 3, 4, 5)); PENS <- c(1e-2, 1e-4, 1e-6); PEN_LAB <- c("0.01", "0.0001", "0.000001")
MODELS <- rbind(data.frame(model = "additive", base = "additive", variant = "additive", penalty = NA_real_, stringsAsFactors = FALSE),
  do.call(rbind, lapply(c("all", "TP"), function(v) data.frame(model = paste0("interaction_", v, "_", PEN_LAB), base = "interaction", variant = paste0("interaction_", v), penalty = PENS, stringsAsFactors = FALSE))))
tp_only <- function(d) { for (v in setdiff(d$vars, "TP")) { j <- match(v, d$vars); ii <- 6L * (j - 1L) + 1:6; blk <- d$dir[, ii[1:3], drop = FALSE] + d$dir[, ii[4:6], drop = FALSE]; d$dir[, ii[1:3]] <- blk / 2; d$dir[, ii[4:6]] <- blk / 2 }; d }
POSITIONS <- c(Low = 10, Intermediate = 35, High = 100)   # target TP, ug/L (as in the existing Figure 6 surfaces)
jobs <- expand.grid(g = names(D), stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame(); if (nzchar(ONLY)) jobs <- jobs[as.integer(strsplit(ONLY, ",")[[1]]), , drop = FALSE]
one <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; e <- get_env(d); xy <- get_xy(d); k <- j$o; n <- nrow(e)
  idx <- sample_sets(n, k, 6000, seed = 7000 + k); y <- ugdm_component_response(d$x, idx); S <- list(); Cf <- list(); L <- list()
  zg <- seq(quantile(e$TP, .02), quantile(e$TP, .98), length.out = 60); z20 <- quantile(e$TP, .2); z80 <- quantile(e$TP, .8)
  for (mi in seq_len(nrow(MODELS))) { mo <- MODELS[mi, ]
    r <- tryCatch({ dd <- ugdm_component_design(e, idx, mo$base, coordinates = xy, lonlat = TRUE); if (identical(mo$variant, "interaction_TP")) dd <- tp_only(dd)
      fit <- if (is.na(mo$penalty)) ugdm_fit_components(dd, y, architecture = "coherent") else ugdm_fit_components(dd, y, architecture = "coherent", interaction_shrink = mo$penalty)
      lay <- zetaGDUM:::.interaction_layout(dd, !is.null(dd$geo)); nm <- colnames(lay$X); Xk <- zetaGDUM:::.component_prediction_matrix(fit, dd)
      tb <- dd$basis; ctr <- tb$center[["TP"]]; scl <- tb$scale[["TP"]]; spl <- tb$splines[[match("TP", dd$vars)]]; bnd <- attr(spl, "Boundary.knots")
      Bz <- function(x) { zz <- pmin(pmax((x - ctr) / scl, bnd[1]), bnd[2]); unclass(predict(spl, zz))[, 1:3, drop = FALSE] }
      uof <- function(lx) { mz <- (lx - ctr) / scl; min(1, max(0, (mz - bnd[1]) / (bnd[2] - bnd[1]))) }
      meanTP <- rowMeans(matrix(e$TP[idx], nrow(idx))); rows <- list(); cf <- list(); lk <- list()
      for (st in names(fit$fits)) { f <- fit$fits[[st]]; b <- setNames(drop(fit$transform %*% f$coef), nm); eta <- drop(Xk %*% f$coef)
        cf[[st]] <- data.frame(group = j$g, pool = POOL, order = k, model = mo$model, variant = mo$variant, penalty = mo$penalty, stage = st, term = nm, coef = b, converged = isTRUE(f$converged), stage_constant = isTRUE(f$constant))
        bL <- if (mo$base == "additive") b[paste0("TP_f", 1:3)] else b[paste0("TP_f", 1:3, "_L")]; bU <- if (mo$base == "additive") bL else b[paste0("TP_f", 1:3, "_U")]
        for (pl in c(names(POSITIONS), "u0_lower_endpoint", "u1_upper_endpoint")) { lx <- if (pl %in% names(POSITIONS)) log10(POSITIONS[[pl]]) else NA_real_; u <- if (pl == "u0_lower_endpoint") 0 else if (pl == "u1_upper_endpoint") 1 else uof(lx)
          bb <- (1 - u) * bL + u * bU; fz <- drop(Bz(zg) %*% bb); rows[[length(rows) + 1]] <- data.frame(group = j$g, pool = POOL, order = k, model = mo$model, variant = mo$variant, penalty = mo$penalty, stage = st, position = pl, TP_position_ugL = if (is.na(lx)) NA_real_ else 10^lx, u = u, TP_log10 = zg, TP_ugL = 10^zg, f = fz, f_centred = fz - fz[1], beta1 = bb[1], beta2 = bb[2], beta3 = bb[3])
          if (pl %in% names(POSITIONS)) { near <- order(abs(meanTP - lx))[1:200]; fd <- drop(Bz(z80) %*% bb) - drop(Bz(z20) %*% bb); lk[[length(lk) + 1]] <- data.frame(stage = st, position = pl, TP_position_ugL = 10^lx, u = u, eta0 = median(eta[near]), delta_model_scale = fd, window_mean_TP_ugL = 10^mean(meanTP[near])) } } }
      ld <- do.call(rbind, lk); w <- ld |> select(stage, position, eta0, delta_model_scale, TP_position_ugL, u, window_mean_TP_ugL) |> pivot_wider(names_from = stage, values_from = c(eta0, delta_model_scale))
      if (all(c("eta0_total", "eta0_share") %in% names(w))) { w <- w |> mutate(total0 = plogis(eta0_total), het0 = total0 * plogis(eta0_share), total1 = plogis(eta0_total + delta_model_scale_total), het1 = total1 * plogis(eta0_share + delta_model_scale_share),
          d_total_pp = 100 * (total1 - total0), d_het_pp = 100 * (het1 - het0), d_occ_pp = d_total_pp - d_het_pp, closure_err = d_total_pp - d_het_pp - d_occ_pp) }
      w <- cbind(group = j$g, pool = POOL, order = k, model = mo$model, variant = mo$variant, penalty = mo$penalty, converged = isTRUE(fit$converged), w)
      list(s = bind_rows(rows), c = bind_rows(cf), l = w) }, error = function(e) { message("ERR ", j$g, " k", k, " ", mo$model, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) { S[[length(S) + 1]] <- r$s; Cf[[length(Cf) + 1]] <- r$c; L[[length(L) + 1]] <- r$l } }
  list(s = bind_rows(S), c = bind_rows(Cf), l = bind_rows(L)) }
t0 <- Sys.time(); R <- mclapply(seq_len(nrow(jobs)), one, mc.cores = NC, mc.preschedule = FALSE); cat("elapsed", format(Sys.time() - t0), "jobs", nrow(jobs), "\n")
write.csv(bind_rows(lapply(R, function(x) if (is.list(x)) x$s)), file.path(OUT, sprintf("interaction_spline_values_%s.csv", POOL)), row.names = FALSE)
write.csv(bind_rows(lapply(R, function(x) if (is.list(x)) x$c)), file.path(OUT, sprintf("interaction_coefficients_allsites_%s.csv", POOL)), row.names = FALSE)
write.csv(bind_rows(lapply(R, function(x) if (is.list(x)) x$l)), file.path(OUT, sprintf("interaction_link_decomposition_%s.csv", POOL)), row.names = FALSE); cat("done\n")
