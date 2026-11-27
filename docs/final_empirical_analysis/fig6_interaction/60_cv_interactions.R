# 60_cv_interactions.R: Figure 6 interaction analysis across ALL informative zeta orders (additive vs explicit interaction).
# Models (all coherent component MS-GDUM, zetaGDUM f7e115f, unmodified, primary predictors log10 TP, TN, Chl-a + geographic distance):
#   additive; explicit interaction for ALL predictors (package model "interaction") at endpoint penalties 0.01, 1e-4, 1e-6;
#   explicit interaction for TP ONLY (TN and Chl-a endpoint blocks tied to be identical: both endpoint columns set to half the additive block;
#   design-level construction applied identically to training and evaluation designs; package code unmodified) at the same three penalties.
# Identical across models: informative orders, 5 complete-linkage spatial blocks, sampled combinations (training seed 10*seed+fold, evaluation seed 100*seed+fold, M = 3000 each), seeds 1:3, ridge 0.001.
# For every fit: in-sample (training design/response) and spatially held-out D2 for total, occupancy-associated and abundance-heterogeneity dissimilarity, using
#   ugdm_assess_components(fit, design, response, training_response = y) for both splits; full fitted coefficient vectors (both stages) are exported.
# Usage (run from final/, ZG_LIB6 set):  POOL=shared|full [ONLY_JOBS=1,2] Rscript fig6_interaction/60_cv_interactions.R
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr); library(parallel)})
POOL <- Sys.getenv("POOL", "shared"); NC <- as.integer(Sys.getenv("NCORES", "9")); OUT <- "fig6_interaction/tables"; ONLY <- Sys.getenv("ONLY_JOBS", "")
D <- readRDS(sprintf("fig6_%s/out/prepared.rds", POOL))
get_env <- function(d) data.frame(TP = log10(d$env$TP), TN = log10(d$env$TN), Chla = log10(d$env$Chla)); get_xy <- function(d) as.matrix(d$env[, c("lon", "lat")])
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10, 12), ZP = c(2, 3, 4, 5)); NBLK <- 5; M_TR <- 3000L; M_EV <- 3000L; SEEDS <- 1:3; PENS <- c(1e-2, 1e-4, 1e-6)
hav <- function(xy) { r <- pi / 180; lo <- xy[, 1] * r; la <- xy[, 2] * r; n <- nrow(xy); m <- matrix(0, n, n); for (i in seq_len(n)) { dl <- lo - lo[i]; a <- sin((la - la[i]) / 2)^2 + cos(la[i]) * cos(la) * sin(dl / 2)^2; m[i, ] <- 2 * 6371 * asin(pmin(1, sqrt(a))) }; m }
fold_of <- lapply(D, function(d) cutree(hclust(as.dist(hav(get_xy(d))), "complete"), NBLK))
MODELS <- rbind(data.frame(model = "additive", base = "additive", variant = "additive", penalty = NA_real_, stringsAsFactors = FALSE),
  do.call(rbind, lapply(c("all", "TP"), function(v) data.frame(model = paste0("interaction_", v, "_", c("0.01", "0.0001", "0.000001")), base = "interaction", variant = paste0("interaction_", v), penalty = PENS, stringsAsFactors = FALSE))))
# TP-only interaction: for non-TP predictors both endpoint columns equal half the full additive block, so the endpoint-average ridge acts exactly as for the additive spread and the
# spread effect is position independent (the interaction design with identical endpoint curves).
tp_only <- function(d) { for (v in setdiff(d$vars, "TP")) { j <- match(v, d$vars); ii <- 6L * (j - 1L) + 1:6; blk <- d$dir[, ii[1:3], drop = FALSE] + d$dir[, ii[4:6], drop = FALSE]; d$dir[, ii[1:3]] <- blk / 2; d$dir[, ii[4:6]] <- blk / 2 }; d }
jobs <- expand.grid(g = names(D), seed = SEEDS, stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame(); if (nzchar(ONLY)) jobs <- jobs[as.integer(strsplit(ONLY, ",")[[1]]), , drop = FALSE]
fullnames <- function(dd) colnames(zetaGDUM:::.interaction_layout(dd, !is.null(dd$geo))$X)
coef_rows <- function(fit, dd, id) { nm <- fullnames(dd); do.call(rbind, lapply(names(fit$fits), function(s) { f <- fit$fits[[s]]; b <- drop(fit$transform %*% f$coef); data.frame(id, stage = s, term = nm, coef = b, stage_constant = isTRUE(f$constant), stringsAsFactors = FALSE) })) }
metric_rows <- function(a, id, split) { m <- a$metrics; cbind(id, split = split, m, evaluated_fraction = a$evaluated_fraction, closure_RMSE = a$closure_RMSE, row.names = NULL) }
one <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; e <- get_env(d); xy <- get_xy(d); fo <- fold_of[[j$g]]; k <- j$o; M <- list(); C <- list()
  for (f in sort(unique(fo))) { tr <- which(fo != f); te <- which(fo == f); if (length(te) < k + 2) next
    il <- sample_sets(length(tr), k, M_TR, seed = 10 * j$seed + f); el <- sample_sets(length(te), k, M_EV, seed = 100 * j$seed + f)
    y <- ugdm_component_response(d$x[tr, , drop = FALSE], il); v <- ugdm_component_response(d$x[te, , drop = FALSE], el)
    for (mi in seq_len(nrow(MODELS))) { mo <- MODELS[mi, ]
      r <- tryCatch({ dd <- ugdm_component_design(e[tr, , drop = FALSE], il, mo$base, coordinates = xy[tr, , drop = FALSE], lonlat = TRUE)
        zz <- ugdm_component_design(e[te, , drop = FALSE], el, reference = dd, coordinates = xy[te, , drop = FALSE], lonlat = TRUE)
        if (identical(mo$variant, "interaction_TP")) { dd <- tp_only(dd); zz <- tp_only(zz) }
        fit <- if (is.na(mo$penalty)) ugdm_fit_components(dd, y, architecture = "coherent") else ugdm_fit_components(dd, y, architecture = "coherent", interaction_shrink = mo$penalty)
        id <- data.frame(group = j$g, pool = POOL, order = k, seed = j$seed, fold = f, model = mo$model, variant = mo$variant, penalty = mo$penalty, converged = isTRUE(fit$converged), share_identified = isTRUE(fit$share_identified), n_train_sites = length(tr), n_test_sites = length(te))
        tr_a <- ugdm_assess_components(fit, dd, y, training_response = y); te_a <- ugdm_assess_components(fit, zz, v, training_response = y)
        list(m = rbind(metric_rows(tr_a, id, "train"), metric_rows(te_a, id, "test")), c = coef_rows(fit, dd, id[, c("group", "pool", "order", "seed", "fold", "model", "variant", "penalty")])) },
        error = function(e) { message("ERR ", j$g, " k", k, " s", j$seed, " f", f, " ", mo$model, ": ", conditionMessage(e)); NULL })
      if (!is.null(r)) { M[[length(M) + 1]] <- r$m; C[[length(C) + 1]] <- r$c } } }
  list(m = bind_rows(M), c = bind_rows(C)) }
t0 <- Sys.time(); R <- mclapply(seq_len(nrow(jobs)), one, mc.cores = NC, mc.preschedule = FALSE); cat("elapsed", format(Sys.time() - t0), "jobs", nrow(jobs), "bad", sum(!vapply(R, is.list, logical(1))), "\n")
Mx <- bind_rows(lapply(R, function(x) if (is.list(x)) x$m)); Cx <- bind_rows(lapply(R, function(x) if (is.list(x)) x$c))
write.csv(Mx, file.path(OUT, sprintf("interaction_metrics_%s.csv", POOL)), row.names = FALSE); write.csv(Cx, file.path(OUT, sprintf("interaction_coefficients_cv_%s.csv", POOL)), row.names = FALSE)
cat("rows", nrow(Mx), nrow(Cx), "non-converged fits", sum(!Mx$converged[!duplicated(Mx[c("group", "order", "seed", "fold", "model")])]), "\n")
