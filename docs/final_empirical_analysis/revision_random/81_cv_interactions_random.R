# 81_cv_interactions_random.R: Figure 6 interaction analysis under RANDOM SITE-LEVEL validation (30 independent 70/30 splits, revision_random/tables/random_splits.csv).
# Same models as fig6_interaction/60_cv_interactions.R (coherent component MS-GDUM, zetaGDUM f7e115f unmodified; log10 TP, TN, Chl-a + geographic distance):
#   additive; explicit interaction for all predictors or TP only at endpoint penalties 0.01, 1e-4, 1e-6. Training combinations are sampled only from training sites and evaluation combinations only from
#   test sites; environmental bases, centring/scaling and I-spline bases are built from the training design (evaluation design uses reference = training design). Seeds depend on split and order only,
#   so every model sees identical sets. Orders: the informative orders (BM 2,3,4,5,6,8,10,12; ZP 2,3,4,5). ridge 0.001, M = 3000 training / 3000 evaluation sets per split.
# Usage (final/, ZG_LIB6 set): POOL=shared|full [ONLY_JOBS=1,2] Rscript revision_random/81_cv_interactions_random.R
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr); library(parallel)})
POOL <- Sys.getenv("POOL", "shared"); NC <- as.integer(Sys.getenv("NCORES", "9")); OUT <- "revision_random/tables"; ONLY <- Sys.getenv("ONLY_JOBS", "")
D <- readRDS(sprintf("fig6_%s/out/prepared.rds", POOL))
get_env <- function(d) data.frame(TP = log10(d$env$TP), TN = log10(d$env$TN), Chla = log10(d$env$Chla)); get_xy <- function(d) as.matrix(d$env[, c("lon", "lat")])
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10, 12), ZP = c(2, 3, 4, 5)); NS <- 30L; M_TR <- 3000L; M_EV <- 3000L; PENS <- c(1e-2, 1e-4, 1e-6)
hav <- function(xy) { r <- pi / 180; lo <- xy[, 1] * r; la <- xy[, 2] * r; n <- nrow(xy); m <- matrix(0, n, n); for (i in seq_len(n)) { dl <- lo - lo[i]; a <- sin((la - la[i]) / 2)^2 + cos(la[i]) * cos(la) * sin(dl / 2)^2; m[i, ] <- 2 * 6371 * asin(pmin(1, sqrt(a))) }; m }
SPL <- read.csv(file.path(OUT, "random_splits.csv"), stringsAsFactors = FALSE); LAB <- c(BM = "EPA macroinvertebrates", ZP = "EPA zooplankton")
role_of <- lapply(names(D), function(g) { z <- SPL[SPL$pool == POOL & SPL$community == LAB[[g]], ]; split(setNames(z$role, z$site_id), z$split) }); names(role_of) <- names(D)
MODELS <- rbind(data.frame(model = "additive", base = "additive", variant = "additive", penalty = NA_real_, stringsAsFactors = FALSE),
  do.call(rbind, lapply(c("all", "TP"), function(v) data.frame(model = paste0("interaction_", v, "_", c("0.01", "0.0001", "0.000001")), base = "interaction", variant = paste0("interaction_", v), penalty = PENS, stringsAsFactors = FALSE))))
# TP-only interaction: for non-TP predictors both endpoint columns equal half the full additive block, so the endpoint-average ridge acts exactly as for the additive spread and the
# spread effect is position independent (the interaction design with identical endpoint curves).
tp_only <- function(d) { for (v in setdiff(d$vars, "TP")) { j <- match(v, d$vars); ii <- 6L * (j - 1L) + 1:6; blk <- d$dir[, ii[1:3], drop = FALSE] + d$dir[, ii[4:6], drop = FALSE]; d$dir[, ii[1:3]] <- blk / 2; d$dir[, ii[4:6]] <- blk / 2 }; d }
jobs <- expand.grid(g = names(D), split = seq_len(NS), stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame(); if (nzchar(ONLY)) jobs <- jobs[as.integer(strsplit(ONLY, ",")[[1]]), , drop = FALSE]
fullnames <- function(dd) colnames(zetaGDUM:::.interaction_layout(dd, !is.null(dd$geo))$X)
coef_rows <- function(fit, dd, id) { nm <- fullnames(dd); do.call(rbind, lapply(names(fit$fits), function(s) { f <- fit$fits[[s]]; b <- drop(fit$transform %*% f$coef); data.frame(id, stage = s, term = nm, coef = b, stage_constant = isTRUE(f$constant), stringsAsFactors = FALSE) })) }
metric_rows <- function(a, id, split) { m <- a$metrics; cbind(id, split = split, m, evaluated_fraction = a$evaluated_fraction, closure_RMSE = a$closure_RMSE, row.names = NULL) }
one <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; e <- get_env(d); xy <- get_xy(d); k <- j$o; s <- j$split; M <- list(); C <- list()
  ro <- role_of[[j$g]][[as.character(s)]][d$env$id]; tr <- which(ro == "train"); te <- which(ro == "test"); stopifnot(length(intersect(tr, te)) == 0, length(tr) + length(te) == nrow(e))
  il <- sample_sets(length(tr), k, M_TR, seed = 100000L + 1000L * s + k); el <- sample_sets(length(te), k, M_EV, seed = 200000L + 1000L * s + k)
  y <- ugdm_component_response(d$x[tr, , drop = FALSE], il); v <- ugdm_component_response(d$x[te, , drop = FALSE], el)
  for (mi in seq_len(nrow(MODELS))) { mo <- MODELS[mi, ]
    r <- tryCatch({ dd <- ugdm_component_design(e[tr, , drop = FALSE], il, mo$base, coordinates = xy[tr, , drop = FALSE], lonlat = TRUE)
      zz <- ugdm_component_design(e[te, , drop = FALSE], el, reference = dd, coordinates = xy[te, , drop = FALSE], lonlat = TRUE)
      if (identical(mo$variant, "interaction_TP")) { dd <- tp_only(dd); zz <- tp_only(zz) }
      fit <- if (is.na(mo$penalty)) ugdm_fit_components(dd, y, architecture = "coherent") else ugdm_fit_components(dd, y, architecture = "coherent", interaction_shrink = mo$penalty)
      id <- data.frame(group = j$g, pool = POOL, order = k, rsplit = s, model = mo$model, variant = mo$variant, penalty = mo$penalty, converged = isTRUE(fit$converged), share_identified = isTRUE(fit$share_identified), n_train_sites = length(tr), n_test_sites = length(te))
      tr_a <- ugdm_assess_components(fit, dd, y, training_response = y); te_a <- ugdm_assess_components(fit, zz, v, training_response = y)
      list(m = rbind(metric_rows(tr_a, id, "train"), metric_rows(te_a, id, "test")), c = coef_rows(fit, dd, id[, c("group", "pool", "order", "rsplit", "model", "variant", "penalty")])) },
      error = function(e) { message("ERR ", j$g, " k", k, " s", s, " ", mo$model, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) { M[[length(M) + 1]] <- r$m; C[[length(C) + 1]] <- r$c } }
  list(m = bind_rows(M), c = bind_rows(C)) }
t0 <- Sys.time(); R <- mclapply(seq_len(nrow(jobs)), one, mc.cores = NC, mc.preschedule = FALSE); cat("elapsed", format(Sys.time() - t0), "jobs", nrow(jobs), "bad", sum(!vapply(R, is.list, logical(1))), "\n")
Mx <- bind_rows(lapply(R, function(x) if (is.list(x)) x$m)); Cx <- bind_rows(lapply(R, function(x) if (is.list(x)) x$c))
write.csv(Mx, file.path(OUT, sprintf("rv6_metrics_%s%s.csv", POOL, if (nzchar(ONLY)) "_test" else "")), row.names = FALSE); write.csv(Cx, file.path(OUT, sprintf("rv6_coefficients_splits_%s%s.csv", POOL, if (nzchar(ONLY)) "_test" else "")), row.names = FALSE)
cat("rows", nrow(Mx), nrow(Cx), "non-converged fits", sum(!Mx$converged[!duplicated(Mx[c("group", "order", "split", "model")])]), "\n")
