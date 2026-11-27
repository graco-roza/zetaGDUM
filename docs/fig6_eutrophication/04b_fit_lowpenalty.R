# 04b_fit_lowpenalty.R (copy of 04_fit.R with interaction_shrink = 1e-4, primary predictors only; sensitivity of the explicit interaction to its endpoint penalty): coherent component MS-GDUM fits per group and informative order.
#  (a) predictor importance (Shapley hierarchy: spread / position / explicit interaction per predictor + geographic distance), fitted on training spatial blocks
#      and attributed on the held-out block (leave-one-block-out), repeated over combination-sampling seeds
#  (b) all-site additive and interaction fits -> conditional prediction surfaces over TP position x spread (supported configurations only) + diagnostics
source("lib6.R"); library(parallel); NC <- as.integer(Sys.getenv("NCORES", "4"))
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10), ZP = c(2, 3, 4, 5, 6), PP = c(2, 3, 4, 5)); NBLK <- 5; M_TR <- 3000L; M_EV <- 3000L; SEEDS <- 1:3
hav <- function(xy) { r <- pi / 180; lo <- xy[, 1] * r; la <- xy[, 2] * r; n <- nrow(xy); m <- matrix(0, n, n); for (i in seq_len(n)) { dl <- lo - lo[i]; a <- sin((la - la[i]) / 2)^2 + cos(la[i]) * cos(la) * sin(dl / 2)^2; m[i, ] <- 2 * 6371 * asin(pmin(1, sqrt(a))) }; m }
folds <- function(xy) cutree(hclust(as.dist(hav(xy)), "complete"), NBLK)
jobs <- expand.grid(g = names(D), chla = TRUE, seed = SEEDS, stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame()
fold_of <- lapply(D, function(d) folds(get_xy(d))); saveRDS(fold_of, "out/folds_unused.rds")
imp1 <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; e <- get_env(d, j$chla); xy <- get_xy(d); fo <- fold_of[[j$g]]; k <- j$o; out <- list()
  for (f in sort(unique(fo))) { tr <- which(fo != f); te <- which(fo == f); if (length(te) < k + 2) next
    r <- tryCatch({ il <- sample_sets(length(tr), k, M_TR, seed = 10 * j$seed + f); el <- sample_sets(length(te), k, M_EV, seed = 100 * j$seed + f)   # local (within-partition) set indices
      y <- ugdm_component_response(d$x[tr, , drop = FALSE], il); v <- ugdm_component_response(d$x[te, , drop = FALSE], el)
      dd <- ugdm_component_design(e[tr, , drop = FALSE], il, "interaction", coordinates = xy[tr, ], lonlat = TRUE)
      zz <- ugdm_component_design(e[te, , drop = FALSE], el, reference = dd, coordinates = xy[te, ], lonlat = TRUE)
      fit <- ugdm_fit_components(dd, y, architecture = "coherent", interaction_shrink = 1e-4); im <- ugdm_component_importance(fit, dd, y, evaluation = zz, evaluation_response = v, seed = j$seed)
      h <- im$hierarchy; h$fold <- f; h }, error = function(e) { message("ERR imp ", i, " fold ", f, ": ", conditionMessage(e)); NULL }); if (!is.null(r)) out[[length(out) + 1]] <- r }
  if (length(out)) cbind(group = j$g, chla = j$chla, order = k, seed = j$seed, bind_rows(out)) else NULL }
t0 <- Sys.time(); IMP <- bind_rows(mclapply(seq_len(nrow(jobs)), imp1, mc.cores = NC, mc.preschedule = FALSE)); cat("importance elapsed", format(Sys.time() - t0), "\n"); saveRDS(IMP, "out/importance_folds_pen1e-4.rds")
# ---- (b) surfaces from all-site fits
sj <- expand.grid(g = names(D), stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame()
surf1 <- function(i) { j <- sj[i, ]; d <- D[[j$g]]; e <- get_env(d); xy <- get_xy(d); k <- j$o; n <- nrow(d$x)
  idx <- sample_sets(n, k, 6000, seed = 7000 + k); ev <- sample_sets(n, k, 8000, seed = 9000 + k); y <- ugdm_component_response(d$x, idx); v <- ugdm_component_response(d$x, ev); res <- list(); diag <- list()
  for (m in c("interaction")) { dd <- ugdm_component_design(e, idx, m, coordinates = xy, lonlat = TRUE); zz <- ugdm_component_design(e, ev, reference = dd, coordinates = xy, lonlat = TRUE)
    fit <- ugdm_fit_components(dd, y, architecture = "coherent", interaction_shrink = 1e-4); dg <- ugdm_component_diagnostics(fit, dd, y)
    for (bn in c(8, 12)) for (pr in c("TP")) { s <- ugdm_component_surface(fit, zz, pr, bins = bn, min_sets = 8, min_sites = 12); s$model <- m; s$bins <- bn; res[[length(res) + 1]] <- s }
    diag[[m]] <- list(stages = dg$stages, interaction_residual = dg$interaction_residual, share_identified = dg$share_identified) }
  list(group = j$g, order = k, surfaces = bind_rows(res), diag = diag) }
t0 <- Sys.time(); SF <- mclapply(seq_len(nrow(sj)), surf1, mc.cores = NC, mc.preschedule = FALSE); cat("surfaces elapsed", format(Sys.time() - t0), "\n"); saveRDS(SF, "out/surfaces_pen1e-4.rds")
