# 71m_fit_random_matched.R (supersedes 71_fit_random.R for Figure 5): Figure 5 under random site-level validation (30 independent 70/30 splits, revision_random/tables/random_splits.csv).
# Models (original MS-GDUM structure, log10 TP, TN, Chl-a + geographic distance; incidence ("total") and relative-abundance responses):
#   spread-only (NNLS comparator of the original Figure 5), spread + linear position, spread + nonlinear position.
# Identical across models, responses and (where possible) orders: the split. Sampled set combinations are drawn ONLY from training sites (training sets) or ONLY from test sites (evaluation sets),
# so no site occurs in both; seeds depend on split and order only. Environmental transformations, standardisation and I-spline bases are fitted on the training sites (ugdm_design on training sets;
# the test design is built with ref = training design). Held-out D2 is evaluated against the training-set mean (heldout_D2). Within-order Shapley/Owen allocations (absolute D2 units) evaluated on test sets.
# Usage (final/, ZG_LIB6 set): POOL=shared|full [ONLY_JOBS=1,2] Rscript revision_random/71_fit_random.R
source("common.R"); source("lib_ugdm.R"); library(parallel)
# Matched-penalty version: uses the principal package (ugdm_position_design / ugdm_fit_position). The spread-only comparator is fitted by ugdm_fit_position(spread_only = TRUE) with the SAME bases, link, constraints and penalties as the position models (the earlier comparator was an unpenalised NNLS fit).
build <- function(e, il, dmat, model, ref = NULL) { m <- MODELS[MODELS$model == model, ]; d0 <- if (is.null(ref)) ugdm_design(e, il, dist_matrix = dmat) else ugdm_design(e, il, dist_matrix = dmat, ref = ref$base); list(base = d0, int = ugdm_position_design(d0, e, position = m$position)) }
fit1 <- function(dsg, r, pen = PEN, spread_only = FALSE) { k <- r$keep; ugdm_fit_position(dsg$int, r$y[k], r$trials[k], spread_smooth = pen$spread_smooth, position_smooth = pen$position_smooth, ridge = pen$ridge, spread_only = spread_only) }
NC <- as.integer(Sys.getenv("NCORES", "9")); POOL <- Sys.getenv("POOL", "shared"); M_TR <- 3000L; M_EV <- 3000L; ONLY <- Sys.getenv("ONLY_JOBS", ""); OUT <- "revision_random/tables"
SPL <- read.csv(file.path(OUT, "random_splits.csv"), stringsAsFactors = FALSE); NS <- max(SPL$split)
O <- read.csv("tables/evaluated_orders.csv"); ords <- function(cc) as.integer(O[O$community == COMM[[cc]] & O$pool == POOL & O$informative, "order"]); RESP <- c(incidence = "total", relative = "relative")
jobs <- list(); for (cc in c("BM", "ZP")) for (k in ords(cc)) for (s in seq_len(NS)) jobs[[length(jobs) + 1]] <- list(cc = cc, k = k, split = s); if (nzchar(ONLY)) jobs <- jobs[as.integer(strsplit(ONLY, ",")[[1]])]
DATA <- lapply(c("BM", "ZP"), function(cc) { d <- POOLS[[POOL]][[cc]]; xy <- as.matrix(d$env[, c("lon", "lat")]); list(x = d$x, e = get_env(d$env), xy = xy, DM = hav(xy), id = d$env$id) }); names(DATA) <- c("BM", "ZP")
ROLE <- lapply(c("BM", "ZP"), function(cc) { z <- SPL[SPL$pool == POOL & SPL$community == COMM[[cc]], ]; split(setNames(z$role, z$site_id), z$split) }); names(ROLE) <- c("BM", "ZP")
POSM <- c(linear = "linear_const", nonlinear = "nonlinear_const")
with_warn <- function(expr) { w <- character(); v <- withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }); list(v = v, w = w) }
one <- function(j) { D <- DATA[[j$cc]]; k <- j$k; cm <- COMM[[j$cc]]; V <- list(); IM <- list(); s <- j$split
  ro <- ROLE[[j$cc]][[as.character(s)]][D$id]; tr <- which(ro == "train"); te <- which(ro == "test"); stopifnot(length(intersect(tr, te)) == 0, length(tr) + length(te) == length(D$id))
  il <- sample_sets(length(tr), k, M_TR, seed = 100000L + 1000L * s + k); el <- sample_sets(length(te), k, M_EV, seed = 200000L + 1000L * s + k); ET <- D$e[tr, ]; EE <- D$e[te, ]; DT <- D$DM[tr, tr]; DE <- D$DM[te, te]
  for (rn in names(RESP)) { r <- ugdm_response(D$x[tr, , drop = FALSE], il, RESP[[rn]]); rt <- ugdm_response(D$x[te, , drop = FALSE], el, RESP[[rn]]); if (sum(r$keep) < 50 || sum(rt$keep) < 50) next
    row0 <- function(model, d2h, d2t, conv, nw, msg) data.frame(community = cm, pool = POOL, order = k, split = s, response = rn, model = model, D2_heldout = d2h, D2_train = d2t, converged = conv, n_warnings = nw, warning_text = msg, n_train_sites = length(tr), n_test_sites = length(te), n_sets_train = sum(r$keep), n_sets_test = sum(rt$keep))
    so <- with_warn(tryCatch({ dsg <- build(ET, il, DT, "linear_const"); fit <- fit1(dsg, r, spread_only = TRUE); dte <- build(EE, il = el, DE, "linear_const", ref = dsg); p <- ugdm_predict(fit, dte$int); list(d2h = heldout_D2(r, rt, p), d2t = fit$D2) }, error = function(e) list(d2h = NA_real_, d2t = NA_real_)))
    V[[length(V) + 1]] <- row0("spread_only", so$v$d2h, so$v$d2t, TRUE, length(so$w), paste(unique(so$w), collapse = " | "))
    for (pm in names(POSM)) { m <- POSM[[pm]]
      res <- with_warn(tryCatch({ dsg <- build(ET, il, DT, m); fit <- fit1(dsg, r); dte <- build(EE, el, DE, m, ref = dsg); p <- ugdm_predict(fit, dte$int); d2h <- heldout_D2(r, rt, p)
        im <- tryCatch(ugdm_within_order_importance(fit, dsg$int, r$y[r$keep], r$trials[r$keep], evaluation = dte$int, evaluation_y = rt$y[rt$keep], evaluation_trials = rt$trials[rt$keep], seed = s, failure = "warn")$hierarchy, error = function(e) NULL)
        list(d2h = d2h, d2t = fit$D2, conv = isTRUE(fit$converged), im = im) }, error = function(e) { message("ERR ", cm, " k", k, " s", s, " ", m, ": ", conditionMessage(e)); NULL }))
      if (is.null(res$v)) { V[[length(V) + 1]] <- row0(m, NA_real_, NA_real_, FALSE, length(res$w), paste(unique(res$w), collapse = " | ")); next }
      V[[length(V) + 1]] <- row0(m, res$v$d2h, res$v$d2t, res$v$conv, length(res$w), paste(unique(res$w), collapse = " | "))
      if (!is.null(res$v$im)) { h <- res$v$im; h$community <- cm; h$pool <- POOL; h$model <- m; h$order <- k; h$split <- s; h$response <- rn; h$D2_heldout <- res$v$d2h; IM[[length(IM) + 1]] <- h } } }
  list(v = bind_rows(V), im = bind_rows(IM)) }
t0 <- Sys.time(); R <- mclapply(jobs, one, mc.cores = NC, mc.preschedule = FALSE); cat("elapsed", format(Sys.time() - t0), "jobs", length(jobs), "bad", sum(!vapply(R, is.list, logical(1))), "\n")
V <- bind_rows(lapply(R, function(x) if (is.list(x)) x$v)); IM <- bind_rows(lapply(R, function(x) if (is.list(x)) x$im))
sfx <- if (nzchar(ONLY)) "_test" else ""; write.csv(V, file.path(OUT, sprintf("rv5m_validation_splits_%s%s.csv", POOL, sfx)), row.names = FALSE); write.csv(IM, file.path(OUT, sprintf("rv5m_importance_splits_%s%s.csv", POOL, sfx)), row.names = FALSE)
cat("rows", nrow(V), nrow(IM), "nonconverged", sum(V$converged %in% FALSE), "NA D2", sum(!is.finite(V$D2_heldout)), "\n")
