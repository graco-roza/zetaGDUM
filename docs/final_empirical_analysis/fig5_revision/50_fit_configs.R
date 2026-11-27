# 50_fit_configs.R: Figure 5 revision, model sensitivity to the environmental transformation.
# Four configurations per community/pool: {log10, raw} nutrient predictors (TP, TN, Chl-a; geographic distance unchanged) x {linear, nonlinear} position, each with its own
# spread-only baseline (same predictors/transformation, same NNLS comparator as the original Figure 5). Identical spatial folds (5 complete-linkage blocks), sampled training/evaluation
# sets (seed = f(seed, fold, order) only, never f(configuration)), seeds 1:3, informative orders and penalties (PEN) in every configuration. Package functions only (zetaGDUM f7e115f, unmodified).
# Held-out D2, training D2, convergence, warnings and within-order Shapley/Owen absolute contributions are stored for every configuration.
# Usage (run from final/):  POOL=shared|full TRIM=0|1 Rscript fig5_revision/50_fit_configs.R
#   TRIM=1: extreme-nutrient sensitivity; lakes above the 99th percentile of TP, TN or Chl-a (raw units, computed within pool/taxon table) are dropped from training AND evaluation; original fold labels are retained.
source("common.R"); source("lib_ugdm.R"); library(parallel)
NC <- as.integer(Sys.getenv("NCORES", "9")); POOL <- Sys.getenv("POOL", "shared"); TRIM <- Sys.getenv("TRIM", "0") == "1"; SEEDS <- 1:3; M_TR <- 3000L; M_EV <- 3000L
ONLY <- Sys.getenv("ONLY_JOBS", ""); OUT <- "fig5_revision/tables"; SUF <- paste0(POOL, if (TRIM) "_trim99" else "")
O <- read.csv("tables/evaluated_orders.csv"); ords <- function(cc) as.integer(O[O$community == COMM[[cc]] & O$pool == POOL & O$informative, "order"]); RESP <- c(incidence = "total", relative = "relative")
jobs <- list(); for (cc in c("BM", "ZP")) for (k in ords(cc)) for (s in SEEDS) jobs[[length(jobs) + 1]] <- list(cc = cc, k = k, seed = s)
if (nzchar(ONLY)) jobs <- jobs[as.integer(strsplit(ONLY, ",")[[1]])]
env_set <- function(env, tf) if (tf == "log10") get_env(env) else data.frame(TP = env$TP, TN = env$TN, Chla = env$Chla)
DATA <- lapply(c("BM", "ZP"), function(cc) { d <- POOLS[[POOL]][[cc]]; xy <- as.matrix(d$env[, c("lon", "lat")]); fold <- spatial_folds(xy); keep <- rep(TRUE, nrow(xy))
  if (TRIM) keep <- with(d$env, TP <= quantile(TP, .99) & TN <= quantile(TN, .99) & Chla <= quantile(Chla, .99))
  list(x = d$x[keep, , drop = FALSE], e = list(log10 = env_set(d$env, "log10")[keep, ], raw = env_set(d$env, "raw")[keep, ]), xy = xy[keep, ], DM = hav(xy[keep, ]), fold = fold[keep], n_dropped = sum(!keep)) }); names(DATA) <- c("BM", "ZP")
POSM <- c(linear = "linear_const", nonlinear = "nonlinear_const")
with_warn <- function(expr) { w <- character(); v <- withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }); list(v = v, w = w) }
one <- function(j) { D <- DATA[[j$cc]]; k <- j$k; V <- list(); IM <- list(); cm <- COMM[[j$cc]]
  for (f in sort(unique(D$fold))) { tr <- which(D$fold != f); te <- which(D$fold == f); if (length(te) < k + 2) next
    il <- sample_sets(length(tr), k, M_TR, seed = 10 * j$seed + f); el <- sample_sets(length(te), k, M_EV, seed = 100 * j$seed + f)
    for (rn in names(RESP)) { r <- ugdm_response(D$x[tr, , drop = FALSE], il, RESP[[rn]]); rt <- ugdm_response(D$x[te, , drop = FALSE], el, RESP[[rn]]); if (sum(r$keep) < 50 || sum(rt$keep) < 50) next
      for (tf in c("log10", "raw")) { E <- D$e[[tf]]; ET <- E[tr, ]; EE <- E[te, ]; DT <- D$DM[tr, tr]; DE <- D$DM[te, te]
        base_row <- function(model, d2h, d2t, conv, nw, msg) data.frame(community = cm, pool = POOL, variant = if (TRIM) "trim99" else "all", transform = tf, order = k, seed = j$seed, fold = f, response = rn, model = model, D2_heldout = d2h, D2_train = d2t, converged = conv, n_warnings = nw, warning_text = msg, n_train_sites = length(tr), n_test_sites = length(te))
        so <- with_warn(tryCatch({ d0 <- ugdm_design(ET, il, dist_matrix = DT); d1 <- ugdm_design(EE, el, dist_matrix = DE, ref = d0); xa <- cbind(d0$Xdir, d0$Xdist); xb <- cbind(d1$Xdir, d1$Xdist)
          fit <- zetaGDUM:::.fit_logit_nnls(-xa[r$keep, , drop = FALSE], r$y[r$keep], r$trials[r$keep]); p <- plogis(fit$intercept - drop(xb %*% fit$coef)); heldout_D2(r, rt, p) }, error = function(e) NA_real_))
        V[[length(V) + 1]] <- base_row("spread_only", so$v, NA_real_, NA, length(so$w), paste(unique(so$w), collapse = " | "))
        for (pm in names(POSM)) { m <- POSM[[pm]]
          res <- with_warn(tryCatch({ dsg <- build(ET, il, DT, m); fit <- fit1(dsg, r); dte <- build(EE, el, DE, m, ref = dsg); p <- ugdm_predict(fit, dte$int); d2h <- heldout_D2(r, rt, p)
            im <- tryCatch(ugdm_within_order_importance(fit, dsg$int, r$y[r$keep], r$trials[r$keep], evaluation = dte$int, evaluation_y = rt$y[rt$keep], evaluation_trials = rt$trials[rt$keep], seed = j$seed, failure = "warn")$hierarchy, error = function(e) NULL)
            list(d2h = d2h, d2t = fit$D2, conv = isTRUE(fit$converged), im = im) }, error = function(e) { message("ERR ", cm, " k", k, " f", f, " ", tf, " ", m, ": ", conditionMessage(e)); NULL }))
          if (is.null(res$v)) { V[[length(V) + 1]] <- base_row(m, NA_real_, NA_real_, FALSE, length(res$w), paste(unique(res$w), collapse = " | ")); next }
          V[[length(V) + 1]] <- base_row(m, res$v$d2h, res$v$d2t, res$v$conv, length(res$w), paste(unique(res$w), collapse = " | "))
          if (!is.null(res$v$im)) { h <- res$v$im; h$community <- cm; h$pool <- POOL; h$variant <- if (TRIM) "trim99" else "all"; h$transform <- tf; h$model <- m; h$order <- k; h$seed <- j$seed; h$fold <- f; h$response <- rn; h$D2_heldout <- res$v$d2h; IM[[length(IM) + 1]] <- h } } } } }
  list(v = bind_rows(V), im = bind_rows(IM)) }
t0 <- Sys.time(); R <- mclapply(jobs, one, mc.cores = NC, mc.preschedule = FALSE); cat("elapsed", format(Sys.time() - t0), "jobs", length(jobs), "errors", sum(sapply(R, inherits, "try-error")), "\n")
V <- bind_rows(lapply(R, function(x) if (is.list(x)) x$v)); IM <- bind_rows(lapply(R, function(x) if (is.list(x)) x$im))
write.csv(V, file.path(OUT, sprintf("sens_validation_folds_%s.csv", SUF)), row.names = FALSE); write.csv(IM, file.path(OUT, sprintf("sens_importance_folds_%s.csv", SUF)), row.names = FALSE)
cat("dropped lakes (trim):", sapply(DATA, function(d) d$n_dropped), "\n")
