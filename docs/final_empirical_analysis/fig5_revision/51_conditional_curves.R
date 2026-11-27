# 51_conditional_curves.R: Figure 5 revision, CONDITIONAL model predictions (replace the marginal equal-count bin summaries of the original Figure 5).
# For each community, pool, configuration ({log10, raw} x {linear, nonlinear position}), informative zeta order and response representation, the MS-GDUM is fitted to all lakes
# (training sets: 6000 sampled combinations, seed 7000 + k; evaluation sets: 8000, seed 9000 + k; identical to the original curve fits; penalties PEN).
# A conditional curve for predictor v and axis a (spread or position) is the model linear predictor evaluated on a synthetic design row in which
#   * the design columns of (v, a) take the mean of the 200 evaluation sets whose set statistic (mean of v for position; mean pairwise |difference| of v for spread) is closest to the grid value x
#     (so that the columns of the nonlinear position block, or of the I-spline spread block, stay mutually consistent), and
#   * every other design column (spread, position and geographic distance terms of all predictors, intercept) is held at its median over the evaluation sets (spread of v held constant when position varies and vice versa).
# Predictions are given on the link (logit) and response (dissimilarity) scales. Environmental feasibility: the grid value, together with the other set statistics held at their medians
# (position and spread of TP, TN, Chl-a and mean distance; standardised), must lie within the 95th percentile of the leave-one-out nearest-neighbour distance among real evaluation sets (flag `feasible`).
# Two grids: range_type = "common" (intersection of the 5-95% ranges of the set statistic across all informative orders of that community, so that orders are comparable) and "own" (5-95% range of that order).
# Usage (run from final/):  POOL=shared|full Rscript fig5_revision/51_conditional_curves.R
source("common.R"); source("lib_ugdm.R"); library(parallel)
NC <- as.integer(Sys.getenv("NCORES", "9")); POOL <- Sys.getenv("POOL", "shared"); G <- 25L; NWIN <- 200L; OUT <- "fig5_revision/tables"
O <- read.csv("tables/evaluated_orders.csv"); ords <- function(cc) as.integer(O[O$community == COMM[[cc]] & O$pool == POOL & O$informative, "order"]); RESP <- c(incidence = "total", relative = "relative")
env_set <- function(env, tf) if (tf == "log10") get_env(env) else data.frame(TP = env$TP, TN = env$TN, Chla = env$Chla)
DATA <- lapply(c("BM", "ZP"), function(cc) { d <- POOLS[[POOL]][[cc]]; xy <- as.matrix(d$env[, c("lon", "lat")]); list(x = d$x, e = list(log10 = env_set(d$env, "log10"), raw = env_set(d$env, "raw")), DM = hav(xy)) }); names(DATA) <- c("BM", "ZP")
VARS <- c("TP", "TN", "Chla")
set_stats <- function(D, tf, il) { E <- D$e[[tf]]; k <- ncol(il); out <- list()
  for (v in VARS) { vv <- E[[v]]; out[[paste0("pos_", v)]] <- rowMeans(matrix(vv[il], nrow(il))); s <- rep(0, nrow(il)); for (a in 1:(k - 1)) for (b in (a + 1):k) s <- s + abs(vv[il[, a]] - vv[il[, b]]); out[[paste0("spr_", v)]] <- s / choose(k, 2) }
  s <- rep(0, nrow(il)); for (a in 1:(k - 1)) for (b in (a + 1):k) s <- s + D$DM[cbind(il[, a], il[, b])]; out$dist <- s / choose(k, 2); as.data.frame(out) }
el_sets <- function(D, k) sample_sets(nrow(D$x), k, 8000, seed = 9000 + k)
# pass 1: common ranges (no model fitting)
rng <- list(); for (cc in c("BM", "ZP")) for (tf in c("log10", "raw")) { st <- lapply(ords(cc), function(k) set_stats(DATA[[cc]], tf, el_sets(DATA[[cc]], k))); names(st) <- ords(cc)
  for (nm in setdiff(names(st[[1]]), "dist")) { q <- sapply(st, function(s) quantile(s[[nm]], c(.05, .95))); rng[[paste(cc, tf, nm)]] <- c(lo = max(q[1, ]), hi = min(q[2, ]), per_order_lo = list(q[1, ]), per_order_hi = list(q[2, ])) } }
cfgs <- expand.grid(tf = c("log10", "raw"), pm = c("linear", "nonlinear"), stringsAsFactors = FALSE); POSM <- c(linear = "linear_const", nonlinear = "nonlinear_const")
jobs <- list(); for (cc in c("BM", "ZP")) for (k in ords(cc)) jobs[[length(jobs) + 1]] <- list(cc = cc, k = k)
run <- function(j) { D <- DATA[[j$cc]]; k <- j$k; cm <- COMM[[j$cc]]; out <- list(); il <- sample_sets(nrow(D$x), k, 6000, seed = 7000 + k); el <- el_sets(D, k)
  for (tf in c("log10", "raw")) { E <- D$e[[tf]]; st <- set_stats(D, tf, el); S <- scale(st); Sc <- attr(S, "scaled:center"); Ss <- attr(S, "scaled:scale")
    for (rn in names(RESP)) { r <- ugdm_response(D$x, il, RESP[[rn]]); rt <- ugdm_response(D$x, el, RESP[[rn]]); kk <- which(rt$keep)
      for (pm in c("linear", "nonlinear")) { m <- POSM[[pm]]
        res <- tryCatch(suppressWarnings({ dsg <- build(E, il, D$DM, m); fit <- fit1(dsg, r); dte <- build(E, el, D$DM, m, ref = dsg); X <- zetaGDUM:::.interaction_layout(dte$int, TRUE)$X; list(fit = fit, X = X) }), error = function(e) { message("ERR ", cm, " k", k, " ", tf, " ", pm, ": ", conditionMessage(e)); NULL }); if (is.null(res)) next
        fit <- res$fit; X <- res$X[kk, , drop = FALSE]; cf <- fit$coef; stopifnot(identical(colnames(X), names(cf))); stk <- st[kk, , drop = FALSE]; Sk <- S[kk, , drop = FALSE]
        Xref <- apply(X, 2, median); sref <- apply(stk, 2, median); eta_ref <- sum(Xref * cf)
        sub <- if (nrow(Sk) > 2500) sample.int(nrow(Sk), 2500) else seq_len(nrow(Sk)); Ssub <- Sk[sub, , drop = FALSE]
        loo <- sapply(seq_len(min(800, nrow(Ssub))), function(i) { d <- sqrt(colSums((t(Ssub[-i, , drop = FALSE]) - Ssub[i, ])^2)); min(d) }); thr <- unname(quantile(loo, .95))
        for (v in VARS) for (ax in c("position", "spread")) { stat <- paste0(if (ax == "position") "pos_" else "spr_", v); xs <- stk[[stat]]
          cols <- if (ax == "spread") grep(paste0("^", v, "_f[123]$"), colnames(X)) else which(colnames(X) == v | grepl(paste0("^", v, "_(quadratic|cubic)$"), colnames(X))); stopifnot(length(cols) >= 3 - 2 * (ax == "position" && pm == "linear"))
          rg <- rng[[paste(j$cc, tf, stat)]]; own <- quantile(xs, c(.05, .95));
          for (rt_ in c("common", "own")) { lims <- if (rt_ == "common") c(rg$lo, rg$hi) else unname(own); if (!is.finite(lims[1]) || lims[2] <= lims[1]) next
            xg <- seq(lims[1], lims[2], length.out = G)
            for (x in xg) { w <- order(abs(xs - x))[seq_len(min(NWIN, length(xs)))]; row <- Xref; row[cols] <- colMeans(X[w, cols, drop = FALSE]); eta <- sum(row * cf)
              sv <- sref; sv[stat] <- x; z <- (sv - Sc) / Ss; nn <- min(sqrt(colSums((t(Ssub) - z[colnames(Ssub)])^2)))
              out[[length(out) + 1]] <- data.frame(community = cm, pool = POOL, transform = tf, position_model = pm, order = k, response = rn, predictor = v, axis = ax, range_type = rt_, x = x, link = eta, dissimilarity = plogis(eta), feasible = nn <= thr, nn_distance = nn, nn_threshold = thr, window_halfwidth = max(abs(xs[w] - x)), in_common_range = x >= rg$lo - 1e-12 && x <= rg$hi + 1e-12, link_reference_row = eta_ref, D2_train = fit$D2, converged = isTRUE(fit$converged), n_eval_sets = length(kk)) } } } } } }
  bind_rows(out) }
t0 <- Sys.time(); R <- mclapply(jobs, run, mc.cores = NC, mc.preschedule = FALSE); cat("elapsed", format(Sys.time() - t0), "\n")
bad <- which(!sapply(R, is.data.frame)); if (length(bad)) { print(R[bad]); stop("failed jobs") }
Gd <- bind_rows(R); write.csv(Gd, file.path(OUT, sprintf("sens_curves_grid_%s.csv", POOL)), row.names = FALSE)
# common/own range table
rr <- bind_rows(lapply(names(rng), function(n) { p <- strsplit(n, " ")[[1]]; cc <- p[1]; data.frame(community = COMM[[cc]], pool = POOL, transform = p[2], statistic = p[3], common_lo = rng[[n]]$lo, common_hi = rng[[n]]$hi, order = ords(cc), own_lo = unlist(rng[[n]][grepl("per_order_lo", names(rng[[n]]))]), own_hi = unlist(rng[[n]][grepl("per_order_hi", names(rng[[n]]))])) }))
write.csv(rr, file.path(OUT, sprintf("sens_curves_ranges_%s.csv", POOL)), row.names = FALSE)
# summary: link-scale and response-scale change across the common feasible range, plus closeness to the ceiling
Sm <- Gd |> filter(range_type == "common", feasible) |> group_by(community, pool, transform, position_model, order, response, predictor, axis) |> arrange(x, .by_group = TRUE) |>
  summarise(n_feasible = n(), x_from = first(x), x_to = last(x), link_from = first(link), link_to = last(link), link_change = last(link) - first(link), link_slope_per_unit = (last(link) - first(link)) / max(last(x) - first(x), 1e-12),
            response_from = first(dissimilarity), response_to = last(dissimilarity), response_change = last(dissimilarity) - first(dissimilarity), mean_dissimilarity = mean(dissimilarity), .groups = "drop")
write.csv(Sm, file.path(OUT, sprintf("sens_curves_summary_%s.csv", POOL)), row.names = FALSE); cat("done", nrow(Gd), "grid rows\n")
