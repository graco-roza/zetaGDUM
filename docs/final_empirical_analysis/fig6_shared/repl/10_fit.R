# 10_fit.R: replacement-oriented components x models x orders x spatial blocks. Held-out D2 (three models; interaction penalty sensitivity) and
# leave-one-block-out within-order importance (nonlinear additive and nonlinear varying-spread models). Same data, predictors, folds and orders as Figure 6.
source("repl/lib_repl.R"); library(parallel); NC <- as.integer(Sys.getenv("NCORES", "9"))
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10, 12), ZP = c(2, 3, 4, 5)); M_TR <- 3000L; M_EV <- 3000L; SEEDS <- 1:3; fold_of <- readRDS("out/folds.rds")
DMs <- lapply(D, function(d) hav(get_xy(d)))
jobs <- expand.grid(g = names(D), seed = SEEDS, stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame()
pen_alt <- modifyList(PEN, list(interaction_shrink = 1e-4))
run1 <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; e <- get_env(d); fo <- fold_of[[j$g]]; DM <- DMs[[j$g]]; k <- j$o; V <- list(); IM <- list()
  for (f in sort(unique(fo))) { tr <- which(fo != f); te <- which(fo == f); if (length(te) < k + 2) next
    il <- sample_sets(length(tr), k, M_TR, seed = 10 * j$seed + f); el <- sample_sets(length(te), k, M_EV, seed = 100 * j$seed + f)
    for (comp in names(COMPS)) { r <- ugdm_response(d$x[tr, , drop = FALSE], il, comp); rt <- ugdm_response(d$x[te, , drop = FALSE], el, comp); kk <- rt$keep; if (sum(r$keep) < 50 || sum(kk) < 50) next
      mu0 <- sum(r$y[r$keep] * r$trials[r$keep]) / sum(r$trials[r$keep]); nd <- zetaGDUM:::.binom_dev(rt$y[kk], rep(mu0, sum(kk)), rt$trials[kk])
      for (m in c(MODELS$model, "nonlinear_varying_pen1e-4")) { mm <- sub("_pen1e-4", "", m); pn <- if (m == "nonlinear_varying_pen1e-4") pen_alt else PEN
        res <- tryCatch({ dsg <- build(e[tr, , drop = FALSE], il, DM[tr, tr], mm); fit <- fit1(dsg, r, pn); dte <- build(e[te, , drop = FALSE], el, DM[te, te], mm, ref = dsg); p <- ugdm_predict(fit, dte$int)
            data.frame(group = j$g, order = k, seed = j$seed, fold = f, component = comp, model = m, D2_train = fit$D2, D2_heldout = 1 - zetaGDUM:::.binom_dev(rt$y[kk], p[kk], rt$trials[kk]) / nd, converged = isTRUE(fit$converged), mean_y = mean(rt$y[kk]), sd_y = sd(rt$y[kk]), n_sets = sum(kk))
          }, error = function(e) NULL); if (!is.null(res)) V[[length(V) + 1]] <- res
        if (m %in% c("nonlinear_const", "nonlinear_varying") && !is.null(res)) { im <- tryCatch(ugdm_within_order_importance(fit, dsg$int, r$y[r$keep], r$trials[r$keep], evaluation = dte$int, evaluation_y = rt$y[kk], evaluation_trials = rt$trials[kk], seed = j$seed, failure = "warn"), error = function(e) { message("ERR imp ", conditionMessage(e)); NULL })
          if (!is.null(im)) { h <- im$hierarchy; h$group <- j$g; h$order <- k; h$seed <- j$seed; h$fold <- f; h$comp <- comp; h$model <- m; h$D2_heldout <- res$D2_heldout; IM[[length(IM) + 1]] <- h } } } } }
  list(v = bind_rows(V), im = bind_rows(IM)) }
t0 <- Sys.time(); R <- mclapply(seq_len(nrow(jobs)), run1, mc.cores = NC, mc.preschedule = FALSE); cat("elapsed", format(Sys.time() - t0), "jobs", nrow(jobs), "errors", sum(sapply(R, function(x) inherits(x, "try-error"))), "\n")
saveRDS(bind_rows(lapply(R, `[[`, "v")), "out/repl_validation.rds"); saveRDS(bind_rows(lapply(R, `[[`, "im")), "out/repl_importance.rds")
