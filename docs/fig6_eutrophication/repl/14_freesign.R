# 14_freesign.R: does the non-negative directional constraint of the extended fit cost predictive power for nestedness and gradient?
# Original MS-GDUM (ugdm_design + ugdm_fit, linear position + I-spline spread) with free_sign = FALSE vs TRUE, same spatial blocks and sets.
source("repl/lib_repl.R"); library(parallel); ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10), ZP = c(2, 3, 4, 5, 6), PP = c(2, 3, 4, 5)); fold_of <- readRDS("out/folds.rds"); DMs <- lapply(D, function(d) hav(get_xy(d)))
jobs <- expand.grid(g = names(D), comp = c("nestedness", "gradient", "turnover", "balanced"), stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame()
run1 <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; e <- get_env(d); fo <- fold_of[[j$g]]; DM <- DMs[[j$g]]; k <- j$o; out <- list()
  for (f in sort(unique(fo))) { tr <- which(fo != f); te <- which(fo == f); if (length(te) < k + 2) next; il <- sample_sets(length(tr), k, 3000, seed = 10 + f); el <- sample_sets(length(te), k, 3000, seed = 100 + f)
    r <- ugdm_response(d$x[tr, ], il, j$comp); rt <- ugdm_response(d$x[te, ], el, j$comp); kk <- rt$keep; mu0 <- sum(r$y[r$keep] * r$trials[r$keep]) / sum(r$trials[r$keep]); nd <- zetaGDUM:::.binom_dev(rt$y[kk], rep(mu0, sum(kk)), rt$trials[kk])
    d0 <- ugdm_design(e[tr, ], il, dist_matrix = DM[tr, tr]); d1 <- ugdm_design(e[te, ], el, dist_matrix = DM[te, te], ref = d0)
    for (fs in c(FALSE, TRUE)) { res <- tryCatch({ fit <- ugdm_fit(d0, r$y[r$keep], r$trials[r$keep], free_sign = fs); p <- ugdm_predict(fit, d1); data.frame(group = j$g, order = k, component = j$comp, fold = f, free_sign = fs, D2_heldout = 1 - zetaGDUM:::.binom_dev(rt$y[kk], p[kk], rt$trials[kk]) / nd) }, error = function(e) NULL); if (!is.null(res)) out[[length(out) + 1]] <- res } }
  bind_rows(out) }
R <- bind_rows(mclapply(seq_len(nrow(jobs)), run1, mc.cores = 9, mc.preschedule = FALSE)); saveRDS(R, "out/repl_freesign.rds")
s <- R |> group_by(group, order, component, free_sign) |> summarise(D2 = median(D2_heldout), .groups = "drop") |> pivot_wider(names_from = free_sign, values_from = D2, names_prefix = "free_") |> mutate(gain = free_TRUE - free_FALSE); write.csv(s, "tables_repl/free_sign_sensitivity.csv", row.names = FALSE)
options(width = 200); print(s |> group_by(group, component) |> summarise(median_gain = median(gain), max_gain = max(gain), .groups = "drop") |> as.data.frame(), digits = 3)
print(s |> filter(component %in% c("nestedness", "gradient")) |> mutate(across(where(is.numeric), ~ round(., 3))) |> as.data.frame())
