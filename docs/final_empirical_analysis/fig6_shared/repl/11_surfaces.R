# 11_surfaces.R: all-site fits (nonlinear additive and nonlinear varying-spread) for each component and informative order, evaluated on independent
# sets of all sites. Predictions and observed means are averaged in cells of mean TP position x TP spread (mean pairwise log10 TP difference);
# a cell is displayed only if >= 8 training sets, >= 12 training sites, >= 8 evaluation sets are in it. Trophic positions: TP 10, 35, 100 ug/L (+/-0.15 log10).
source("repl/lib_repl.R"); library(parallel); NC <- as.integer(Sys.getenv("NCORES", "9"))
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10, 12), ZP = c(2, 3, 4, 5)); NB <- 12; NS <- 8
cellstat <- function(il, tp) { pos <- rowMeans(matrix(tp[il], nrow(il))); k <- ncol(il); sp <- rep(0, nrow(il)); for (a in 1:(k - 1)) for (b in (a + 1):k) sp <- sp + abs(tp[il[, a]] - tp[il[, b]]); list(pos = pos, spread = sp / choose(k, 2)) }
jobs <- expand.grid(g = names(D), stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame()
run1 <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; e <- get_env(d); n <- nrow(d$x); k <- j$o; DM <- hav(get_xy(d)); tp <- e$TP
  il <- sample_sets(n, k, 6000, seed = 7000 + k); el <- sample_sets(n, k, 8000, seed = 9000 + k); ct <- cellstat(il, tp); ce <- cellstat(el, tp)
  pb <- seq(min(ce$pos), max(ce$pos), length.out = NB + 1); sb <- seq(min(ce$spread), max(ce$spread), length.out = NS + 1); cut2 <- function(x, b) pmin(pmax(findInterval(x, b, rightmost.closed = TRUE), 1), length(b) - 1)
  ct$cp <- cut2(ct$pos, pb); ct$cs <- cut2(ct$spread, sb); ce$cp <- cut2(ce$pos, pb); ce$cs <- cut2(ce$spread, sb); out <- list()
  for (comp in names(COMPS)) { r <- ugdm_response(d$x, il, comp); rt <- ugdm_response(d$x, el, comp)
    for (m in c("nonlinear_const", "nonlinear_varying")) { res <- tryCatch({ dsg <- build(e, il, DM, m); fit <- fit1(dsg, r); dte <- build(e, el, DM, m, ref = dsg); p <- ugdm_predict(fit, dte$int); kk <- rt$keep
        cells <- expand.grid(cp = 1:NB, cs = 1:NS); cells$pos <- (pb[cells$cp] + pb[cells$cp + 1]) / 2; cells$spread <- (sb[cells$cs] + sb[cells$cs + 1]) / 2
        st <- t(sapply(seq_len(nrow(cells)), function(q) { it <- which(ct$cp == cells$cp[q] & ct$cs == cells$cs[q] & r$keep); ie <- which(ce$cp == cells$cp[q] & ce$cs == cells$cs[q] & kk); c(length(it), length(unique(il[it, ])), length(ie), if (length(ie)) mean(p[ie]) else NA, if (length(ie)) mean(rt$y[ie]) else NA) }))
        cells$train_sets <- st[, 1]; cells$train_sites <- st[, 2]; cells$eval_sets <- st[, 3]; cells$pred <- st[, 4]; cells$obs <- st[, 5]; cells$supported <- cells$train_sets >= 8 & cells$train_sites >= 12 & cells$eval_sets >= 8
        cells$group <- j$g; cells$order <- k; cells$component <- comp; cells$model <- m; cells$D2_train <- fit$D2; cells }, error = function(e) { message("ERR ", conditionMessage(e)); NULL }); if (!is.null(res)) out[[length(out) + 1]] <- res } }
  bind_rows(out) }
S <- bind_rows(mclapply(seq_len(nrow(jobs)), run1, mc.cores = NC, mc.preschedule = FALSE)); saveRDS(S, "out/repl_surfaces.rds"); cat("cells", nrow(S), "supported", sum(S$supported), "\n")
