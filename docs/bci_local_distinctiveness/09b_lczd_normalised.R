# 09b_lczd_normalised.R: richness-normalised incidence LCZD (zeta_lczd(type = "normalised")), the incidence counterpart that is standardised by mean set richness like the
# occurrence LCBD and occurrence T. It is estimated from sampled sets (exact at order 2: all 780,625 pairs). Two independent replicate samples per order give a precision check;
# large sample sizes are processed in chunks of 1e6 sets whose expected-dissimilarity estimates are averaged. Usage: CENSUS=1985_tree Rscript 09b_lczd_normalised.R
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(parallel)}); nm <- Sys.getenv("CENSUS", "1985_tree"); NC <- as.integer(Sys.getenv("NCORES", "4"))
B <- readRDS("out/bci_matrices.rds"); x <- B$mat[[nm]]; x <- x[, colSums(x) > 0]; n <- nrow(x); CH <- 1e6
ORD <- c(2, 3, 4, 5, 8, 12, 20, 40, 80); MS <- c(`2` = choose(n, 2), `3` = 1e6, `4` = 1e6, `5` = 2e6, `8` = 2e6, `12` = 3e6, `20` = 4e6, `40` = 8e6, `80` = 2e7)
ed1 <- function(k, rep) { M <- MS[[as.character(k)]]; if (k == 2) { z <- zeta_lczd(x, 2, "normalised", idx = t(combn(n, 2))); return(z$expected_dissimilarity) }
  nch <- ceiling(M / CH); acc <- numeric(n); for (c in seq_len(nch)) { idx <- sample_sets(n, k, M = CH, seed = 1000 * rep + 10 * c + k); acc <- acc + zeta_lczd(x, k, "normalised", idx = idx)$expected_dissimilarity }; acc / nch }
jobs <- expand.grid(k = ORD, rep = 1:2); jobs <- jobs[order(-jobs$k), ]; t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(i) { if (jobs$k[i] == 2 && jobs$rep[i] == 2) return(NULL); ed1(jobs$k[i], jobs$rep[i]) }, mc.cores = NC, mc.preschedule = FALSE); cat(nm, "elapsed", format(Sys.time() - t0), "\n")
out <- list(); for (k in ORD) { a <- res[[which(jobs$k == k & jobs$rep == 1)]]; b <- if (k == 2) a else res[[which(jobs$k == k & jobs$rep == 2)]]; ed <- (a + b) / 2
  out[[as.character(k)]] <- data.frame(set = nm, order = k, quad = rownames(x), ED_norm = ed, LCZD_norm = ed / sum(ed), ED_rep1 = a, ED_rep2 = b, sets_per_replicate = MS[[as.character(k)]], seed_rho = cor(a, b, method = "spearman")) }
O <- bind_rows(out); write.csv(O, sprintf("out/lczd_normalised_%s.csv", nm), row.names = FALSE); print(O |> group_by(order) |> summarise(seed_rho = first(seed_rho), sets = first(sets_per_replicate)))
