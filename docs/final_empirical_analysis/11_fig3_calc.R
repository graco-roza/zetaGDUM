# 11_fig3_calc.R: Figure 3, species and taxon representation. EXACT analytical calculations (expected_species_contributions(), uniform sets of k distinct sites) for incidence and relative abundance:
#   zeta profiles; taxon contributions to similarity (shared composition) and dissimilarity; contribution-weighted occupancy and local abundance; effective number of contributing taxa;
#   occupancy-associated (D_O) and abundance-heterogeneity (D_H) components, 1 = zeta^A_k + D_O + D_H.
# Uncertainty: 100 site-resampling (bootstrap) draws of the exact quantities (seed 20261201 + draw). Null: 499 within-site abundance permutations (abundances shuffled among the occupied taxa of each site; occupancy fixed; seed 20261126 + b).
# Orders: all grid orders with sharing probability >= 0.10 (tables/evaluated_orders.csv). Usage: POOL=shared|full Rscript 11_fig3_calc.R
source("common.R"); library(parallel); NC <- as.integer(Sys.getenv("NCORES", "9")); POOL <- Sys.getenv("POOL", "shared"); NB <- 100L; NP <- 499L
O <- read.csv("tables/evaluated_orders.csv"); ord_of <- function(cc) { o <- O[O$community == COMM[[cc]] & O$pool == (if (cc == "BCI") "bci" else POOL) & O$informative, "order"]; as.integer(o) }
tabs <- function(x, k) { p <- colMeans(x > 0); R <- expected_species_contributions(x, k, "relative"); I <- expected_species_contributions(x, k, "presence"); a <- ifelse(p > 0, R$marginal_mean / p, 0)
  data.frame(taxon = R$species, occupancy = p, local_abundance = a, P_share = R$sharing_probability, sim_rel = R$expected_minimum, dis_rel = R$dissimilarity, sim_inc = I$expected_minimum, dis_inc = I$dissimilarity, marg_rel = R$marginal_mean) }
summ <- function(t, k) { w <- function(v) v / sum(v); Dh <- sum(t$local_abundance * t$P_share - t$sim_rel); Do <- sum(t$marg_rel - t$local_abundance * t$P_share)
  out <- data.frame(order = k, zeta_inc = sum(t$sim_inc), zeta_rel = sum(t$sim_rel), D_O = Do, D_H = Dh, closure_error = 1 - sum(t$sim_rel) - Do - Dh)
  for (q in c("sim_inc", "dis_inc", "sim_rel", "dis_rel")) { ww <- w(t[[q]]); out[[paste0("neff_", q)]] <- 1 / sum(ww^2); out[[paste0("occ_w_", q)]] <- sum(ww * t$occupancy); out[[paste0("abund_w_", q)]] <- sum(ww * t$local_abundance) }; out }
run_all <- function(x, ords) bind_rows(lapply(ords, function(k) summ(tabs(x, k), k)))
res <- list(); taxa_tab <- list(); boot <- list(); null <- list()
for (cc in names(COMM)) { x <- get_comm(cc, POOL); x <- x[, colSums(x > 0) > 0, drop = FALSE]; ords <- ord_of(cc); cat(cc, nrow(x), "x", ncol(x), "orders", paste(ords, collapse = ","), "\n")
  obs <- run_all(x, ords); obs$community <- COMM[[cc]]; res[[cc]] <- obs; stopifnot(max(abs(obs$closure_error)) < 1e-10)
  tt <- bind_rows(lapply(ords, function(k) { t <- tabs(x, k); t$order <- k; t$community <- COMM[[cc]]; t })); taxa_tab[[cc]] <- tt
  set.seed(20261201); bs <- mclapply(seq_len(NB), function(b) { set.seed(20261201 + b); i <- sample.int(nrow(x), replace = TRUE); xb <- x[i, , drop = FALSE]; r <- run_all(xb, ords); r$draw <- b; r }, mc.cores = NC); boot[[cc]] <- bind_rows(bs) |> mutate(community = COMM[[cc]])
  q <- x / rowSums(x); occ <- lapply(seq_len(nrow(x)), function(i) which(x[i, ] > 0))
  nl <- mclapply(0:NP, function(b) { y <- q; if (b > 0) { set.seed(20261126 + b); for (i in seq_len(nrow(y))) { j <- occ[[i]]; if (length(j) > 1) y[i, j] <- q[i, j[sample.int(length(j))]] } }
    r <- bind_rows(lapply(ords, function(k) { R <- expected_species_contributions(y, k, "relative"); p <- colMeans(y > 0); a <- ifelse(p > 0, R$marginal_mean / p, 0); t <- data.frame(taxon = R$species, occupancy = p, local_abundance = a, P_share = R$sharing_probability, sim_rel = R$expected_minimum, dis_rel = R$dissimilarity, sim_inc = 0, dis_inc = 0, marg_rel = R$marginal_mean); s <- summ(t, k); s[, c("order", "zeta_rel", "D_O", "D_H", grep("_rel$", names(s), value = TRUE))] }))
    r$permutation <- b; r }, mc.cores = NC); null[[cc]] <- bind_rows(nl) |> mutate(community = COMM[[cc]]) }
write.csv(bind_rows(res), sprintf("tables/fig3_observed_%s.csv", POOL), row.names = FALSE); write.csv(bind_rows(taxa_tab), sprintf("tables/fig3_taxon_contributions_%s.csv", POOL), row.names = FALSE)
write.csv(bind_rows(boot), sprintf("tables/fig3_bootstrap_draws_%s.csv", POOL), row.names = FALSE); write.csv(bind_rows(null), sprintf("tables/fig3_permutation_draws_%s.csv", POOL), row.names = FALSE); cat("done\n")
