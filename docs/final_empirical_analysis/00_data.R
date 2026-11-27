# 00_data.R: build the analysis pools. Shared EPA lake pool = lakes present (with complete TP, TN, Chl-a, coordinates and non-empty samples) in BOTH the macroinvertebrate (LUSBM_2) and zooplankton (LUSZP_2) tables; full pools = all valid lakes of each taxon.
# Inputs: ../fig6/out/prepared.rds (from ../fig6/01_prep.R: US NLA 2017 workbooks of the Liu Guohao eutrophication synthesis), ../bci/out/bci_matrices.rds (BCI 1985 tree records, ../bci/01_prepare.R).
suppressMessages(library(dplyr)); P <- readRDS("../fig6/out/prepared.rds"); B <- readRDS("../bci/out/bci_matrices.rds")
ids <- intersect(P$BM$env$id, P$ZP$env$id); cat("shared lakes:", length(ids), "\n")
mk <- function(d, keep = NULL) { e <- d$env; x <- d$x; if (!is.null(keep)) { i <- match(keep, e$id); e <- e[i, ]; x <- x[i, , drop = FALSE] }; x <- x[, colSums(x > 0) > 0, drop = FALSE]; stopifnot(all(rowSums(x) > 0), identical(rownames(x), e$id)); list(x = x, env = e) }
shared <- list(BM = mk(P$BM, ids), ZP = mk(P$ZP, ids)); stopifnot(identical(shared$BM$env$id, shared$ZP$env$id), all(shared$BM$env[, c("TP", "TN", "Chla", "lon", "lat")] == shared$ZP$env[, c("TP", "TN", "Chla", "lon", "lat")]))
full <- list(BM = mk(P$BM), ZP = mk(P$ZP)); bx <- B$mat[["1985_tree"]]; bx <- bx[, colSums(bx) > 0]
saveRDS(list(shared = shared, full = full, bci = list(x = bx, habitat = B$habitat)), "out/pools.rds")
# dataset characteristics
rows <- list(); add <- function(name, pool, x, env = NULL) rows[[length(rows) + 1]] <<- data.frame(community = name, pool = pool, sites = nrow(x), taxa = ncol(x), mean_richness = mean(rowSums(x > 0)), median_total_abundance = median(rowSums(x)), proportion_singleton_taxa = mean(colSums(x > 0) == 1),
  TP_ug_L_min = if (is.null(env)) NA else min(env$TP), TP_ug_L_median = if (is.null(env)) NA else median(env$TP), TP_ug_L_max = if (is.null(env)) NA else max(env$TP), log10TP_range = if (is.null(env)) NA else diff(range(log10(env$TP))), cor_logTP_logTN = if (is.null(env)) NA else cor(log10(env$TP), log10(env$TN)), cor_logTP_logChla = if (is.null(env)) NA else cor(log10(env$TP), log10(env$Chla)), cor_logTN_logChla = if (is.null(env)) NA else cor(log10(env$TN), log10(env$Chla)))
add("BCI trees (1985)", "1250 20-m quadrats", bx); for (p in c("shared", "full")) for (cc in c("BM", "ZP")) add(c(BM = "EPA macroinvertebrates", ZP = "EPA zooplankton")[[cc]], p, get(p)[[cc]]$x, get(p)[[cc]]$env)
D <- bind_rows(rows); dir.create("tables", FALSE); write.csv(D, "tables/dataset_characteristics.csv", row.names = FALSE); print(t(D))
