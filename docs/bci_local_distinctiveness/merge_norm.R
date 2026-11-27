# merge_norm.R: add the richness-normalised incidence LCZD (09b_lczd_normalised.R output) to the site table as LCZD_norm_<k> / ED_norm_<k> (tree-record sets only)
add_norm <- function(S) { fs <- list.files("out", "^lczd_normalised_.*csv$", full.names = TRUE); if (!length(fs)) return(S); N <- dplyr::bind_rows(lapply(fs, read.csv)); key <- paste(S$set, S$quad)
  for (k in sort(unique(N$order))) { Nk <- N[N$order == k, ]; m <- match(key, paste(Nk$set, Nk$quad)); S[[paste0("LCZD_norm_", k)]] <- Nk$LCZD_norm[m]; S[[paste0("ED_norm_", k)]] <- Nk$ED_norm[m] }; S }
NORM_ORDERS <- c(2, 3, 4, 5, 8, 12, 20, 40, 80)
