# 01_prep.R: load the three selected lake communities, harmonise, write analysis-ready RDS + summary table
suppressMessages({library(readxl); library(dplyr)})
root <- "/Users/graco-roza/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofHelsinki(2)/Liu, Guohao - Nature_R01/raw_data_new/By_groups"
SEL <- c(BM = "BM_Benthic Macroinvertebrates/LUSBM_2", ZP = "ZP_Zooplankton/LUSZP_2", PP = "PP_Phytoplankton/LPOPP")
LAB <- c(BM = "Macroinvertebrates (US lakes)", ZP = "Zooplankton (US lakes)", PP = "Phytoplankton (Polish lakes)")
load1 <- function(code) {
  f <- file.path(root, paste0(SEL[[code]], ".xlsx")); sp <- suppressMessages(read_excel(f, "species")); en <- suppressMessages(read_excel(f, "environment")); co <- suppressMessages(read_excel(f, "coordinates"))
  ids <- as.character(sp[[1]]); x <- as.matrix(sp[, -1]); suppressWarnings(storage.mode(x) <- "numeric"); x[is.na(x)] <- 0; rownames(x) <- ids
  e <- data.frame(id = as.character(en[[1]]), TP = as.numeric(en$TP), TN = as.numeric(en$TN), Chla = as.numeric(en$Chla)); xy <- data.frame(id = as.character(co[[1]]), lon = as.numeric(co$X), lat = as.numeric(co$Y))
  d <- data.frame(id = ids) |> left_join(e, by = "id") |> left_join(xy, by = "id")
  ok <- complete.cases(d[, c("TP", "TN", "Chla", "lon", "lat")]) & d$TP > 0 & d$TN > 0 & d$Chla > 0 & rowSums(x) > 0 & !duplicated(d$id)
  x <- x[ok, colSums(x[ok, ] > 0) > 0, drop = FALSE]; d <- d[ok, ]
  list(code = code, label = LAB[[code]], x = x, env = d, n_raw = length(ids), n_dup = sum(duplicated(ids)), n_dropped = sum(!ok))
}
D <- lapply(names(SEL), load1); names(D) <- names(SEL); dir.create("out", FALSE); saveRDS(D, "out/prepared.rds")
s <- bind_rows(lapply(D, function(d) with(d, data.frame(group = code, label = label, file = SEL[[code]], sites_raw = n_raw, sites_used = nrow(x), dropped = n_dropped, species = ncol(x), mean_richness = mean(rowSums(x > 0)),
  median_total = median(rowSums(x)), prop_species_singletons = mean(colSums(x > 0) == 1), TP_min = min(env$TP), TP_q10 = quantile(env$TP, .1), TP_median = median(env$TP), TP_q90 = quantile(env$TP, .9), TP_max = max(env$TP),
  logTP_range = diff(range(log10(env$TP))), Chla_median = median(env$Chla), TN_median = median(env$TN), lon_range = diff(range(env$lon)), lat_range = diff(range(env$lat)),
  cor_logTP_logTN = cor(log10(env$TP), log10(env$TN)), cor_logTP_logChla = cor(log10(env$TP), log10(env$Chla)), cor_logTN_logChla = cor(log10(env$TN), log10(env$Chla))))))
write.csv(s, "out/selected_summary.csv", row.names = FALSE); options(width = 220); print(t(s), quote = FALSE)
cat("\nshared site IDs BM & ZP:", length(intersect(D$BM$env$id, D$ZP$env$id)), "\n"); cat("head ids:", head(D$BM$env$id, 3), "|", head(D$ZP$env$id, 3), "|", head(D$PP$env$id, 3), "\n")
cat("coords ranges\n"); print(lapply(D, function(d) range(d$env$lon))); print(lapply(D, function(d) range(d$env$lat)))
cat("abundance sample (BM, ZP, PP):\n"); for (d in D) print(summary(d$x[d$x > 0]))
