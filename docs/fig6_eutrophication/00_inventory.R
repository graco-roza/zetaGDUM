# 00_inventory.R: profile every lake workbook (no model fitting; selection is by data attributes only)
suppressMessages(library(readxl)); suppressMessages(library(dplyr))
root <- "/Users/graco-roza/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofHelsinki(2)/Liu, Guohao - Nature_R01/raw_data_new/By_groups"
fs <- list.files(root, "^L.*xlsx$", recursive = TRUE, full.names = TRUE)
out <- lapply(fs, function(f) {
  sp <- suppressMessages(read_excel(f, "species")); en <- suppressMessages(read_excel(f, "environment")); co <- suppressMessages(read_excel(f, "coordinates"))
  nm <- names(en); x <- as.matrix(sp[, -1]); suppressWarnings(storage.mode(x) <- "numeric")
  ids <- intersect(intersect(as.character(sp[[1]]), as.character(en[[1]])), as.character(co[[1]]))
  g <- function(v) if (v %in% nm) { z <- suppressWarnings(as.numeric(en[[v]])); c(n = sum(is.finite(z)), lo = min(z, na.rm = TRUE), hi = max(z, na.rm = TRUE)) } else c(n = NA, lo = NA, hi = NA)
  tp <- g("TP"); tpr <- g("TP_raw"); tn <- g("TN"); ch <- g("Chla")
  xs <- x[as.character(sp[[1]]) %in% ids, , drop = FALSE]; xs[is.na(xs)] <- 0
  data.frame(file = sub(".xlsx", "", basename(f)), sites_species = nrow(sp), sites_matched = length(ids), species = ncol(x), nonbinary = mean(xs[xs > 0] != 1), mean_rich = mean(rowSums(xs > 0)),
             prop_zero_total = mean(rowSums(xs) == 0), n_TP = tp["n"], TP_lo = tp["lo"], TP_hi = tp["hi"], TPraw_lo = tpr["lo"], TPraw_hi = tpr["hi"], has_TN = !is.na(tn["n"]), n_TN = tn["n"], has_Chla = !is.na(ch["n"]), n_Chla = ch["n"], env_cols = paste(nm[-1], collapse = ","))
}) |> bind_rows()
write.csv(out, "out_inventory_all_lake_files.csv", row.names = FALSE); options(width = 250); print(out |> select(-env_cols), digits = 3)
