# 03_validate.R: independent spatial-block validation of the three models (spread only / additive / explicit interaction), both architectures,
# repeated combination sampling (R seeds), per group and order. Primary predictors: log10 TP, TN, Chla (+ geographic distance); sensitivity: without Chla.
source("lib6.R"); library(parallel); NC <- as.integer(Sys.getenv("NCORES", "4"))
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10), ZP = c(2, 3, 4, 5, 6), PP = c(2, 3, 4, 5)); REPS <- 8; NBLK <- 5; COMB <- 1500L
jobs <- expand.grid(g = names(D), chla = c(TRUE, FALSE), rep = seq_len(REPS), stringsAsFactors = FALSE) |> rowwise() |> mutate(orders = list(ORD[[g]])) |> tidyr::unnest(orders) |> rename(order = orders) |> as.data.frame()
run1 <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; t0 <- Sys.time()
  r <- tryCatch(ugdm_validate_components(d$x, get_env(d, j$chla), get_xy(d), order = j$order, blocks = NBLK, combinations = COMB, geographic = TRUE, lonlat = TRUE, seed = 1000 * j$rep + j$order), error = function(e) { message("ERR ", i, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL); r$group <- j$g; r$chla <- j$chla; r$rep <- j$rep; r$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")); r }
t0 <- Sys.time(); res <- bind_rows(mclapply(seq_len(nrow(jobs)), run1, mc.cores = NC, mc.preschedule = FALSE)); cat("elapsed", format(Sys.time() - t0), "jobs", nrow(jobs), "\n")
saveRDS(res, "out/validation_folds.rds"); write.csv(res, "out/validation_folds.csv", row.names = FALSE)
print(res |> filter(chla, architecture == "coherent") |> group_by(group, order, component, model) |> summarise(D2 = median(D2, na.rm = TRUE), frac = median(supported_fraction), conv = mean(converged), .groups = "drop") |> pivot_wider(names_from = model, values_from = D2) |> as.data.frame(), digits = 3)
