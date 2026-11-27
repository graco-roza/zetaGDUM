# 03b_penalty.R: sensitivity of the explicit interaction to the endpoint penalty (default 0.01; package docs report 1e-4 sensitivity). Scored on held-out spatial blocks (matched sets).
source("lib6.R"); library(parallel); NC <- as.integer(Sys.getenv("NCORES", "4"))
ORD <- list(BM = c(2, 3, 4, 5, 6, 8, 10, 12), ZP = c(2, 3, 4, 5)); PEN <- c(1e-2, 1e-4, 1e-6); REPS <- 8
jobs <- expand.grid(g = names(D), pen = PEN, rep = seq_len(REPS), stringsAsFactors = FALSE) |> rowwise() |> mutate(o = list(ORD[[g]])) |> tidyr::unnest(o) |> as.data.frame()
run1 <- function(i) { j <- jobs[i, ]; d <- D[[j$g]]; r <- tryCatch(ugdm_validate_components(d$x, get_env(d), get_xy(d), order = j$o, blocks = 5, models = c("additive", "interaction"), architectures = "coherent", combinations = 1500L, interaction_shrink = j$pen, geographic = TRUE, lonlat = TRUE, seed = 1000 * j$rep + j$o), error = function(e) NULL)
  if (is.null(r)) return(NULL); r$group <- j$g; r$pen <- j$pen; r$rep <- j$rep; r }
res <- bind_rows(mclapply(seq_len(nrow(jobs)), run1, mc.cores = NC, mc.preschedule = FALSE)); saveRDS(res, "out/penalty_folds.rds")
w <- res |> select(group, order, component, pen, rep, fold, model, D2, supported_D2) |> pivot_wider(names_from = model, values_from = c(D2, supported_D2)) |> mutate(gain = D2_interaction - D2_additive, gain_sup = supported_D2_interaction - supported_D2_additive)
print(w |> group_by(group, order, component, pen) |> summarise(gain = median(gain, na.rm = TRUE), gain_supported = median(gain_sup, na.rm = TRUE), pos = mean(gain > 0, na.rm = TRUE), .groups = "drop") |> as.data.frame(), digits = 3)
