# 62_tables.R: paired interaction gains and shape summaries from 60/61 (no refitting). Run from final/: Rscript fig6_interaction/62_tables.R
#   delta_D2_train = D2_train_interaction - D2_train_additive ; delta_D2_test = D2_test_interaction - D2_test_additive (same sampled sets, folds, seeds)
#   Seed means within fold first, then median / min / max across the five spatial folds. All fits are kept (non-converged fits flagged and counted, none removed).
#   TP-only vs all-predictor interaction contrast at the same penalty (delta of the two interaction models).
suppressMessages({library(dplyr); library(tidyr)}); OUT <- "fig6_interaction/tables"
for (pool in c("shared", "full")) {
  M <- read.csv(file.path(OUT, sprintf("interaction_metrics_%s.csv", pool))); key <- c("group", "pool", "order", "seed", "fold", "component", "split")
  stopifnot(!anyDuplicated(M[c(key, "model")]))
  A <- M |> filter(model == "additive") |> select(all_of(key), D2_add = D2, conv_add = converged)
  I <- M |> filter(model != "additive") |> select(all_of(key), model, variant, penalty, D2_int = D2, conv_int = converged) |> inner_join(A, by = key)
  stopifnot(nrow(I) == sum(M$model != "additive"))
  P <- I |> mutate(delta_D2 = D2_int - D2_add, any_nonconverged = !(conv_add & conv_int))
  # seed mean within fold
  Pf <- P |> group_by(group, pool, order, fold, component, split, model, variant, penalty) |> summarise(delta_D2 = mean(delta_D2), D2_int = mean(D2_int), D2_add = mean(D2_add), n_seeds = n(), n_nonconverged = sum(any_nonconverged), .groups = "drop")
  write.csv(Pf, file.path(OUT, sprintf("interaction_paired_folds_%s.csv", pool)), row.names = FALSE)
  Ps <- Pf |> group_by(group, pool, order, component, split, model, variant, penalty) |> summarise(median = median(delta_D2), min = min(delta_D2), max = max(delta_D2), mean = mean(delta_D2), n_folds = n(), n_folds_positive = sum(delta_D2 > 0), n_nonconverged_pairs = sum(n_nonconverged), .groups = "drop")
  write.csv(Ps, file.path(OUT, sprintf("interaction_paired_summary_%s.csv", pool)), row.names = FALSE)
  # TP-only vs all-predictor interaction at matched penalty
  T <- M |> filter(variant %in% c("interaction_all", "interaction_TP")) |> select(all_of(key), variant, penalty, D2) |> pivot_wider(names_from = variant, values_from = D2) |> mutate(delta_TPonly_minus_all = interaction_TP - interaction_all) |>
    group_by(group, pool, order, fold, component, split, penalty) |> summarise(delta = mean(delta_TPonly_minus_all), .groups = "drop") |> group_by(group, pool, order, component, split, penalty) |> summarise(median = median(delta), min = min(delta), max = max(delta), n_folds = n(), n_folds_positive = sum(delta > 0), .groups = "drop")
  write.csv(T, file.path(OUT, sprintf("interaction_tp_only_vs_all_%s.csv", pool)), row.names = FALSE)
  Cv <- M |> distinct(group, pool, order, seed, fold, model, variant, penalty, converged, share_identified) |> group_by(group, pool, order, model, variant, penalty) |> summarise(n_fits = n(), n_converged = sum(converged), n_nonconverged = sum(!converged), n_share_not_identified = sum(!share_identified), .groups = "drop")
  write.csv(Cv, file.path(OUT, sprintf("interaction_convergence_%s.csv", pool)), row.names = FALSE)
  # shape magnitude: explicit interaction contrast (High - Low position) on the model scale
  S <- read.csv(file.path(OUT, sprintf("interaction_spline_values_%s.csv", pool))) |> filter(position %in% c("Low", "High"))
  W <- S |> select(group, pool, order, model, variant, penalty, stage, position, TP_log10, f_centred) |> pivot_wider(names_from = position, values_from = f_centred) |> mutate(contrast_high_minus_low = High - Low) |>
    group_by(group, pool, order, model, variant, penalty, stage) |> summarise(max_abs_contrast = max(abs(contrast_high_minus_low)), end_contrast = last(contrast_high_minus_low), range_low = diff(range(Low)), range_high = diff(range(High)), .groups = "drop") |> mutate(relative_contrast = max_abs_contrast / pmax(range_low, range_high, 1e-12))
  write.csv(W, file.path(OUT, sprintf("interaction_shape_magnitude_%s.csv", pool)), row.names = FALSE)
  cat("== ", pool, "\n"); print(Ps |> filter(split == "test", component == "total", variant == "interaction_all") |> group_by(group, penalty) |> summarise(median_of_medians = median(median), orders_positive = sum(median > 0), n = n(), .groups = "drop") |> as.data.frame(), digits = 3)
  print(Ps |> filter(split == "train", component == "total") |> group_by(group, variant, penalty) |> summarise(median_train_gain = median(median), .groups = "drop") |> as.data.frame(), digits = 3)
}
