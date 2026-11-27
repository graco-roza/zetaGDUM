# 52_tables.R: Figure 5 revision, numerical results derived from the fold-level fits of 50_fit_configs.R (no refitting).
# Outputs (fig5_revision/tables/, one file per pool/variant where indicated):
#   sens_paired_folds_<pool>.csv     within-fold paired differences in held-out D2 (seed means first), all configurations
#   sens_paired_summary_<pool>.csv   median, min, max, mean and n of fold differences across folds
#   sens_performance_<pool>.csv      held-out D2 (median, range) and training D2 (mean) for every configuration, order and representation
#   sens_convergence_<pool>.csv      convergence and warning counts per configuration
#   sens_importance_<pool>.csv       absolute contributions (mean across seeds, then folds; fold min/max), relative share without clipping, flags
#   sens_importance_fold_values_<pool>.csv   fold-level absolute contributions (seed means)
#   sens_extremes_distribution.csv   skewness and tail diagnostics of the nutrient predictors (raw vs log10), per pool and taxon table
#   sens_extremes_set_exposure.csv   share of sampled sets containing at least one lake above the 99th percentile of TP, TN or Chl-a
#   sens_extremes_trim_comparison_<pool>.csv  raw/log10 held-out D2 and contributions with and without the extreme lakes (trim99 variant)
source("common.R"); suppressMessages(library(dplyr)); OUT <- "fig5_revision/tables"
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE); med_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
mn <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE); mx <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
MODEL_LAB <- c(spread_only = "Spread only", linear_const = "Spread + linear position", nonlinear_const = "Spread + nonlinear position")
for (suf in c("shared", "full", "shared_trim99", "full_trim99")) { vf <- file.path(OUT, sprintf("sens_validation_folds_%s.csv", suf)); if (!file.exists(vf)) next
  V <- read.csv(vf, stringsAsFactors = FALSE); I <- read.csv(file.path(OUT, sprintf("sens_importance_folds_%s.csv", suf)), stringsAsFactors = FALSE)
  key <- c("community", "pool", "variant", "transform", "order", "response")
  # seed means within fold
  Vs <- V |> group_by(across(all_of(c(key, "fold", "model")))) |> summarise(D2h = mean_na(D2_heldout), D2t = mean_na(D2_train), n_seeds_valid = sum(is.finite(D2_heldout)), .groups = "drop")
  W <- Vs |> select(all_of(key), fold, model, D2h) |> tidyr::pivot_wider(names_from = model, values_from = D2h)
  P <- W |> transmute(across(all_of(key)), fold, delta_linear_vs_spread = linear_const - spread_only, delta_nonlinear_vs_spread = nonlinear_const - spread_only, delta_nonlinear_vs_linear = nonlinear_const - linear_const, D2_spread_only = spread_only, D2_linear = linear_const, D2_nonlinear = nonlinear_const)
  Pl <- P |> tidyr::pivot_longer(starts_with("delta_"), names_to = "contrast", values_to = "delta_D2")
  write.csv(Pl, file.path(OUT, sprintf("sens_paired_folds_%s.csv", suf)), row.names = FALSE)
  Ps <- Pl |> group_by(across(all_of(key)), contrast) |> summarise(median = med_na(delta_D2), min = mn(delta_D2), max = mx(delta_D2), mean = mean_na(delta_D2), n_folds = sum(is.finite(delta_D2)), n_folds_positive = sum(delta_D2 > 0, na.rm = TRUE), .groups = "drop")
  write.csv(Ps, file.path(OUT, sprintf("sens_paired_summary_%s.csv", suf)), row.names = FALSE)
  Pf <- Vs |> group_by(across(all_of(key)), model) |> summarise(heldout_median = med_na(D2h), heldout_min = mn(D2h), heldout_max = mx(D2h), heldout_mean = mean_na(D2h), train_mean = mean_na(D2t), n_folds = sum(is.finite(D2h)), .groups = "drop") |> mutate(model_label = MODEL_LAB[model])
  write.csv(Pf, file.path(OUT, sprintf("sens_performance_%s.csv", suf)), row.names = FALSE)
  Cv <- V |> filter(model != "spread_only") |> group_by(across(all_of(key)), model) |> summarise(n_fits = n(), n_converged = sum(converged %in% TRUE), prop_converged = mean(converged %in% TRUE), n_with_warnings = sum(n_warnings > 0), skew_warning = sum(grepl("skewed", warning_text)), collinearity_warning = sum(grepl("collinearity", warning_text)), n_failed = sum(!is.finite(D2_heldout)), .groups = "drop")
  write.csv(Cv, file.path(OUT, sprintf("sens_convergence_%s.csv", suf)), row.names = FALSE)
  # importance: seed mean within fold, then fold summaries; relative share not clipped
  Ik <- c(key, "model"); If <- I |> group_by(across(all_of(c(Ik, "fold", "predictor", "component")))) |> summarise(absolute = mean_na(absolute), D2h = mean(D2_heldout), .groups = "drop")
  write.csv(If, file.path(OUT, sprintf("sens_importance_fold_values_%s.csv", suf)), row.names = FALSE)
  Ia <- If |> group_by(across(all_of(c(Ik, "predictor", "component")))) |> summarise(abs_mean = mean_na(absolute), abs_median = med_na(absolute), abs_min = mn(absolute), abs_max = mx(absolute), heldout_D2 = mean(D2h), n_folds = sum(!is.na(absolute)), .groups = "drop") |>
    mutate(relative_share = ifelse(heldout_D2 > 0, abs_mean / heldout_D2, NA_real_), relative_flag = case_when(heldout_D2 <= 0 ~ "undefined (held-out D2 <= 0)", heldout_D2 < .02 ~ "unstable (held-out D2 < 0.02)", TRUE ~ "ok"))
  write.csv(Ia, file.path(OUT, sprintf("sens_importance_%s.csv", suf)), row.names = FALSE) }
# extreme-value diagnostics
sk <- function(x) { z <- (x - mean(x)) / sd(x); mean(z^3) }
rows <- list(); rows2 <- list(); for (p in c("shared", "full")) for (cc in c("BM", "ZP")) { e <- POOLS[[p]][[cc]]$env
  for (v in c("TP", "TN", "Chla")) { x <- e[[v]]; q99 <- quantile(x, .99); rows[[length(rows) + 1]] <- data.frame(pool = p, community = COMM[[cc]], predictor = v, n_lakes = length(x), min = min(x), median = median(x), q99 = unname(q99), max = max(x), max_over_median = max(x) / median(x), skewness_raw = sk(x), skewness_log10 = sk(log10(x)), prop_above_q99 = mean(x > q99), share_of_sum_top1pct = sum(x[x >= q99]) / sum(x), prop_z_above_3_raw = mean(abs((x - mean(x)) / sd(x)) > 3), prop_z_above_3_log10 = mean(abs((log10(x) - mean(log10(x))) / sd(log10(x))) > 3)) }
  ext <- with(e, TP > quantile(TP, .99) | TN > quantile(TN, .99) | Chla > quantile(Chla, .99)); O <- read.csv("tables/evaluated_orders.csv"); ks <- O[O$community == COMM[[cc]] & O$pool == p & O$informative, "order"]
  for (k in ks) rows2[[length(rows2) + 1]] <- data.frame(pool = p, community = COMM[[cc]], order = k, n_lakes = nrow(e), n_extreme_lakes = sum(ext), prop_lakes_extreme = mean(ext), expected_share_sets_with_extreme_lake = 1 - (1 - mean(ext))^k) }
write.csv(bind_rows(rows), file.path(OUT, "sens_extremes_distribution.csv"), row.names = FALSE); write.csv(bind_rows(rows2), file.path(OUT, "sens_extremes_set_exposure.csv"), row.names = FALSE)
for (p in c("shared", "full")) { fa <- file.path(OUT, sprintf("sens_performance_%s.csv", p)); ft <- file.path(OUT, sprintf("sens_performance_%s_trim99.csv", p)); if (!all(file.exists(c(fa, ft)))) next
  A <- read.csv(fa) |> mutate(variant = "all"); T <- read.csv(ft) |> mutate(variant = "trim99"); B <- bind_rows(A, T)
  Ba <- read.csv(file.path(OUT, sprintf("sens_importance_%s.csv", p))) |> mutate(variant = "all"); Bt <- read.csv(file.path(OUT, sprintf("sens_importance_%s_trim99.csv", p))) |> mutate(variant = "trim99")
  perf <- B |> select(community, transform, order, response, model, variant, heldout_median, heldout_min, heldout_max, train_mean) |> tidyr::pivot_wider(names_from = variant, values_from = c(heldout_median, heldout_min, heldout_max, train_mean)) |> mutate(change_heldout_median_trim = heldout_median_trim99 - heldout_median_all, change_train_mean_trim = train_mean_trim99 - train_mean_all)
  imp <- bind_rows(Ba, Bt) |> select(community, transform, model, order, response, predictor, component, variant, abs_mean) |> tidyr::pivot_wider(names_from = variant, values_from = abs_mean) |> mutate(change_abs_trim = trim99 - all) |> rename(abs_mean_all = all, abs_mean_trim99 = trim99)
  write.csv(perf, file.path(OUT, sprintf("sens_extremes_trim_comparison_%s.csv", p)), row.names = FALSE); write.csv(imp, file.path(OUT, sprintf("sens_extremes_trim_importance_%s.csv", p)), row.names = FALSE) }
cat("tables written\n")
