# figure2_summary.R: replicate-derived plotting table for Figure 2 (rarity and species representation), order 2.
# Source: results/rarity_forms/<scenario>_<rep>.rds (sim_rarity_forms.R; 20 simulated communities per scenario; weights from zeta_decompose()).
# Per scenario: w_occ / w_ab = share of the shared composition (occurrence / abundance weight) carried by the scenario's dominant form; ratio = effective species number of abundance / of occurrence.
# Median and 25th-75th percentiles across the replicates. Balanced scenario: ratio only; selected forms for panel C. No data are transcribed or rounded.
suppressMessages(library(dplyr)); D <- "../results/rarity_forms"; fs <- list.files(D, full.names = TRUE); L <- lapply(fs, readRDS); cat("replicates:", length(L), "\n")
forms <- c("Lr-Bn-La", "Lr-Bn-Sa", "Lr-Nn-La", "Lr-Nn-Sa", "Sr-Bn-La", "Sr-Bn-Sa", "Sr-Nn-La", "Sr-Nn-Sa", "balanced"); q <- function(v, p) unname(quantile(v, p, na.rm = TRUE))
W <- bind_rows(lapply(L, function(x) { w <- x$wts; w$scenario <- x$scenario; w$rep <- x$rep; w$ratio <- w$eff_abu / w$eff_occ; w[w$order == min(w$order), ] }))   # ratio: effective species number of abundance weights / of occurrence weights; stopifnot(all(W$order == 2))
n_rep <- W |> distinct(scenario, rep) |> count(scenario); print(n_rep)
dom <- W |> filter(scenario != "balanced", form == scenario) |> select(scenario, rep, w_occ, w_abu, ratio); bal <- W |> filter(scenario == "balanced") |> distinct(rep, ratio)
S <- bind_rows(dom |> group_by(scenario) |> summarise(w_occ_lo = q(w_occ, .25), w_occ = median(w_occ), w_occ_hi = NA_real_, w_ab_lo = q(w_abu, .25), w_ab = median(w_abu), w_ab_hi = NA_real_, ratio_lo = q(ratio, .25), ratio = median(ratio), ratio_hi = NA_real_, n_replicates = n(), .groups = "drop"),
  bal |> summarise(scenario = "balanced", w_occ_lo = NA_real_, w_occ = NA_real_, w_occ_hi = NA_real_, w_ab_lo = NA_real_, w_ab = NA_real_, w_ab_hi = NA_real_, ratio_lo = q(ratio, .25), ratio = median(ratio), ratio_hi = NA_real_, n_replicates = n()))
hi <- bind_rows(dom |> group_by(scenario) |> summarise(w_occ_hi = q(w_occ, .75), w_ab_hi = q(w_abu, .75), ratio_hi = q(ratio, .75)), bal |> summarise(scenario = "balanced", ratio_hi = q(ratio, .75)))
S <- S |> select(-w_occ_hi, -w_ab_hi, -ratio_hi) |> left_join(hi, by = "scenario"); S <- S[match(forms, S$scenario), ]; S$scenario[S$scenario == "balanced"] <- "Balanced"
S <- S[, c("scenario", "w_occ", "w_ab", "ratio", "w_occ_lo", "w_occ_hi", "w_ab_lo", "w_ab_hi", "ratio_lo", "ratio_hi", "n_replicates")]; write.csv(S, "figure2_summary.csv", row.names = FALSE)
# quartiles are computed from the replicate values before the median is stored under the same name (summarise() evaluates arguments in order; the median must not overwrite w_occ before w_occ_hi)
B <- W |> filter(scenario == "balanced", form %in% c("Lr-Bn-La", "Lr-Bn-Sa")) |> group_by(form) |> summarise(w_occ_lo = q(w_occ, .25), w_occ_hi = q(w_occ, .75), w_occ_median = median(w_occ), w_ab_lo = q(w_abu, .25), w_ab_hi = q(w_abu, .75), w_ab_median = median(w_abu), n_replicates = n()) |> rename(w_occ = w_occ_median, w_ab = w_ab_median) |> select(form, w_occ_lo, w_occ, w_occ_hi, w_ab_lo, w_ab, w_ab_hi, n_replicates); write.csv(B, "figure2_balanced_forms.csv", row.names = FALSE)
options(width = 200); print(as.data.frame(S), digits = 3); print(as.data.frame(B), digits = 3)
