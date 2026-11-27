# Figure 6: context-dependent eutrophication effects with component-specific MS-GDUM

Package: zetaGDUM at commit f7e115f (installed from GitHub, not modified). Manuscript and package untouched. Everything is descriptive and predictive; nothing here is a causal effect.

## 1. Dataset inventory and selection
All 25 lake workbooks in `raw_data_new/By_groups` (Liu, Guohao - Nature_R01) were profiled by `00_inventory.R` (output `out_inventory_all_lake_files.csv`); selection used only data attributes (sites with TP, TN, Chl-a and coordinates; abundance data; TP range; richness), before any model was fitted.

| Group | File | Sites | Species | TP range (µg/L) | Why |
|---|---|---|---|---|---|
| Macroinvertebrates | LUSBM_2 (US NLA 2017) | 1081 | 558 | 1.4–11,239 | counts, largest n, widest TP gradient |
| Zooplankton | LUSZP_2 (same lakes, 1077 shared) | 1093 | 170 | 1.4–11,239 | abundance, different taxon, same environment |
| Phytoplankton | LPOPP (Poland) | 117 | 574 | 10–1,030 | abundance, only large phytoplankton set with TP+TN+Chl-a; independent region |

Not selected: other US NLA years (same lakes, repeats), Finnish sets (n=100, TP 22–116, narrow), LCAZP/LCAPP (n=69–137, TP < 45), fish/plants (presence only or no Chl-a). Chl-a, TN, TP log10-transformed; geographic distance (great-circle) as fourth predictor. Dropped for missing/zero values: 15 (BM), 16 (ZP), 0 (PP).

## 2. Informative orders (`zeta_order_information`, `tables/order_information.csv`)
All orders are flagged "candidate" by the package. I additionally required P(shared species) ≳ 0.1: macroinvertebrates k = 2–10, zooplankton k = 2–6, phytoplankton k = 2–5. Relative variance falls steeply (BM relative variance at k=3 is 0.23 of k=2, k=10 is 0.001); at k ≥ 12–16 shared composition vanishes.

## 3. Models and validation
Coherent architecture, three models (spread only, additive nonlinear position + spread, explicit interaction), all with geographic distance. Five spatial blocks (complete-linkage on great-circle distance), 8 repeats of 1500 combinations per fold (`03_validate.R`; 23,040 fits, 3 not converged). Importance by leave-one-block-out Shapley (`04_fit.R`).
Collinearity (log scale, VIF TP/TN/Chl-a): US sets 3.3/2.9/2.3 (r TP–TN 0.79, TP–Chl-a 0.73); Polish set ≈1.0 (r < 0.1).

## 4. Main findings
1. **Nonlinear position adds predictive value for total dissimilarity in the US lakes** (held-out ΔD² median +0.01 to +0.03 for macroinvertebrates, k = 2–6; zooplankton components +0.01–0.05 at k = 3–6), but not for phytoplankton (negative, n = 117, folds as small as 13 lakes).
2. **Explicit spread × position interaction adds essentially nothing.** At the default penalty (0.01) its gain is ~1e-5 D²; with penalty 1e-4 ≤ 0.001 (median); with 1e-6 it ranges −0.02 to +0.02 with inconsistent sign across folds (`tables/interaction_penalty_sensitivity.csv`). Interaction attributions are ~1e-5 (default) to ~1e-3 (1e-4). Context dependence seen in Panel C therefore comes from nonlinear I-splines of position, not from the explicit interaction. The hypothesis that interactions become more informative at multisite orders is **not supported** in these data, and not rejected either, because the interaction is strongly penalised/aliased with additive splines (as documented by Codex).
3. **Context dependence (Panel C, order 3, spread range supported at all positions).** Macroinvertebrates: per unit TP spread, dissimilarity rises about 1.7× faster among oligotrophic (TP≈10 µg/L) than eutrophic (≈110 µg/L) lakes for total dissimilarity (slope ratio eutrophic/oligotrophic 0.59), occupancy 0.77, heterogeneity 0.83; at k = 5 ratios fall to 0.27–0.35. Zooplankton: opposite direction at k = 3 (ratios 1.2–1.8) and strongly unstable at k = 5 (heterogeneity 7.7, total −5.8, total slope near zero). Phytoplankton: oligotrophic not supported (min TP 10); eutrophic slopes ≈ 0 against mesotrophic. So context dependence exists for US benthos, is inconsistent across taxa, and is small in absolute terms (total dissimilarity changes < 0.03 across the supported spread range at k=3).
4. **Component contrast.** With coherent closure, occupancy and heterogeneity predictions move in opposite directions with spread (occupancy up, heterogeneity down, e.g. benthos k=3: +0.099 and −0.071 per unit spread at low TP, total +0.028). This partly reflects that heterogeneity is defined only among shared species: wider environmental spread reduces sharing, shifting dissimilarity from H to O. Occupancy and heterogeneity effects are thus not independent evidence; total dissimilarity hides an almost-cancelling pair. This is the ecological information pairwise total dissimilarity cannot show.
5. **Importance (Panel B).** Geographic distance is the leading attribution for benthos at k = 2–3 and phytoplankton at all orders; TP leads (and persists to k=10) for benthos at k ≥ 4; Chl-a leads for zooplankton occupancy/heterogeneity at k ≥ 4 (up to 0.07 ΔD²), TN contributes little. Chl-a is partly a response to enrichment (r with TP 0.73), so this should not be read as an independent effect. Attributions are small (< 0.16 ΔD² in total) and fold ranges are wide (`tables/importance_predictor.csv`).

## 5. Negative findings and limits
- Phytoplankton: negative held-out D² for components at all orders, TP and TN attribution ≤ 0; five spatial blocks with 13–31 lakes. Not informative for the context question.
- Explicit interaction: no reliable gain (see above).
- Zooplankton surface at k = 5 is unstable (slope ratios of ±6–8 where slopes are ≈0).
- Failed coalitions (predictor spline/sample unsupported in a fold): benthos 21–30 of 90 per predictor in primary run, phytoplankton TP 12, zooplankton TP 6; averages use valid folds only and counts are in the tables.
- 3 of 23,040 validation fits did not converge.
- Supported fraction of evaluation sets falls from 0.92 (k=2) to 0.40 (k=10) for benthos (`validation_summary_coherent.csv`).
- No interval calibration exists for importance or interactions (Codex doc); fold ranges are descriptive, only five folds.
- Prediction surfaces are conditional averages over actual site configurations, not rectangular grids; oligotrophic position at k ≥ 8 is only 1 supported spread bin (common range <3 bins).

## 6. Issues for Codex
1. Phytoplankton Polish set: TP, TN, Chl-a are essentially uncorrelated (|r| < 0.1), unusual for lake data; check source units/site matching (Guohao).
2. `ugdm_component_importance` returns NA for a predictor in some folds without a stated reason; a diagnostic field naming the failed coalition (and the 4-distinct-values rule) would help.
3. Default `interaction_shrink = 0.01` makes the explicit interaction numerically equal to additive; recommend documenting which penalty to use for empirical hypotheses, selected on training folds.
4. Surface supports: `ugdm_component_surface` position bins are fixed grid centres, not user-selected TP values; I picked the nearest supported bin within 0.15 log10.

## 7. Files
Figure candidates: `figures/Figure6_main_order3.{pdf,svg,png}` (Panel C at order 3), `figures/Figure6_alt_order5.*` (Panel C at order 5). Scripts: `00_inventory.R, 01_prep.R, 02_orders.R, 03_validate.R, 03b_penalty.R, 04_fit.R, 04b_fit_lowpenalty.R, 05_tables.R, 06_surface_tables.R, 07_figure6.R` (run with env `ZG_LIB6` pointing to the library holding the f7e115f build). Tables in `tables/`; replicate-level results in `out/*.rds`.
