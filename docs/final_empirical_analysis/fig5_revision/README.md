# Figure 5 revision: transformation sensitivity, paired position effects, conditional responses

Separate from the original Figure 5 (`../figures/Figure5_*`, `../figure5_environmental_redesign.R`), which is untouched. Candidate figures for comparison before a final selection. zetaGDUM f7e115f, unmodified. No manuscript text.

## Scripts (run from `final/`, `ZG_LIB6` set)
| script | purpose |
|---|---|
| `50_fit_configs.R` (`run_fits.sh`) | 4 configurations {log10, raw} x {linear, nonlinear position} + spread-only baseline per transformation; held-out/training D2, convergence, warnings, Shapley/Owen absolute contributions. `TRIM=1`: extreme-lake sensitivity |
| `51_conditional_curves.R` | conditional model predictions (link and response scale), all informative orders, TP/TN/Chl-a x spread/position |
| `52_tables.R` | paired within-fold contrasts, performance, convergence, importance, extreme-value diagnostics |
| `53_figures.R` | candidate main figure and supplements (PDF/SVG/600 dpi PNG, 183 mm) + plotted-value CSVs |

## Design
* Identical spatial folds (5 complete-linkage blocks), sampled training/evaluation sets (seeds depend on seed, fold and order only), seeds 1:3, informative orders, penalties and representations in every configuration. Geographic distance is unchanged. Spread-only baselines use the same predictors and transformation.
* Paired contrasts: seed means within fold first; then fold values, median, min, max (`sens_paired_*`). Contrasts: linear vs spread-only, nonlinear vs spread-only, nonlinear vs linear.
* Conditional curves: other design columns (spread, position and distance terms of all predictors) held at their medians over evaluation sets; the focal columns come from the 200 sets nearest to the grid value (keeps the I-spline/polynomial blocks consistent). Feasibility flag: the focal statistic plus the other set statistics at their medians must be within the 95th percentile of leave-one-out nearest-neighbour distance among real sets. `range_type = common` uses the intersection of the 5-95% ranges of the set statistic across all informative orders of a taxon, `own` the order-specific range.
* Extreme nutrients: lakes above the 99th percentile of TP, TN or Chl-a (24 lakes per pool) removed from training and evaluation (`*_trim99`, same fold labels; test sets therefore differ, so the comparison with the full data is not strictly paired).
* The main-figure configuration (log10 predictors, nonlinear position for B and C) was fixed a priori as the most general model, not chosen from results. Panels A and S1 show both position models.

## What the results show (shared pool; see tables)
* Convergence: 100% of fits converged in all four configurations and both pools (1,440/1,440 in the full pool). Raw predictors trigger the package's skewness warning in every fit (|skewness| TP 20.4, TN 9.7, Chl-a 22.9 vs 0.6, 0.03, 0.17 after log10).
* Position adds held-out skill for macroinvertebrates (median paired delta D2 about +0.02 to +0.03, positive in 61-63 of 80 fold-order combinations on the log scale; spread + nonlinear position slightly ahead of linear), but the fold spread is large, with negative values at orders 10-12 and for relative abundance. For zooplankton the gains are small and inconsistent (positive in 24-26 of 40 fold-order combinations; at k = 2 median about zero).
* Raw scale: lower training D2 everywhere; held-out D2 equal or lower than log10, and for zooplankton the median paired gain of position is slightly negative (-0.003). Removing the 24 extreme lakes changes raw-scale results in both directions: held-out D2 of the raw macroinvertebrate incidence model rises at orders 2-6 (about +0.03 to +0.05), whereas zooplankton incidence raw-scale D2 falls from about 0 to about -0.2 at k = 5; log10 models change little. Raw-scale contributions of TP position, TN position and geographic distance shift by up to about 0.05 (S7B). The raw scale is therefore sensitive to the extreme lakes, but not in a single direction. The extremes also enter many sampled sets (about 4% of sets at k = 2, 11% at k = 5, 24% at k = 12 contain an extreme lake).
* Conditional TP curves: on the link scale the TP position and spread effects do not weaken with order: the link-scale change across the common range generally increases (position, macroinvertebrate incidence 0.22 at k = 2 to 0.88 at k = 12), while the change on the response scale falls towards zero because the mean fitted dissimilarity approaches 1 (0.63 at k = 2, 0.995 at k = 12). The flattening of the original panel C is therefore mostly ceiling compression of the response scale, not a weaker environmental response (S3, S3b, S6). The link-scale slopes depend on the fitted model and are not independent evidence of the mechanism; fold-wise uncertainty of the curves was not computed.
* Negative contributions are retained: nonlinear position of TN and geographic distance are negative at higher orders.

## Outputs
Tables (`tables/`): `sens_validation_folds_*` (fold/seed fits), `sens_paired_folds_*`, `sens_paired_summary_*`, `sens_performance_*`, `sens_convergence_*`, `sens_importance_*`, `sens_importance_folds_*`, `sens_importance_fold_values_*`, `sens_curves_grid_*` (prediction grids, both scales, all predictors and configurations), `sens_curves_ranges_*`, `sens_curves_summary_*`, `sens_extremes_*`. Suffixes: `shared`, `full`, `shared_trim99`, `full_trim99`.

Figures (`figures/`, + `*_<panel>.csv` plotted values): `Figure5_revised_main` (shared pool), `FigureS5r_1_performance_four_configurations`, `_2_importance_raw_vs_log`, `_3_TP_conditional_response_scale`, `_3b_TP_conditional_centred_link`, `_4_TN_Chla_conditional`, `_5_relative_importance_unclipped`, `_6_saturation_link_vs_response`, `_7_extreme_nutrient_sensitivity`, and the same set with `full_pool` in the name for the full pools. Relative share = absolute contribution / mean held-out D2, undefined where held-out D2 <= 0 and drawn hollow where < 0.02.
