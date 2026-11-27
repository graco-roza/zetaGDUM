# Final empirical analysis batch (Figures 3-6)

Communities: **BCI trees, 1985 census** (tree records, status alive, dbh >= 10 mm; 1250 20 m x 20 m quadrats; documented in `../bci/README.md`), **EPA National Lakes Assessment 2017 macroinvertebrates** and **zooplankton** (workbooks LUSBM_2 and LUSZP_2 of the shared Liu Guohao eutrophication synthesis folder; `../fig6/01_prep.R`). Distance Decay 2.0 data and Polish phytoplankton are not used anywhere. Package: zetaGDUM at commit f7e115f, unmodified; no new statistical procedures. No manuscript text is included.

**EPA pools.** The *shared pool* (principal) is the 1,077 lakes with valid TP, TN, Chl-a, coordinates and non-empty samples in both taxon tables; environmental values are identical in the two tables (checked). The *full pools* (1,081 macroinvertebrate and 1,093 zooplankton lakes) are sensitivity analyses. Taxa (not species) are the EPA units; macroinvertebrates are mostly genus-level.

## Layout and run order (set `ZG_LIB6` to the library containing the f7e115f build)
| Step | Script | Output |
|---|---|---|
| Data pools | `00_data.R` | `out/pools.rds` (not committed), `tables/dataset_characteristics.csv` |
| Evaluated orders | `10_orders.R` | `tables/evaluated_orders.csv`, `evaluated_orders_summary.csv` |
| Fig. 3 exact calculations | `11_fig3_calc.R` (POOL=shared and full) | `tables/fig3_*` (observed, taxon contributions, bootstrap and permutation draws) |
| Fig. 3 figures | `12_fig3_figure.R` | `figures/Figure3_main`, `FigureS3_*`, `tables/fig3_plotted_values_*.csv` |
| Fig. 4 | `40_fig4_copy.sh` (analysis in `../bci`) | `figures/Figure4_main`, `FigureS4_*`, `tables/fig4/` |
| Fig. 5 fits | `20_fig5_fit.R` (POOL=shared and full) | `tables/fig5_validation_folds_*`, `fig5_importance_folds_*`, `fig5_response_curves_*` |
| Fig. 5 figures | `21_fig5_figure.R` | `figures/Figure5_main`, `FigureS5_*`, `tables/fig5_plotted_*` |
| Fig. 6 and supplement | `30_setup_fig6.R`, then `fig6_shared/run_fig6.sh` and `fig6_full/run_fig6.sh` | `fig6_<pool>/figures`, `figures_repl`, `tables`, `tables_repl` |
| Documentation | `41_documentation_tables.R` | `tables/model_settings.csv` |

All stochastic steps have fixed seeds (listed in `tables/model_settings.csv` and the script headers). Large intermediate `.rds` files are not committed.

## Evaluated orders
Orders at which a random set of k sites shares at least one taxon with probability >= 0.10 (`zeta_order_information`, grid 2-640): macroinvertebrates 2-12 (2,3,4,5,6,8,10,12), zooplankton 2-5, BCI 2-480 (19 orders). The same criterion is used in Figures 3, 5 and 6. The earlier eutrophication analysis used macroinvertebrates 2-10 and zooplankton 2-6 (a looser reading of the same criterion); orders are harmonised here.

## Deviations from the earlier analyses (stated for transparency)
* Figure 3 (formerly five-community letter analysis): the exact analytical contributions of `expected_species_contributions()` replace the sampled-set estimates (4,000 sets) and the paired Bayesian site-weight draws; uncertainty is shown with 100 site-resampling draws of the exact quantities. Within-site abundance permutations (499) are retained. Hold-out and position-share analyses of the earlier letter now appear in Figure 5.
* Figure 5 uses the original MS-GDUM structure (spread terms, linear position, geographic distance) for incidence ("total") and relative-abundance dissimilarity and is distinct from the coherent component analysis of Figure 6 (total, occupancy-associated and heterogeneity components, nonlinear position, optional interaction).
* The held-out D2 is computed against the training-set mean on sets drawn only from held-out blocks; it can exceed the training D2 when block compositions differ.
* Importance allocations can be negative (held-out gains); relative importance is shown only where the held-out D2 >= 0.02.

## Figure 6 (closing figure) — added
`figure6_closing.R` (A: centred TP-spread response curves at low/intermediate/high TP by component; B: diverging occupancy-associated vs abundance-heterogeneity endpoint bars with their sum) and `figure6_conditional.R` (matched synthetic site-set conditional predictions from the coherent additive model; zetaGDUM f7e115f built from git into a separate library; support screen: package component support + joint nearest-neighbour support + TP inside training range; common profiles/spread cells across all positions). Outputs: `figures/figure6_closing/` (`..._binned_*` = supported cell means, NOT isolated TP effects; `..._conditional_*` = matched-profile predictions) and `figures/figure6_conditional/` (predictions, per-profile diagnostics, `conditional_support_summary.csv`). Conditional figures exist only where common support was established (shared k=3, full k=3); k=5 conditional support is insufficient for the shared pool and for full-pool macroinvertebrates (see summary), so the labelled binned figures are the k=5 record. Deletion-fit (leave-one-spatial-block-out) conditional sensitivity: no deletion fit retained full common support, so no deletion ranges are available.

## Figure 6 TP spline candidates
`figure6_eutrophication_splines.R` (run with `FIG6_RAW_SPREAD=1 ZG_LIB6=<f7e115f lib> Rscript figure6_eutrophication_splines.R fig6_interaction/tables figures/figure6_splines`): fitted TP-specific interaction spline shapes (penalty 1e-6 main; 1e-2, 1e-4 sensitivities; shared and full pools; total and heterogeneity-share stages) for predeclared low/high orders (BM k = 2, 10; ZP k = 2, 5) at TP positions 10/35/100 ug/L, plus raw-spread profiles evaluated only within TP spread-position support (`Figure6_raw_spread_support.csv`; unsupported panels are labelled, never extrapolated). `*_values.csv` and `*_axis_checks.csv` accompany every figure. Only change to the supplied script: the micro sign in two axis labels is a plotmath expression (PNG device rendered it as dots).

## Figure 6 synthesis
`figure6_synthesis.R` (`ZG_LIB6=<f7e115f lib> Rscript figure6_synthesis.R fig6_interaction/tables figures/figure6_synthesis`): matched synthetic TP profiles (observed non-TP configurations held fixed; TP position and spread set by construction; no refit; TP-only interaction penalty 1e-6 vs additive). Low order k = 2; high order = highest k > 2 with P(shared) >= 0.50 and a common supported TP-spread interval (selected on the shared pool, from support and sharing only): macroinvertebrates k = 5 (12 common profiles, spread 0.45-0.70), zooplankton k = 3 (18 profiles, spread 0.05-1.10). Full-pool sensitivity: no common supported contrast for the higher order (panels labelled "Insufficient common environmental support"). Support diagnostics: `support_*.csv`, `Figure6_panel_support.csv`, `Figure6_candidate_order_support.csv`, `Figure6_axis_checks_*.csv`.

## Figure 6 GDM-style synthesis
`figure6_gdm_synthesis.R` (`ZG_LIB6=<f7e115f lib> Rscript figure6_gdm_synthesis.R fig6_interaction/tables figures/figure6_gdm_synthesis`): predeclared orders (macroinvertebrates k = 2, 5; zooplankton k = 2, 4). A: coherent component predictions (nine bars per supported panel); B: independently fitted I-spline transformations f(TP | position) for total, occupancy-associated and heterogeneity stages (nested facets; linear TP axis in the main figure, `*_logaxis` versions for the supplement; `Figure6_GDM_splines_*` = standalone twelve-panel figure). Checks: independent total fit reproduces the coherent total (max prediction/coefficient error ~1e-15), all nine bars rendered in every supported barplot (`Figure6_A_bar_checks.csv`), zooplankton k = 4 support (`Figure6_panel_support.csv`, `support_*`). Cosmetic edits to the supplied script: row-strip text size 5.2 (labels were clipped) and an ASCII hyphen in the spread-interval note (en dash rendered as dots in PNG).
Update: every panel of Panel B (and the standalone twelve-panel figure) now has its own free x and y scale (`facet_nested(..., scales = 'free', independent = 'all')`); the shared-limit anchors were removed.
