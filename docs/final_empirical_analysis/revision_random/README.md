# revision_random: random site-level validation, model equivalence, Figure 6 versions, order-selection diagnostics

See `REPORT.md` for results, discrepancies and limitations. Scripts (run from `final/`, order in `run_all.sh`):

| script | purpose | main outputs |
|---|---|---|
| `70_splits.R` | 30 random 70/30 site splits per community pool | `tables/random_splits.csv` |
| `71_fit_random.R`, `72_tables_random.R`, `73_fig5_random_figure.R` | Figure 5 under random validation | `tables/rv5_*`, `figures/Figure5_random_validation*`, `FigureS5_random_*` |
| `81_cv_interactions_random.R`, `82_tables_fig6_random.R` | interaction CV under random validation | `tables/rv6_*` (coefficients as `.csv.gz` in the repository) |
| `83_figure6_versions.R` | two Figure 6 versions (`FIG6_VERSION=additive` or `k2_additive_higher_TP_interaction`) | `figures/figure6_*` |
| `90_equivalence.R`, `91_equivalence_figure.R` | gdmmTMB / gdm / zetadiv equivalence | `tables/eq_*`, `figures/FigureS_model_equivalence*` |
| `95_order_diagnostics.R` (+ `95a_support_machinery.R`), `96_order_figure.R` | three axes of order selection | `tables/order_diagnostics_*`, `figures/FigureS_order_selection_diagnostics*` |

`figures/FigureS5_spatial_block_validation.*` is the previous (spatial-block) Figure 5, unchanged.
