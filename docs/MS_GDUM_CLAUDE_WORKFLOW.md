# Component MS-GDUM workflow for Claude

Run the example from the repository root after installing the development package in an R library.

```sh
R CMD INSTALL .
Rscript docs/final_gdum/workflow.R
```

The example simulates 168 sites on a transect. Phosphorus is in simulated concentration units and is spatially autocorrelated. Species occurrence and conditional abundance respond to it, with opposing position-dependent niche widths. The simulation is for demonstrating the workflow. Its fitted curves provide no evidence about a real freshwater community.

## Prepare and screen

Map the rows of the community, environmental predictors, coordinates and fold assignments to the same site IDs before fitting. Supply finite observations and positive community totals. Choose the sampling design and spatial blocks before comparing models. A cluster or block needs at least as many validation sites as the order being evaluated. Continuous spline predictors need four distinct values in each training sample.

`zeta_order_information(community, environment, response = "all")` profiles incidence, total relative-abundance, occupancy and heterogeneity variation, as well as design support and integration error. Inspect both absolute and relative variance. Record exact degeneracy, missing pairwise denominators and gaps in candidate orders. The profile supplies context for choosing which models to assess. Its candidates still need validation.

## Compare independent spatial performance

Use `ugdm_validate_components()` with predetermined `fold` labels. Alternatively supply `coordinates` for complete-linkage spatial blocks and set `lonlat = TRUE` for geographic coordinates. Coordinates only enter the model when `geographic = TRUE`. The function draws training and evaluation sets solely from their respective site partitions. It re-estimates the basis in each training fold and evaluates competing models on matched sets. Examine convergence, calibration bias, RMSE, explained deviance relative to the training null and `supported_fraction`. A fold with little overlap should be visible even if its prediction score is positive.

Use the additive coherent model when joint compositional predictions are required. Fit the independent architecture when comparing marginal predictions or auditing whether the coherent link loses accuracy. The two fit the same total model under matched settings. Consider an interaction only when an explicit spread by position hypothesis has been defined and independent, supported predictive performance justifies the extra coefficients. In an analysis with penalty tuning, do the tuning on training sites and verify the chosen penalty on separate validation sites.

## Fit, attribute and visualise

1. Draw site sets with `sample_sets()`. Calculate `ugdm_component_response()` from those same sets.
2. Build `ugdm_component_design()` with `model = "additive"`, `"spread_only"` or `"interaction"`. Call `ugdm_fit_components()` with the desired architecture and declared penalties.
3. Build a new design with `reference = fit` to reuse the same environmental and geographic basis. A prediction design must use the fitted order.
4. Use `ugdm_component_diagnostics()` for unpenalised information rank, condition, active bounds and the residual interaction geometry.
5. Use `ugdm_component_importance()` on independent assessment sets where possible. Attributions are conditional on the model, penalty, null, score and assessment distribution. For many predictors use `resolution = "predictor"` to reduce coalition cost. Record failed coalitions and the returned numerical error for sampled permutations.
6. Use `ugdm_component_surface()` to display predicted dissimilarities on actual configurations. The returned grid masks bins with inadequate distinct training sets or training sites and bins without supported evaluation configurations. The default support counts are display settings. Inspect the effect of binning and vary them when warranted. The function describes conditional predictive averages over the supplied sets. Different points can represent different underlying configurations or other predictor values.
7. For a spread comparison at low and high phosphorus, choose represented position bins and only compare spread bins supported at both positions. Compare the corresponding additive surface before attributing the contrast to an explicit interaction. Nonlinear I-spline distances can generate context-specific effects in the additive model.

The example writes a CSV of the masked surface and a plot of cross-sections with common spread support to `ms_gdum_workflow_output/` by default. Set `MS_GDUM_WORKFLOW_OUT` to choose another directory. The example fits an interaction to illustrate extraction, with the choice left separate from the empirical comparison.

## Confirm numerical precision and interpretation

Repeat independent set integration with the same sites, model specification, basis and evaluation sets. Pass a matrix with one row per independent integration and one column per predictor to `ugdm_within_order_diagnostic()`. Use the component name and the explained deviance from independent site validation. This estimates numerical sensitivity on fixed sites. It does not represent ecological uncertainty. Keep `calibrated = FALSE` for the workflow. The completed study has not established interval calibration for empirical predictor attribution or explicit interactions.

Report descriptive importance and surfaces only when relevant environmental support and prediction are defensible. Keep unsupported predictions and failed coalition results visible. Explain whether the intended use is modelling total, occupancy-associated or heterogeneity dissimilarity.

[Executable workflow](final_gdum/workflow.R) · [Example surface](final_gdum/results/phosphorus_surface.csv) · [Example cross-sections](final_gdum/results/phosphorus_spread_sections.pdf) · [Statistical validation](MS_GDUM_FINAL_VALIDATION.md)
