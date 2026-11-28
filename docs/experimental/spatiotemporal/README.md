# Experimental spatiotemporal MS-GDUM

This prototype operates on repeatedly sampled sites and separates spatial occupancy, persistence across surveys, and positive abundance. It is sourced from this directory and uses the existing `zetaGDUM` package for the spatial response and environmental basis. The [technical report](REPORT.md) gives the definitions, likelihood, simulation results and limits of inference.

## Reproduce

Run from the repository root with R ≥ 4.0. Install the repository package with `R CMD INSTALL .` or use a compatible installed build. The results in this directory used the committed development baseline `a49359d`. Tests require `testthat` and compare the multisite beta implementation to `betapart` when available. SVG generation requires an R build with Cairo support.

```sh
Rscript docs/experimental/spatiotemporal/test.R
Rscript docs/experimental/spatiotemporal/run_simulations.R
Rscript docs/experimental/spatiotemporal/make_figures.R
```

Set `ST_REPLICATES` for a smaller pilot and `ST_OUTPUT` for an alternative output directory. Use the same output directory for the simulations and figures. The default is 32 replicates per scenario and writes the tables and figures in this directory. Timing depends on the platform.

## Example

```r
source("docs/experimental/spatiotemporal/load.R")
st_load()

sim <- simulate_st_gdum(c("interaction", "losses_only"), seed = 8)
panel <- sim$latent
surface <- st_zeta_retention(st_zeta_decline(panel, spatial_max = 5, temporal_max = 5))
plot_st_surface(surface, "zeta", "Mean shared richness")
events <- st_decompose(panel, spatial_order = 2)
species <- st_occupancy(panel, before = 1, after = 5)

design <- st_design(panel, sim$env, spatial_order = 1)
terms <- c("position_initial.gradient", "position_change.gradient",
           "position_by_change.gradient")
fit <- st_gdum(panel, design, gain_terms = terms, loss_terms = terms[1:2],
               blocks = seq_len(dim(panel$x)[1]))
prediction <- st_predict(fit)
```

Supply a numeric sites × times × species array, with values of zero for recorded absences and `NA` for missing observations. Use `st_panel(x, times = survey_times)` to set the ordered numeric survey times and units used for lag. Site and species names are retained. Environmental data use a sites × times × predictors array. Named environmental axes are aligned to the panel. An explicit `missing = "complete_sites"` selects sites with every survey and species observation available across the full declared panel and reports the dropped sites.

- `st_zeta()` returns site and time subsets, shared-species counts and species labels per rectangle. `st_zeta_decline()` computes exhaustive spatial means without enumerating site sets. `st_zeta_retention()` adds ratios to the preceding order along each axis.
- `st_decompose()` returns endpoint species events and full sharing histories for selected survey times. `st_occupancy()` returns species-specific spatial and numerical changes plus exact mean zeta contributions across orders.
- `st_beta()` returns both native pairwise partitions. `st_beta_multisite()` returns Baselga’s multisite Sørensen partition for a specified set of sites. `st_beta_change()` accounts for pairwise beta change using observed gain and loss masks.
- `st_design()` constructs environmental predictors for specified transitions and spatial combinations. `st_gdum()` estimates conditional gain and loss probabilities. `st_spatial_gdum()` fits native MS-GDUM snapshots on a common environmental basis. `st_abundance_gdum()` fits positive abundance or paired log-abundance changes at spatial order one.
- `st_order_diagnostic()` shows support, zero fractions and the information available for a common position contrast across orders. It reports Monte Carlo error for sampled zero fractions. `plot_st_gdum()` displays predictions on user-supplied design rows, which should follow a specified environmental contrast.

## Outputs

| Figure | Diagnostic |
| --- | --- |
| [01](figures/01-decline-retention.svg) | Joint decline and retention along each axis |
| [02](figures/02-order-information.svg) | Zero sharing, effective species and relative position information |
| [03](figures/03-directional-sharing.svg) | Gains, losses and net shared richness across spatial orders |
| [04](figures/04-occupancy-abundance.svg) | Occupancy and conditional-abundance changes for individual species |
| [05](figures/05-effect-recovery.svg) | Recovery of known environmental coefficients |
| [06](figures/06-environmental-functions.svg) | Fitted spatial spread, gains, losses and positive abundance |
| [07](figures/07-detection-and-prediction.svg) | Spurious events under imperfect detection and held-site prediction |
| [08](figures/08-overlap.svg) | Working and site-jackknife uncertainty under overlapping sets |

The [results directory](results/) includes raw replicate-level tables, order diagnostics, scenario-level species examples and run metadata.
