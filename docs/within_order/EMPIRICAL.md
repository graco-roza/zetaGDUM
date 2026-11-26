# Within-order profiles for the five manuscript communities

The validation uses the existing DD2 data and documented preprocessing for Finnish birds, Forest Hemiptera, Pond diatoms, Arctic plants and Atlantic marine fish. It fits the nonlinear additive MS-GDUM with a fixed slope ridge of .001. Each community is evaluated at orders 2, 4, 8 and 12 and at its previously analysed middle and high orders, with duplicates removed. Each order must fit inside the smallest spatial validation block. The responses are incidence dissimilarity normalised by average set richness and relative-abundance dissimilarity, both with constant working weights.

Two independent uniform combination pools estimate importance on the same observed sites. Their discrepancy on a single, independent, common assessment pool measures sensitivity to sampled training sets. For predictor collections above eight, 48 common permutations are used for both fits, with the corresponding conditional permutation Monte Carlo standard errors. The paired two-pool gap describes numerical variability and does not estimate ecological uncertainty. Both Shapley predictor contributions and Owen allocations of spread, position and distance are reported. Relative contributions are defined when total in-sample explained deviance exceeds .0001.

Four spatially contiguous held-out blocks assess site transfer. The established DD2 predictor preprocessing is fixed before cross-validation. Outside the Finnish bird dataset this includes a skewness-based transform chosen from the full covariate table. The transfer diagnostic is conditional on that choice. Each fold then estimates the standardisation, spline basis and distance range on the other three blocks. The pooled held-out explained deviance is computed from the ratio of the summed fold mean deviances to the summed fold mean null deviances. Each fold null predicts the mean training dissimilarity on the held-out sites. The figure ticks identify orders with positive pooled held-out explained deviance. This is descriptive evidence for predictive transfer. Replicated score integration and ecological interval calibration remain necessary to declare predictor separation.

The community profiles can be regenerated from the repository root with the installed development package. Supply the original DD2 data directory as `DD2`.

```sh
for i in 1 2 3 4 5; do
  Rscript docs/within_order/empirical_profile.R "$DD2" "$i" 1200
done
Rscript docs/within_order/summarise_empirical.R
Rscript docs/within_order/plot_empirical.R
```

The output includes response variation, nonzero sharing, design rank and conditioning, convergence, importance estimates, independent-pool differences, and spatial validation results. The plots show descriptive ranking trajectories. [Absolute importance](results/empirical_absolute.pdf) · [Relative importance](results/empirical_relative.pdf) · [Spread, position and distance](results/empirical_components.pdf) · [Compositional variation](results/empirical_composition.pdf) · [Combined diagnostics](results/empirical_overview.csv)
