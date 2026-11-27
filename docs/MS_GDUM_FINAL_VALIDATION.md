# Final statistical validation of component MS-GDUM

The package now fits total relative-abundance dissimilarity, occupancy-associated dissimilarity, and abundance heterogeneity among shared species. The coherent architecture guarantees joint compositional closure. The additive and spread-only models provide interpretable predictive comparators. The optional interaction is usable for targeted, supported predictive comparisons. The validation does not establish general ecological uncertainty calibration for predictor importance.

[Figure 3 component validation PDF](final_gdum/results/figure3_component_validation.pdf) · [Figure 3 SVG](final_gdum/results/figure3_component_validation.svg) · [R workflow for Claude](MS_GDUM_CLAUDE_WORKFLOW.md)

## Models and estimands

Let $T$ be total relative-abundance dissimilarity, $O$ the occupancy-associated component and $H$ abundance heterogeneity among taxa shared by all sites in the set. The existing response partition gives $T=O+H$. The independent architecture estimates each marginal by a fractional-logit model. It returns three bounded predictions, with independently estimated $O+H$ generally different from predicted $T$. The signed error in that equality is returned for diagnosis.

The coherent architecture models $μ=E[T\mid X]$. The second fractional-logit criterion uses response $H/T$ at positive total and working weight $T$. At a fixed $X$, its expected negative log score is

$$
-E[H\mid X]\log\pi-\{E[T\mid X]-E[H\mid X]\}\log(1-\pi).
$$

The unrestricted optimum is the ratio $E[H\mid X]/E[T\mid X]$ when the denominator is positive. The fitted fractional regression approximates that ratio. Zero-total rows contribute to the total model alone. If the training total is identically zero, predicted shared composition is one and the share stage is flagged unidentified. Constant boundary stages return exact zero or one. The reported composition is $(1-μ,\,μ(1-\pi),\,μ\pi)$.

The three designs use directional I-spline spread alone, additive nonlinear position, or two directional endpoint curves interpolated by the clipped mean set position. Only the total model constrains directional coefficients by default. Occupancy, heterogeneity and the coherent share allow unrestricted signs. A ridge acts on the endpoint average and a separate penalty acts on differences. Tied endpoint coefficients reproduce the additive design and its ridge penalty. The model can therefore distinguish explicit context dependence from the nonlinear position dependence already induced by a directional I-spline. At order two, raw spread and midpoint determine the pair endpoints. A raw spread contrast between positions can therefore change the I-spline difference even with identical endpoint coefficients.

The package computes unpenalised working-information rank and condition for each stage. A penalised converged solution by itself supplies no rank evidence. It also projects the interaction columns off the additive design to measure residual environmental support. For a display, it counts distinct training combinations and sites in local spread–position bins, retaining only represented evaluation configurations. The binning checks each predictor's bivariate geometry. Simultaneous bivariate support is weaker than full multivariate overlap. Positive cross-validation scores can accompany sparse support. The surface therefore averages predictions on represented configurations rather than evaluating a rectangular grid of arbitrary site sets.

## Predictive attribution

Refits for each coalition use the same training response, design coordinates, sampling weights and penalty restrictions. The null uses training component means. Between predictors the allocation is Shapley. Within a predictor it uses Owen spread and position allocations, with explicit interaction entering after both main effects. Each interaction allocation is the gain from relaxing that predictor's equal-endpoint restriction in the corresponding coalition. These gains telescope to the predictor allocation and the sum of predictor allocations is total explained deviance for that component. The coherent architecture refits both stages and scores the final $O$ and $H$ predictions. Coefficients or allocations from the share stage alone do not measure importance of either component. Held-out gains may be negative. Relative shares are reported only when the total gain is meaningfully positive.

## Simulation and reference design

Seven scenarios alter the prevalence distribution, niche strength, position-dependent niche width, predictor correlation and spatial dependence. Occurrence-only, abundance-only, opposing context, and null scenarios each used **200 independent 110-site communities** at orders 2, 4 and 8. The same community and two independent 900-combination pools are used across models for each order. All seven scenarios had separate 12-community pilots with 600 combinations. Additional pilots examine weak environmental effects, strong predictor correlation and spatial dependence. Each trial uses unit working weights for total and independent marginals. The coherent share uses dissimilarity working weights.

Two independent 2,100-site communities and 16,000 sets per reference fit give population-model approximations. A separate set of 480 evaluation sites and 900 fixed configurations provides both an independent realised community and response means estimated over 140 simulated assemblages at that geometry. Predictive transfer uses the realised holdout. Recovery of the conditional surfaces and associated Shapley game uses the estimated conditional mean responses. These two normalisations of explained deviance differ. The conditional-mean null can be especially small in the null scenario. The observed-response null audit therefore uses realised held-out responses and the theoretical zero predictor contributions of the independence mechanism. Reference-range, error and availability are retained. No causal effect is inferred from either predictive allocation.

The site-assembly simulation deliberately allows an occurrence driver to alter abundance heterogeneity and a conditional-abundance driver to alter occupancy dissimilarity. Relative standardisation changes means and minima even when only one underlying process responds to environment. The null scenario assesses chance predictive gains. A distinct noise-free positive control on the fitted basis assesses implementation-level recovery of an explicit interaction.

[Simulation mechanisms](final_gdum/design.R) · [Model results](final_gdum/results/model_summary.csv) · [Reference contributions](final_gdum/results/reference_importance.csv) · [Reference response surfaces](final_gdum/results/reference_surfaces.csv.gz)

## Architecture and interaction results

Coherent predictions satisfied closure to floating-point precision. Total predictions from the two architectures were identical under matched settings. Across the 12 principal additive scenario–order cells, the median difference in oracle RMSE (coherent minus independent) was +.00046 for occupancy and −.00017 for heterogeneity. The ranges were −.00211 to +.00553 and −.00989 to +.00042, respectively. The coherent architecture's main advantage is valid joint predictions. These differences do not justify a general accuracy advantage. Observed mean calibration bias and slope are reported for each component, and mean null or closure error should be checked directly in new datasets.

In the contrasting-context order-2 scenario, the coherent spread-only, additive and interaction designs gave heterogeneity oracle RMSE of .07399, .05314 and .05274. Their median deviance explained on the realised evaluation response was .271, .414 and .416. Occupancy oracle RMSE changed from .10818 to .05112 and .05116 across those designs. The large gain came from adding nonlinear position. The default endpoint penalty of .01 left predictions close to the additive model. A matched 32-community sensitivity study lowered the interaction penalty to .0001 and improved median oracle RMSE in the opposing-context scenario by .0096 at order 2 and .0067 at order 4, averaging comparisons over the three responses by their median. The corresponding change at order 8 was below .00001, and the occurrence-only scenario showed little benefit. The penalty was compared against the same evaluation sets. Those results are sensitivity evidence, not an independently selected tuning rule. The noise-free identification control at order 4 had total RMSE .0120 under the additive design and .0000114 under the explicit interaction.

At order 2 in the interior of the training domain, an I-spline span with linear and quadratic polynomials contains both the endpoint difference $b-a$ and $b^2-a^2$. Multiplying the linear difference by the raw midpoint gives $(b^2-a^2)/2$. A linear combination of the explicit varying-coefficient columns can therefore duplicate additive spread. Penalisation can recover a useful prediction in the resulting span without identifying all of its coefficients. The median interaction condition in the contrasting order-2 cell was infinite, compared with approximately 963 for the additive design. At higher orders the same equality need not hold because the set mean and an individual pair midpoint differ. This is a reason to compare the nested designs and focus on supported predictions rather than individual endpoint coefficients.

The principal cell convergence rate was .990–1.000, with valid importance estimates in .945–1.000 of replicates per cell. In the contrasting scenario the coherent additive model recovered the leading predictor in .995–1.000 of replicates when the reference difference exceeded .02. Its contrast RMSE nevertheless ranged from .012 at order 2 for total to .115 at order 8 for occupancy. The incidence of correct ordering alone does not establish precision. In the null scenario the coherent additive median held-out explained deviance was negative across all evaluated orders and components, ranging from −.023 to −.008. Occurrence-only heterogeneity still had median deviance explained .111 at order 2, and abundance-only occupancy .415. These results demonstrate cross-component signal under single-mechanism simulations.

[Penalty sensitivity](final_gdum/results/interaction_penalties.csv) · [Interaction positive control](final_gdum/results/interaction_positive_control.csv) · [Matched contrast recovery](final_gdum/results/contrast_summary.csv) · [Pre-fit simulation profiles](final_gdum/results/simulation_profiles.csv.gz)

## Low-order uncertainty

The focused study used the coherent additive fit and exact predictor coalitions. At orders 2 and 4 it drew 59 conditional-Poisson Bayesian bootstrap samples in each of 200 occurrence-only communities with 70 sites. Uncertainty was calculated on the same held-out configurations and the same predictor/response for every draw. The 850-set fit was repeated with an independent pool. For the finite-community target the pairwise calculation enumerated all 2,415 pairs. The order-4 counterpart used 16,000 combinations and retains integration uncertainty. The population-model reference used two independent large site samples. The study compared percentile intervals with their basic reflections centred at the original estimate. It also evaluated the matched bootstrap mean bias correction to the point estimator. The assessment covers 80% and 95% coverage, bias, RMSE, width and practical false declarations.

The 64-community observed-response sensitivity makes the attribution scale comparable to held-out empirical deviance. For order-2 heterogeneity secondary importance, point bias was −.0237 and the exhaustive finite-community target was already biased by −.0193 relative to the population reference. The point-to-exhaustive integration bias was −.0043. Finite-site predictive attribution accounted for more of this mean bias than the incomplete integration. The analogous occupancy biases were −.0066, −.0048 and −.0018. The bootstrap itself moved the corresponding order-2 average heterogeneity estimate a further −.0104 relative to the fitted estimate.

| Observed-response secondary importance | Percentile 95% coverage | Basic 95% coverage |
| --- | ---: | ---: |
| Order 2 total | .922 | .750 |
| Order 2 occupancy | .797 | .844 |
| Order 2 heterogeneity | .672 | .953 |
| Order 4 total | .844 | .906 |
| Order 4 occupancy | .828 | .938 |
| Order 4 heterogeneity | .828 | .953 |

Reflection improved several coverages but worsened pairwise total coverage. The correction did not uniformly improve point RMSE. In the 64-community realised null audit, percentile intervals declared a practical advantage in 1 of 64 occupancy contrasts and none of the total or heterogeneity contrasts. The basic reflection declared 4 of 64 heterogeneity advantages versus 0 for the percentile procedure. The bootstrap results do not support switching to a single correction.

For the matched 32-community order-2 conditional-surface audit, removing monotonicity did not improve heterogeneity percentile coverage (.531 in both configurations) or remove its mean bias (approximately −.097 in both). Stronger ridge reduced that bias to −.0475 and raised observed coverage to .906 in those same communities. The change also changed the population regression estimand. This is a mechanism check of shrinkage and predictive overfit, with 39 draws per community and a small repeated sample. It does not establish transferable interval calibration or an empirical penalty choice. The normalised explained deviance of the near-flat null conditional surface was unstable, which is why practical false declarations under independence are based on the realised holdout scale. Predictive attribution remains descriptive.

[Interval calibration and bias decomposition](final_gdum/results/interval_summary.csv)

## Verification on the five manuscript communities

The verification used the existing DD2 data, published preprocessing settings and four spatial blocks. For each community it assessed order 2 and a later feasible order, using matching combinations across architectures and designs. Spline and scaling parameters were estimated on training sites only in each fold. Existing non-bird preprocessing includes choosing an environmental skewness transform on the initial covariate table. Spatial scores are conditional on that choice. The archived empirical results contain diagnostics and predictions aggregated over set combinations. The underlying sites and community tables remain in their existing data location.

All 240 fold fits completed. For coherent additive Forest Hemiptera at order 2, pooled spatial explained deviance was .418 for total, .338 for occupancy and −.154 for heterogeneity. At order 6 it was .388, .229 and .138. Pond diatoms had negative transfer at order 2 and positive scores for all three responses at order 8. Arctic plants had negative pooled scores at both verified orders. Finnish birds showed positive pooled scores for all responses at order 15, with zero evaluation sets passing the configured joint bivariate support checks. This is a contrast between predictive scores and support for interpreting a local response surface. Configuration-specific support should be checked separately before using any surface. The data were used to verify package behaviour, without ecological or inferential conclusions.

Across the empirical networks, additive fixed-site prediction sensitivity between independent pools ranged from .00021 to .0320 RMSE. Scaled interaction condition numbers exceeded 134,000 and included rank-deficient designs. Coherent predictor-only allocations used 24 matched paths above three predictors, with conditional permutation errors reported. The full models and those sampled coalition fits converged in this audit. Such convergence says little about ecological uncertainty. The larger model is an option for hypotheses about filtering, with additive and spread-only alternatives retained for comparison.

[Spatial validation summary](final_gdum/results/empirical_summary.csv) · [Fold results and numerical sensitivity](final_gdum/results)

## Computational cost

The standalone benchmark used 36 matched scenarios of size 70, 110 or 250 across orders 2, 4 and 8, with 900 training and 900 evaluation combinations. It times design preparation, fitting and full hierarchical attribution. Median elapsed seconds for spread-only, additive and interaction were .068, .173 and .341 with the coherent model, and .073, .188 and .457 with independent marginals. Timings include shared pre-processing per comparator and vary with hardware and concurrent workloads. Attribution can dominate fitting cost, and these times are from two-predictor designs. Repeated bootstrap refits need appreciably more time. [Benchmark data](final_gdum/results/benchmark.csv)

## Status and scope

| Status | Result |
| --- | --- |
| Implemented and validated | Component partition, coherent closure and boundary handling, both fractional-logit architectures, nested design construction, spatial validation, independent integration diagnostics, and attribution additivity. |
| Implemented and experimental | Optional spread–position interactions, endpoint regularisation and supported surface displays. Predictive gains were scenario- and penalty-dependent. The diagnostics expose aliasing and limited support. |
| Restricted to descriptive interpretation | Predictor, spread, position and interaction allocations, their rankings and conditional response surfaces. No ecological importance interval is recommended from this study. |
| Unsupported | Universal order or penalty thresholds, joint ecological inference across orders, causal separation of drivers from component models, inference from unsupported surface cells, and interpretation of aliased endpoint coefficients. |

Use the retained `zeta_order_information()` for ecological and design context. Use independent spatial assessment and the extended `ugdm_within_order_diagnostic()` for each response. Coherent occupancy and heterogeneity require meaningful working information in both stages. The diagnostic requires explicit calibration before identifying inferential support. Positive explained deviance by itself does not establish calibrated attribution. The empirical task remains separate from this verification.

## Package checks

Full source-package `R CMD check` passed on arm64 macOS with R 4.6.1. It recorded 1,713 passing expectations, with one pre-existing CRAN-guarded skip. The workflow was also run from the checked package installation. The manuscript files were untouched.

## Reproduction

Install the package and run the scripts from the repository root. The first loop writes a reproducible cache of held-out targets. Set `MS_GDUM_REFERENCE_DIR` to choose the cache location. It defaults to `docs/final_gdum/results/reference_cache`. The CSVs record the reference summaries and recovered trajectories. The transient RDS cache is unnecessary for reading them. On Unix `MS_GDUM_WORKERS` controls parallel simulation workers.

```sh
for scenario in occurrence abundance contrasting null weak correlated spatial; do
  Rscript docs/final_gdum/reference.R "$scenario" 2,4,8
  Rscript docs/final_gdum/simulate.R "$scenario" 110 12 2,4,8 600 pilot
done
for scenario in occurrence abundance contrasting null; do
  Rscript docs/final_gdum/simulate.R "$scenario" 110 200 2,4,8 900 principal
done
for order in 2 4; do
  Rscript docs/final_gdum/uncertainty.R occurrence "$order" 200 59 baseline
  Rscript docs/final_gdum/uncertainty.R occurrence "$order" 64 39 baseline observed
done
for mode in matched free ridge; do
  Rscript docs/final_gdum/uncertainty.R occurrence 2 32 39 "$mode"
done
Rscript docs/final_gdum/uncertainty.R null 2 64 39 baseline
Rscript docs/final_gdum/uncertainty.R null 2 64 39 baseline observed
Rscript docs/final_gdum/profiles.R
Rscript docs/final_gdum/sensitivity.R
Rscript docs/final_gdum/benchmark.R
Rscript docs/final_gdum/reference_surface.R

for i in 1 2 3 4 5; do
  Rscript docs/final_gdum/empirical.R "$DD2" "$i"
done
Rscript docs/final_gdum/summarise.R
Rscript docs/final_gdum/figures.R
MS_GDUM_WORKFLOW_OUT=docs/final_gdum/results Rscript docs/final_gdum/workflow.R
```

The main figure shows architecture differences, supported contextual response, low-order coverage and spatial transfer. The simulated phosphorus workflow uses public package functions for fit, prediction, validation, masking, attribution and confirmation. [Workflow guide](MS_GDUM_CLAUDE_WORKFLOW.md)
