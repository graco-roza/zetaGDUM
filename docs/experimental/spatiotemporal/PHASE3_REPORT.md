# Phase III report

The community-level partition recovered substantial opposing contributions beneath an exactly zero abundance-dissimilarity trend. Across 200 order-two panels, occupancy expansion contributed **−0.1503** and numerical change among persistent occurrences contributed **+0.1503**. The largest absolute panel mean net change was $5.9\times10^{-17}$. The framework also distinguished zero change in complete-sampling incidence zeta under balanced redistribution from changes in local composition.

The environmental model recovered the direction and typical magnitude of the independent-panel contrasts. Its gains in held-out accuracy were clearest for numerical change. Uncertainty varied with the dependence design, and high-order shared-vector distances could be undefined even in nonempty panels. The method is suitable for descriptive accounting and conditional community-level prediction on the supported designs. Ecological process attribution remains conditional on the block allocation and the observed sampling panel.

## Reproduce and use

All outputs are confined to the experimental directory. [Phase III methods](PHASE3_METHODS.md) gives the derivations, population contrasts and estimands. [Phase II](PHASE2_REPORT.md) and [Phase I](REPORT.md) remain available as separate archived experiments.

Run from the repository root with the development `zetaGDUM` package installed. The tests use `testthat`, and the optional published comparison uses `ecopart`. The committed benchmark used R 4.6.1, native zetaGDUM 0.1.0 and ecopart 0.2.0. Cairo is needed for the SVG exports.

```sh
Rscript docs/experimental/spatiotemporal/test_phase3.R
Rscript docs/experimental/spatiotemporal/phase3/run_validation.R
Rscript docs/experimental/spatiotemporal/phase3/run_orders.R
Rscript docs/experimental/spatiotemporal/phase3/run_holdout.R
Rscript docs/experimental/spatiotemporal/phase3/run_sensitivities.R
Rscript docs/experimental/spatiotemporal/phase3/make_figures.R
```

`ST3_OUTPUT` and `ST3_FIGURES` redirect outputs. Set `ST3_PILOT=1` and reduce `ST3_REPS`, `ST3_TRUTH_REPS`, `ST3_GRID_REPS`, `ST3_HOLD_REPS` and `ST3_BOOT` for shorter checks. `ST3_WORKERS` controls parallel simulation. Use the default complete runs before making all eight figures.

```r
source("docs/experimental/spatiotemporal/phase3_load.R")
st3_load()
sim <- st3_simulate("filtering", seed = 7)
windows <- st3_pair_windows(sim$panel$times)
sets <- st3_sample_sets(n = 20L, k = 2L, budget = 36L)
d <- st3_environment_design(sim$panel, sim$env, sets, windows)
y <- st3_collect(sim$panel, sets, windows, crossed = TRUE)
fit <- st3_fit_process_gdum(d, y$y)
prediction <- st3_predict_process_gdum(fit, d, initial = y$before)
```

A single paired contrast can be inspected with `st3_plan()`, `st3_shapley()`, `st3_crossed()` and `st3_site_reorganisation()`. Species identities and process labels are available in the accounting ledger. Models accept a matrix with one row per combination and contrast, in the design's ordering.

## Mathematical validation and ecological accounting

All 227 Phase III assertions passed, including the comparison to native incidence, abundance and dominance zeta and the defined pairwise Baselga components. The largest efficiency error over the primary simulations was $3.2\times10^{-15}$, including the fitted total identity. Finite subset values assure common support within each decomposition. An undefined quantity in an intermediate counterfactual can invalidate all four contributions even when observed endpoints are both defined.

The four pure cases assigned the raw change to their own block. At $(2,1)$ the mean abundance-dissimilarity changes were −0.0174 for expansion, +0.0150 for contraction, +0.0004 for balanced redistribution, and −0.0052 for numerical change. The respective mean absolute combination contributions were 0.0411, 0.0408, 0.0595 and 0.0699. The magnitude differences show why averaging signed changes first can conceal reorganisation. In the pure redistribution case, the exact complete-combination mean incidence-zeta change is zero by the occupancy identity. The separate grid's mean local temporal Sørensen and relative Bray–Curtis distances were 0.0619 and 0.1010.

In the critical opposing case, mean gross attributed magnitude was 0.3007. Mean site-resolved temporal relative Bray–Curtis was 0.383 in the order-grid experiments, despite the unchanged spatial abundance dissimilarity. The model estimated the high-minus-low climate-change contrast as −0.1026 for expansion and +0.1026 for numerical change. Independent Monte Carlo estimates were −0.1097 and +0.1097. The total was structurally zero; the interval tables suppress inference on machine-roundoff values and on structurally inactive components. [Primary panels](phase3/results/main_panels.csv) · [Site-resolved changes](phase3/results/site_reorganisation.csv) · [Independent contrasts](phase3/results/independent_truth.csv)

The abundance decomposition also distinguishes the triggering process from the dimension affected. In the filtering example at $(2,1)$, occupancy-related dissimilarity averaged 0.317 initially and 0.269 finally. Numerical dissimilarity averaged 0.174 initially and 0.208 finally. The numerical process reduced the occupancy-related component by about 0.037 and increased the numerical component by about 0.033. Relative-abundance normalization is one route for these crossed effects. The two entries are retained even when their row sum is small. [Order and crossed results](phase3/results/order_processes.csv)

### Allocation and pairwise comparators

In 24 mixed-process panels with exactly enumerable assignments, changing the allocation altered mean expansion attributions over complete triples by a mean range of 0.057. The mean of the largest combination-specific ranges was 0.476. For expansion, magnitudes of Monte Carlo error against the exact mean fell from 0.00784 with four draws to 0.00355 with 16 and 0.00241 with 64. Their mean estimated MC SEs were 0.01178, 0.00551 and 0.00263. Sampling shares assignment masks across combinations. Assess aggregated assignment uncertainty by aggregating each draw before calculating its SE.

Tatsumi's published pairwise sums and the observed change agreed to within $10^{-15}$ on their defined support. For relative Bray–Curtis the published mean loss homogenisation and differentiation magnitudes were 0.078 and 0.129, and gain magnitudes were 0.110 and 0.087. The four-process Shapley positive and negative magnitudes averaged 0.101 and 0.073. These allocations use different reference paths and process groups. The signed totals are the meaningful equality benchmark. [Allocation sensitivity](phase3/results/allocation_sensitivity.csv) · [Pairwise comparisons](phase3/results/pairwise_comparator.csv)

## Environmental validation

The primary study fit 200 panels for each of 12 ecological and null scenarios plus 200 shared-shock panels. The independent reference used 240 further panels per scenario. Across 22 supported nonzero process contrasts, the correlation between the independent Monte Carlo truths and mean model estimates was 0.993, median absolute bias was 0.0030 and median replicate-level RMSE was 0.0127. This is calibration to ecological contrasts across simulated panels on a fixed environmental design. Individual fits and prediction on novel environments remain less precise.

The table shows mean holdout RMSE for relevant $(2,1)$ abundance outcomes, averaged across the 13 scenarios with 24 held-out panels per scenario. Complete site histories use separate accounting panels for training and testing. Forward scoring uses the final four surveys.

| Response | Held sites, additive | Held sites, intercept | Forward, additive | Forward, intercept |
| --- | ---: | ---: | ---: | ---: |
| Expansion | 0.0533 | 0.0546 | 0.0504 | 0.0510 |
| Contraction | 0.0320 | 0.0308 | 0.0270 | 0.0262 |
| Redistribution | 0.0122 | 0.0117 | 0.0166 | 0.0162 |
| Numerical | 0.0589 | 0.0666 | 0.0575 | 0.0658 |
| Net change | 0.0908 | 0.0924 | 0.0894 | 0.0926 |

The optional position-by-change interaction worsened mean held-site net RMSE from 0.0908 to 0.0917 and forward RMSE from 0.0894 to 0.0898, so the additive model remains primary. About 21% of these evaluation rows were outside at least one training marginal range, and approximately 0.003% of held-site rows produced an incompatible normalized endpoint. These flags identify bounds or support problems without altering forecasts. More stringent joint-support checks and new environmental and sampling designs remain to be tested. [Holdout panels](phase3/results/holdout_panels.csv) · [Signed response curves](phase3/results/environmental_curves.csv)

### Dependence and interval coverage

| Uncertainty method | Clean null net coverage | Shared shocks net coverage |
| --- | ---: | ---: |
| Spatial-block deletion | 95.5% | 78.5% |
| Moving temporal blocks | 87.5% | 82.0% |
| Crossed resampling | 100.0% | 95.5% |

These are 200-panel results per regime. The shared-shock net sampling SD was 0.0523, compared with mean SE 0.0343, 0.0391 and 0.0590 for the three corrections respectively. Crossed resampling was conservative in the clean null and under the tested alternative contrasts, where coverage ranged from 98% to 100%. Under shared shocks, crossed rejection rates for zero expansion, numerical change and net change were 2.0%, 1.5% and 4.5%. The MC SE for the 95.5% coverage estimate is about 1.5 percentage points. These results are conditional on the tested block structure. Eight paired contrasts and five spatial blocks give limited support for inference under different shock or spatial dependence processes. [Interval summaries](phase3/results/interval_summary.csv) · [Panel-level uncertainty](phase3/results/interval_panels.csv)

## Temporal alignment and higher orders

The grid evaluated $k=2$–6 and $m=1$–4 on 12 panels per scenario. In the filtering example, the mean effective model degrees of freedom dropped from 12.06 at $(2,1)$ to 6.02 at $(6,4)$ although the design rank stayed 17. The respective grid evaluations took about 0.10 and 0.41 wall-clock seconds per panel in this run. Retained model objects occupied about 0.14 MB and 0.08 MB, respectively. The high-order fit used fewer windows, and these object sizes exclude peak working memory. The larger temporal orders have only two paired windows under the chosen schedule, so the report omits their predictive contrast uncertainty.

The site-specific anchor stress test retained positive abundance in every original and counterfactual assemblage. Temporal Sørensen and temporal shared-relative Bray–Curtis still lost common support, from 96.3% at $(2,1)$ to 1.6% at $(6,4)$, because the initial rectangle often had zero jointly shared species. Spatial/joint abundance dissimilarity remained defined. A high-order interval with too few windows or too few shared species cannot be recovered by adding combinations of the same observations.

Shifting the aligned persistence windows at $(3,3)$ changed the mean redistribution contribution to raw joint zeta from −0.063 to −0.036 species. Increasing the gap gave +0.378. Irregular time labels changed mean spans and lags but left identical survey matches unchanged, as expected. These contrasts combine occupancy and time information. Design selection should explicitly consider alignment, duration, lag, common support and held-out uncertainty. [Order support and cost](phase3/results/order_support.csv) · [Window sensitivity](phase3/results/window_sensitivity.csv)

The results support process-based ST-MS-GDUM as a community-level accounting and prediction framework. Analysis of empirical panels should report allocation sensitivity, keep the temporal distance estimands beside signed spatial change, and calibrate uncertainty for the actual site and survey design. The experiment assumes matched sites, comparable effort and complete detection. Missing surveys, inconsistent detection, inferred movement and prediction beyond the observed regional panel require separate evidence.

## Figures

| Figure | Diagnostic |
| --- | --- |
| [1](phase3/figures/01-framework.svg) | Four intervention blocks and three compositional summaries |
| [2](phase3/figures/02-abundance-components.svg) | Occupancy-related and numerical dissimilarity over both orders |
| [3](phase3/figures/03-signed-processes.svg) | Process contributions and separate net change |
| [4](phase3/figures/04-crossed-processes.svg) | Crossed abundance contributions at representative orders |
| [5](phase3/figures/05-opposing-processes.svg) | Opposing ecological processes with zero net change |
| [6](phase3/figures/06-environmental-curves.svg) | Signed GDUM response curves |
| [7](phase3/figures/07-validation.svg) | Contrast recovery, holdout accuracy and coverage |
| [8](phase3/figures/08-order-selection.svg) | Order support, information and predictive uncertainty |
