# Pre-fit information and modelling reliability by zeta order

`zeta_order_information()` profiles the variation available to MS-GDUM before fitting any environmental model. It identifies exact degeneracy, shows how much variation remains, and reports numerical precision. Across the simulation families, no tested threshold made reliable inferential recommendations on independent communities. The package therefore leaves candidate orders for post-fit confirmation using independent predictions and locally calibrated uncertainty. [Package function](../man/zeta_order_information.Rd) · [Figure 3 PDF](order_information/results/figure3_order_selection.pdf) · [Figure 3 SVG](order_information/results/figure3_order_selection.svg)

## Using the diagnostic

```r
profile <- zeta_order_information(community, environment, coordinates,
                                  orders = 2:30, response = "all", seed = 71)
plot(profile)
attr(profile, "environment_detail")
attr(profile, "candidate_orders")
```

Rug marks on the relative-information plots indicate candidate orders. The profile keeps the original order grid, so isolated candidates and gaps remain visible. Inspect sampled variances with their standard errors and increase the numerical budget where the candidate interpretation is sensitive to integration. `candidate` calls for post-fit checking. `verify_variation` identifies a numerically flat estimate that has not been proved constant. In the plotted condition panel, crosses mark rank-deficient designs. Use `ugdm_within_order_diagnostic()` to check the fitted model against independent predictions, repeated integrations and a calibrated uncertainty procedure.

## Compositional information

The target is the distribution over uniformly chosen sets of $k$ distinct sites. Let $Z(S)$ be shared richness, $r_i$ site richness, $q_{is}$ the site-relative abundance of species $s$, and $a_s(S)=\min_{i\in S}q_{is}$. The four response variances use the following definitions.

| Response | Definition |
| --- | --- |
| Incidence | $1-Z(S)/\bar r_S$ |
| Relative abundance | $1-\sum_s a_s(S)$ |
| Occupancy-associated | $\sum_{s:a_s(S)=0}\bar q_{s,S}$ |
| Abundance heterogeneity among shared species | $\sum_{s:a_s(S)>0}(\bar q_{s,S}-a_s(S))$ |

The last two responses sum to relative-abundance dissimilarity for each set. Their variances add with twice their covariance. The partition describes dissimilarity attributable to species absence in at least one site and to unequal relative abundance among species present in every site. A species missing from part of a set contributes its full mean abundance to the occupancy component. This is a descriptive partition and does not identify ecological mechanisms.

For each response report absolute variance $V_k$ and $I_k=V_k/V_2$. An exactly constant pairwise response leaves $I_k$ undefined, even when later orders vary. The index can increase or exceed one. A zero estimated from sampled pairs also leaves the ratio undefined and is flagged by the pairwise calculation method. A small denominator can make the ratio unstable. The complete reporting therefore keeps both absolute and relative variation.

The exact occurrence moments, overlap bounds for sharing, mean abundance zeta, retention and relative contributions reuse `zeta_order_support()`. Incidence-based zeta variance was established by McGeoch et al. (2019). Here it gives the variance of incidence dissimilarity when site richness is constant, scaled by the squared common richness. Different site richness makes the incidence dissimilarity a ratio with a variable denominator. Exhaustive calculations either enumerate small populations or use analytical cases where available. [McGeoch et al. (2019)](https://doi.org/10.1002/ecy.2832)

To calculate relative-abundance variance when the number of distinct eligible thresholds is small, write each abundance minimum as a weighted sum of indicators. A threshold contributes its increment when the entire set lies in its occupied threshold set $A$. For threshold sets $A$ and $B$, the indicator covariance is

$$
\frac{\binom{|A\cap B|}{k}}{\binom{n}{k}}
-\frac{\binom{|A|}{k}\binom{|B|}{k}}{\binom{n}{k}^{2}}.
$$

Thresholds present at fewer than two sites cannot contribute at the evaluated orders. The algorithm caps retained thresholds because the cross-product grows quadratically. Near cancellation in a variance formula, the implementation integrates the centred response directly. Differences at the floating-point roundoff scale are numerically unresolved. The exact mean shared-species abundance before subtracting the minima is

$$
\sum_s \bar q_s\,\frac{\binom{o_s-1}{k-1}}{\binom{n-1}{k-1}}.
$$

This supplies exact means for both abundance components. Zero or constant cases also give exact variances. Generic component covariance and richness-normalised incidence variance use enumeration or independent uniform pools. Pool summaries use unbiased sample variances and their across-pool standard error. For approximate ratios, the delta-method uncertainty combines separately sampled numerator and denominator. Sharing has exact bounds, and sampled sharing also has a Wilson interval so that an observed zero does not imply certainty. The uniform target requires positive finite totals at all sites for abundance diagnostics.

## Environmental information

For any fixed site-level feature $f$, including linear position and each nonlinear MS-GDUM basis function,

$$
\operatorname{Var}(\bar f_S)=s_f^2(1-k/n)/k.
$$

The public profile reports exact raw position variances per predictor and summary variances of the standardised and nonlinear design columns. It uses one basis across orders. The unchanged site values are the population in these formulas.

For a pairwise spread kernel $h_{ij}$, centre edges by their finite-population mean. Write $A=\sum_{i<j}h_{ij}^2$, $B=\sum_i(\sum_{j\ne i}h_{ij})^2$, and $p_m=\binom{k}{m}/\binom{n}{m}$. Then

$$
\operatorname{Var}\left[\binom{k}{2}^{-1}\sum_{i<j\in S}h_{ij}\right]
=\frac{(p_2-2p_3+p_4)A+(p_3-p_4)B}{\binom{k}{2}^2}.
$$

An impossible edge pattern has probability zero. This counts disjoint and incident edges without enumerating sets. Absolute pair differences use sorted prefix sums, giving exact variance in $O(n\log n)$ work per predictor or I-spline column. Geographic distances use the same edge expression. Joint rank and scaled condition are calculated on the actual joint design, using enumeration or pooled independent combinations. The diagnostics distinguish shrinking column variance from collinearity among scaled columns. Constant or sparsely valued predictors are flagged because the package spline requires four distinct values.

## Calibration design

The same ecological realization was fitted across a continuous sequence of orders in the pilot. The six niche-process scenarios changed regional richness, occupancy distributions, abundance–occupancy coupling, within-species abundance heterogeneity, environmental niche breadth and effect strength. The withheld scenarios included strong predictor correlation, spatial dependence, and a null environment–composition relationship. Site counts were 70, 150 and 250. The continuous pilot ran at orders 2–30 and extended to 40 for broad and dominant communities. Both incidence and total relative-abundance fits used constant working weights. The occupancy and heterogeneity profiles were calculated for every simulated community. The modelling reliability study calibrates only the incidence and total relative-abundance responses.

Three principal scenario-size combinations each used **200 independent ecological communities**, with matching seeds across orders and responses. The order grid concentrated additional inference near the early precision changes and bracketed later losses. Each combination used 900 training sets, a second independent pool of 900, 79 Gaussian score draws, and 1,500 evaluation sets on separate simulated sites. The validation held back correlated, spatial and null scenarios from rule fitting. Besides the continuous pilot on those scenarios, a targeted validation sample independently simulated 64 ecological communities for each withheld scenario-size combination. Spatial score draws used contiguous blocks. Block results depend on their grouping and provide limited evidence for spatial calibration.

Three independently simulated 2,200-site populations and 16,000 combinations per pool defined each scenario, response and order reference. A model fitted on one simulated community was evaluated on a different one, cycling through the three. Reference spread and usable counts are archived. The estimand is held-out absolute explained-deviance Shapley attribution in the additive constrained nonlinear MS-GDUM. It is conditional on the specified predictor basis and ridge penalty. It is not a direct ecological effect parameter. With two predictors, the four predictor-coalition refits give exact Shapley values for that game. Paired site perturbations for the four coalitions used the package's projected Gaussian score approximation. The optimized calibration code agreed with the public Shapley implementation to $1.5\times10^{-6}$ in 18 benchmark comparisons. [Scenario definitions](order_information/design.R) · [Population targets](order_information/results/population_targets.csv) · [Model reliability](order_information/results/reliability.csv)

### Independent modelling criteria

The modelling criteria never use $I_k$, retention, predictor variance or any pre-fit classifier. A reference predictor difference must exceed .02 absolute explained deviance to count as meaningful. With a negligible reference difference, recovery of a named leader is not required. The summary reports false declarations of a practically meaningful difference.

- **Unsupported** requires failure of predictive transfer or fit availability. Prediction requires median held-out explained deviance at least .025, positive held-out deviance in at least 75% of replicates, and at least 95% converged fits.
- **Inferentially supported** additionally requires usable intervals in 90% of replicates. Mean absolute bias must be at most .03, RMSE at most .06, numerical RMSE at most .04, and mean 95% interval width at most .20 for both allocations and their contrast. Observed coverage must reach .70 and .85 for nominal 80% and 95% intervals. Correct leadership must reach 85% when the reference gap exceeds .02. Practical superiority power must reach 80% for a reference gap of at least .04. False or reversed declarations must stay at or below 10%.
- **Prediction only** identifies the cells that meet the prediction criteria with unresolved allocation or contrast inference. The archive separately flags cells meeting meaningful-contrast requirements with poor coverage for individually near-null allocations. This distinction exposes the low-order calibration problem.

These cutoffs are operational definitions for this simulation study. They are sufficient to compare the screens on an identical reference, with appreciable classification noise for the eight-community pilot. The code retains coverage, width, bias and other component metrics so users can inspect another scientific tolerance. A stringent threshold or a reference close to the decision boundary can change a case's class.

## Evidence from matched orders

Of 48 principal order-response cells, seven met the strict inference criteria, 38 supported prediction only, and three were unsupported. Limiting orders depended on the response and community. At 250 sites the dominant scenario had supported relative-abundance orders 4, 6, 8, 12 and 28 on the evaluated grid. Order 18 between them failed precision or calibration, and order 2 undercovered the strong allocation and contrast. The same scenario supported incidence inference at order 12 on the evaluated grid, while held-out incidence prediction persisted at order 40. Untested intermediate orders cannot be inferred from these lists. [Per-cell estimates](order_information/results/reliability.csv) · [Supported and screened order lists](order_information/results/selection_maxima.csv)

At the same site count of 150, broad communities at order 40 had abundance retention .995 and sharing probability one, yet median abundance zeta was .015 and only about 1.05 species contributed effectively. Relative-abundance variance was approximately $1.9\times10^{-6}$, or .00041 of pairwise variance. The dominant scenario at order 40 retained mean abundance zeta .979, but only 1.38 effective species and absolute response variance $2.8\times10^{-6}$. This was .084 of an already small pairwise variance. The linear design variance was .018 for both scenarios on the common standardised scale, with pairwise spread variance about .0055. Retention and absolute zeta cannot substitute for compositional variation. The abundance-relative index alone can also obscure small starting variance. [Figure 3](order_information/results/figure3_order_selection.pdf) · [Absolute variation and sensitivity](order_information/results/figureS1_order_sensitivity.pdf)

A matched negative control shuffled the two environmental predictors together across sites and made an independent shuffle on the evaluation sites. This preserves the response and all exact unsupervised environmental variation. The distribution of the joint environmental design under uniform sets is also invariant. It eliminates the pairing of environmental and compositional differences. Median held-out explained deviance at orders 2 and 24 fell from .628 and .090 to −.006 and −.003. Rank and condition estimated from finite numerical pools can differ slightly after shuffling. This experiment gives a direct reason the proposed marginal screens cannot establish prediction or predictor importance. [Matched permutation control](order_information/results/permutation_control.csv)

### Low-order uncertainty audit

At order 2 in the broad 70-site relative-abundance cell, held-out explained deviance was .413 and nominal 95% coverage of the meaningful contrast was .985. Coverage of the near-null secondary allocation was .180. Its point bias was −.0217 explained deviance. With a matched 32-community subset, order-2 secondary coverage was .188 under direct conditional-Poisson Bayesian bootstrap, .375 with larger score pools and .438 with a larger independent evaluation pool, all below the nominal level. Larger combination budgets reduced numerical variability but did not restore calibration. The direct bootstrap confirmed that this failure was not specific to linearised score updates. The narrow 150-site relative-abundance cell had .780 and .680 coverage for the meaningful contrast at orders 2 and 4, respectively, despite held-out explained deviance .652 and .613. Low-order statistical failure is reported separately from information loss. [Matched audits](order_information/results/reliability.csv)

In the withheld null communities the median held-out deviance was negative across the evaluated orders. The practical false-superiority rate was zero for relative abundance in this validation run. Those figures do not imply calibration for other null mechanisms or nominal one-predictor intervals.

## Order-selection rules on independent data

Thresholds were trained on medians from the continuous training pilot in the narrow, broad and dominant scenarios. Single-feature rules searched simple cutoffs in both directions. The combined rule searched conjunctions of absolute response variance, environmental variation and scaled design reciprocal condition. The retention breakpoint rule used orders at or below the package's ecological breakpoint. Balanced errors on training orders determined the threshold. Independent validation aggregated profiles by scenario-size-response before comparison to the separately estimated reliability labels. The archive additionally applies the frozen thresholds to each individual held-out community. Community-level evaluations share the scenario-specific reliability reference.

| Screening rule | New communities overclaimed | Informative orders excluded | Highest-order MAE | Unseen scenarios overclaimed | Informative orders excluded | Highest-order MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Zeta decline | 61.5% | 28.6% | 4.33 | 50.0% | 70.0% | 2.33 |
| Retention rate | 83.7% | 0% | 21.33 | 79.2% | 0% | 26.33 |
| Retention breakpoint | 91.7% | 71.4% | 15.00 | 70.6% | 0% | 7.00 |
| Absolute response variance | 64.3% | 28.6% | 22.00 | 100% | 100% | 21.50 |
| Relative response variance | 78.6% | 14.3% | 7.50 | 69.2% | 20.0% | 7.33 |
| Environmental variation | 80.0% | 14.3% | 8.33 | 70.6% | 0% | 7.00 |
| Combined diagnostics | 82.4% | 14.3% | 9.67 | 71.9% | 10.0% | 7.00 |

Overclaiming is the fraction of recommended orders without strict inference support. It includes prediction-only orders. The combined rule recommended 34 principal and 32 withheld cells. Zero of its principal recommendations and 13 of the withheld ones were outright unsupported for prediction. The highest-order error uses only the evaluated grid and encodes no supported order as one. Consequently a smaller error can coexist with systematically missing valid orders, as with zeta decline on the withheld scenarios. The tuned combined inference rule dropped its response variance cutoff, relying on geometry and conditioning. This is evidence against generalising the fitted threshold. [Frozen rules](order_information/results/rules_inference.csv) · [Full accuracy table](order_information/results/selection_accuracy.csv) · [Individual community error rates](order_information/results/selection_community_accuracy.csv)

Even the separate predictive screen transferred poorly. Absolute variance recommended 45 of 48 cells in new communities in familiar scenarios, with one unsupported recommendation and one excluded predictive cell. On completely held-out scenarios it recommended 42 cells, including 18 unsupported ones. The combined predictive screen had similar transfer errors. Differences in sample size and response should be read from the per-scenario files, because the principal comparison jointly varies community structure and site count. The separate continuous pilot crosses all three site counts with every scenario. The empirical use supported by these results is structural screening, numerical budgeting and prioritising which candidates to fit. Calibrated ranking requires independent fits and ecological uncertainty. A universal maximum order or threshold is unsupported.

## Numerical checks and computational cost

Exhaustive enumeration verified the two response partitions, exact and enumerated variance, undefined pairwise normalisation, ratios above one, absolute and nonlinear position variance, directional I-spline spread variance, and geographic spread variance. Independent-pool checks verified repeatability and nonzero numerical standard errors. At a 29-order grid with equal two-pool budgets of 500 sets, median profile time was about .87–.90 seconds across site counts. Fitting the two-predictor game and producing 39 score draws with held-out assessment took 1.55–1.60 seconds, a 1.77–1.81-fold difference. There is no evidence for a larger speed-up on this small model. Potential savings from eliminating fits depend on whether a defensible screen can exclude orders on the actual dataset. [Timing benchmark](order_information/results/benchmark.csv)

`R CMD check` completed with status OK, including examples, vignettes and the PDF manual. The package check recorded 1,678 passing test expectations and one pre-existing CRAN-guarded skip. The manuscript was left untouched.

## Reproduction and outputs

Run from the repository root with the package installed. `ZETA_ORDER_WORKERS` controls parallel R workers for simulations on Unix. The stored CSV and compressed CSV files retain condition-level and replicate-level numerical results, including binary recovery, all interval endpoints and recorded failures. Normal reference results use 16,000 sets. Pilot results use eight independent communities per scenario, site count, response and split, with 19 score draws and 500 combinations. Only the narrow, broad and dominant scenarios contribute a training split.

```sh
for scenario in narrow broad dominant correlated spatial null; do
  Rscript docs/order_information/reference.R "$scenario"
  for n in 70 150 250; do
    for response in total relative; do
      case "$scenario" in
        narrow|broad|dominant)
          Rscript docs/order_information/calibration.R "$scenario" "$n" train "$response" 8 19 500 ;;
      esac
      Rscript docs/order_information/calibration.R "$scenario" "$n" test "$response" 8 19 500
    done
  done
done

for response in total relative; do
  Rscript docs/order_information/calibration.R narrow 150 principal "$response" 200 79 900 2,3,4,6,8,12,20,30
  Rscript docs/order_information/calibration.R broad 70 principal "$response" 200 79 900 2,4,6,8,12,18,28,40
  Rscript docs/order_information/calibration.R dominant 250 principal "$response" 200 79 900 2,4,6,8,12,18,28,40
  Rscript docs/order_information/calibration.R correlated 150 validation "$response" 64 79 900 2,3,4,6,8,12,20,30
  Rscript docs/order_information/calibration.R spatial 250 validation "$response" 64 79 900 2,3,4,6,8,12,20,30
  Rscript docs/order_information/calibration.R null 70 validation "$response" 64 79 900 2,3,4,6,8,12,20,30
done

Rscript docs/order_information/calibration.R narrow 150 budget relative 32 79 3600 2,4,6,8,12
Rscript docs/order_information/calibration.R broad 70 budget relative 32 79 3600 2,4,6,8,12
Rscript docs/order_information/calibration.R dominant 250 budget relative 32 79 3600 2,4,6,8,12
Rscript docs/order_information/calibration.R broad 70 direct relative 32 79 900 2,4,6
Rscript docs/order_information/calibration.R broad 70 assessment relative 32 79 900 2,4,6
Rscript docs/order_information/calibration.R narrow 150 assessment relative 32 79 900 2,4,6

for scenario in narrow broad dominant correlated spatial null; do
  Rscript docs/order_information/profiles.R "$scenario"
done
Rscript docs/order_information/reliability.R
Rscript docs/order_information/selection.R
Rscript docs/order_information/benchmark.R
Rscript docs/order_information/permutation_control.R
Rscript docs/order_information/figures.R
```

Figure 3 panels A through C use independent continuous-pilot trajectories for the same site count of 150. Panel D uses the larger principal and withheld calibration samples at their evaluated orders. Gaps in the evaluated grid remain blank. Supplementary Figure S1 shows estimand-specific coverage, integration budget sensitivity and the absolute scale of the compositional variance. [Supplementary PDF](order_information/results/figureS1_order_sensitivity.pdf) · [Supplementary SVG](order_information/results/figureS1_order_sensitivity.svg)
