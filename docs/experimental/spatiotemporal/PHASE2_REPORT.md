# Phase II feasibility report

M1 improved held-out shared-event predictions over the Phase I pooled model, M0, at spatial orders two and three. Its predictions differed for the same site–species history in different combinations, and some probabilities violated spatial and temporal nesting. The dynamic benchmark, M2, had better held-out scores and produced coherent surfaces when forecast from a common origin. The factorial experiment did not establish recovery of the affected commonness group without its true labels. Shared survey shocks also defeated nominal 95% inference under every uncertainty method tested.

Keep Phase I as a baseline, M1 as a set-specific predictor, and M2 as the generative benchmark. The results support the exact metrics and conditional prediction. They leave ecological attribution and inference under temporal shocks unresolved.

## Scope and reproduction

The [methods](PHASE2_METHODS.md) specify the estimands, likelihood, assumptions and validation design. The first two commits implement exact metrics and transition identities, then the models and community prediction. The third adds the analyses and figures. Phase I's code and outputs remain accessible through its original loader and scripts.

Run the following from the repository root with the development `zetaGDUM` installed. The R packages `testthat`, `lme4`, `ecopart` and `glmmTMB` support tests, dynamic fits and the published beta and count comparisons. The run used R 4.6.1, `zetaGDUM` 0.1.0 and `ecopart` 0.2.0. On this machine `glmmTMB` emitted a TMB build-version warning when loaded. The held-out count comparison and the count prediction test ran successfully. Cairo produces the SVG and PDF figures. Results retain the committed seed protocol. Runtime and memory are platform-dependent.

```sh
Rscript docs/experimental/spatiotemporal/test_phase2.R
Rscript docs/experimental/spatiotemporal/phase2/run_factorial.R
Rscript docs/experimental/spatiotemporal/phase2/summarize_factorial.R
Rscript docs/experimental/spatiotemporal/phase2/run_calibration.R
Rscript docs/experimental/spatiotemporal/phase2/run_predictions.R
Rscript docs/experimental/spatiotemporal/phase2/run_trajectories.R
Rscript docs/experimental/spatiotemporal/phase2/run_diagnostics.R
Rscript docs/experimental/spatiotemporal/phase2/make_figures.R
```

`ST2_WORKERS` controls parallel jobs on supported systems. `ST2_OUTPUT` redirects simulation tables and `ST2_FIGURES` redirects graphics. Pilot settings `ST2_PILOT=1`, `ST2_FACTOR_REPS`, `ST2_PRED_REPS`, `ST2_NULL_REPS` and `ST2_DRAWS` are available for debugging. Use the full default outputs to regenerate the complete figure set.

## Mathematical and order diagnostics

The complete Phase II suite passed 422 assertions. Single-survey products, minima, dominance and the species-level occurrence–abundance decomposition agreed with the native package at every tested spatial order. Exhaustive rectangles confirmed the analytical means for all three responses. Tests cover axis-wise decline, invariance to within-assemblage scaling, abundance-only incidence invariance, the dominance–incidence bound, undefined relative abundance, gain–loss budgets, exact Poisson-binomial mass, rank and window checks, and the additive spatial beta comparisons.

The diagnostic panel illustrates that temporal sampling changes the estimand. At $(k,m)=(5,4)$, unrestricted windows had mean span 7.8 survey intervals and consecutive windows spanned 3.0. Mean incidence zeta was 0.692 and 1.798, abundance zeta 0.00422 and 0.01200, and dominance zeta 0.00410 and 0.01175, respectively. In the consecutive design, the sampled zero-sharing frequency was 16.8% for incidence and abundance, and 98.6% for dominance. Effective species counts were 22.29, 15.62 and 8.25, respectively. Larger curvature for a small-scale continuous response can reflect a smaller working residual variance. The three response scales need separate calibration.

At $(5,4)$, the 64-set diagnostic kept about 7.4 MB of R objects for the paired incidence and dominance fits and took roughly 1.4 seconds on this run, including response diagnostics. This is retained object size, with peak memory and exhaustive large-order transition fits outside the scaling benchmark. For order six and 100 draws, the native species-informed proposal increased positive-sharing sets from 78.7% to 88.4% and reduced the between-draw SD of the uniform-target richness estimate from 0.094 to 0.075, at mean effective set size 85.4. This is Monte Carlo precision conditional on one panel and anchor window. [Order diagnostics](phase2/results/order_diagnostic.csv) · [Sampling precision](phase2/results/sampling_precision.csv)

## Held-out prediction and coherence

In four replicate panels each for the null, additive, interaction and opposing scenarios, the mean shared-event log loss over gains and losses at $(2,1)$ was:

| Holdout | M0 pooled | M1 occupancy-conditioned | M2 dynamic |
| --- | ---: | ---: | ---: |
| Complete site histories | 0.5069 | 0.4937 | 0.4099 |
| Regions | 0.5102 | 0.4970 | 0.4089 |
| Forward transitions | 0.5065 | 0.4924 | 0.4104 |
| Joint space and forward time | 0.5078 | 0.4934 | 0.4099 |

Across the assessed cells through $(3,3)$, mean log loss and Brier score matched at spatial order one and favoured M1 over M0 at the larger orders. Available cells differ across holdouts because windows must fit the training and test horizon. At $(2,1)$ the within-predictor interaction reduced mean held-out log loss by 0.0204 in the planted-interaction scenario and increased it by 0.000015 under the clean null. Exploratory higher-order scoring and unequal-support cases are in the replicate tables. At $(2,1)$ M1 and M2 mean Brier scores were 0.130 and 0.106 for expansion, 0.127 and 0.103 for contraction, and 0.0296 and 0.0271 for pure redistribution. Balanced redistribution is scored against its count; at order two the count equals the pure redistribution indicator. M0 has only shared-event scores at larger spatial orders because its likelihood leaves local outcomes undetermined.

For order-two, one-step community forecasts, M1 and M2 mean absolute errors were 0.137 and 0.064 for abundance zeta, 3.41 and 1.88 species for dominance zeta, 0.074 and 0.057 for native turnover, and 0.043 and 0.030 for native nestedness. The corresponding M2 joint site-and-future forecast trajectory errors at $(2,3)$ were 0.782 species for incidence, 0.00745 for abundance and 0.201 species for dominance. The trajectory scores use a fixed observed initial assemblage. Future environmental values were supplied at evaluation time.

For M1, 99.12% of 3.94 million tested repeated site–species–window contexts showed different probabilities across combinations. Its maximum gap was 0.816. Of the nested predictions, 8.11% of spatial comparisons and 8.83% of temporal comparisons violated ordering, with maximum excesses 0.542 and 0.313. M0 also showed nesting violations. M2 showed zero violations on common-origin simulated trajectories and zero repeated-context differences in its conditional predictions. Temporal comparisons for M0 and M1 can involve different observed conditioning histories, so the latter violations describe incompatibility of the predictions when interpreted as a single unconditional surface. The direct probabilities were kept as fitted. [Scores](phase2/results/prediction_scores.csv) · [Calibration bins](phase2/results/prediction_calibration.csv) · [Coherence](phase2/results/prediction_coherence.csv) · [Trajectory predictions](phase2/results/trajectory_predictions.csv)

## Factorial ecological validation

The full experiment used 34 scenarios with four independent 60-site, 12-survey panels per scenario and four common-random-number forecast draws per panel for simulation-based predictions. This includes all 24 single-target pathways and the non-targeted and stress scenarios. Forecasts used fitted species intercepts only where provided by M2. The true group memberships were used after fitting to score target effects. Results distinguish conditional-window and trajectory-surface effects.

| Targeted pathway | M1 window RMSE / zero forecast | M2 trajectory RMSE / zero forecast | M1 highest-scoring group was target | M2 trajectory highest-scoring group was target |
| --- | ---: | ---: | ---: | ---: |
| Gains, scored by incidence | 0.982 | 0.991 | 4/32 | 2/32 |
| Losses, scored by incidence | 1.001 | 0.976 | 5/32 | 3/32 |
| Positive abundance, scored by relative minima | 0.998 | 0.989 | 6/32 | 4/32 |

Group ranking uses RMS fitted effect over the available orders. Its chance reference under equally likely groups is 4/32. For the abundance pathway, changes in local normalization can move non-target species’ relative minima, so even true response ranking can favour a non-target group. The aggregate relative-abundance surface alone cannot establish which group drove an effect.

Pooled truth surfaces for different target groups were sometimes closer than the typical within-group variation across four panels. This happened in four of 28 group pairs for gains, four for losses and two for abundance. This is evidence of limited finite-panel discrimination, rather than a proof of global structural non-identifiability. The clearest held-out failures were opposing and nonlinear effects, where pooled fitted change contrasts could be near zero despite nonzero true effects. At $(2,2)$, the planted interaction had a pooled true incidence contrast of −5.35 species. The conditional M1 estimate was −0.43 and the conditional M2 estimate −5.39. Strong prediction on some scenarios therefore coexisted with poor recovery of mixed species effects. A more detailed model would need independent evidence of heterogeneity and validation of any inferred species classifications. [Replicate contrasts](phase2/results/factorial_contrasts.csv) · [Group scoring](phase2/results/factorial_group_discrimination.csv) · [Surface similarity](phase2/results/factorial_surface_similarity.csv)

## Dependence, detection and beta change

There were 200 clean null panels and 200 panels with shared temporal and regional shocks. The table shows coverage for gain slopes at order two. Monte Carlo SE for coverage of 0.815 is 0.027 and for 0.470 is 0.035.

| SE method | Clean null coverage | Shared-shock coverage |
| --- | ---: | ---: |
| Independent trials | 89.0% | 8.0% |
| Delete site | 98.0% | 47.0% |
| Delete spatial block | 96.5% | 59.0% |
| Moving temporal blocks | 83.0% | 73.0% |
| Crossed blocks | 100.0% | 81.5% |

Under shocks the true sampling SD for order-two gain coefficients was 0.185, versus average estimated SE 0.128 for crossed blocks. The corresponding loss coverage was 81.0%. The fixed +0.5 SD contrast had 84.5% gain coverage under the crossed procedure. Aliased, boundary and nonconverged primary fits were absent in these runs. The data support neither nominal temporal inference with 12 shared-shock surveys nor the independent-trial approximation. Exact deletion was run on the first four panels per regime and one-step deletion on the rest. Results quantify that approximation and block design jointly. [Null summaries](phase2/results/null_summary.csv) · [Panel-level estimates](phase2/results/null_panels.csv)

In 16 matched detection panels, climate-dependent detection shifted the mean apparent gain slope from −0.006 to +0.539 and the loss slope from −0.005 to −0.785. Effort and detection require separate data or a replicated observation model in later work. In six held-region count panels the zero-truncated NB sensitivity had positive-count MAE 7.41, versus 7.43 for lognormal. The fixed-neighbourhood sensitivity had essentially the same site-gain log loss as main M2 (0.42188 versus 0.42173). These runs were checks on supported variants, with limited replication. [Detection sensitivity](phase2/results/detection_sensitivity.csv) · [Model sensitivities](phase2/results/model_sensitivities.csv)

In the 50-pair beta benchmark, the published Tatsumi Sørensen loss components averaged −0.111 for homogenisation and +0.183 for differentiation. The gain components averaged −0.185 and +0.117. For Bray–Curtis the corresponding means were −0.162, +0.200, −0.176 and +0.183. The separate magnitudes show cancellation that a net change obscures. The native two-event accounting and the published partitions both recover endpoint change where their formulas are defined. Their allocations answer different questions. [Native accounting](phase2/results/native_beta_change.csv) · [Published pairwise comparison](phase2/results/tatsumi_pair.csv)

## Acceptance criteria

| Criterion | Outcome |
| --- | --- |
| Required mathematical identities and native spatial agreement | Passed on the automated suite. |
| Improve or match M0 for shared-event prediction | Passed on mean held-out scores for the tested cells. M0 remains the baseline. |
| Optional interaction recovers planted predictive improvement with a clean-null check | Passed for the pre-specified order-two, one-transition score. |
| Inference calibrated against independent panels | Mixed under the clean null and failed under shared shocks. |
| Quantify cross-order incoherence | Completed. M1 failed joint coherence, while M2 passed the common-origin trajectory checks. |
| Recover target ecological effects without supplied groups | Failed. Effect errors were close to the zero-contrast forecast and group ranking was near chance. |

Further empirical use requires design-specific validation of availability, sampling effort, observation intervals and spatial coverage. The prototype estimates sampled occupancy, persistence and positive abundance on the stated designs. It provides no independent evidence of dispersal or geographical range change. The most immediate statistical work is a comparison with replicated-region and replicated-survey designs, continued scrutiny of joint dependence, and targeted models for species-specific slopes trained and scored without oracle group labels.

## Figures

[Eight-page PDF](phase2/figures/phase2-diagnostics.pdf) · [SVG versions](phase2/figures/) · [Replicate-level and diagnostic CSVs](phase2/results/)
