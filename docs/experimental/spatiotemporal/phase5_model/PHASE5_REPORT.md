# Phase V results

The experimental implementation models community combinations on their observed spatial, local temporal and spatiotemporal support. The coherent model couples a bounded endpoint to signed ecological process predictions. These predictions add exactly to the fitted net change. An independent crossed model predicts occupancy-related dissimilarity and numerical heterogeneity jointly.

The Swedish lake results show that environmental predictors explain more of the **static spatial differentiation** than of the **changes in differentiation**. Held-out process prediction improved most for littoral redistribution. In simulations, a known expansion driver was often recovered, while contraction rankings were unstable even with a larger combination budget. These are predictive associations and compositional accounting results.

## Mathematical and computational validation

All 32 automated checks passed. They cover the three-block incidence versus four-block dummy identity, observed Shapley efficiency, valid denominators, actual irregular lags, native position/spread, bounded endpoint predictions, all crossed row/column sums, fixed-support telescoping, and separated spatial/forward validation folds. The maximum observed Shapley error across the simulations was below \(10^{-14}\). The complementary tests check closure-based relative distances, site-level gain and loss accounts, and use of full regional support even when only one focal site contributes a local response. The examples and validation tables are in `results/`.

The out-of-sample example removes 25% of a simulated irregular network. The full environmental model's forward endpoint RMSE was 0.0354, 0.0295 and 0.0235 in three successive test blocks. The nuisance-only model achieved 0.0302, 0.0267 and 0.0224. The intercept-only model achieved 0.0836, 0.0808 and 0.0725. This example separates algebraic coherence from predictive gain.

## Finnish moths

The full archive contributed 282 trap-years from 14 traps over 24 years. The exact counts were 1,551 within-year pairs, 253 one-calendar-year local comparisons and 1,293 matched spatial pair transitions. Consecutive observed visits gave 268 local comparisons, including 15 longer-lag visits. No suitable environmental covariates were archived for fitting. [BioTIME study 604](https://biotime.st-andrews.ac.uk/selectStudy.php?study=604)

Mean exact spatial Sørensen fell from 0.56447 in 1993 to 0.54137 in 2016, a change of −0.02310 on the all-observation networks. Summing the annual matched-panel changes gave +0.01749. The mean absolute annual discrepancy between the two estimands was 0.01799. On the previous fixed seven-trap panel, the corresponding cumulative change and direct difference were both −0.01759 for 2004–2013. The fixed panel and the 24-year variable network therefore answer different questions.

The annual matched-network cumulative expansion, contraction and redistribution accounts were −0.61763, +0.62531 and +0.00981 Sørensen units. These are sums over annual contrasts with changing matched networks. They expose opposing processes that are largely absent from the small net total. Average one-year local Sørensen was 0.26537 across contrast years and 0.24430 on the earlier seven-trap panel. Averaging at the comparison-row level gives 0.26559 for the full archive. Reorganisation remained substantial on the observed support.

## Swedish lakes

The environmental join retained complete pH and conductivity records for 2,842 of 2,844 phytoplankton visits, all 2,695 littoral visits and all 2,727 profundal visits. Climate keys matched every biological visit used. Chemistry and climate were joined by archive site/year. Temperature and precipitation describe the calendar year, so the forward-block analysis is conditional on the environmental observations being available. [Swedish lake archive](https://doi.org/10.5878/03t1-sv66)

| Assemblage | Observed visits | Exact spatial pairs | One-year local contrasts | Eligible matched pair transitions | Sampled matched pair transitions |
| --- | ---: | ---: | ---: | ---: | ---: |
| Phytoplankton | 2,844 | 138,443 | 2,686 | 128,800 | 6,889 |
| Littoral | 2,695 | 129,664 | 2,530 | 118,574 | 6,480 |
| Profundal | 2,727 | 127,296 | 2,534 | 114,075 | 7,170 |

Over a common 1995–2022 comparison, the exact all-observation spatial changes were −0.01686, +0.02257 and +0.03861 for phytoplankton, littoral and profundal. The Phase IV fixed panels gave −0.01858, −0.00762 and −0.01176. The all-observation first-year phytoplankton series extends back to 1992 and shows a larger rise (+0.18150 to 2022) on the changing network. These differences establish sensitivity to the declared time and sampling support. They have no unique interpretation as effects of network change.

For annual contrasts, the average absolute all-minus-matched discrepancy was 0.00902, 0.00667 and 0.01320 in assemblage order. Matched estimates use uniform combination draws, so the discrepancy includes combination-sampling variation. Local one-year Sørensen averaged 0.373, 0.376 and 0.381 over contrast years. Reclassifying processes with independent training/test geographic panels changed process labels by a mean absolute 0.017–0.021 in phytoplankton, despite identical held-out endpoint observations. The refitted spatial validation for both estimands is in `empirical_spatial_estimands.csv`.

### Blocked prediction and attribution

The table gives mean temporal-block RMSE for the annual spatial incidence change models. The coherent endpoint always lay in [0,1]. Its predictive advantage over the nuisance model was assemblage dependent. The unconstrained comparator had comparable errors and generated some incompatible endpoints, especially in the spatial holdout for profundal (5.0% of predictions by mean fold rate).

| Assemblage | Coherent, intercept | Coherent, nuisance | Coherent, environment | Direct ridge, environment | Joint process improvement over nuisance |
| --- | ---: | ---: | ---: | ---: | ---: |
| Phytoplankton | 0.1266 | 0.0792 | 0.0793 | 0.0792 | −0.24% |
| Littoral | 0.1186 | 0.0841 | 0.0811 | 0.0811 | +4.50% |
| Profundal | 0.1891 | 0.1629 | 0.1612 | 0.1607 | +1.76% |

Static spatial, local temporal and two-survey rectangular dissimilarities had temporal-block RMSE ranges of 0.0986–0.1811, 0.0990–0.1960 and 0.0802–0.1525. The crossed model's total derived from jointly fitted components had endpoint RMSEs of 0.1152, 0.1056 and 0.1785 in assemblage order. Figure 7 shows the contribution of spatial spread to predictive skill in static responses. Constraining static spread coefficients to be nonnegative slightly worsened temporal-block static RMSE in all three assemblages. Switching to equal-row training weights made little difference to the primary change model.

The environmental predictor blocks improved held-out littoral redistribution loss by 6.30%, of which the air temperature block contributed 0.0406 in predictive improvement units. The corresponding improvement for profundal redistribution was 2.58%, with air temperature contributing 0.0165. Phytoplankton redistribution showed only 0.78% overall improvement. Its pH contribution was positive (+0.0164) and other features offset part of it. The full environmental model worsened temporal-block contraction prediction in all three assemblages. The full environmental improvement for the joint response was negative for phytoplankton, so no relative importance percentages are reported. The held-out block ranges in Figure 6 describe sensitivity to the chosen temporal blocks, with only four blocks per assemblage.

## Controlled trajectory simulations

The main, null and independently seeded reference configurations each used 100 independent 60-site × 200-taxon × 30-survey panels. Eight smaller scenario suites used ten panels each at the same dimensions. Pure expansion, contraction, redistribution and numerical cases gave zero contributions from their inactive blocks. Exact mean order-two incidence zeta increased by 7.42 taxa under pure expansion and fell by 6.46 under pure contraction, on average over the simulated endpoints. It was unchanged under pure redistribution and pure numerical change. The latter nonetheless generated mean local relative Bray–Curtis of 0.06421. The stable redistribution scenario had mean annual net spatial change 0.00015 alongside mean local Sørensen 0.03379. Stationary expected spatial beta is the generator's target in this scenario. The realized mean spatial change need not be exactly zero.

Common-random-number interventions gave RMS changes in the realized expansion and contraction contributions of 0.00411 when E1 was disabled and 0.00356 when E2 was disabled, respectively. E2's RMS effect on realized redistribution was 0.01657. These intervention scales are in diversity units. Independently estimated predictive importance uses the deployed model class and a larger, 128-combination-per-period budget. In the 32-combination fits, E1 ranked first for expansion in 88 of 100 panels and E2 first for redistribution in 57. In the higher-budget reference the frequencies were 99 and 94, respectively. E2 ranked first for contraction in 32 panels in each budget, the same count as for E2 under the null. Joint top-rank recovery for E2 was 47% at budget 32 and 80% in the reference. This does not support a general ranking-recovery claim.

Mean main-config temporal endpoint RMSE was 0.01676 with environmental predictors, 0.01635 with nuisance predictors and 0.04292 with an intercept. Joint process squared loss was 0.0002733 versus 0.0002694 for nuisance alone. The null environmental model also deteriorated relative to nuisance alone (RMSE 0.01427 vs 0.01383). All four mean null joint predictive importance values were negative. A positive peak in a replicate or a first-place ranking by itself is weak evidence of an ecological signal.

The 25-panel incomplete-network experiment evaluated 15 combinations of missingness type and rate. Random removal of 10%, 25% and 50% of site-times preserved roughly 65.7%, 31.6% and 6.3% of complete-panel annual pair transitions. At 50% environmentally associated loss, the all-minus-matched RMSE rose to 0.01341, compared with 0.01028 under random loss. The environmentally selected matched changes averaged +0.00102 above the complete reference on comparable contrasts. Endpoint RMSE remained around 0.017 for many scenarios, even while combination support and the estimand changed. Figure 10 and the replicate results keep these outcomes separate.

## Limitations and deliverables

- The primary empirical fits use order two for spatial models and two-survey rectangles for spatiotemporal models. The empirical order screen reaches spatial order six and temporal order five. At that order 94.7% of profundal rectangles had zero shared taxa, versus 10.2% at pairwise spatial order. Only 5.8% of its paired windows supported complete shared-composition attribution, versus 89.0% at pairwise spatial order. Moths still shared taxa at the highest order, with only 17 sampled matched-window combinations. Higher orders can be limited by either saturation or coverage. The simulation order grid is smaller.
- Absolute abundance, population dynamics and detection were unavailable in the empirical archives used here. The abundance response and crossed attribution are relative-composition summaries. Moths were not environmentally enriched.
- The ecological interpretation is conditional on observed support, the allocation convention and the selected predictor set. Modest allocation draws are used for scalable fitting. The reusable response API exposes allocation Monte Carlo SE. Ecological inference from independent overlapping rows is unsupported, and the run does not establish ecological interval coverage or calibrate coefficient significance.
- Predictive coalitions inherit penalties chosen for their full outer training fold. Full, nuisance and intercept scores share their test rows. Inference from a small number of temporal/spatial blocks is limited. Partial curves hold other design columns fixed and are useful as conditional evaluations of the fitted functions. Their points need not form an observed multivariate environmental trajectory.
- Reclassification within regions changes the attribution target even when the endpoint stays the same. The numerical results separate model evaluation with realized regional labels and reclassified train/test regions.
- Figures 1–10 are in `figures/`, with vector PDF and PNG exports. CSVs in `results/` include support diagnostics, blocked results, importance, curves, response summaries, and replicate-level simulations. The methods and reproduction sequence are in [PHASE5_METHODS.md](PHASE5_METHODS.md).
