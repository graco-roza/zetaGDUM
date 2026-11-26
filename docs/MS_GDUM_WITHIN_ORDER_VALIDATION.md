# Within-order predictor importance and calibration

The fitted model is additive on the logit scale. Directional I-spline coefficients are nonnegative and environmental position can be nonlinear. Each zeta order is an independent analysis. The scientific quantities are explained-deviance contributions and paired differences between predictors within an order.

## Estimand and uncertainty

For a coalition $C$, refit the fractional-logit loss with the declared training weights and the corresponding restriction of the penalty. Evaluate the resulting prediction on one fixed set of assessment rows. Use the intercept-only training fit to define the common deviance denominator for all coalitions in that draw. Absolute Shapley importance is the mean marginal gain across predictor orderings. A predictor consists of its spread and position terms, and geographic distance forms its own group. Within each predictor, Owen values average the order of the spread and position terms conditional on the previously entered predictors. They sum to that predictor's absolute contribution. Held-out gains can be negative. Relative shares are suppressed unless the full held-out or in-sample explained deviance exceeds the declared positive threshold. [Derivation and exact-versus-approximate benchmark](within_order/COALITIONS.md)

Small games use exact subset weights. Larger games reuse sampled predictor orderings and cached coalition fits. The reported permutation standard error is conditional on the data, so it should be checked separately from ecological interval width. The exact evaluator shares the training standardisation, basis and penalty coordinates with every coalition. Failed fits are retried cold after a warm start fails. Unresolved failures are surfaced.

The direct Bayesian bootstrap draws Dirichlet weights over sites and samples a set with probability proportional to the product of its site weights. It builds one set pool per draw and refits every coalition using that pool. The ordinary site bootstrap draws site copies with replacement, then samples sets of distinct copies. Spatial sensitivity uses common weights or score multipliers within declared blocks.

The approximate methods fit the coalitions separately at the original data. Let $J_i(C)$ be the centred pseudovalue of the coalition estimating equation. The approximate score perturbation is a weighted sum of $J_i(C)/n$. Each draw reuses the same multipliers for every coalition. At a regular interior optimum the Gaussian update is the flat-prior Gaussian GMM approximation from the square estimating equation, with covariance $H^{-1}\widehat\Omega H^{-T}/n$. In the multiplier version the multipliers are Rademacher variables. At active constraints the one-step update solves a projected quadratic problem including the original fitted gradient in the Hessian metric. Projection applies to the coalition update. It is an approximation to constrained refitting. The refit comparator adds the same random linear term to the nonlinear fractional-logit loss. It measures local approximation error for an identical perturbation.

For a draw with allocations $A_v$, ranking uses the same vector for all predictors. Ties divide leading probability equally. Pairwise probabilities and central 80% and 95% intervals use matched differences $A_v-A_w$. We evaluate practical superiority at an absolute explained-deviance difference of .01. Importance under this model is a predictive allocation. In the correlated-predictor scenario, a predictor omitted from the ecological generating mechanism can receive positive population importance.

## Population reference and simulation assumptions

Niche simulations generate Poisson abundance at independent sites given three continuous environmental predictors. Species optima span the first two predictor axes, with narrower first-axis niches in the baseline. One background species keeps composition defined. Environmental effects on incidence and relative abundance arise from species niches and abundance. The fit uses three basis terms for directional spread and nonlinear local position, with fixed small quadratic penalties of $10^{-4}$. This defines the regression and attribution target. The generating mechanism supplies no Shapley coefficients directly.

The common background taxon guarantees some sharing in the simulated sets. The sparse scenario reduces shared composition beyond that taxon and does not test a zero-sharing threshold.

Three independent 3,600-site reference communities provide large training and assessment pools of 85,000 sets. Cyclic independent assessment across those communities gives three numerical reference estimates. Ecological replicates use 70 or 150 sites and independent assessment sets. The archived reference spread between the strong and weak predictors is small for relative abundance at several orders. We do not interpret a prediction model's high full deviance as evidence that it can resolve such gaps.

The order diagnostic records response variance, shared composition, scaled information, fixed-site integration variability, ecological interval widths and independent-site predictive deviance. It takes user-declared precision and validity requirements. Each order is checked on its own. An order with good held-out performance can still have unstable predictor attribution. Empirical descriptive results are archived separately. [Five-community empirical analysis](within_order/EMPIRICAL.md)

## Main calibration results

Each principal cell used 200 independent ecological communities, 700 training combinations per pool, 79 uncertainty draws per community and 4,500 fixed independent assessment combinations. The score methods used two extra independent training pools. Coverage is conditional on the fixed assessment pool and on obtaining at least 80% finite draws. It treats the numerical population reference as fixed. The usable count exposes the selection. A matched 59-draw run of the score methods used the same ecological seeds and is a Monte Carlo sensitivity check. It does not add independent communities. Reported times are observed elapsed seconds per community and vary with concurrent workloads. All replicate-level distributions of the interval summaries, including failures, are archived.

| Sites | Order | Response | Method | Usable / 200 | Strong − weak 80% coverage | 95% coverage | Mean 95% width | Correct leader |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 70 | 2 | Incidence | Direct Bayesian bootstrap | 200 | .985 | 1.000 | .065 | 1.000 |
| 70 | 2 | Incidence | Gaussian GMM | 199 | .910 | 1.000 | .080 | 1.000 |
| 70 | 2 | Incidence | Multiplier | 199 | .925 | .990 | .080 | 1.000 |
| 150 | 4 | Incidence | Direct Bayesian bootstrap | 200 | .860 | .995 | .063 | 1.000 |
| 150 | 4 | Incidence | Gaussian GMM | 200 | .965 | 1.000 | .113 | 1.000 |
| 150 | 4 | Incidence | Multiplier | 200 | .970 | 1.000 | .113 | 1.000 |
| 70 | 8 | Relative abundance | Direct Bayesian bootstrap | 200 | .925 | 1.000 | .173 | .277 |
| 70 | 8 | Relative abundance | Gaussian GMM | 200 | .980 | 1.000 | .284 | .270 |
| 70 | 8 | Relative abundance | Multiplier | 200 | .995 | 1.000 | .281 | .240 |
| 150 | 12 | Relative abundance | Direct Bayesian bootstrap | 200 | .965 | .995 | .110 | .693 |
| 150 | 12 | Relative abundance | Gaussian GMM | 200 | .995 | 1.000 | .205 | .640 |
| 150 | 12 | Relative abundance | Multiplier | 200 | 1.000 | 1.000 | .212 | .708 |

At order 2 all three methods identified the population leader, with a strong-minus-weak reference gap of .715. The 95% intervals for the irrelevant predictor covered its reference in just .105, .090 and .065 of replicates for the direct bootstrap, Gaussian GMM and multipliers, respectively. The score methods each lost one replicate to a coalition fit that failed both warm and cold starts. The direct bootstrap took 19.3 seconds per replicate, compared with about 3.5–3.6 seconds for the score methods plus .068 seconds for their extra pools under the same run.

At order 4 the population strong-minus-weak gap was .483. The absolute allocation for the irrelevant predictor was approximately zero in the population reference. Its 95% interval covered the reference in only 24% of direct bootstrap and 15.5–16.5% of score-method replicates. The original point estimator had mean bias −.006 for this predictor. The corresponding mean biases of the uncertainty medians were −.006 and approximately −.014. The score procedures showed stronger negative held-out allocation bias, even where their paired contrasts covered. Finite-sample optimism of the fitted training objective and held-out loss of irrelevant variables make this allocation a difficult near-null target. Separate score covariance by order left this bias. The earlier joint-order calibration considered different estimands and data-generating mechanisms. The present results therefore do not establish a matched coverage advantage over joint-order methods. At least one monotonicity bound was active for the irrelevant predictor in 98% of the order 4 fits. When all three were active, nominal 95% coverage of its individual allocation was .312 for the direct bootstrap and .273 and .195 for the score methods. The coverage problem persists within an observed active-face subgroup. [Active-constraint diagnostic](within_order/results/boundary_diagnostic.csv)

Bias of the strong-minus-weak contrast was −.012 for the direct bootstrap and −.006 to −.007 for the score procedures, and RMSE was .013–.015. The direct bootstrap took about 7.5 seconds per replicate for this cell and the two score methods took about 1.6 seconds each before the separate pool preparation.

At order 8 the relative-abundance gap was .0021, well below the practical difference of .01. The direct bootstrap assigned mean leading probability .583 to its top-ranked predictor but that predictor led in only .277 of ecological replicates. For the score procedures those figures were about .52 and .24–.27. The three reference gaps at this order ranged from .0015 to .0029 with the same leader. These are observed ecological reliability summaries for a single generating process, with a common target per scenario. They do not establish calibration of posterior probabilities across processes. At order 12 the reference gap reversed to −.0048 and even the narrowest nominal 95% contrast interval averaged .110. High coverage of these contrasts coincided with poor resolution. Success percentages also depend on the numerical precision of the reference when the gap is very small. [Per-estimand calibration](within_order/results/calibration_summary.csv) · [Replicate-level rank assignments](within_order/results/calibration_rank.csv) · [Coverage figure](within_order/results/calibration_coverage.pdf)

A small conventional full-site-resampling comparison used 12 communities with 70 sites at order 2. For the incidence strong-minus-weak difference its 80% and 95% coverage were .917 and 1.000, mean 95% width .086 and median bias .010. Sites sampled more than once were treated as distinct copies when forming sets. This is a small reference for the direct procedure, with approximately 23 elapsed seconds per community under concurrent workloads.

### Stress cases at order 8

Each stress case used 24 independent communities with 70 sites, incidence and 39 score draws. In the symmetric equal-importance case, neither score procedure falsely declared superiority at a probability of .95, and both covered the null in 24 of 24 replicates at both interval levels. Their 95% widths averaged .445 and .488. Neither procedure declared a winner at .95 in the correlated, broader niche or sparse incidence cases, whose reference differences were between .016 and .027. This limits claims about power. The correlated reference allows positive attribution to a predictor omitted from the biological data mechanism. The spatial case used an AR(1) gradient and one weight per contiguous block. Compared with exchangeable scores on the same sites, it produced 95% widths .642 versus .655 for Gaussian draws and .686 versus .544 for multiplier draws, with full observed 95% coverage in both designs. These small tests do not establish that the selected block grouping makes spatial inference valid. For nonlinear environmental position, the reference difference was .157. Gaussian and multiplier methods had 80% coverage of .708, 95% coverage of .958 and .917, and median bias of −.066 and −.062 respectively.

In a matched 24-community order 8 relative-abundance sensitivity analysis, increasing the fit and two score pools from 700 to 2,800 sets narrowed the average 95% strong-minus-weak width from .272 to .177 for Gaussian GMM and from .301 to .189 for multipliers. RMSE fell from .036 to .025 and from .046 to .020. The pool-generation cost rose from .058 to .215 seconds on average. The contrast remained too small to resolve. This supports an incomplete-integration contribution to the original width, with a remaining ecological, curvature or boundary component that this contrast cannot separate.

## Approximation and numerical integration

For matched score draws at 70 sites and 1,000 combinations per pool, replacing exact perturbed coalition refits with one-step updates introduced paired importance-difference RMSE between .0025 and .0132 across order 4 and 12 incidence and relative abundance. Hybrid updates reduced that range to .00017–.00204. In order 12 incidence, for example, total elapsed time including shared preparation was 5.01 seconds for refits, 3.36 for hybrid and 2.19 for one-step, over 19 draws. Exact refits returned 17 complete draws in that cell. The errors are calculated on matched finite draws. This comparison isolates the coalition update, with Rademacher score perturbations held constant. It cannot establish that the perturbed target is ecologically calibrated. [Resampling derivation](within_order/RESAMPLING.md) · [Larger update benchmark](within_order/results/benchmark_updates.csv) · [Focused final benchmark](within_order/results/benchmark_updates_final.csv)

The independent fixed-site integration study varied the pool size over 700, 2,800 and 11,200 combinations, with 24 replicates in each order-response cell. The first 12 selected the smallest budget for which all three predictor allocations and the strong-minus-weak comparison had estimated numerical SE below .01. The second 12 assessed that threshold independently. Both order 2 responses became exact at 2,800 because all 2,415 pairs could be enumerated. At order 4 incidence, the 2,800-set pilot estimated .0056 but the second half gave .0117. At order 8 incidence only 11,200 qualified, with confirmation SE .0083. At order 8 relative abundance, 2,800 qualified with confirmation SE .0078, whereas 11,200 gave .0043. This fixed tolerance is illustrative. When the scientific contrast of interest is .01, a numerical SE of .008 still obscures the contrast. Recheck budgets on the actual sites, response and estimand using independent pools. [Budget verification](within_order/results/integration_budget.csv)

Exact local zeta control variates gave a reduction in median scaled score and covariance errors from .210 and .044 to .133 and .029 at order 8 for relative abundance. The corresponding elapsed CPU times were .070 and .068 seconds. For order 8 incidence the error rose from .205 and .047 to .226 and .069 at similar cost. Two independent 70,000-combination order 8 reference pools still differed by .090 in scaled score RMS. The order 2 benchmark had exhaustive references and yielded mixed results. This supports a conditional, independently piloted control policy. [Control-variate benchmark](within_order/results/control_benchmark.csv)

## Empirical implications

The five communities have descriptive point estimates and spatial held-out prediction checks for both responses. Finnish birds and Forest Hemiptera had the same leading predictor in both integration pools at 12/12 and 11/12 order-response cells, respectively. The agreement was 9/10 for Pond diatoms, 3/10 for Arctic plants and 9/10 for Atlantic marine fish. Arctic plant predictions had positive pooled spatial explained deviance at incidence orders 4, 8 and 12 and relative-abundance order 12. Atlantic marine fish relative-abundance response SD reached .0008 at high order. A relative attribution in that regime can represent very little variation in dissimilarity. The empirical diagnostic curves give the predictive and numerical ranges separately. No ecological intervals or ranking probabilities are presented for the observed communities because the calibration benchmark did not justify interpreting them. [Community figures and data](within_order/EMPIRICAL.md)

The supported use in this phase is comparative description on independent assessment sets, accompanied by combination-pool sensitivity and response magnitude. In these simulations the direct Bayesian bootstrap is the narrower uncertainty comparator and distinguished the very large incidence gaps at orders 2 and 4, but its near-null individual allocation intervals were poorly calibrated. The score procedures are faster and permit internally paired draws, while projection, incomplete score integration, and weak predictor separation can still make their results unsuitable for inference. Retain the existing bootstrap defaults. Treat positive spatial-transfer explained deviance and stable two-pool rankings as distinct diagnostics. Reserve inferential rankings for a design with ecological validation at the relevant order, adequate contrast precision and reliable coalition convergence. These results do not support a universal informative-order cutoff.

## Reproduction

Install the package from this branch and run the scripts at the repository root. Each ecological replicate has a deterministic seed. Set `MS_GDUM_CALIBRATION_WORKERS` to the available fork worker count on Unix. The scripts use one worker by default. Elapsed times depend on machine and concurrent work.

```sh
for k in 2 4 8 12; do Rscript docs/within_order/reference.R baseline "$k"; done
for scenario in equal correlated position sparse spatial broad; do
  Rscript docs/within_order/reference.R "$scenario" 8
done

Rscript docs/within_order/calibration.R 70 2 incidence baseline 200 79 700
Rscript docs/within_order/calibration.R 150 4 incidence baseline 200 79 700
Rscript docs/within_order/calibration.R 70 8 relative baseline 200 79 700
Rscript docs/within_order/calibration.R 150 12 relative baseline 200 79 700
Rscript docs/within_order/calibration.R 70 2 incidence baseline 12 49 700 site_bootstrap
for scenario in equal correlated position sparse spatial spatial-naive broad; do
  Rscript docs/within_order/calibration.R 70 8 incidence "$scenario" 24 39 700 gaussian_gmm,multiplier
done
for budget in 700 2800; do
  Rscript docs/within_order/calibration.R 70 8 relative baseline 24 39 "$budget" gaussian_gmm,multiplier
done

Rscript docs/within_order/boundary_diagnostic.R
Rscript docs/within_order/benchmark_updates.R
Rscript docs/within_order/benchmark_updates.R --extended
Rscript docs/within_order/coalition_benchmark.R
Rscript docs/within_order/control_benchmark.R
Rscript docs/within_order/integration_study.R
Rscript docs/within_order/summarise_calibration.R
Rscript docs/within_order/plot_calibration.R
R CMD build .
R CMD check zetaGDUM_0.1.0.tar.gz
```

For the score-only draw-count sensitivity, repeat the three principal commands at orders 4, 8 and 12 with 59 draws and `gaussian_gmm,multiplier` as the final argument. The additional baseline pilot settings `(sites, order, response)` were (70, 4, relative), (70, 12, incidence), (150, 2, relative) and (150, 8, incidence), each with 12 communities, 49 draws and 700 combinations. The first principal runs at orders 4, 8 and 12 used an early score setup before cold retry was added for failed coalition warm starts. For the archived runs, the separate retry files replay those failed replicates with the final optimizer. The summary replaces the corresponding original rows without duplicating ecological communities. To rerun this audit on the archived data, repeat the three 79-draw principal commands at orders 4, 8 and 12 with the final arguments `gaussian_gmm,multiplier --retry-failed`. A fresh run already uses cold retries. The cold-retry audit recovered all 52 previously failed ecological replicates at orders 4, 8 and 12 for both score methods. It addresses coalition startup failures. It did not resolve statistical undercoverage or weak ranking recovery. The 59-draw historical comparisons retain their original setup results. The final script records missing coalition fits without terminating an entire run. [Binned ranking reliability](within_order/results/calibration_rank_reliability.csv)

## Package verification

The source tarball passed full `R CMD check` with status OK on arm64 macOS under R 4.6.1, including examples, vignettes and manual PDF generation. The check recorded 1,623 passing test expectations and one pre-existing test skipped by its CRAN guard. A separate local run with the CRAN guard disabled passed all 47 expectations in the existing core GDUM test file, including its comparison with gdmmTMB. The within-order tests cover exact versus sampled allocations, component sums, coordinate validation, unresolved fits, matched draws, constraints and weighted precision diagnostics.
