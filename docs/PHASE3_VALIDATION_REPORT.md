# Phase 3 validation report

Exact conditional zeta and the environmental response moments passed exhaustive numerical checks. The new local summaries describe which sites contribute to multisite dissimilarity and how that description changes with order. The guarded hybrid, balanced designs and overlap-aware regression intervals remain experimental. The existing fitting and uncertainty defaults remain in place.

## Reproducibility and estimands

The scripts and replicate-level CSVs are in [docs/phase3](phase3). From the package root, install the package with `R CMD INSTALL .` in the active R library and run:

```sh
Rscript docs/phase3/fitting.R
Rscript docs/phase3/fixed_kernel.R
Rscript docs/phase3/calibration.R null
Rscript docs/phase3/calibration.R weak
Rscript docs/phase3/calibration.R strong
Rscript docs/phase3/calibration.R null 24 10000 72,241
Rscript docs/phase3/calibration.R strong 24 10000 72,241
Rscript docs/phase3/local.R /path/to/DD2/Data
Rscript docs/phase3/summarise.R
```

The empirical run requires the same DD2 tables and preprocessing as the manuscript-scale validation. Omit the directory argument to run the simulated local comparisons and computational benchmark independently. The repository contains only aggregate empirical summaries. The Rcpp helpers for calibration are the previously checked conditional Poisson sampler and interior logit solver in [calibration_helpers.cpp](manuscript_validation/calibration_helpers.cpp). All other fitting benchmarks use the package fitter. Elapsed times are indicative single-machine measurements on macOS arm64 with R 4.6.1. Benchmarks report R object sizes in KiB. Peak memory and machine-independent time were outside this measurement.

The fitting benchmarks target a uniform distribution of distinct-site sets in a fixed community. Relative abundance uses unit working weights. Incidence uses the legacy mean-richness working weight on richness-normalised total dissimilarity. Both use the same working weights and predictor bases for each method and its reference. Constant-weight normalised incidence remains sampled. The uncertainty benchmarks instead compare intervals to known generating effects across independent site samples. Fit error against a complete finite-community objective and interval coverage of an ecological generating parameter are different quantities.

## Exactness and computational scaling

The tests enumerate all orders from one through six on binary and sparse abundance examples, with ties, ubiquitous species, narrow occupancy, identical communities and communities with no shared species. They compare unconditional, single-site, distinct-pair and arbitrary forced-subset expectations to direct enumeration at a tolerance of $10^{-10}$. Relative-abundance inputs with empty sites are rejected. Zero dissimilarity produces undefined normalised LCZD. A separate test of sampled normalised incidence checks that set-level ratios are averaged directly. Averaging the exact conditional site, pair and forced-subset results recovers the aggregate target under the corresponding uniform conditioning design.

For a 12-site complete regression design, environmental cross-products agreed at $10^{-10}$ for both weighting targets. A geographic design verified that its nonlinear columns stay sampled. On those complete designs the experimental fits agreed with the conventional deviance at tolerance $10^{-7}$ and fitted probabilities within $10^{-4}$. Replicating sets unevenly and supplying the exact inverse proposal ratios recovered the exhaustive weighted fit. Deliberately incompatible moments triggered a diagnostic fallback. Tests also verify the response-scale local Shapley sum, disjoint groups, exact complete affine pair counts and the signed covariance algebra.

Computational measurements use 350 species, 12 calls per cell, and $k=\lfloor .08n\rfloor$. The figures below are medians for abundance at approximate density .85. At density .03 the pair calculation took 1 ms at $n=241$ and 3 ms at $n=600$, where most species have occupancy below the order. The dense input was 659 KiB and 1,641 KiB. At $n=600$ and density .032, a compressed numeric column payload would be about 81 KiB using eight bytes per value, four per row index and four per column pointer. The current code traverses occupied sites for accumulation and materialises the community input as dense. The pair matrix remains dense regardless of that conversion.

| Sites | Site-only time (ms) | Pair-inclusive time (ms) | Site-only result (KiB) | Pair-inclusive result (KiB) |
|---:|---:|---:|---:|---:|
| 241 | 4 | 16 | 2.5 | 456.7 |
| 600 | 11 | 99 | 5.3 | 2,818.2 |

## Local contributions in the five manuscript communities

The standard LCBD comparator is the site share of total squared distance from the centroid of Hellinger-transformed compositions. Pairwise GDUM position is the site-standardised environmental vector multiplied by the position coefficients from the order-two relative-abundance model. The corresponding site-level predicted dissimilarity averages both the spread and position terms after the logistic transformation. At order two the exact abundance LCZD input $U_i$ matched the observed mean pairwise dissimilarity involving $i$, with maximum absolute discrepancy $1.0\times10^{-15}$ across communities.

| Community | Sites | Highest order | SD of $U_i$, order 2 → high | High-order rank correlation with mean pairs | With Hellinger LCBD | With pairwise GDUM position | Raw incidence $U_i$ with richness |
|---|---:|---:|---:|---:|---:|---:|---:|
| Finnish birds | 241 | 28 | .112 → .0038 | .869 | .864 | .692 | .974 |
| Pond diatoms | 102 | 15 | .069 → .0009 | .833 | .871 | .200 | .977 |
| Forest Hemiptera | 147 | 11 | .094 → .0039 | .918 | .884 | .452 | .982 |
| Arctic plants | 72 | 7 | .060 → .0062 | .749 | .790 | .285 | .874 |
| Atlantic marine fish | 202 | 14 | .050 → .000056 | .761 | .771 | .459 | .999 |

At order two, relative-abundance ranks correlated .943–.985 with Hellinger LCBD. At the highest orders the correlations were .771–.884. Absolute relative-abundance differences between sites compressed with order, so relative ranks can remain different even when the site means are nearly equal. Relative-abundance distinctiveness correlated −.513 to .262 with site richness at the highest orders and .302–.622 with an abundance-weighted inverse occupancy index. Raw incidence includes a mean-richness term, which explains why richness can dominate when shared composition has largely disappeared.

In the 72-site simulation at order seven, multiplying all counts at a site by its own effort factor preserved proportional composition and every $U_i$ to numerical tolerance. Uneven individual detection changed the baseline rank correlation to .744 and moved individual $U_i$ by up to .105. Adding variable numbers of singleton individuals changed it to .998 with a maximum movement of .0045. Multiplying taxa by distinct abundance factors gave a rank correlation of .962 and a maximum movement of .051. These are sensitivity experiments on specific transformations. They do not identify sampling effort or establish intrinsic uniqueness in the field datasets.

## Estimating the finite-community GDUM objective

The exhaustive study used one simulated 13-site community with two environmental predictors at orders three and five. All 286 and 1,287 distinct sets respectively define the reference objectives. Each budget used 24 independently sampled designs per order. IID, species mixture ($\alpha=.5$), disjoint site blocks, affine pair cycles, and environmental quartile strata used the conventional fitter. Hybrid variants used the same IID draws. Parameter errors use the total spread height and the position coefficient of the first environmental predictor. Importance and deviance Shapley values are compared on a common reference design, with each variable and its spread/position terms defining one deviance coalition. Local Shapley dissimilarity is a separate response-scale comparison. Convergence, sharing, omission support, rank agreement, prediction RMSE and loss differences are recorded in the CSVs.

| Response and method, 600 sets | Fallback | Spread RMSE | Position RMSE | Importance RMSE | Deviance Shapley RMSE | Prediction RMSE |
|---|---:|---:|---:|---:|---:|---:|
| Incidence, IID | — | .138 | .0107 | .0437 | .0045 | .0073 |
| Incidence, site-balanced | — | .112 | .0131 | .0410 | .0045 | .0077 |
| Incidence, hybrid | 6% | .991 | .1040 | .2715 | .2330 | .0538 |
| Incidence, control variate | 0% | .152 | .0168 | .0401 | .0058 | .0087 |
| Relative, IID | — | .136 | .0147 | .0173 | .0031 | .0057 |
| Relative, site-balanced | — | .126 | .0164 | .0181 | .0029 | .0059 |
| Relative, hybrid | 60% | .748 | .0927 | .1273 | .1368 | .0357 |
| Relative, control variate | 0% | .234 | .0328 | .0393 | .0061 | .0116 |

At 120 sets, the hybrid fell back in 62.5% of weighted-incidence fits and 97.9% of abundance fits. A strict certificate is a feasibility and optimisation safeguard. It admits some objectives that still give large finite-budget errors. At 600 sets, the control variate passed every full-fit certificate but the relative-abundance position RMSE was 2.2 times the IID error, and the median optimisation pipeline took 11 versus 4 ms apart from precomputing the moments. The paired relative-abundance spread change between orders three and five also had greater error with the control variate. Each small-community set shared a species in this simulation, so the high-order simulations supply the more informative comparison for targeted sampling.

### High-order simulations at manuscript dimensions

Five other simulated communities match the manuscript site and species counts before empty-site handling, and the highest orders $(n,k)=(72,7),(102,15),(147,11),(202,14),(241,28)$. Niche breadth, gradient smoothness and abundance vary between cells. Sites with zero simulated counts receive a site-specific singleton to permit the relative-abundance model while leaving sharing at $k\ge2$ unchanged. Ten independent designs of 1,500 sets per cell use a common IID reference of 30,000 sets. The reference retains Monte Carlo error and uses the same community, so reference prediction error is a numerical convergence measure. Pair-balanced high-order results only apply to the prime site count $n=241$.

| Response and method | Fallback | Spread RMSE | Position RMSE | Importance RMSE | Deviance Shapley RMSE | Pipeline median (ms) |
|---|---:|---:|---:|---:|---:|---:|
| Incidence, IID | — | 1.338 | .0335 | .0486 | .0032 | 11.0 |
| Incidence, species mixture | — | .847 | .0252 | .0396 | .0023 | 37.5 |
| Incidence, site-balanced | — | 1.172 | .0283 | .0506 | .0028 | 8.0 |
| Incidence, environment-balanced | — | .899 | .0354 | .0510 | .0030 | 27.0 |
| Incidence, control variate | 68% | 1.339 | .0332 | .0481 | .0017 | 18.0 |
| Relative, IID | — | 1.344 | .0394 | .0456 | .0038 | 16.0 |
| Relative, species mixture | — | .937 | .0278 | .0398 | .0024 | 43.0 |
| Relative, site-balanced | — | 1.176 | .0394 | .0501 | .0039 | 13.0 |
| Relative, environment-balanced | — | 1.142 | .0382 | .0505 | .0043 | 34.5 |
| Relative, control variate | 60% | 1.337 | .0393 | .0433 | .0026 | 28.0 |

All returned fits converged and had numerical Hessian rank 9 of 9. Directional constraints were frequently active, averaging 2.49 of six spread coefficients in IID fits. The exact-term hybrid fell back on every full fit at this budget, giving the IID point estimates, while some reduced weighted-incidence coalitions passed and changed the deviance Shapley estimate. This is why hybrid full-model fallback and its Shapley error can both be high. Timings include set generation, design, response and model fitting. They exclude common-reference scoring and Shapley refits. Cached cross-products took a median of 7 ms to compute and occupied roughly 2 KiB. Retained training designs were 147–281 KiB.

The species mixture increased the mean count of shared sets from 794 to 1,148 while the weighted sharing estimate remained near .53. The strongest increase was at $n=72,k=7$, from 116 to 805 sets with sharing, with a weighted probability of .077 under either sampler. The additional information on sharing came with lower minimum exclusion Kish ESS, from 1,328 to 710 in that cell. At a fixed count of 1,500 sets the species mixture reduced pooled spread and position RMSE for both targets. It was slower than IID by a factor of approximately 2.7–3.4 in this timing run. These data establish a fixed-budget gain in the tested communities, while equal-time superiority remains open. Weighted importance rank agreement was identical to IID in these pooled simulations, .86 for incidence and .92 for abundance.

Site blocks reduced the across-site inclusion count CV from .079 to .018 and increased the mean minimum exclusion ESS from 1,318 to 1,342. The site-block result was mixed across effect and importance metrics. Environmental balance reduced inclusion CV to .029, with similar omission support and slower construction. At $n=241,k=28$, 1,500 rows cover only a small fraction of the 57,840 affine cycle. The sampled pair-count CV was .225 under affine sampling versus .220 under IID. Exact balance was demonstrated by enumeration of a complete prime cycle in the tests. Truncated affine cycles offer no consistent advantage here. ESS uses the proposal ratios before response working weights and measures set support under omission, with ecological replication still determined by the sites. Local response-scale Shapley RMSE, per-cell losses, biases and the paired-order comparisons can be reproduced from the supplied CSVs.

## Overlap-aware uncertainty

### Fixed-kernel validation

An additive kernel $g(S)=k^{-1}\sum_{i\in S}U_i$ used iid standard normal site values and 250 ecological replicates per manuscript dimension. Its complete mean equals the site mean exactly and has variance $1/n$. Each matched draw used about 800 sets. The overlap-resolved estimator used 100 independent union draws at each overlap less than $k$.

| Sites, order | Mean matched variance / truth | Mean overlap variance / truth | Matched SD / truth | Overlap SD / truth | Negative overlap variance |
|---|---:|---:|---:|---:|---:|
| 72, 7 | 1.026 | .978 | .178 | .916 | 15.6% |
| 102, 15 | .981 | .999 | .153 | .935 | 12.8% |
| 147, 11 | .981 | 1.092 | .132 | 1.240 | 17.6% |
| 202, 14 | 1.015 | 1.160 | .106 | 1.668 | 23.6% |
| 241, 28 | 1.009 | .981 | .111 | 1.359 | 23.2% |

The matched means were within 3% of the expected complete and incomplete variances in these runs. The overlap-resolved estimator was much noisier. Its largest mean discrepancy was about 16% of the theoretical variance at $n=202$, where its replicate SD was about 167% of the truth. Negative estimates were retained.

### Fitted-score calibration

The simulation generated a fractional logit response from a prespecified one-variable I-spline design with site-shared, bounded, mean-zero noise independent of the predictors. The unconstrained intercept and spread and position coefficients had known population values. Spread height and position equalled $(0,0)$, $(.20,.007)$ or $(1,.035)$. The same 48 ecological draws per dimension were reused across signal levels, which makes the sensitivity to signal paired. Each point estimate used about 1,600 matched sets. An external deterministic spline reference fixed the predictor scale. Bayesian intervals used the 2.5th and 97.5th percentiles of 79 independent conditional Poisson bootstrap fits on the original sites. Site deletion and the Gaussian deletion-score multiplier used a Student critical value with $n-1$ degrees of freedom. Matched and overlap-specific sandwiches used a normal critical value. The latter used 180 unions per overlap with the incomplete matched-design correction. No method changed the point estimate.

The table pools 240 ecological replicates per signal. Within each cell, conditional coverage and rejection use available intervals with positive contrast variance. For the overlap-aware methods a positive contrast variance does not imply that the full fitted covariance matrix is positive semidefinite. The final column gives strong-signal rejection among available intervals, so it should be read alongside the availability column. The individual site-count cells and their biases and diagnostics are in the reproducible summaries.

| Method | Effect | Null interval availability | Null false-positive rate | Weak coverage / power | Strong coverage / power | Median time at null |
|---|---|---:|---:|---:|---:|---:|
| Conditional Poisson Bayesian bootstrap | Spread | 100% | 18.3% | 81.7% / 19.2% | 81.3% / 42.9% | 204 ms |
| Conditional Poisson Bayesian bootstrap | Position | 100% | 10.4% | 89.6% / 15.8% | 87.9% / 57.1% | 204 ms |
| Delete-one-site jackknife | Spread | 100% | 0.4% | 99.6% / 0.4% | 99.6% / 2.9% | 33 ms |
| Delete-one-site jackknife | Position | 100% | 1.3% | 98.8% / 1.3% | 98.8% / 10.0% | 33 ms |
| Deletion-score multiplier | Spread | 100% | 3.3% | 96.7% / 3.3% | 96.7% / 9.2% | 1 ms |
| Deletion-score multiplier | Position | 100% | 2.5% | 97.5% / 1.7% | 97.5% / 15.8% | 1 ms |
| Matched-score sandwich | Spread | 33.3% | 17.5% | 82.5% / 13.8% | 81.5% / 30.9% | 3 ms |
| Matched-score sandwich | Position | 43.3% | 10.6% | 90.3% / 8.7% | 90.4% / 37.5% | 3 ms |
| Overlap-resolved sandwich | Spread | 47.1% | 2.7% | 97.3% / 2.7% | 99.1% / 1.8% | 17 ms |
| Overlap-resolved sandwich | Position | 43.8% | 2.9% | 97.1% / 2.9% | 98.1% / 5.8% | 17 ms |

Pooled biases in the point estimates were −.011, −.011, −.007 for spread and +.002, +.002, +.002 for position for the null, weak and strong settings. Their across-replicate SDs were 1.230–1.342 for spread and .031–.034 for position. Near-collinear spread basis columns and scarce independent information make power low even when point bias is small. In the fitted-score experiment the matched covariance was indefinite in every replicate. The overlap covariance was positive semidefinite in one of 240 replicates at each signal level. Conditional coverage after excluding nonpositive contrast variances therefore provides little evidence for using these sandwiches in practice.

A sensitivity run at exactly 10,000 sets used the first 24 independent ecological replicates per signal at each of $(72,7)$ and $(241,28)$. The full covariance matrices of both overlap-aware methods were indefinite in every run. At the null, 9 of 48 matched spread variances and 16 of 48 overlap-resolved spread variances were positive. Null spread rejection was 20.8% for the Bayesian bootstrap, 10.4% for the multiplier, and 2.1% for the jackknife. Strong spread coverage and power were 77.1% / 52.1%, 87.5% / 27.1%, and 97.9% / 14.6%, respectively. The multiplier's calibration degraded for spread in this larger-budget sensitivity despite its better pooled calibration at 1,600 sets. Only two dimensions and 48 draws per signal contribute to the sensitivity table. The median calculation times were 1.53 s for the Bayesian bootstrap, 0.25 s for site deletion, and 0.004 s for the multiplier at the null. The larger-budget Hessian's median condition number was about 22,000 at $n=72$ and 392,000 at $n=241$.

## Implications and unresolved issues

- Use the exact site statistics for relative abundance and raw incidence when those response definitions match the ecological question. Compare absolute distinctiveness as well as rank. Retain sampled set ratios for constant-weight incidence.
- Keep conventional fitting as the default. The exact cross-products are suitable as diagnostics and experimental inputs. Exact response terms by themselves produced poor finite-budget estimates. The intercept-centred control variate improved moment compatibility but did not reliably improve numerical precision or computation time.
- Site-balanced sampling is worth a numerical sensitivity check where site-omission support matters. Species-informed sampling increased sharing and improved fixed-count recovery in these high-order simulations, with an exclusion-ESS and time cost. The experiment supports neither a universal sampler choice nor a statistical-order cutoff.
- The fixed-kernel comparison supports the matched estimator in its stated sampling design. Fitted-score sandwiches are unsuitable as the primary regression interval in these experiments. Among the existing alternatives, deletion-score intervals are a fast regular-model sensitivity analysis and explicit site refits expose deletion support and provide a conservative comparator here. The results do not justify replacing the package's Bayesian-bootstrap default or relying on that bootstrap alone to substantiate a spread effect.
- Calibration used one independent-site, conditional-mean-correct generating mechanism with unconstrained coefficients. Spatial dependence, active constraints, irregular inference for nearly unidentified splines, and uncertainty in local ranks remain open. No interval claims are made for the five empirical communities from this phase. The empirical results are descriptive and cannot determine whether manuscript hypotheses remain supported under alternative intervals.
- Package verification used a source tarball and a full `R CMD check`, including the examples, vignettes and new exhaustive tests. Result: **Status OK**.
