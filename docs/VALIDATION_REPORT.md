# Validation of species-informed sampling

Phase 2 investigates order-specific inference, ecological sampling uncertainty,
and adaptive set sampling. [Phase 2 validation report](PHASE2_VALIDATION_REPORT.md)

The mixture targets the same expected weighted model loss as conventional sampling. It is useful for drawing more shared sets, particularly when sharing is rare. That gain comes with lower weight-based effective sample size and sometimes a higher computational cost. The tests below evaluate numerical exactness, estimator behaviour, and sensitivity of selected manuscript results at the default targeted fraction of 0.5.

## Exhaustive verification

Enumerated all 20 order-three sets of a six-site community, including sets sharing zero, one, or two species. Tested a uniform target, a product-of-site-weights target with an ineligible site, and a product target with site weights scaled by $10^{270}$. For each target, the test compared the log-domain normalising calculations, the per-species selection probabilities, the set probability ratios, and 30,000 sampled sets at `alpha = 0.6`.

| Base target | Exact mean shared species | Absolute error in mean shared | Maximum absolute weight error | Maximum absolute sample frequency error |
| --- | ---: | ---: | ---: | ---: |
| Uniform | 0.500000 | 3.89e-16 | 4.44e-16 | 0.00490 |
| Weighted, with an ineligible site | 0.759271 | 1.11e-16 | 4.44e-16 | 0.00394 |
| Weighted, with extreme common scaling | 0.646057 | 7.66e-15 | 6.22e-15 | 0.00143 |

The tests include unshared sets, an order above every species' eligible occupancy, and a single supported set. A separate log-sum check used site weights spanning 510 orders of magnitude. On exhaustive fits with uniform and product-weighted targets, applying the proposal correction recovered the model, variable importance, and Shapley results to the test tolerance of $10^{-8}$.

Across every order of the small reference community, the largest absolute errors in the exact contributions were $2.78\times10^{-16}$ for unnormalised incidence and $8.33\times10^{-17}$ for relative abundance. The tests also verify that averaging the set-level richness-normalised incidence contributions can give different results from the analytical terms.

## Parameter recovery and interval coverage

Used 20 independently simulated communities per structure, with up to 18 sites after removing empty sites. The broad-niche communities used median breadth 1.50 and the narrow-niche communities used 0.53. The patchy communities used a mixture of geographic ranges, niche breadths and local abundances, with 75% assigned to the dominant small-range, narrow-niche, high-abundance form. At order 10 the median exact fraction of shared sets was 1.00, .187, and 1.00 in the three structures, respectively. Patchy communities varied considerably, down to .069 at that order. Each order used 220 sets and 49 paired site-weight draws, and the reference was the exhaustive fit to the same community under the same response working weights. The study covered 360 community–order–response targets.

Values are conventional / species-informed. Error columns give the median absolute error pooled across the two gradients. The last column is the median signed training D² error relative to the exhaustive fit.

| Response | Order | Intercept error | Spread-height error | Position-coefficient error | D² error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Total incidence | 3 | .041 / .041 | .102 / .106 | .0154 / .0149 | +.0077 / +.0117 |
| Total incidence | 6 | .083 / .079 | .169 / .158 | .0107 / .0118 | +.0103 / +.0116 |
| Total incidence | 10 | .137 / .114 | .283 / .278 | .0099 / .0122 | +.0045 / +.0060 |
| Relative abundance | 3 | .040 / .046 | .099 / .090 | .0154 / .0168 | +.0129 / +.0126 |
| Relative abundance | 6 | .097 / .075 | .194 / .160 | .0122 / .0114 | +.0123 / +.0067 |
| Relative abundance | 10 | .164 / .138 | .314 / .272 | .0104 / .0098 | +.0165 / +.0101 |

In narrow-niche communities at order 10, the median shared draw count increased from 37 to 133 in point samples and from 30.5 to 127 in bootstrap samples. Point ESS went from 220 to 130 and bootstrap ESS from 220 to 127. Median absolute spread-height errors fell from 1.451 to .891 for incidence and from 1.978 to .940 for relative abundance. Across all structures at that order, median signed spread-height errors were +.044 and 0 for incidence and +.034 and 0 for relative abundance. The mixed position-coefficient errors and the small positive D² errors persisted after correction.

The nominal 95% site-bootstrap intervals contained the exhaustive finite-community spread and position estimates in all tested cases, under both methods. The intervals include site-weight variation around a fixed community. At order 10 the median incidence spread-height interval widths were 3.23 and 3.48 and the position widths were .102 and .108. This finite-community reference exercise measures whether intervals bracket the exhaustive fit.

For repeated-experiment coverage, an additional simulation used 24 independent realisations of each community structure and an artificial bounded response. Its conditional mean was a logistic GDUM with intercept .4, spread heights 1.15 and .79, and position coefficients .07 and −.04. Zero-mean site disturbances contributed to every set containing a disturbed site. The benchmark used 350 sets and 59 bootstrap draws at orders 4 and 14. The responses in this experiment were generated by the known model. The separate exhaustive-fit experiment above used ecological dissimilarities.

| Order | Parameter | Mean signed bias | Mean absolute error | 95% interval coverage |
| --- | --- | ---: | ---: | ---: |
| 4 | Intercept | −.0073 / −.0014 | .0992 / .1005 | 86.1% / 83.3% |
| 4 | Spread heights | +.0974 / +.0877 | .2661 / .2676 | 80.6% / 84.7% |
| 4 | Position coefficients | +.0036 / +.0038 | .0223 / .0209 | 89.6% / 88.2% |
| 14 | Intercept | −.0074 / −.0054 | .1399 / .1342 | 63.9% / 65.3% |
| 14 | Spread heights | +.0868 / +.0801 | .3244 / .3155 | 62.5% / 68.1% |
| 14 | Position coefficients | +.0011 / +.0011 | .0068 / .0069 | 74.3% / 74.3% |

Both methods undercovered the generating parameters, especially at the higher order. These coverage results apply to this bounded shared-site error mechanism and the stated draw budgets. They provide no basis for treating high-order percentile intervals as nominally calibrated in the manuscript.

## Higher orders and computational budgets

Used six independent communities per structure with up to 100 sites. The across-community medians of median species occupancy were .745, .207, and .080 for broad, narrow, and patchy structure. The corresponding medians of maximum occupancy were 1.000, .859, and .945. Orders were 3, 8, 16, and 24. Timing covers the point fit plus nine Bayesian-bootstrap fits with two-group Shapley decomposition. Held-out results use two identical 60/40 site splits per community and the same 2,000 conventional test sets per split for every training comparison. The tables give medians of the six community results. Each paired cell is conventional / species-informed.

### Equal draw budget

Both methods drew 500 sets for each point and bootstrap fit. The baseline fraction sharing is from a separate 4,000-set conventional sample per community.

| Structure | Order | Baseline sharing | Shared point draws | Point ESS | Point D² | Held-out D² | Held-out RMSE | Seconds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Broad | 3 | 1.000 | 500 / 500 | 500 / 493 | 0.672 / 0.660 | 0.536 / 0.537 | 0.0628 / 0.0628 | 0.059 / 0.090 |
| Broad | 8 | 1.000 | 500 / 500 | 500 / 491 | 0.639 / 0.614 | 0.520 / 0.509 | 0.0575 / 0.0585 | 0.071 / 0.099 |
| Broad | 16 | 1.000 | 500 / 500 | 500 / 491 | 0.549 / 0.517 | 0.166 / 0.157 | 0.0485 / 0.0470 | 0.088 / 0.116 |
| Broad | 24 | 1.000 | 500 / 500 | 500 / 489 | 0.483 / 0.460 | -0.205 / -0.122 | 0.0367 / 0.0351 | 0.104 / 0.134 |
| Narrow | 3 | 0.925 | 458 / 482 | 500 / 435 | 0.678 / 0.684 | 0.624 / 0.618 | 0.0714 / 0.0726 | 0.067 / 0.087 |
| Narrow | 8 | 0.479 | 238 / 372 | 500 / 368 | 0.331 / 0.322 | 0.201 / 0.216 | 0.0365 / 0.0346 | 0.075 / 0.093 |
| Narrow | 16 | 0.117 | 61 / 281 | 500 / 287 | 0.266 / 0.306 | -0.106 / -0.132 | 0.0200 / 0.0195 | 0.095 / 0.116 |
| Narrow | 24 | 0.027 | 10 / 258 | 500 / 258 | 0.284 / 0.315 | -0.661 / -0.216 | 0.0110 / 0.0072 | 0.120 / 0.136 |
| Patchy | 3 | 0.995 | 499 / 500 | 500 / 481 | 0.141 / 0.144 | 0.072 / 0.081 | 0.0867 / 0.0871 | 0.056 / 0.071 |
| Patchy | 8 | 0.869 | 433 / 470 | 500 / 451 | 0.093 / 0.141 | -0.003 / -0.003 | 0.0665 / 0.0679 | 0.067 / 0.084 |
| Patchy | 16 | 0.532 | 270 / 371 | 500 / 380 | 0.104 / 0.123 | 0.053 / 0.016 | 0.0566 / 0.0574 | 0.085 / 0.099 |
| Patchy | 24 | 0.288 | 151 / 328 | 500 / 319 | 0.116 / 0.097 | 0.083 / -0.009 | 0.0471 / 0.0475 | 0.103 / 0.125 |

For narrow-niche communities at order 24, the median bootstrap shared count increased from 8.5 to 252, with median bootstrap ESS falling from 500 to 257. Conventional sampling yielded one undefined point D² in six communities at each of orders 16 and 24 because the response was constant in the sampled sets. The mixture yielded defined point D² in all six communities at each order. At order 24 the mixture reduced median held-out RMSE by about 34%, while both held-out D² medians were below zero. At order 16 the narrow-niche comparison gave lower RMSE and lower held-out D² under the mixture. RMSE measures squared prediction error. Held-out D² measures deviance relative to each method's working-weighted training mean.

### Matched elapsed time

Separate timings at 200 and 800 sets calibrated each sampler's draw budget to a common target elapsed time per structure and order. Realised medians at order 24 were within about 8% between methods. The table shows orders 16 and 24.

| Structure | Order | Draw budget | Seconds | Shared point draws | Point ESS | Held-out D² | Held-out RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Broad | 16 | 620 / 395 | 0.103 / 0.103 | 620 / 395 | 620 / 388 | 0.128 / 0.159 | 0.0482 / 0.0457 |
| Broad | 24 | 593 / 413 | 0.122 / 0.121 | 593 / 413 | 593 / 404 | -0.208 / -0.165 | 0.0364 / 0.0375 |
| Narrow | 16 | 573 / 437 | 0.108 / 0.103 | 68 / 240 | 573 / 251 | -0.063 / -0.140 | 0.0199 / 0.0193 |
| Narrow | 24 | 522 / 483 | 0.123 / 0.132 | 11 / 249 | 522 / 250 | -0.647 / -0.256 | 0.0110 / 0.0069 |
| Patchy | 16 | 554 / 456 | 0.091 / 0.094 | 298 / 338 | 554 / 347 | 0.056 / 0.026 | 0.0550 / 0.0572 |
| Patchy | 24 | 523 / 481 | 0.112 / 0.119 | 157 / 314 | 523 / 309 | 0.090 / 0.064 | 0.0465 / 0.0458 |

### Identical CPU limit

An independent point-fit benchmark used the same 0.120-second CPU allowance for both samplers on each of six communities per structure and order. Each method added batches of 256 sets and refitted on the accumulated sets until the next fit exceeded the allowance. The last fit completed inside the allowance was scored using 3,000 common conventional held-out sets. These timings include repeated point fitting and exclude bootstrap and held-out scoring.

| Structure | Order | Draws completed | CPU seconds used | Shared draws | ESS | Point D² | Held-out RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Narrow | 16 | 3,200 / 2,176 | .110 / .112 | 567 / 1,287 | 3,200 / 1,275 | .390 / .391 | .0158 / .0152 |
| Narrow | 24 | 2,944 / 2,048 | .108 / .107 | 179 / 1,024 | 2,944 / 1,039 | .301 / .454 | .0037 / .0066 |
| Patchy | 16 | 3,584 / 2,560 | .115 / .111 | 3,154 / 2,417 | 3,584 / 2,332 | .127 / .127 | .0520 / .0516 |
| Patchy | 24 | 3,072 / 2,304 | .111 / .109 | 2,478 / 2,083 | 3,072 / 2,046 | .119 / .113 | .0385 / .0385 |

Under this fixed CPU limit, the mixture produced more shared draws for narrow-niche communities and fewer for patchy communities. In this order-24 narrow-niche benchmark, the mixture had higher training D² and higher held-out RMSE. The result differs from the small equal-draw experiment at order 24, which used a different set of simulated communities and a smaller training budget.

## Empirical sensitivity analysis

Re-ran whole-community incidence and abundance models from the Finnish breeding-bird, Southwest Australian plant, and Swiss butterfly case studies. The comparison kept the site filters, predictors, response families, and predictor groups from the existing analyses. Orders were 2, 8, 15, and 28 for birds, 2, 6, and 9 for plants, and 2, 8, and 13 for butterflies. Fits used 2,000 sets and 39 paired bootstrap draws per method. Held-out results are medians over 10 identical 70/30 site splits per method, using common conventional test sets. These are targeted sensitivity runs at a reduced budget relative to the original 10,000-set, 100-draw model runs.

The table gives the highest tested order, with conventional followed by species-informed sampling in each cell. Shared counts and ESS refer to the point sample. RMSE uses the working-weighted response scale.

| Case and response | Order | Shared draws | ESS | Point D² | Held-out D² | Held-out RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Birds, turnover | 28 | 1,599 / 1,781 | 2,000 / 1,692 | .224 / .223 | .079 / .119 | .0438 / .0434 |
| Birds, balanced | 28 | 1,599 / 1,781 | 2,000 / 1,692 | .320 / .312 | .326 / .338 | .0514 / .0518 |
| Birds, relative | 28 | 1,599 / 1,781 | 2,000 / 1,692 | .336 / .316 | .299 / .281 | .0248 / .0253 |
| Plants, total | 9 | 311 / 1,148 | 2,000 / 1,160 | .258 / .281 | .234 / .237 | .0050 / .0048 |
| Plants, turnover | 9 | 311 / 1,148 | 2,000 / 1,160 | .209 / .234 | .166 / .181 | .0103 / .0100 |
| Butterflies, turnover | 13 | 1,655 / 1,820 | 2,000 / 1,733 | .114 / .111 | .103 / .115 | .0726 / .0725 |
| Butterflies, balanced | 13 | 1,655 / 1,820 | 2,000 / 1,733 | .150 / .150 | .122 / .117 | .0149 / .0149 |

Median bootstrap Shapley rankings at the highest tested order agreed between methods in all seven comparisons. In plant total dissimilarity at order 9, for example, the median precipitation contribution was .150/.142 and the temperature contribution was .081/.084. The birds' relative-abundance held-out D² decreased while their incidence-turnover D² increased. The three case studies show mixed out-of-sample changes overall. Reruns cover the whole-community sampler comparisons. The full HIATES archive and all occupancy-class, rarity-form, and manuscript figure workflows were outside this sensitivity run.

## Limitations and manuscript implications

- A large shared draw count indicates greater representation of shared combinations in the proposal. Correction weights can still give those draws little influence under the target. When sharing is already common, extra sampling work can reduce the number of completed shared draws under a fixed time limit.
- The mixture is a targeted alternative for studying orders where conventional draws rarely share species. The comparisons support higher shared-draw counts in some sparse regimes and smaller coefficient errors in the narrow-niche exhaustive comparison. ESS declined, and results for parameter accuracy, bootstrap coverage, and held-out prediction depended on the scenario and metric.
- High-order weighting can produce a small effective set sample and conditional Poisson plus per-species normalisation has computational overhead. Training and held-out deviance are unstable when the response is almost constant. At orders beyond possible sharing, the model can only estimate an intercept for total or turnover dissimilarity on nonempty sites.
- Randomly splitting sites evaluates the stated held-out task. Transfer to geographically independent regions and performance across a broader empirical collection remain unmeasured. Timing results are specific to the sampled communities, implementation, and hardware.
- The manuscript can present the analytical incidence and relative-abundance contributions with their distinct units and include the mixture as a sensitivity analysis. Report the order, original response working weights, proposal, retained shared draws, and ESS alongside high-order model summaries. The geographic-distance basis is fixed across proposals in these comparisons, so fits with coordinates may differ from older fits that used sample-derived distance bounds.

## Reproduction and package checks

The script [validate_sampling.R](validate_sampling.R) generates the numerical results with base R and the development package on macOS arm64, R 4.6.1. It reads the bird inputs from the sibling manuscript workspace and the plant inputs from the installed `gdm` package. For butterflies, set `ZETAGDUM_BUTTERFLIES` to the cleaned `bdm_butterflies_clean.rds` input.

```sh
R CMD INSTALL .
Rscript docs/validate_sampling.R /tmp/zeta-validation ../paper_mee
```

The final source-package `R CMD check` completed with **0 errors, 0 warnings, and 0 notes**. It used RStudio's bundled Pandoc and the installed TeX fonts via `R_RD4PDF=times,hyper`. The check omitted a local startup override that forced obsolete C++14. The complete development test run passed 258 expectations, including the optional installed `gdmmTMB` comparison. Standard `R CMD check` passed 257 expectations and skipped that comparison by its existing CRAN guard.
