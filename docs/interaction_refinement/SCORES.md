# Conditional score integration

For fixed parameters and an integer order $2\le k<n$, let $H_k(S;\theta)=c(S)X(S)[Y(S)-\mu(S;\theta)]$. Uniform global sets include a particular site with probability $k/n$. The full score $G_k$, conditional inclusion score $G_{ki}$ and deletion score $G_{k,-i}$ therefore satisfy

$$G_k=\frac{k}{n}G_{ki}+\left(1-\frac{k}{n}\right)G_{k,-i}.$$

Writing $a_k=k(n-1)/(n-k)$ gives $J_{ki}=G_k+a_k(G_{ki}-G_k)$. Every set appears in exactly $k$ conditional strata, so $G_k=n^{-1}\sum_i G_{ki}$. A focal path independently takes a random distinct prefix from the $n-1$ other sites. Giving every site weight $1/n$ retains the unconditional target even when production path counts are unequal. Equal weights here apply to the site means. The working weights remain part of $H_k$.

One prefix can be evaluated at several orders. For a cross-order linear contrast with block weights $L_k$, first evaluate $\sum_k L_ka_k H_k$ along each matched path and average within the focal site. Centre that average across sites and add $\sum_k L_kG_k$. This is algebraically the projection of the joint pseudovalues and retains within-path numerical covariance. Independently sampled order paths estimate the same contrast with a different numerical covariance. Each order evaluation of every path is included in the cost, even if several focal sites draw the same unordered combination. The software preserves the original global-pool alternatives.

## Exact controls

For any vector $Z(S)$ with exact conditional expectation $u_{ki}$ use $H_k(S;\theta) - A_i^T[Z(S)-u_{ki}]$. Estimate the coefficients $A_i$ from independently sampled focal paths, separately by site and order. The production estimator is conditionally unbiased given the pilot. The score plan freezes these coefficients, the reference basis and sampled paths before repeated evaluations at different parameters. Its exactness refers to control expectations. Numerical variance reduction depends on the correlation between each control and the evaluated score.

For relative abundance, $Z$ can include $Y=1-\sum_s\min_{j\in S}q_{js}$, whose exact mean comes from the local zeta threshold formula. For incidence, the model's set-normalised response is $Y=1-Z_k/r(S)$ where $r(S)$ is the mean richness. The default incidence response control evaluates $r(S)Y=r(S)-Z_k$. Its expectation follows by subtracting the exact shared richness mean from exact mean site richness. This variable remains an exact control for constant or richness working weights. The optional focal control uses $1-Z_k/r_i$, with the denominator fixed by the identity of the focal site. We define the control to be zero for an empty focal site. A richness-weighted model similarly gives all-empty sets zero score. The constant-weight incidence model requires sets of defined dissimilarity. The expectation of the model's mean-richness-normalised ratio is never approximated by a ratio of expectations.

For a site covariate $z$, the conditional mean of the set mean is

$$\frac{z_i+(k-1)\sum_{j\ne i}z_j/(n-1)}{k}.$$

For a spread basis $B$ and edge values $d_{jl}=|B(z_j)-B(z_l)|$, put $R_i=\sum_{j\ne i}d_{ij}$ and $T=\sum_{j<l}d_{jl}$. The exact conditional spread expectation is

$$\frac{\frac{k-1}{n-1}R_i + \frac{(k-1)(k-2)}{(n-1)(n-2)}(T-R_i)}{\binom{k}{2}}.$$

The implementation evaluates the $R_i$ with sorted prefix sums, requiring no site-by-site distance matrix. These exact predictor controls apply to the original environmental spread and mean position. They can also act as auxiliary variables for custom score matrices.

## Adaptive allocation

Let $s_i$ be a pilot estimate of the within-site standard deviation of the amplified path-level contrast score. For several contrasts, use the square root of the sum of their variance estimates in the user-supplied contrast units. With independent strata, the average estimation variance of the centred pseudovalues is proportional to $\sum_i s_i^2/M_i$. In the idealised equal-cost case this favours $M_i\propto s_i$. Divide the remaining budget after guaranteeing a minimum count to every site, then allocate the integer remainder by fractional priorities. The method keeps $1/n$ weights and uses independent production draws. When multiple effects enter the allocation criterion, its finite-budget performance depends on their relative scales and the stability of the pilot SDs.

## Numerical validation

Tests verify nested prefixes, focal-site sampling, exact predictor and response means, distinct incidence controls, exact matched deletion pseudovalues under exhaustive evaluation and the equality of direct contrasts to the joint projection. The numerical study uses one spatially structured niche community at each of 24 and 70 sites, with orders 2/4 and 2/8, respectively. The 24-site references enumerate both orders. At 70 sites, the order-two reference is exact and the order-eight reference averages two fresh independent 250,000-path runs. The accompanying tables report their discrepancy. At equal evaluation counts, 24 independent combination samples hold ecology, fixed coefficients and the information matrix fixed. The fresh timing-matched comparison uses 16 repetitions. Its budgets use linear interpolation or extrapolation between independent timing pilot runs at the two path counts and target the slowest median CPU cost at the larger budget. Score moments are projected through the inverse working information to express errors in linearised integrated spread and position coefficient units. Their global mean measures numerical one-step sensitivity of the point estimate. Centred site pseudovalue RMSE and covariance error measure the numerical uncertainty approximation. The latter compares the fixed-reference cross-order covariance with the sampled covariance, scaled by the product of the two reference marginal SDs.

Run from the installed package's repository root.

```sh
Rscript docs/interaction_refinement/score_validation.R 24 24
Rscript docs/interaction_refinement/score_validation.R 70 24
```

### Exact-reference comparison at 24 sites

| Method | Paired spread / position site RMSE, 11,520 evaluations | Median CPU for 11,520 evaluations | Paired spread / position site RMSE, matched CPU | Median matched CPU |
|---|---:|---:|---:|---:|
| Global deletion | .297 / .161 | .116 s | .290 / .158 | .081 s |
| Global pooled conditional ratios | .317 / .165 | .116 s | .300 / .163 | .073 s |
| Focal independent | .679 / .351 | .122 s | .474 / .257 | .091 s |
| Focal matched | .665 / .340 | .106 s | .467 / .239 | .093 s |
| Focal matched, response control | .695 / .348 | .148 s | .542 / .254 | .086 s |
| Focal matched, response and predictor controls | .692 / .330 | .219 s | .513 / .239 | .084 s |
| Focal matched, controls and adaptive allocation | .690 / .328 | .234 s | .543 / .256 | .086 s |

At median CPU costs of .073–.093 seconds, matching focal paths reduced spread and position change errors slightly relative to independent focal paths. Global pools had smaller site errors. A single global set contributes to all $k$ sites, whereas the focal estimator uses that draw for its focal site alone. The corresponding point sensitivity did little to reveal that difference. At 11,520 evaluations, paired spread one-step RMSE was .022 for global deletion and .033 for matched focal paths. At similar CPU the errors were .030 and .023. For paired position, cross-order covariance RMSE at similar CPU was .075 for pooled global conditional ratios and .101 for matched focal paths. At equal evaluations, adaptation lowered paired position site error slightly relative to unadjusted focal paths. Its pilot and computation costs outweighed this gain in the timed comparison.

### Manuscript-scale comparison at 70 sites

| Method | Paired spread / position site RMSE, 33,600 evaluations | Median CPU for 33,600 evaluations | Paired spread / position site RMSE, matched CPU | Median matched CPU |
|---|---:|---:|---:|---:|
| Global deletion | 1.086 / .391 | .391 s | 1.098 / .375 | .245 s |
| Global pooled conditional ratios | 1.085 / .377 | .381 s | 1.083 / .367 | .239 s |
| Focal independent | 3.249 / 1.141 | .362 s | 2.649 / .959 | .237 s |
| Focal matched | 3.277 / 1.147 | .315 s | 2.513 / .900 | .245 s |
| Focal matched, response control | 3.541 / 1.215 | .365 s | 2.822 / .978 | .242 s |
| Focal matched, response and predictor controls | 3.691 / 1.256 | .440 s | 3.007 / 1.047 | .243 s |
| Focal matched, controls and adaptive allocation | 4.112 / 1.362 | .434 s | 3.244 / 1.101 | .231 s |

The two high-budget references disagreed by .271 for spread and .104 for position on the scale of a single reference's estimated pseudovalue RMSE. That residual uncertainty is small relative to the comparison errors. At 33,600 evaluations, point spread one-step RMSE was .053 both for matched focal paths and for global deletion. Paired site error was about three times larger under the matched focal design. At similar CPU it remained about twice as large even though focal matching evaluated more combinations. Matching versus independent focal paths modestly reduced the site-level paired RMSE at similar CPU. The comparison did not show a consistent improvement in covariance recovery, which depends on the numerical precision of both levels. At equal CPU the pooled global cross-order covariance RMSE was .048 for spread and .042 for position, compared with .109 and .106 under the matched focal estimator.

The exact controls preserve the estimating equation. With the sites, models and pilot sizes used here the observed gains were insufficient to recover the precision lost by reducing production draws and assigning each path to a single site. Neither directly matched focal scores nor adaptive allocation provides a reason to change the global-pool default. These fixed-community experiments address numerical precision only. They do not establish ecological coverage for spatially dependent sites or the optimum control allocation in other communities.

The [24-site](results/score_24_summary.csv) and [70-site](results/score_70_summary.csv) summaries also retain the lower and higher order results, additional global score variants, the smaller evaluation budget and both point and site RMSE. The raw repeats retain every method and achieved CPU cost. Times can vary across calibration runs and machine load, so the timed comparisons use the measured fresh CPU costs.
