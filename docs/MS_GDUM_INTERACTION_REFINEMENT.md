# MS-GDUM refinement

This development branch extends the original nonlinear MS-GDUM in three places. The regression can represent nonlinear position and position-dependent monotone spread. Focal-site sampling and exact local zeta are available for evaluating fixed-parameter scores. A frozen-covariance Bayesian estimating-equation method provides a regular-interior sensitivity analysis. The original estimators and defaults remain available.

## Reproducing the results

Install this branch and run from the repository root. The study scripts write replicate-level results to [the results directory](interaction_refinement/results). For the ecological regression and focal-site score experiments, the linked reports document scenario generation, computational budgets and summary calculations.

```sh
R CMD INSTALL .
Rscript docs/interaction_refinement/model_validation.R recovery 24
Rscript docs/interaction_refinement/model_validation.R allocation 12
Rscript docs/interaction_refinement/model_validation.R curves 12
Rscript docs/interaction_refinement/summarise_models.R
Rscript docs/interaction_refinement/score_validation.R 24 24
Rscript docs/interaction_refinement/score_validation.R 70 24
Rscript docs/interaction_refinement/retention_validation.R
Rscript docs/interaction_refinement/conditioned_bayes_validation.R 12
Rscript docs/phase5_refinement/posterior.R long_reference null 1 1 900
Rscript docs/phase5_refinement/posterior.R long_reference null 2 2 900
Rscript docs/phase5_refinement/posterior.R long_reference weak 2 2 900
Rscript docs/phase5_refinement/posterior.R long_reference weak 3 3 900
Rscript docs/interaction_refinement/collect_jel.R
for signal in null weak; do
  Rscript docs/interaction_refinement/bayes_validation.R calibration "$signal" 20
  Rscript docs/interaction_refinement/bayes_validation.R nonlinear "$signal" 10
  Rscript docs/interaction_refinement/bayes_validation.R joint "$signal"
done
```

The simulation commands can run independently except where a summary consumes their results. The joint comparisons also use the two independently checked Phase 5 references already in this branch.

## Extended fractional-logit model

Let $z_{iv}$ be a site predictor standardised on training sites. For a set $S$ and predictor $v$, define $m_v=|S|^{-1}\sum_{i\in S}z_{iv}$ and the directional design $D_{vb}=\operatorname{mean}_{i<j\in S}|B_b(z_{iv})-B_b(z_{jv})|$ using the original I-spline basis. Predictions of nonlinear position use $|S|^{-1}\sum_{i\in S}p_v(z_{iv})$. In particular, $\overline{z^2}=m_v^2+\overline{(z-m_v)^2}$. A quadratic position function therefore responds to the within-set variance as well as to the squared mean. Comparisons that vary only the mean of transformed predictors can mistakenly assign this structure to directional spread.

For fixed bounds $L_v,U_v$ and $u_v=(m_v-L_v)/(U_v-L_v)$, use the columns $(1-u_v)D_{vb}$ and $u_vD_{vb}$. Both sets of coefficients are nonnegative, so the interpolated transform is monotone at each fixed position. For new sets outside the training mean range, interpolation is clipped and reported as extrapolation. This is a statement about the fitted directional transform at fixed position. The combined fitted prediction can move differently when mean, spread and site-level position functions change together.

For a monotone quadratic $f(z)=az+bz^2$, $f(z_i)-f(z_j)=(z_i-z_j)[a+b(z_i+z_j)]$. At order two the pair midpoint is exactly $m_v$. The constant-plus-midpoint slope supplied by linear varying coefficients can consequently duplicate a nonlinear I-spline transform. At higher orders the pair midpoints vary within the set. They can still track the site-set mean closely, and the resolution depends on observed joint support. Penalised estimates can be stable even with weak unpenalised information. The package reports working-information rank and scaled condition, active coefficients, interaction columns residualised on the main terms, and projections of declared coefficient contrasts into weak information directions.

The four fitted variants are linear position with constant spread, nonlinear position with constant spread, linear position with varying spread and nonlinear position with varying spread. The nonlinear position basis is cubic on the observed training interval with linear tails. Its coefficients are unrestricted. When position curvature is sufficiently penalised, its nonlinear terms shrink towards zero. The spread curvature matrices are evaluated for each endpoint and averaged. A separate penalty shrinks endpoint differences towards zero, and a fourth penalty shrinks non-intercept coefficients. Shrinkage and shape constraints affect fitted contrasts, so reported fit diagnostics use the unpenalised working information. Site-block smoothing fits standardisation and basis bounds on training sites only. The large number of sets within a block supplies numerical integration precision. The number of blocks supplies the validation units. This implementation uses held-out block selection and exposes the choice and its variation across folds. [Hastie & Tibshirani (1993)](https://doi.org/10.1111/j.2517-6161.1993.tb01939.x), [Pya & Wood (2015)](https://doi.org/10.1007/s11222-013-9448-7), [Marra & Wood (2011)](https://doi.org/10.1016/j.csda.2011.02.004)

Give difference-in-differences contrasts on both the link and response scales. The inverse logit alone can produce a nonzero response-scale difference-in-differences under an additive link. Provide matched configurations drawn from observed support before interpreting ecological mechanisms. For hierarchical deviance attribution, fit the additive submodel under the same restrictions and derive spread and position Shapley allocations there. Attribute the change in explained deviance from that additive fit to the full model as the incremental interaction. Penalty-selected fits can make the increment negative. Ecological interpretation should also consider the mean–variance relation in the community counts that generated the dissimilarities. [Warton et al. (2012)](https://doi.org/10.1111/j.2041-210X.2011.00127.x)

In 24 simulations with curved turnover and zero explicit interaction, the unpenalised joint model reported a fitted explicit contrast with absolute magnitude above .05 in 11 cases at order two and zero at orders six and twelve. The .05 value is a descriptive contrast threshold. Every order-two varying-spread information matrix was rank deficient. At order six with genuine positive/negative interactions, the unpenalised joint contrast RMSE was .0063/.0079, compared with .0291/.0362 under the combined penalty. With collinear predictors the median joint information condition number exceeded two million at order six. Site-block tuning in spatial simulations chose the same penalty after every fold omission in only eight of 24 selections. These results support reporting interaction contrasts and the support and penalty sensitivity together. [Extended model derivation, simulations and results](interaction_refinement/MODEL_VALIDATION.md)

## Focal-site scores

A uniform random permutation of the $n-1$ other sites gives the conditional uniform set law at every order when combined with the focal site. Prefixes give matched paths across orders. Each focal site's conditional mean receives weight $1/n$ regardless of draw allocation. Average these means to obtain the global mean. The complete-set deletion decomposition is $G=(k/n)G_i+(1-k/n)G_{-i}$, so its exact pseudovalue is $G+k(n-1)(G_i-G)/(n-k)$. This identity requires $k<n$. For a paired contrast the same draws can estimate its conditional influence directly. A fixed control variate with an analytically known conditional mean retains the target; a pilot independent of the final draws estimates the regression adjustment. Overlapping site sets motivate distinguishing combination-sampling error from ecological sampling uncertainty. [Chen & Kato (2019)](https://arxiv.org/abs/1712.00771)

At matched CPU cost in the 24-site exhaustive comparison, pairing focal paths reduced paired spread/position score RMSE from .474/.257 to .467/.239, while pooled global scores reached .300/.163. Among 70 sites at 33,600 total evaluations, the corresponding RMSE was 3.277/1.147 for matched focal paths and 1.085/.377 for pooled global scores. At .231–.245 CPU seconds in the 70-site study, paired score errors were 2.513/.900 for matched focal paths and 1.083/.367 for pooled global scores. At the matched evaluation budget, both the focal and global deletion point spread one-step RMSEs were .053. The targeted focal draws improve the direct pairing relative to independent focal draws, while a global path supplies information for several focal strata at once. With these communities and budgets, that reuse outweighed the conditional control and allocation gains. [Conditional-score derivations and equal-budget comparisons](interaction_refinement/SCORES.md)

## Local retention and environmental interpretation

The full conditional zeta curve is the survival function of retained focal composition along a uniformly sampled path. Occurrence gives each focal species weight $1/r_i$. Relative abundance weights thresholds by the corresponding increments of relative abundance. A threshold occupied at $m$ sites including the focal site survives order $k$ with probability $\binom{m-1}{k-1}/\binom{n-1}{k-1}$. Summing this curve gives $n/(n-m+1)$, so the contribution to the full index follows the requested closed forms. The index conditional on retention through order two along the path, $T_i^+$, is computed using the directly summed tail. The loss distribution reports remaining mass through the final evaluated order as right-censored.

Across 40 enumerated communities and both response types, maximum profile and full-index errors were $2.22\times10^{-16}$ and $1.78\times10^{-15}$. In the simulated saturation scenario, median occurrence retention reached 70 at every site while median relative-abundance retention was 22.17. Across 12 such communities, the median within-community correlation of relative retention with fitted order-eight dissimilarity was −.630. It was −.927 with unsaturated occupancy. These are descriptive associations over the same sites. The detailed report gives spread and position allocations and a temporal comparison using common stations or an explicitly truncated order. [Retention derivation and validation](interaction_refinement/RETENTION.md)

## Frozen-covariance Bayesian estimating equations

Stack the unconditional fractional-logit estimating equations $G(\theta)$ across orders. Each order averages its working-weighted score under the target set measure, with the package sampling correction applied separately. At a regular interior pilot root, compute the joint site pseudovalues $J_i$ under independent set integration when available. Their sample covariance $\Omega$ preserves covariance between orders. Evaluate $-nG(\theta)^T\Omega^{-1}G(\theta)/2$ with $\Omega$ fixed and add the log normal prior. This estimates the ecological covariance only insofar as the site pseudovalues resolve their site signal. Higher-order dependence or spatial dependence can invalidate the simple site scaling. [Yin (2009)](https://saasresearch.hku.hk/~gyin/materials/2009YinBA.pdf), [Lyddon et al. (2019)](https://doi.org/10.1093/biomet/asz006)

Let $H=-\partial G/\partial\theta$ at the unrestricted root $\widehat\theta$. The local loss precision is $Q=nH^T\Omega^{-1}H$, with prior-free sandwich covariance $H^{-1}\Omega H^{-T}/n$. For zero-mean prior precision $P$, the Gaussian reference has covariance $(Q+P)^{-1}$ and mean $(Q+P)^{-1}Q\widehat\theta$. The nonlinear integration uses the same prior and exact moment loss with frozen covariance. Singular matrices stop under a declared scaled eigenvalue tolerance. A custom unpenalised model matrix can be used for regular interior sensitivity analysis, with the extended-model basis held fixed and the same layout supplied to site-score integration.

The optional experimental boundary sensitivity projects unrestricted draws onto the environmental spread orthant by minimising their displacement in metric $Q$. Profiling out unconstrained coefficients leaves a nonnegative least-squares problem, after which the free-coordinate adjustment is restored. This produces draws with boundary mass. An active constrained fit instead satisfies a KKT condition, and the unrestricted zero-score quasi-posterior does not supply its sampling distribution.

### Repeated calibration

Twenty independent 70-site simulations per signal used orders 2 and 8, 5,000 combinations for the global score and an independent 15,000 for the frozen covariance. The pilot used centred deletion scores. The computational pseudovalue comparator bootstrapped linearised coefficient influences with common Dirichlet site weights. Its timing covers bootstrap generation after the influences have been prepared. Point identification in this experiment uses three unrestricted coefficients per order and a controlled conditional model with shared site noise. Reported rates average spread and position within each layer while the archives preserve each effect. There are 20 independent ecological samples per signal and within-sample outcomes are correlated.

| Signal and method | 50% coverage, lower / higher / paired | 95% coverage, lower / higher / paired | One-sided positive decision rate, lower / higher / paired |
|---|---:|---:|---:|
| Null, pseudovalue | .475 / .475 / .450 | .875 / .850 / .925 | .075 / .075 / .000 |
| Null, Gaussian GMM | .500 / .475 / .450 | .850 / .850 / .925 | .075 / .075 / .000 |
| Weak, pseudovalue | .475 / .475 / .525 | .900 / .850 / .925 | .375 / .350 / .475 |
| Weak, Gaussian GMM | .500 / .500 / .525 | .850 / .825 / .925 | .375 / .300 / .475 |

A positive decision required posterior probability at least .95. The weak-signal paired effects have opposite signs. For Gaussian GMM, the positive decision rates were .950 for the positive position change and .000 for the negative spread change. The negative spread change was detected in .600 at the corresponding negative decision threshold. Gaussian GMM took a median .14 seconds including the fit and covariance, versus .04 seconds for pseudovalue draws after score preparation. Some 95% coverages were below nominal with both methods. The first four replicates per scenario were also used to test projection. In the four matched null replicates, one-sided positive decisions were .125 at both level orders with either method. In the weak subset, lower-order decisions were .500 before projection and .625 after it. This subset is too small for calibration.

### Nonlinear integration and JEL

The comparison used six independently verified JEL references on identical 900-combination pools. The two Phase 5 references were augmented with four independent long-chain runs. All six references met the joint criteria of split $\widehat R\leq1.05$ and minimum ESS 200. The additional runs took 78–81 seconds each. Nonlinear GMM took a median 1.18 seconds under the null and 1.24 under weak effects. All its minimum ESS values exceeded 534 and all split $\widehat R$ values were below 1.012.

| Signal and comparison | Largest absolute gap in positive probability, lower / higher / paired | Largest absolute gap in mean link effect, lower / higher / paired |
|---|---:|---:|
| Null, Gaussian versus nonlinear GMM | .033 / .028 / .025 | .011 / .027 / .016 |
| Weak, Gaussian versus nonlinear GMM | .012 / .037 / .045 | .014 / .022 / .017 |
| Null, nonlinear GMM versus JEL | .087 / .404 / .286 | .038 / .362 / .375 |
| Weak, nonlinear GMM versus JEL | .043 / .418 / .384 | .033 / .567 / .587 |

The differences between GMM and JEL reflect differences between update rules on fixed data, and cannot establish which method is better calibrated. As a separate check, nonlinear GMM was run on the first ten independent 5,000-set simulation replicates per signal, with the same frozen 15,000-set pilot. All twenty chain sets passed the numerical check. The null 95% coverage was .750/.800/.950, compared with .800/.800/.950 for Gaussian GMM on the same replicates. Weak-signal rates were .750/.750/.950 with both integrations. On these samples approximately five seconds of additional nonlinear integration did little to change coverage. This points to uncertainty in the site-score covariance and update rule as issues to investigate alongside posterior computation. Repeated JEL calibration, spatial calibration and applications with active constraints remain open.

### Covariance from conditioned scores

Twelve 70-site niche communities used orders 2 and 8. The covariance pilot compared 14,000 global coupled paths with 14,000 paths divided equally among focal sites, either directly or with 45 pilot paths and 155 production paths per site. Both conditional plans count all focal paths and both orders. The median relative Frobenius distance to the global moment covariance was .078 without controls and .055 with response and predictor controls. The median order-eight position SE was .218 using the global pilot and .314 or .312 using the conditional pilots. For the paired position change the corresponding SEs were .196, .298 and .297. None of these three pilot covariances is an exact ecological reference for order eight. This comparison establishes integration of the new scores into the Gaussian inference and illustrates sensitivity of parameter uncertainty to the covariance pilot. It does not establish which covariance is better calibrated.


## Interpretation and open questions

- Nonlinear position can improve recovery when position itself is curved. Use the rank and support diagnostics before interpreting a fitted interaction. At order two the explicit varying coefficient can absorb ordinary nonlinear turnover, and at higher orders the estimated effect can be sensitive to the curvature and interaction penalties.
- Use held-out site or spatial blocks to choose candidate penalties. Report the sensitivity of declared contrasts to candidate strengths and fold omissions. Post-selection interval calibration under spatial dependence remains open.
- The new focal-site estimator has an exact unconditional target and closes the connection with analytical local zeta. The matched global-pool estimators performed better for site-score uncertainty in the tested communities. Results could depend on site count, cost per kernel evaluation, variance of the control variables and the particular contrast.
- Gaussian GMM is a fast sensitivity analysis for regular interior fits with independently verified site-score integration. Its agreement with nonlinear GMM addressed numerical integration in the tested cases. Undercoverage for some level effects and discrepancies with JEL leave ecological calibration unresolved. Projection remains experimental.
- The original nonlinear MS-GDUM is the primary ecological model. The selected results do not establish interaction mechanisms in the manuscript communities. The validation uses simulated communities and fixed-posterior audits.

The source package returned `R CMD check` status OK. The package suite passed 1,496 expectations with one existing CRAN-only skip. Model tables were regenerated from the versioned replicate files and matched the committed summaries.
