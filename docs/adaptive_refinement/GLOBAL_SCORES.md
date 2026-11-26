# Global score control variates and independent-pool covariance

## Score estimand

Use global coupled paths throughout. A uniform size-$k$ set contains any particular site with probability $p_k=k/n$. Write the vector kernel as $g_k(S;\theta)=c(S)X(S)\{Y(S)-\mu(S;\theta)\}$ and denote the target/proposal ratio by $w(S)$. For any fixed $\theta$ and basis, the sample mean of $w(S) I_i(S)g_k(S;\theta)/p_k$ is unbiased for the conditional score $G_{ki}$ under the original uniform-set target. This uses raw ratios, without dividing by their sample mean or the observed inclusion frequency. Each path contributes to all its included sites.

For each order and site, let $Z(S)=\{1,Y^*(S)\}$ and $u_i=E_P[Z(S)\mid i\in S]$. Estimate $A_i$ by regressing $wI_i g_k$ on $wI_i Z$ in the independent pilot, including the zeros from sets that omit the site. The adjusted production estimator is
$$
\widehat G_{ki}=\frac{1}{p_k}\Big[\overline{w I_i g_k}
-A_i^T\big(\overline{w I_i Z}-p_k u_i\big)\Big].
$$
The centred control expectation is exactly zero under the proposal distribution. The variance reduction depends on the pilot and the evaluated parameter. The estimator remains conditionally unbiased at other fixed parameter values because the control coefficients and predictor basis are frozen. Substituting the same sampled moments into the pilot and production steps invalidates this argument.

For relative abundance, $Y^*=Y$ with exact conditional expectation from LCZD. For incidence, the control is $Y^*=r(S)-Z_k(S)$, where $r(S)$ is mean richness and $Z_k$ is shared richness. This works with either constant working weights or richness working weights because the score kernel retains its declared $c(S)$ and $Y(S)$. No calculation treats the expected ratio as a ratio of conditional expectations. The exact unconditional score is the equal mean across focal sites. The software estimates it by the equal mean of the adjusted conditional scores and uses the complete-set deletion identity to form each pseudovalue. Both the unadjusted and adjusted score estimates use the same global coupled allocation.

The implementation uses sparse site accumulations. It stores path-level moments without creating a combination-by-site incidence matrix. Standardising the pilot control covariance separates weak scale from rank. Unidentified control directions receive zero coefficients.

## Independent-pool covariance

Let $A_r$ contain the joint site pseudovalues, centred within each pool across the same sites. The average covariance is $\Omega_{avg}=\{\sum_r A_r^T A_r\}/[R(n-1)]$. The pooled covariance uses the average site vector. All distinct cross-products are estimated through
$$
\Omega_{cross}=\frac{R\Omega_{pooled}-\Omega_{avg}}{R-1},\qquad
\Omega_{MC}=\Omega_{avg}-\Omega_{cross}.
$$
At $R=2$ this gives $\Omega_{MC}=(A_1-A_2)^T(A_1-A_2)/[2(n-1)]$. Tests also compare the $R=4$ formula with averaging all six pairwise cross-products. For a common fixed parameter, identical basis, and independent conditionally unbiased pool estimates, the expected cross-product equals the finite-community covariance of exact pseudovalues conditional on the sites. This is a statement about numerical error. The exact pseudovalue covariance can itself differ from the ecological sampling covariance of an incomplete U-statistic or a fitted nonlinear parameter.

Scale the selected covariance by the average within-pool score SDs, which remain available even when a cross-covariance diagonal is negative. The diagonal target loses cross-order dependence. Block shrinkage instead retains the observed within-block correlations of the average covariance. Blocks in the study join the same coefficient type across orders. Where there are at least three pools, the approximate shrinkage intensity is a ratio of the standardised leave-one-pool-out numerical variance to the squared distance from the target, capped at one. This is a conditional numerical adaptation inspired by linear shrinkage. It does not use the number of combinations as an ecological effective sample size and does not assert the optimality of the large-dimensional estimator for these dependent score rows. [Ledoit and Wolf (2004)](https://doi.org/10.1016/S0047-259X(03)00096-4)

The PSD projection replaces negative eigenvalues with zero in standardised score coordinates. This can remain singular. An explicit positive eigenvalue floor is a distinct sensitivity analysis and is reported along with the before and after spectra and amount of adjustment. The unregularised covariance and requested contrast covariances are always available. The leave-one-pool-out SE quantifies conditional numerical variability in entries and declared contrasts. It becomes noisy with few pools and does not capture ecological sampling variability.

## Reproduction

Run these scripts from an installed development package at the repository root.

```sh
Rscript docs/adaptive_refinement/score_validation.R
Rscript docs/adaptive_refinement/covariance_validation.R
```

The score study compares independent single-community replicates of coupled set-sampling errors at orders 2/4, 2/8, and 3/12. It uses 16 numerical replicates per budget and model, with both response types at 70 sites. Four independent production pools contribute to each estimate, with 20% of the evaluation budget reserved for a shared independent control-coefficient pilot when controls are used. We report the realised median CPU time including the pilot, path generation, environmental design, response, and score evaluation. Covariance postprocessing is common to the sampling methods and is small relative to building the path matrices. Reference site scores are exact at 24 sites. At 70 and 150 sites the comparison uses the mean of two fresh 80,000-path references. The archived reference gap measures their discrepancy. All effects use inverse working information at the same fixed reference parameter to place errors in coefficient units.

The covariance study uses 240 independent ecological replicates per site count (70, 150 and 250) and per cross-order dependence scenario. It has six Gaussian score coordinates with matched intercept, spread and position effects across orders. The matched cross-order correlation is .65 or .99, and four independent pools have additive, correlated numerical error. We compare directly with the exact finite-community covariance and the known population covariance. The reported probability comparison uses a fixed contrast centre and an oracle SD to isolate covariance differences. It does not assess parameter bias or full posterior interval calibration.

## Global score validation results

| Sites, response | Median CPU, no controls / both (s) | Paired spread site RMSE, no controls / both | Cross-order covariance RMSE, no controls / both | Paired spread variance ratio, pooled / cross without controls |
|---|---:|---:|---:|---:|
| 24, relative | .072 / .072 | .109 / .126 | .042 / .036 | 1.031 / .970 |
| 70, relative | .074 / .081 | .654 / .646 | .045 / .033 | 1.147 / .992 |
| 70, incidence | .076 / .079 | .555 / .637 | .050 / .043 | 1.182 / 1.035 |
| 150, relative | .090 / .100 | .533 / .605 | .063 / .061 | 1.224 / 1.022 |

Site RMSE compares the pooled mean pseudovalues with the finite-community reference. Covariance error is the standardised cross-order block RMSE of the cross-pool estimator. Held-out pools preserved the declared target for both response types. The inclusion-plus-response controls made cross-order covariance more accurate by this measure at similar times in these communities, while their paired spread scores were usually slightly noisier. For 70 relative the .008 score RMSE difference is small compared with the .188 discrepancy between the two numerical reference scores. At 150 sites only .003 of the order-12 reference sets shared any species apart from the background species. The pooled variance remained inflated despite the site-score gains from recycling each global set. Cross-pool paired variance was closer to the reference in all four cases shown. [Global score summary](results/score_summary.csv)

## Known covariance validation results

| 70-site method | Ordinary finite covariance RMSE | Cancellation finite covariance RMSE | Cancellation paired spread variance / truth | Cancellation probability available | Cancellation probability MAE |
|---|---:|---:|---:|---:|---:|
| One raw pool | .235 | .238 | 29.790 | 1.000 | .243 |
| Average within pools | .202 | .205 | 30.559 | 1.000 | .244 |
| Pooled mean site score | .078 | .080 | 8.415 | 1.000 | .170 |
| Cross-pool | .062 | .064 | 1.092 | .729 | .054 |
| Cross-pool with PSD projection | .062 | .064 | 1.400 | 1.000 | .052 |
| Cross-pool with block shrinkage | .062 | .065 | 2.702 | .925 | .083 |
| Cross-pool with diagonal shrinkage | .062 | .064 | 3.581 | .979 | .100 |
| Cross-pool with explicit .03 eigen floor | .062 | .065 | 4.729 | 1.000 | .128 |

Probability available refers only to a positive variance for the declared scalar contrast. A full joint GMM still requires an identified, positive-definite moment covariance. Under cancellation the raw cross-pool matrix was indefinite in .950, .817 and .738 of samples at 70, 150 and 250 sites. At 70 sites the PSD projection still left at least one eigenvalue at zero to numerical precision in .875 of samples. The regularisation preserving matched coefficients in blocks inflated the paired variance less than diagonal shrinkage or a uniform positive eigenvalue floor. Neither shrinkage guaranteed positive definiteness with the estimated strengths. Full matrix risk changed little between the different regularisations, while paired contrast variance and probability changed appreciably. Report both and treat scalar variance estimates separately when joint scores are too weak. Even a near-unbiased cross-pool estimate can fail to resolve a small ecological difference at a given integration budget. [Known covariance summary](results/covariance_summary.csv)
