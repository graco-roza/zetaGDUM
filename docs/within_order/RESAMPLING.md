# Within-order resampling of coalition importance

For each coalition \(C\), let \(X_C\) contain the intercept and eligible basis columns. Estimate \(\widehat\theta_C\) by minimising the average working loss
\[
\ell_C(\theta)=\sum_m w_m\{\log(1+\exp(X_{mC}\theta))-y_m X_{mC}\theta\}
 +\theta^T P_C\theta/2,
\]
where the effective training weights sum to one and \(P_C\) is the relevant principal submatrix of the fitted penalty. Directional coefficients satisfy nonnegativity bounds. Evaluate all coalition and intercept-only predictions on identical evaluation rows within a draw.

## Direct resampling

In the Bayesian bootstrap draw one vector \(w_i\) of Dirichlet weights for the sites. Conditional Poisson sampling draws sets with probability proportional to \(\prod_{i\in S}w_i\). Every coalition shares these sampled sets, the ecological response, working weights and evaluation rows. With no separate evaluation design, recompute the response and explained deviance on the resampled training rows. With an explicit evaluation design, keep that design fixed across draws. The ordinary site bootstrap resamples site IDs with replacement and subsequently samples sets of distinct site copies. Distinct copies may refer to the same original site. Both methods require a data builder with a fixed training basis.

## Score resampling

At each fitted coalition's coefficient vector, compute the full and delete-one scores. Use the existing native kernel to tally the score summed over sets containing each site. Subtract those tallies from the full score, divide by the retained working-weight mass, and centre the resulting score pseudovalues. An intercept-only fit receives the same treatment. Average pseudovalues over independent pools before estimating covariance. This pools their numerical errors while retaining each site's ecological score. A user-supplied site-score function can use independently piloted global control variates. It must return score pseudovalues in the units of the average working loss, with the fixed penalty excluded.

For site or block \(b\), let \(J_{bC}\) be the sum of centred site scores in the group. Use the same iid normal or Rademacher multipliers \(\xi_b\) across coalitions, with
\[
s_C=\sqrt{G/(G-1)}\,n^{-1}\sum_{b=1}^{G}\xi_b J_{bC}.
\]
At an interior fit with square estimating equations, the normal update has the flat-prior Gaussian GMM sandwich covariance. This representation preserves cross-coalition covariance without inverting the stacked covariance, which can be singular through shared predictors and sites. Rademacher multipliers have the same second moments and a different finite-sample distribution. Spatial block versions require defensibly independent groups of sites and use one multiplier or Bayesian-bootstrap weight per block. Unequal block sizes change the finite-sample behaviour.

For monotone coefficients, include the gradient \(g_C=\nabla\ell_C(\widehat\theta_C)\) in addition to the curvature \(H_C=\nabla^2\ell_C(\widehat\theta_C)\). Approximate a draw by
\[
\arg\min_{\theta_C\in\mathcal K_C}
\left[g_C^T(\theta_C-\widehat\theta_C)
+\tfrac12(\theta_C-\widehat\theta_C)^TH_C(\theta_C-\widehat\theta_C)-s_C^T\theta_C\right].
\]
The baseline gradient keeps the KKT information at coefficients on a constraint boundary. Refit mode minimises the original perturbed objective \(\ell_C(\theta)-s_C^T\theta\). Hybrid mode checks the projected score residual in the optimiser's scaled coordinates before deciding whether to refit. Gaussian GMM and its projection are approximations at active constraints. A linear score perturbation can leave an exact objective unbounded if a sufficiently large perturbation overwhelms the logistic curvature. Such draws are reported as failures when the optimiser fails.

Score methods evaluate predictions on fixed rows. An independent evaluation design defines a conditional predictive importance target. Reusing training rows gives fixed-design assessment and does not replicate the variation of a new in-sample empirical distribution. Exact resampling and score approximations can therefore differ when no separate evaluation is supplied.

## Ranking summaries

From each draw, calculate predictor Shapley importance, aggregate spread and position Owen allocations, and all intra-order predictor contrasts. Split tied leading ranks equally. Report strictly positive superiority and differences above the declared practical threshold, with 80% and 95% central intervals. Superiority probabilities are conditional on finite returned draws. Retain draw-level failures and report their frequency whenever summarising an analysis. Relative importance is available only above the declared positive total explained-deviance threshold.
