# Species-informed sampling for multisite GDUM

## Set distributions

Let $S$ be an unordered set of $k$ distinct sites from $n$ eligible sites. Under the uniform target $P(S)=1/{n\choose k}$. For a Bayesian-bootstrap draw with site weights $\omega_i$, the target is

```math
P_\omega(S)=\frac{\prod_{i\in S}\omega_i}{e_k(\omega)},\qquad
e_k(\omega)=\sum_{|A|=k}\prod_{i\in A}\omega_i.
```

Zero-weight sites are outside the weighted support. Let $O_s$ be the sites occupied by species $s$, $o_s=|O_s|$, and $z(S)=\sum_s 1\{S\subseteq O_s\}$. The sharing probability of species $s$ under the two targets is

```math
P_s(k)=\frac{{o_s\choose k}}{{n\choose k}}
\quad\text{or}\quad
P_{\omega,s}(k)=\frac{e_k(\omega_{O_s})}{e_k(\omega)}.
```

Write $Z=E_P[z(S)]=\sum_s P(S\subseteq O_s)$. With probability $1-\alpha$, draw from $P$. With probability $\alpha$, choose $s$ with probability $P(S\subseteq O_s)/Z$, then draw from $P(\cdot\mid S\subseteq O_s)$. The latter step samples uniformly within occupied sites for the uniform target and uses conditional Poisson sampling within occupied sites for a bootstrap target. The same set may be reached through any of its shared species. Summing over those paths gives

```math
Q(S)=P(S)\left[1-\alpha+\alpha\frac{z(S)}{Z}\right],\qquad
w_{\mathrm{sample}}(S)=\frac{P(S)}{Q(S)}.
```

The default targeted fraction is $\alpha=0.5$. Values in $[0,1)$ retain the full support of $P$. Zero-sharing sets enter through the conventional branch and receive weight $1/(1-\alpha)$. If no eligible species has occupancy at least $k$, the sampler uses $P$ and assigns unit corrections. Both model pipelines use exhaustive pairs at order two when their budget covers all pairs. The public `sample_sets_species()` always returns $M$ independent draws and permits repeated sets.

## Numerical implementation

The uniform sampler uses log binomial coefficients and the weighted sampler uses log elementary symmetric sums. Starting at $L_0=0$ and $L_j=-\infty$ for $j>0$, for each positive site weight update the eligible entries in descending $j$,

```math
L_j\leftarrow\operatorname{logaddexp}(L_j,\log\omega_i+L_{j-1}).
```

The descending update adds each site at most once. The dynamic program uses $O(k)$ memory and $O(km)$ time on a pool of $m$ eligible sites. Normalising log masses uses their maximum. The correction uses the log expected sharing, including when the numeric value of that expectation underflows. Conditional Poisson uses a common shift on the log odds to bring the expected sample size to $k$, then accepts Bernoulli draws of exactly $k$ sites. Every size-$k$ outcome has probability proportional to its product of weights because the common shift cancels after conditioning. A pool of exactly $k$ eligible sites returns its only set directly. Targeted draws are batched by selected species before conditional Poisson sampling.

## Fitting and validation

`zeta_orders(sampling = "species", alpha = 0.5)` applies the mixture to the uniform point fit and each direct Bayesian-bootstrap draw. With `weighting = "product"`, it reuses the point sample and combines its uniform proposal correction with the product of that draw's site weights. Exhaustive pairs receive the product weight under either bootstrap workflow. Conventional sampling remains the default. Default order selection still uses the conventional estimated share of sets with positive sharing, so increased sharing under the proposal cannot move the order-selection threshold. Use explicit `orders` to study higher orders.

Each model loss term is weighted by the product of its response working weight, sampling correction, and any additional weight. This identity holds pointwise for the expected loss because $Q(S)w_{\mathrm{sample}}(S)=P(S)$. Finite-sample fits, deviance ratios, and normalised summaries can vary between proposals. The same effective weights enter the fitted and null deviances, importance covariances, and every Shapley subset fit. By default, incidence retains its original richness working weights and relative abundance uses units. `working_weights = "constant"` assigns both unit working weights and leaves the definition of each response intact. Raw-count working weights remain controlled by `abundance_trials` in `ugdm_response()`. Sets with undefined responses are dropped for either proposal. Constant responses return an intercept-only fit with undefined explained deviance and importance.

To use a common model basis across proposals and bootstrap draws, the pipelines scale mean pairwise geographic distance using the range of eligible site-pair distances. Validation obtains that range from the training sites and reuses it for held-out predictions, clamping values outside the training bounds. This changes the distance basis for workflows whose earlier bounds came from the realised set sample. Direct `ugdm_design()` calls can still use sample bounds, or specify `dist_range` or a reference design.

`ugdm_validate()` can draw either proposal on each side of a site split. Its held-out deviance and RMSE use the product of sampling corrections and working weights. Its squared Pearson correlation and slope through the origin use the sampling correction alone, preserving the original set-level metric definitions. `test_sampling = "conventional"` with the same `M_test` and seed evaluates alternative training proposals on the same held-out sets. The independent test seed keeps those sets fixed even when the training sampler consumes a different number of random draws.

Point and bootstrap diagnostics report the number of retained draw rows, the number sharing species, the effective mixture fraction and Kish effective sample size $(\sum_j w_j)^2/\sum_j w_j^2$. Here $w_j$ combines sampling and bootstrap weights and excludes response working weights. Repeated sets count as draws. Shared sites induce dependence, so this ESS measures weight concentration only.

## Exact species contributions under uniform site sets

Write $I_{is}=1\{x_{is}>0\}$, $p_s=o_s/n$, and $q_{is}=x_{is}/\sum_t x_{it}$ for a nonempty site. Uniform set sampling implies that a set mean has the same expectation as a site mean. Incidence uses the set-level difference between mean presence and shared presence, giving expected contribution

```math
E[\overline I_s(S)-1\{S\subseteq O_s\}]=p_s-P_s(k).
```

This is in richness units. The set-level richness-normalised contribution divides that same numerator by mean site richness in the set. Its expectation depends on the joint distribution of the numerator and denominator. It generally differs from the unnormalised expectation divided by overall mean richness. Use `species_contributions(..., type = "presence")` for the set-level terms.

For relative abundance, the expected contribution is

```math
E[\overline q_s(S)-\min_{i\in S}q_{is}]
=\operatorname{mean}_i(q_{is})-P_s(k)m_s(k).
```

Here $m_s(k)$ is the expected minimum conditional on sharing species $s$. Sort the occupied-site relative abundances in ascending order. The value at rank $j$ is the minimum when that site and $k-1$ of the $o_s-j$ later sites are selected. Its conditional probability is ${o_s-j\choose k-1}/{o_s\choose k}$, through rank $o_s-k+1$. Rank-based counting also counts tied values correctly. A species with $o_s<k$ has zero expected minimum and undefined conditional mean. The species contributions sum to expected relative-abundance dissimilarity across uniform sets. `expected_species_contributions()` implements these uniform expectations. Its relative-abundance calculation requires every site to have positive total abundance.
