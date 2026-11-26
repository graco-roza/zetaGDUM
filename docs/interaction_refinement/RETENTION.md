# Exact local retention

## Definition and computation

Fix a focal site \(i\) among \(n\) observed sites. Couple the orders by a uniform permutation of the other \(n-1\) sites, retaining its first \(k-1\) elements at order \(k\). A species or abundance threshold occupied at \(m\) sites, including the focal site, then survives through order \(k\) with probability
\[
h_m(k)=\frac{\binom{m-1}{k-1}}{\binom{n-1}{k-1}}.
\]
The recurrence \(h_m(1)=1\) and
\(h_m(k+1)=h_m(k)\max(0,m-k)/(n-k)\) evaluates the full profiles without large binomial coefficients. Writing the ratio as \(\binom{n-k}{n-m}/\binom{n-1}{n-m}\) and summing gives
\[
\sum_{k=1}^n h_m(k) = \frac{\binom{n}{n-m+1}}{\binom{n-1}{n-m}} = \frac{n}{n-m+1}.
\]

For occurrence, a species locally present has weight \(1/r_i\), where \(r_i\) is focal-site richness. For relative abundance, species \(s\)'s distinct positive thresholds are \(t_{s\ell}\) and threshold increments are \(\Delta t_{s\ell}\). If the focal site belongs to \(O_{s\ell}=\{j:q_{js}\geq t_{s\ell}\}\), its threshold contributes \(\Delta t_{s\ell} h_{|O_{s\ell}|}(k)\). The threshold decomposition follows from telescoping the minimum abundance of each species across a set. Both profile types are weighted survival curves. The weights have unit total at a site of positive richness. Occurrence uses a fixed focal denominator, so the retention index corresponds to that focal estimand. Relative-abundance retention also uses the same shared composition as the exact relative-abundance LCZD.

\(T_i(K)=\sum_{k=1}^K\zeta_{k|i}\) is the expected number of retained orders through \(K\). The package evaluates full \(T_i(n)\) using the closed form and also returns \(T_i(K)\) at any declared truncation. The conditional index
\[
T_i^+(K)=1+\sum_{k=2}^K\zeta_{k|i}/\zeta_{2|i}
\]
is undefined if order two is unavailable or has zero survival. Computing the tail sum directly preserves the index when \(\zeta_{2|i}\) is too small to register in \(T_i(K)-1\) at machine precision. A profile returns the full retention-loss distribution through \(K\), with mass \(\zeta_{k|i}-\zeta_{k+1|i}\) for first loss at order \(k+1\) and a separate right-censored mass \(\zeta_{K|i}\). A value near \(n\) with appreciable censored mass describes persistence across the observed sites.

## Numerical validation

Run the following after installing the development version.

```sh
Rscript docs/interaction_refinement/retention_validation.R
```

Forty random communities of three to eight sites were enumerated for occurrence and relative abundance, including repeated abundance values. The worst absolute difference between exact enumeration and the conditional profile was \(2.22\times 10^{-16}\). The largest closed-form full-index error was \(1.78\times 10^{-15}\), and the largest order-three truncated-index error was \(8.88\times 10^{-16}\). The retention-loss distribution summed to one within \(2.22\times 10^{-16}\). Further package checks exercise the empty-site occurrence case, undefined conditional indices, and right-censoring when every station shares all of its composition. [Replicate-level exact checks](results/retention_exact.csv)

## Relation to environmental components

Twelve 70-site metacommunities per scenario used Gaussian species niches along a spatially autocorrelated gradient. The occupancy-saturated scenario added positive baseline counts to every species-site cell before standardisation. For each simulated community, the original GDUM used 12,000 order-eight sets with relative-abundance dissimilarity and constant working weights. Site-level fitted components used 50,000 independently sampled validation sets on the same sites, holding the training basis fixed. The site-level fit uses the original additive spread and position Shapley allocation on the dissimilarity scale. All 24 fits converged.

| Scenario | Median occurrence \(T_i\) | Median relative \(T_i\) | Relative right-censored mass | Relative \(T_i\) vs exact uniqueness | Relative \(T_i\) vs fitted uniqueness | Relative \(T_i\) vs spread / position allocation |
|---|---:|---:|---:|---:|---:|---:|
| Niche gradient | 6.68 | 2.68 | .011 | −.971 | −.927 | −.757 / −.095 |
| Occupancy saturated | 70.00 | 22.17 | .277 | −.965 | −.630 | −.607 / −.286 |

The reported associations are medians of within-community Spearman correlations. Occurrence retention reached the maximum at every site under saturated occupancy, while relative retention was lower and still distinguished stations. The exact relative uniqueness and fitted environmental components describe why sites with long retention may differ in dissimilarity at the fitted order. Association with observed relative uniqueness partly reflects the mathematical overlap with the order-eight point on the retention curve. Association with fitted uniqueness was weaker under saturated occupancy. These simulations use the same sites for retention and regression, and the correlations are descriptive summaries. [Replicate-level environmental results](results/retention_environment.csv)

## Repeated sampling occasions

For comparable communities, intersect station IDs before computing profiles separately by occasion. To keep observed station sets when visits differ, select a common maximum order no larger than the smallest occasion. This matches the range of orders, and the resulting estimates still depend on the composition of the observed station set at each occasion. Full-order \(T_i(n)\) is returned with individual occasion fits. Seasonal averages use the declared comparable \(T_i(K)\), conditional \(T_i^+(K)\), the full profile itself and the loss distribution, with separate counts of available estimates for each. An unvisited station is absent from that occasion's matrix. Actual zero-richness observations leave occurrence retention undefined and are invalid for relative abundance.

A code example removes 6, 8 and 14 observations from three 70-station occasions, then compares two early occasions and one late occasion through order eight. The common-station analysis retained 44 stations on each occasion. The truncated analysis retained 69 distinct early stations (57 observed twice) and 56 late stations. A seasonal SD requires at least two observed values, so the single late occasion has missing SDs throughout. The example illustrates missing-station mechanics using independently simulated occasions with a common station index. [Seasonal output](results/retention_temporal.csv)
