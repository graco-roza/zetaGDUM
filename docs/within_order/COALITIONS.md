# Within-order coalition attribution

For a coalition of fitted terms (C), define (v(C)=1-D_C/D_{arnothing}), where (D_C) is the fractional-binomial deviance on a common assessment set pool. Coefficients and the intercept-only reference are refitted on the training pool using the same design coordinates, monotonicity constraints and principal submatrix of the training penalty. The intercept-only value is zero by definition. If its assessment deviance is zero, the allocations are undefined. The held-out game allows negative values and contributions. Relative allocations are returned only when (v(mathrm{full})) exceeds the declared positive threshold.

Let (G) be the number of predictor groups, (g) the group containing cell (j), and (m_g) the number of cells in that group. A predictor group contains its spread and position cells. Distance is a single-cell group. For (Asubseteq [G]setminus{g}) and (Bsubseteq gsetminus{j}), let (C(A,B)) consist of all cells in (A) and those in (B). The reported cell allocation is the Owen value

$$
phi_j = sum_{A,B}
rac{v(C(A,B)cup{j})-v(C(A,B))}
{Ginom{G-1}{|A|},m_ginom{m_g-1}{|B|}}.
$$

Summing cells within a predictor gives the Shapley value of that predictor in the grouped game. Summing cells by their spread, position or geographic category gives coherent component totals under the same hierarchy. This grouping assigns shared predictive information according to the declared coalition game. It cannot by itself distinguish a predictor carrying correlated information from an independent ecological mechanism.

At a small group count, subset weights evaluate the exact allocation without enumerating factorial predictor orderings. At a larger group count, randomly order the groups and independently randomise cell order within each group. Adjacent marginal contributions telescope in every ordering, preserving additivity. The conditional Monte Carlo standard error is the sample standard deviation of the contributions divided by the square root of the number of sampled orderings. For aggregated values, the standard error uses the sum within each ordering. The resampling study concerns additional uncertainty from finite observed sites and combination integration.

A shared cache fits each required coalition at most once during an allocation, apart from one retry after a failed warm start. The starting vector reuses converged shared terms from a nested fit where available. A failed coalition stops by default. Under warning mode it makes affected values missing and records its columns and reason.

The reproducible held-out benchmark uses three and five predictors, 48 training sites, 22 disjoint assessment sites, order four, 500 sampled sets per pool and a fixed training basis. Five independently seeded approximations run at 40, 160 and 640 orderings. [Replicate-level comparisons](results/coalition_benchmark.csv)

```sh
Rscript docs/within_order/coalition_benchmark.R
```

| Predictors | Allocation | Cell allocation RMSE against exact | Mean reported Monte Carlo SE | Median coalition fits | Median elapsed seconds |
|---:|---|---:|---:|---:|---:|
| 3 | Exact | 0 | 0 | 32 | 0.19 |
| 3 | 40 orderings | .0154 | .0117 | 32 | 0.17 |
| 3 | 160 orderings | .0065 | .0059 | 32 | 0.18 |
| 3 | 640 orderings | .0033 | .0029 | 32 | 0.50 |
| 5 | Exact | 0 | 0 | 192 | 0.87 |
| 5 | 40 orderings | .0112 | .0063 | 131 | 0.57 |
| 5 | 160 orderings | .0065 | .0031 | 185 | 1.19 |
| 5 | 640 orderings | .0026 | .0016 | 192 | 1.80 |

The approximation error fell as the number of orderings increased. By 640 orderings it had visited all exact coalitions in these examples. At that point exact subset attribution required fewer seconds. The five replicates quantify conditional permutation error for these fixed fitted datasets.
