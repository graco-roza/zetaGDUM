# zetaGDUM

Zeta diversity with relative abundance, exact species contributions to shared composition and
total dissimilarity, the retention-rate breakpoint, a shared-species partition for sets of sites,
and a multi-site extension of generalised
dissimilarity uniqueness models (GDUM) with paired-bootstrap importance and Shapley decomposition of deviance explained.

For sets of more than two sites, the partition uses species shared by every site; it does not
capture all replacement among subsets. The simulator includes environmental niches, optional
geographic range windows and site-level thinning, but no dispersal dynamics or competition.

```r
# install.packages("remotes"); remotes::install_github("graco-roza/zetaGDUM")
library(zetaGDUM)

sim <- simulate_metacommunity(100, 100, seed = 1)

# 1. zeta curve, retention rate and breakpoint (presence, abundance, dominance)
cur <- zeta_curve(sim$comm); zeta_break(cur)
zeta_curve(sim$comm, type = "abundance")
zeta_curve(sim$comm, type = "dominance", p = 0.8)

# 2. species weights at order k (who is represented)
zeta_weights(sim$comm, k = 5)
zeta_weights(sim$comm, k = 5, abundance = TRUE)

# 3. components for sets of sites
idx <- sample_sets(100, 4, M = 2000, seed = 1)
head(zeta_sets(sim$comm, idx))
head(pair_partition(sim$comm))          # Baselga and Podani-Legendre partitions for pairs

# Species contributions sum exactly to total set dissimilarity.
inc <- species_contributions(sim$comm, idx, "presence")
rel <- species_contributions(sim$comm, idx, "relative")
head(inc$response)
head(rel$response)

# Raw-count and relative-abundance responses address different questions.
head(ugdm_response(sim$comm, idx, "bray_curtis")$y)
head(ugdm_response(sim$comm, idx, "relative")$y)

# 4. multi-site GDUM across orders, paired bootstrap, Shapley groups
res <- zeta_orders(sim$comm, sim$env, sim$coords, M = 5000, ndraw = 50,
                   groups = list(env1 = "env1", env2 = "env2", distance = "distance"))
compare_orders(res, "importance")
compare_orders(res, "shapley")
```

Main functions

| Function | Purpose |
|---|---|
| `zeta_curve()`, `zeta_break()` | exact zeta curve, retention rate, 97% breakpoint |
| `zeta_weights()` | exact share of each species at order `k` (presence or abundance) |
| `species_contributions()` | exact species contributions to total incidence, relative-abundance or raw-count dissimilarity |
| `dominance_matrix()` | locally dominant species (cumulative abundance share `p`) |
| `sample_sets()`, `zeta_sets()` | sets of sites and their turnover, total, nestedness (and abundance) components |
| `pair_partition()` | Baselga and Podani-Legendre partitions for pairs |
| `ugdm_design()`, `ugdm_fit()` | multi-site GDUM (spread = directional, position = non-directional, optional distance) |
| `ugdm_importance()`, `ugdm_shapley()` | variance-partition importance and Shapley contributions to D2 |
| `site_weights()`, `zeta_orders()`, `compare_orders()` | paired Bayesian bootstrap across orders |
| `simulate_metacommunity()` | niche-based simulation of a metacommunity |

Relative-abundance zeta requires non-empty sites. Raw-count abundance responses also depend on
total abundance and therefore require comparable sampling effort across sites. The GDUM fits
nonlinear set responses; similarity weights do not partition its fitted coefficients.

### Species-informed site-set sampling

Set `sampling = "species"` in `zeta_orders()` or `ugdm_validate()` to draw
more sets with shared species. The default mixture uses the targeted branch
half the time. Sampling weights correct the fit and summaries back to the
conventional target, including zero-sharing sets. The default
`sampling = "conventional"` retains conventional sampling. For unit working
weights on incidence or relative-abundance dissimilarity, set
`working_weights = "constant"`.

```r
fit <- zeta_orders(sim$comm, sim$env, orders = c(2, 4, 8), M = 2000,
                   sampling = "species", alpha = 0.5, working_weights = "constant")
fit$point_diagnostics
expected_species_contributions(sim$comm, k = 4, type = "relative")
```

See the [methods specification](docs/METHODS_SPEC.md) and
[validation report](docs/VALIDATION_REPORT.md).
