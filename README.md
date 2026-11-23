# zetaGDUM

Zeta diversity with abundance, exact species weights at any order, the retention-rate breakpoint,
turnover and nestedness components for sets of sites, and a multi-site extension of generalised
dissimilarity uniqueness models (GDUM) with paired-bootstrap importance and Shapley decomposition of deviance explained.

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
| `dominance_matrix()` | locally dominant species (cumulative abundance share `p`) |
| `sample_sets()`, `zeta_sets()` | sets of sites and their turnover, total, nestedness (and abundance) components |
| `pair_partition()` | Baselga and Podani-Legendre partitions for pairs |
| `ugdm_design()`, `ugdm_fit()` | multi-site GDUM (spread = directional, position = non-directional, optional distance) |
| `ugdm_importance()`, `ugdm_shapley()` | variance-partition importance and Shapley contributions to D2 |
| `site_weights()`, `zeta_orders()`, `compare_orders()` | paired Bayesian bootstrap across orders |
| `simulate_metacommunity()` | niche-based simulation of a metacommunity |
