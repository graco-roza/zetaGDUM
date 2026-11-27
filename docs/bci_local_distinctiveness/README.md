# BCI demonstration of local compositional distinctiveness and retention (Figure 4)

Package: zetaGDUM at commit f7e115f (unmodified; `zeta_lczd()`, `zeta_local_retention()`, internal `conditional_zeta_cpp` used by `conditional_zeta()` for per-species local contributions). R packages: adespatial (beta.div), dplyr, tidyr, ggplot2, patchwork. Run order: `01_prepare.R`, `02_metrics.R`, `03_lcbd_tests.R`, `04_analysis.R`, `08_extra_tables.R`, `05_figure4.R`, `06_si_figures.R`, `07_si_figures2.R` (set `ZG_LIB6` to the library holding the f7e115f build; `00_api_check.R` is a smoke test on `vegan::BCI`).

## Data acquisition (not redistributed here)
| Item | Source | Use |
|---|---|---|
| Census tables `bci.tree1-8.rdata`, `bci.stem1-8.rdata`, `bci.spptable.rdata` | Condit, Pérez, Aguilar, Lao, Foster & Hubbell (2019), Dryad, doi:10.15146/5xcp-0d46 (2019 version, files `bci.tree.zip`, `bci.stem.zip`, `bci.spptable.rdata`). Dryad blocks scripted downloads; the files were downloaded manually in a browser and extracted to `raw/dryad/`. | 1985 = census 2, 2015 = census 8 |
| `habitat.txt` (x, y, habitat of each 20 x 20 m quadrat) | https://raw.githubusercontent.com/forestgeo/bciex/master/data-raw/from_suzanne_lao/bcifiles/habitat.txt (Harms et al. 2001 classification, via S. Lao) | habitat of the 1250 quadrats |
| Paper being reproduced | Legendre & Condit (2019) Forest Ecosystems 6:7, doi:10.1186/s40663-019-0164-4 | spatial resolution, LCBD method, published counts |

## Preprocessing (`01_prepare.R`)
Quadrats: 1250 cells of 20 m x 20 m (50 x 25), assigned from tree coordinates (gx, gy). Individuals: **tree records with status "A" (alive) and dbh >= 10 mm, one record per tree** (`bci.tree*` tables). Species by quadrat counts form the abundance matrix; 1985 has 242,062 individuals, 307 species; 2015 has 207,719 individuals, 297 species; no empty quadrat.

Verification against Legendre & Condit (2019): this rule reproduces exactly (i) 244,015 individuals in 1990 (their Results), (ii) a net swamp gain of 1,354 individuals 1985-2015 (2,798 to 4,152), (iii) a non-swamp loss of 35,697 individuals. The stem tables (all stems of each tree) do **not** reproduce these totals (273,037 in 1985), so tree records, not stems, are the units used in the paper; stem-based analyses are provided as a sensitivity (`FigureS6`; Spearman rho between tree-based and stem-based LCBD 0.98 in 1985, 0.97 in 2015; occurrence results identical). The union of species in 1985 and 2015 is 322 here versus 328 reported (a small discrepancy, probably in the treatment of morphospecies).

Habitat table: 1250 unique quadrats, seven categories (high plateau 170, low plateau 620, slope 284, stream 32, young 48, mixed 66, swamp 30). The 30 swamp quadrats (x 240-460 m, y 200-300 m, the central area) match the swamp (group 4, 30 quadrats) of Legendre & Condit. Their six groups cannot be fully reconciled with these seven categories: they report 500 quadrats in group 1, and no combination of categories sums to 500. Analyses therefore use the seven categories as supplied, with swamp versus the rest as the main contrast.

## Metrics (`02_metrics.R`)
* Conventional LCBD: Hellinger transformation, LCBD = row sum of squares of the centred table / total sum of squares (Legendre & De Cáceres 2013); own computation equals `adespatial::beta.div(method = "hellinger")` to machine precision (max abs difference 0; BDtotal 0.3538 in 1985, 0.3732 in 2015). Incidence LCBD uses the Sørensen coefficient.
* LCZD, exact analytical (`zeta_lczd`, relative abundance and unnormalised incidence) at orders 2, 3, 4, 5, 8, 12, 20, 40, 80, 160, 320, 640.
* Retention (`zeta_local_retention(max_order = NULL)`, all 1250 orders, occurrence and relative abundance): **T** (characteristic retention order) is the sum over orders of the conditional zeta profile; **T_plus** = 1 + (T - 1) / zeta_2|i (conditional retention order; retained order count conditional on sharing at order 2). `T == T_full` exactly (max abs difference 0) for both representations in both censuses.
* Species contributions to local distinctiveness (`species_local()` in `bci_lib.R`): for site i and order k, c_is = q_is - zeta_{s|i}(k), computed from the compiled per-species conditional zeta; sums to the exact expected dissimilarity of `zeta_lczd()` (max abs error 1e-15). Split into an occupancy part O_is = q_is (1 - P_s) (P_s = probability that the other k - 1 sites all contain species s) and an abundance-heterogeneity part H_is = c_is - O_is. For incidence the contribution is the expected set-mean occurrence minus the shared occurrence (all occupancy).
* Permutation tests of LCBD (`03_lcbd_tests.R`): 99,999 permutations, Holm correction across the 1250 quadrats (needed because the minimum attainable p must be below 0.05 / 1250).

Large intermediate objects (`out/*.rds`, raw data) are not committed; all tables used in the figures are in `tables/` and `out/site_metrics.csv`.
