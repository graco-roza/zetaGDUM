# Phase II methods

Phase I remains available through `load.R` and commit `7ac8fd6`. Phase II uses `phase2_load.R`. Both are experimental modules within the repository.

## Rectangular zeta

Let $i,t,s,v$ index sites, chronological surveys, species and environmental predictors. $I,J$ are independently identifiable sets of sites and surveys, with $k=|I|$ and $m=|J|$. A balanced panel uses the same sites and species pool, comparable recording effort and perfect detection. For observed abundance $A_{its}$, incidence is $X_{its}=1(A_{its}>0)$ and relative abundance $Q_{its}=A_{its}/\sum_s A_{its}$. Local dominance $D_{its}$ is the native `dominance_matrix(A, p=0.8)` indicator computed on each original assemblage. It follows the native ordering of ties. [Hui and McGeoch (2014)](https://doi.org/10.1086/678125) · [Latombe et al. (2017)](https://doi.org/10.1111/2041-210X.12756) · [McGeoch et al. (2019)](https://doi.org/10.1002/ecy.2832)

For a rectangle, the presence sum uses $H_s=\prod_{i\in I,t\in J} X_{its}$, abundance uses $M_s=\min_{i\in I,t\in J} Q_{its}$, and dominance uses the product of $D_{its}$. Average over equally weighted site sets and either all survey sets or all windows of consecutive surveys. The output reports the range and mean of chronological spans at each order. Equal numbers of surveys can cover different spans in irregular designs. Retention divides by the previous order along the selected axis when that denominator is positive.

For incidence and dominance, first take the temporal product at each site. If $o_s(J)$ sites retain a species, the exact spatial expectation is $\binom{o_s(J)}{k}/\binom{n}{k}$. For abundance, take the minimum of original relative abundances within each site over $J$, and sort these minima across sites. The expected minimum is
$$
a_s(k\mid J)=\sum_{j=1}^{n-k+1} B_{(j)s}(J)\exp\{\operatorname{lchoose}(n-j,k-1)-\operatorname{lchoose}(n,k)\}.
$$
The weights count selections having a given smallest rank, including when several ranks tie. Temporal minima receive no second normalization. Average $a_s(k\mid J)$ and $P_s(k\mid J)$ over the same temporal sets. Then $\mu_s=a_s/P_s$ measures conditional minimum abundance for a shared species. For zero sharing $\mu_s$ is undefined and the abundance contribution is zero when relative abundance is defined. Normalise each contribution vector by its sum to obtain species weights, using zero weights when the sum is zero. The default errors when any included assemblage has total abundance zero, like the native abundance curve. The explicit `empty="undefined"` option gives an undefined abundance mean and retains the incidence and dominance outputs on the same panel.

## Directional and occupancy accounting

For fixed $I$ and endpoints $t<u$, use endpoint products $H_0,H_1$. The species ledger records $P=H_0H_1$, $G=(1-H_0)H_1$, $L=H_0(1-H_1)$. Summing gives $Z_0=P+L$, $Z_1=P+G$, $\Delta Z=G-L$. Divide gains and losses by $2P+G+L$ for directional Sørensen components, with undefined components when the denominator is zero. The same state ledger applies to dominance. Each dominance event is split according to whether species were shared by incidence at both endpoints. [Legendre (2019)](https://doi.org/10.1002/ece3.4984)

For abundance use endpoint relative minima $M_0,M_1$. The persistent amount is $\min(M_0,M_1)$, gain is $[M_1-M_0]_+$ and loss is $[M_0-M_1]_+$. Occurrence gain and loss are $G M_1$ and $L M_0$. Numerical gain and loss multiply the signed positive and negative parts by $P$. Their net sum gives the change in shared relative abundance. These are accounting quantities conditional on the observed communities.

For each species, $r_1=r_0+g-l$, with $g,l$ counting sites gained and lost. Expansion is $[g-l]_+$, contraction is $[l-g]_+$, and balanced redistribution is $\min(g,l)$. Therefore $g+l$ equals expansion plus contraction plus twice balanced redistribution. Pure redistribution requires $g=l>0$. This describes sampled locations and does not identify movements between them. Regional appearance/disappearance is determined on the complete panel. For complete spatial combination sampling, before, after and persistent contribution are $q_k(o_0),q_k(o_1),q_k(h)$ with $q_k(o)=\binom{o}{k}/\binom{n}{k}$. Subtract $q_k(h)$ from the endpoint contribution to recover species gains and losses of shared occupancy.
