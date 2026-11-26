# Functional penalties and multisite context

## Penalty coordinates

For the existing monotone I-spline basis \(B_b\), a fitted endpoint curve is \(f(z)=\sum_b \beta_b B_b(z)\). An interaction penalty can be imposed directly in function space,
\[
P_I = (\beta_U-\beta_L)^T K(\beta_U-\beta_L), \qquad
K_{bc} = \int_L^U B_b'(z)B_c'(z)\,dz.
\]
The derivative Gram matrix uses three-node Gaussian quadrature on each knot interval. The integrand is a polynomial of degree at most four for the package's degree-two I-spline derivatives, so this rule evaluates it exactly up to rounding. Applying the same function shift to both endpoints leaves the penalty unchanged. The existing coefficient-distance penalty is available by default. The optimiser uses the average fractional-logit loss plus one half of the weighted quadratic penalty.

The curvature matrices already supplied by the package integrate the squared second derivative. In unit-range coordinates \(t=(z-L)/(U-L)\), their multipliers are \(R^3\) for curvature and \(R\) for the first-derivative Gram, where \(R=U-L\). Hence all such quadratic forms express integrated link-scale squared changes on the same unit range. Scientific hypotheses can still use different strengths for directional and position curvature, interaction and ridge. Raw coefficient magnitude remains a separate sensitivity penalty. In the contextual model its directional part uses canonical coefficients so changing how a common spline compensates between the two endpoint blocks cannot alter the objective.

The optional slope shrinkage uses the orthogonal projection of \(f\) onto the affine functions under uniform measure on \(t\). After removing the overall level, the amplitude of that projection is
\[
s = 12\int_0^1(t-1/2)f(L+Rt)\,dt.
\]
Penalise \(s^2\) in range units or \((s/R)^2\) in native standardised-predictor units. For a position curve with coefficients on \(z,t_*^2,t_*^3\) and \(t_*=2t-1\), the amplitude is \(R\beta_z + (6/5)\beta_{t_*^3}\). Directional slope shrinkage averages the squared projection for the two endpoints, or the four physical corners of the contextual model. Curvature alone leaves a common linear response unchanged. Shrinking the separate affine projection and curvature together can shrink the entire environmental effect towards zero, apart from its level absorbed by the intercept. As one illustration, on the unit range the coefficient penalty of a linear difference with coefficients \((1,1,1)\) is \(3\) and its derivative penalty is \(9\). A zero-endpoint profile with coefficients \((1,-2,1)\) is charged \(6\) and \(1.8\), respectively. Unit-range derivative penalties were invariant over widths 1, 3 and 7, and their native integrals matched independent numerical integration within \(3.6\times10^{-15}\). [Numerical penalty checks](results/functional_penalty_units.csv)

## Design before fitting

For a set \(S\), let \(D_b\) be the mean of \(d_{ij,b}\) and \(T_b=\overline{u_Sd_{ij,b}}\), with \(u_S\) the set mean mapped to the fixed unit range. Put \(P_b=\overline{u_{ij}d_{ij,b}}\) using the similarly mapped pair midpoint. The contextual contribution in unit-range coordinates is \(T_b-P_b=C_b/R\) whenever both interpolators are inside the training range. Outside it, both are clipped separately for monotone interpolation and the design retains the exact raw \(C_b\) for diagnostics. At order two \(u_S=u_{ij}\), so the fitted and raw contextual columns are identically zero.

A first exact alias persists above order two. When \(B\) is a linear ramp, the pair midpoint contribution \(\overline{m_{ij}|z_i-z_j|}\) is proportional to the mean difference of the corresponding monotone quadratic transform, already in the original spread span on the training range. With the package's single polynomial segment, the canonical pair-midpoint block therefore loses at least this direction. An arbitrary contextual direction scales \(d_{ij,b}\) by the difference between the set mean and pair midpoint. Under reflection-symmetric environmental configurations, if the pair difference is invariant under reflection, the contributions cancel in opposite pairs. Linear-ramp context vanishes there. Some other bases can have a nonzero contextual column in these sets while its independent information remains zero.

The constrained fit uses four nonnegative slope blocks \(L_p,U_p,L_S,U_S\) per predictor, each a coefficient vector for the spread basis. The columns are \(D-P,P,D-T,T\). For each basis, the derivative coefficient at fixed pair midpoint and set mean is an affine function on the unit square,
\[
(1-u_{ij}) L_p + u_{ij} U_p + (1-u_S)L_S + u_S U_S .
\]
It is nonnegative everywhere on the square. Conversely, an affine coefficient surface nonnegative on the whole square has a representation of this form. For intercept \(a\) and two slopes \(b,c\), allocate the minimum \(a+\min(0,b)+\min(0,c)\) between the two blocks, then adjust each lower endpoint by \(\max(0,-b)\) or \(\max(0,-c)\). Because the allocation of the minimum is arbitrary, the four endpoint blocks differ by one redundant common shift per basis. The canonical coefficient form is
\[
(L_p+L_S)D + [(U_p-L_p)+(U_S-L_S)]P + (U_S-L_S)(T-P).
\]
Penalties and predictions are invariant to the redundant shift. The fit retains the nonnegative coordinates and diagnoses rank in the canonical form. Weighted SVD residualisation compares each interaction block with the preceding main-effect span. Fitting to these residual columns with the usual lower bounds would change the allowed functions and is avoided. Interpret joint fitted contrasts when individual components lack support.

At 70 sites and 2,500 combinations per order, uniform and skewed environmental draws each gave ranks 7 (additive), 9 (pair midpoint) and 9 (with context) at order two, using one predictor and nonlinear position. At orders six and twelve the first two ranks were unchanged and the contextual rank was 12. This still leaves an exact alias after removing the endpoint redundancy, and the largest residual fraction of any context column decreased from .430 to .212 under a uniform environment, and from .308 to .119 under a skewed environment. Restricting every set to symmetric pairs reduced the rank sequence to 4, 5, 5 at orders six and twelve. There the raw linear-ramp contextual RMS was below \(1.1\times10^{-16}\). The conditional context can therefore reveal new directions for supported environmental configurations. Order alone cannot establish its estimability. [Design audit](results/context_geometry.csv)

## Site-block selection

Repeated partitions accept independently permuted site blocks or spatially coherent bands in projected coordinates. Every candidate uses the same training and evaluation combinations within a fold and the same training-derived basis. The fit travels from stronger to weaker combinations of penalties by warm starts, retrying nonconverged fits from scratch. Folds receive equal weight within a partition and partitions receive equal weight. The selector can use predictive deviance, a held-out ecological contrast loss or the best ecological contrast loss among candidates within a predeclared deviance tolerance. A loss with a known ecological target applies to simulations. For empirical data the criterion or contrast sensitivity must be declared in advance. Fold omission documents sensitivity across repeated partitions and does not quantify sampling uncertainty after tuning.

Reproduce the geometry and penalty audit with:

```sh
Rscript docs/adaptive_refinement/context_geometry_validation.R
```
