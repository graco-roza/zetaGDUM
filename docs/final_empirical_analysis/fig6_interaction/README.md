# Figure 6 interaction analysis across all informative orders (diagnostic; Figure 6 not finalised)

zetaGDUM f7e115f (unmodified, built into a separate library). Primary predictors log10 TP, TN, Chl-a + geographic distance; coherent component model; shared 1,077-lake pool principal, full pools as sensitivity.

| script | purpose |
|---|---|
| `60_cv_interactions.R` | additive vs explicit interaction (all predictors, or TP only) at endpoint penalties 0.01, 1e-4, 1e-6; identical orders, 5 spatial blocks, sampled combinations (3000 train/3000 test, seeds 10*seed+fold / 100*seed+fold), seeds 1:3; in-sample and spatially held-out D2 (total, occupancy-associated, heterogeneity) via `ugdm_assess_components(fit, design, response, training_response = y)` for both splits; all fitted coefficient vectors |
| `61_spline_shapes.R` | all-site fits: TP-spread spline coefficients/shapes f(z) = sum_j beta_j(u) I_j(z) at TP positions ~10/35/100 ug/L (+ endpoints u = 0, 1) for the total and heterogeneity-share stages, every order; link-scale decomposition |
| `62_tables.R` | paired gains (train/test), TP-only vs all-predictor contrasts, convergence, shape magnitude |
| `63_figures.R` | diagnostic figures (PDF/SVG/600 dpi PNG, 183 mm) |

TP-only interaction: package code has no per-predictor option, so for TN and Chl-a both endpoint columns are set to half the additive block (identical endpoint curves; the endpoint-average ridge acts as for the additive spread), identically in training and evaluation designs.
`delta_D2 = D2_interaction - D2_additive`, separately for the training and the spatially held-out split; seed means within fold, then median/min/max over folds. No fit removed (non-converged fits are flagged: 2 of 1,260 per pool).
Link decomposition: a fixed model-scale TP-spread increment f(z80) - f(z20) is applied at the fitted baseline linear predictors of the sets nearest each TP position and converted to the response scale; in the additive model position differences are purely link-scale.
