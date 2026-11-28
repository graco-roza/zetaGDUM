# zetaGDUM (development)

* Explicit spread-by-position interactions are withdrawn from the public interface. The experimental implementation (design, fitting, diagnostics, contrasts, penalty selection with interaction shrinkage, contextual spread, and the `interaction` component model) is preserved unchanged on the `experimental/interactions` branch.
* The additive position models keep their numerical behaviour under new names: `ugdm_position_design()` (was `ugdm_interaction_design(spread = "constant")`), `ugdm_fit_position()` (was `ugdm_fit_interactions()` without interaction shrinkage), `ugdm_position_contrast()`, `ugdm_position_diagnostics()` (spread and position deviance allocations only) and `ugdm_select_penalties()`.
* `ugdm_component_design()` and `ugdm_validate_components()` accept `spread_only` and `additive` models only; `ugdm_fit_components()` no longer has `interaction_shrink`; `ugdm_component_diagnostics()` no longer returns interaction residuals.
* `ugdm_site_partitions()` is retained (moved to `validate.R`).
