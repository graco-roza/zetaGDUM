#!/bin/sh
# 42_collect.sh: gather the Figure 6 outputs of both pools into figures/ and tables/ under final names (run after fig6_shared/run_fig6.sh and fig6_full/run_fig6.sh and 40_fig4_copy.sh)
for e in pdf svg png; do
  cp fig6_shared/figures/Figure6_main_order3.$e figures/Figure6_main.$e; cp fig6_shared/figures/Figure6_alt_order5.$e figures/FigureS6_order5.$e
  cp fig6_full/figures/Figure6_main_order3.$e figures/FigureS6_full_pool_main.$e; cp fig6_full/figures/Figure6_alt_order5.$e figures/FigureS6_full_pool_order5.$e
  for r in R1_magnitude_and_predictability R2_absolute_importance R3_relative_importance R4_magnitude_vs_importance_TP R5_context_dependence_order3 R5_context_dependence_order5 R6_added_value_position_interaction; do
    cp fig6_shared/figures_repl/FigS_$r.$e figures/FigureS6_components_$r.$e; cp fig6_full/figures_repl/FigS_$r.$e figures/FigureS6_components_full_pool_$r.$e; done
done
for p in shared full; do mkdir -p tables/fig6_$p tables/fig6_${p}_components; cp fig6_$p/tables/*.csv tables/fig6_$p/; cp fig6_$p/tables_repl/*.csv tables/fig6_${p}_components/; cp fig6_$p/out/selected_summary.csv tables/fig6_$p/ 2>/dev/null; done
