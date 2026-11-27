#!/bin/sh
# 40_fig4_copy.sh: Figure 4 (local community analysis, BCI) is the latest completed analysis in ../bci (see ../bci/README.md); figures and tables are copied unchanged under final names.
for e in pdf svg png; do cp ../bci/figures/Figure4_main.$e figures/Figure4_main.$e; done
for f in S1_incidence_vs_abundance S2_conditional_retention S3_census_years S4_species_contributions S5_matched_LCBD_contrasting_retention S6_trees_vs_stems; do n=$(echo $f | sed 's/^S\([0-9]\)_/FigureS4_\1_/'); for e in pdf svg png; do cp ../bci/figures/Figure$f.$e figures/$n.$e; done; done
mkdir -p tables/fig4; cp ../bci/tables/*.csv tables/fig4/; cp ../bci/out/site_metrics.csv ../bci/out/verification_counts.csv ../bci/out/metric_checks.csv ../bci/out/lcbd_tests.csv ../bci/out/lcbd_tests_occ_1985_tree.csv ../bci/out/lcbd_tests_occ_2015_tree.csv ../bci/out/lczd_normalised_1985_tree.csv ../bci/out/lczd_normalised_2015_tree.csv tables/fig4/
