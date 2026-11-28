#!/bin/sh
# run_all.sh: revised final analyses (run from final/ with ZG_LIB6 set to a library containing zetaGDUM f7e115f; needs gdm, gdmmTMB, zetadiv, RANN, ggh4x, svglite, patchwork)
Rscript revision_random/70_splits.R
for p in shared full; do POOL=$p Rscript revision_random/71_fit_random.R > revision_random/logs/71_$p.log 2>&1; done
Rscript revision_random/72_tables_random.R; Rscript revision_random/73_fig5_random_figure.R
for p in shared full; do POOL=$p Rscript revision_random/81_cv_interactions_random.R > revision_random/logs/81_$p.log 2>&1; done
Rscript revision_random/82_tables_fig6_random.R
FIG6_VERSION=additive Rscript revision_random/83_figure6_versions.R
FIG6_VERSION=k2_additive_higher_TP_interaction Rscript revision_random/83_figure6_versions.R
Rscript revision_random/90_equivalence.R; Rscript revision_random/91_equivalence_figure.R
Rscript revision_random/95_order_diagnostics.R; Rscript revision_random/96_order_figure.R
