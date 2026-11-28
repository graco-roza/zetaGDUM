#!/bin/sh
# run_all.sh: redesigned Figure 6 (run from final/; ZG_LIB6 = library built from branch principal/additive-msgdum)
Rscript revision_fig6/10_all_sites.R
NCORES=9 Rscript revision_fig6/20_splits.R
Rscript revision_fig6/30_tables.R
Rscript revision_fig6/40_figures.R
# Figure 5 with matched penalties and retention diagnostics (revision_random)
for p in shared full; do POOL=$p Rscript revision_random/71m_fit_random_matched.R; done
Rscript revision_random/72m_tables_random_matched.R; Rscript revision_random/73m_fig5_random_figure_matched.R; Rscript revision_random/95_order_diagnostics.R; Rscript revision_random/96_order_figure.R
