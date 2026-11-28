#!/bin/sh
# run_all.sh: interaction simulations (run from final/ with ZG_LIB6 = library containing zetaGDUM f7e115f)
NCORES=9 Rscript revision_sim/10_simA.R
TRIALS=1000 NCORES=9 Rscript revision_sim/10_simA.R          # low-noise sensitivity
Rscript revision_sim/11a_simB_truth.R                          # species pool and Monte Carlo truth
NCORES=9 Rscript revision_sim/11b_simB_fit.R
Rscript revision_sim/12_tables.R
Rscript revision_sim/14_projection_test.R
Rscript revision_sim/13_figure.R
