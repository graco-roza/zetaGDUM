#!/bin/sh
# run_all.sh: all-orders interaction analysis for both pools (run from final/ with ZG_LIB6 set).
for p in shared full; do POOL=$p Rscript fig6_interaction/60_cv_interactions.R > fig6_interaction/cv_$p.log 2>&1; done
