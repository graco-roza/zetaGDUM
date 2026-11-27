#!/bin/sh
# run_fits.sh: all configuration fits (run from final/ with ZG_LIB6 set). Shared pool = principal; full pools and the extreme-nutrient trim = sensitivity.
for p in shared full; do for t in 0 1; do POOL=$p TRIM=$t Rscript fig5_revision/50_fit_configs.R > fig5_revision/fit_$p$t.log 2>&1; done; done
