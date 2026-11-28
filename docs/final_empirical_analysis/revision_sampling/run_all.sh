#!/bin/sh
# run_all.sh: sampling-requirement simulation (run from final/; ZG_LIB6 = library built from branch principal/additive-msgdum).
# Single-threaded BLAS is essential (parallel workers otherwise oversubscribe the machine).
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SIM=known NCORES=5 Rscript revision_sampling/20_run.R      # known-parameter experiment (945 unique units, resumable: unit results in revision_sampling/logs/units_known)
SIM=meta  NCORES=5 Rscript revision_sampling/20_run.R      # site-level metacommunity experiment (945 unique units)
Rscript revision_sampling/25_test_null.R                   # variance of the held-out response (reproduced from seeds)
Rscript revision_sampling/26_binding.R                     # diagnostic: which support criterion limits availability
Rscript revision_sampling/30_tables.R
Rscript revision_sampling/40_figures.R
