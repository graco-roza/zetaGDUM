#!/bin/sh
# run_all.sh: decisive interaction validation (run from final/, ZG_LIB6 = library with zetaGDUM f7e115f; needs ../revision_sim/00_sim_common.R)
Rscript revision_ident/20_identifiability_A.R
NCORES=9 Rscript revision_ident/21_balanced_known.R
Rscript revision_ident/31_meta_truth.R
NCORES=9 Rscript revision_ident/32_meta_fit.R
Rscript revision_ident/40_tables.R
Rscript revision_ident/41_figure.R
