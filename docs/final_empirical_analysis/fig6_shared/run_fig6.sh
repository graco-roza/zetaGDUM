#!/bin/sh
# Run inside this directory with ZG_LIB6 set to the zetaGDUM f7e115f library
Rscript 02_orders.R
Rscript 03_validate.R
Rscript 03b_penalty.R
Rscript 04_fit.R
Rscript 04b_fit_lowpenalty.R
Rscript 05_tables.R
Rscript 06_surface_tables.R
Rscript 07_figure6.R
Rscript repl/10_fit.R
Rscript repl/11_surfaces.R
Rscript repl/12_tables.R
Rscript repl/13_figures.R
Rscript repl/14_freesign.R
