#!/bin/sh
set -eu
scripts=docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts
for step in 01_download 02_test 03_audit 04_annual 05_simulation 06_abundance 07_uncertainty 08_figures; do
  Rscript "$scripts/$step.R"
done
