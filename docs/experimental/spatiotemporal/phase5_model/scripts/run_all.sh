#!/usr/bin/env bash
set -euo pipefail
phase=docs/experimental/spatiotemporal/phase5_model
for stage in 01_tests 02_empirical 03_simulation 04_example 05_figures 06_summaries 07_spatial_estimands 08_orders; do
  Rscript "$phase/scripts/$stage.R"
done
