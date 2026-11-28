# Run from repository root after installing zetaGDUM.
source("docs/experimental/spatiotemporal/load.R")
st_load()
testthat::test_dir("docs/experimental/spatiotemporal/tests", stop_on_failure = TRUE)
