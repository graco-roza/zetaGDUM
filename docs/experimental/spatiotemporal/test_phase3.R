root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase3_load.R"))
st3_load(root)
testthat::test_dir(file.path(root, "tests_phase3"), stop_on_failure = TRUE)
