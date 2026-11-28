root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase2_load.R"))
st2_load(root)
testthat::test_dir(file.path(root, "tests_phase2"), stop_on_failure = TRUE)
