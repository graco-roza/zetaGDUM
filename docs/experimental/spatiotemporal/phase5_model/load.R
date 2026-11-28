st5_load <- function(path = "docs/experimental/spatiotemporal", envir = parent.frame(),
                     cache = tempdir()) {
  sys.source(file.path(path, "phase3_load.R"), envir)
  envir$st3_load(path, envir)
  for (name in c("support", "responses", "design", "fit", "validation", "simulate"))
    sys.source(file.path(path, "phase5_model", "R", paste0(name, ".R")), envir)
  for (kernel in c("phase4b_annual_metacommunities/scripts/incidence.cpp",
                   "phase4_swedish_lakes/scripts/rectangle_kernel.cpp"))
    Rcpp::sourceCpp(file.path(path, kernel), cacheDir = cache, env = envir)
  invisible(path)
}
