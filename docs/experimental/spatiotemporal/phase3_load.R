st3_load <- function(path = "docs/experimental/spatiotemporal", envir = parent.frame()) {
  sys.source(file.path(path, "load.R"), envir)
  envir$st_load(path, envir)
  modules <- c("metrics_partition", "temporal_reorganisation", "process_allocation", "process_shapley", "model_process_gdum", "simulate_process_gdum")
  for (name in modules) sys.source(file.path(path, "R", paste0(name, ".R")), envir)
  invisible(path)
}
