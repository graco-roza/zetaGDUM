st2_load <- function(path = "docs/experimental/spatiotemporal", envir = parent.frame()) {
  sys.source(file.path(path, "load.R"), envir)
  envir$st_load(path, envir)
  modules <- c("metrics_abundance", "transitions")
  for (name in modules) sys.source(file.path(path, "R", paste0(name, ".R")), envir)
  invisible(path)
}
