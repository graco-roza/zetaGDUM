# Run from the repository root, or pass the experimental directory to st_load().
st_load <- function(path = "docs/experimental/spatiotemporal", envir = parent.frame()) {
  for (name in c("metrics", "model", "diagnostics", "simulate", "plot")) {
    sys.source(file.path(path, "R", paste0(name, ".R")), envir = envir)
  }
  invisible(path)
}
