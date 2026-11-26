source("docs/phase4/common.R")
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
meta <- read_manifest()
set.seed(8922)
rows <- list()
for (i in seq_len(nrow(meta))) {
  a <- meta[i, ]
  d <- load_dd2(a, args[1L])
  p <- phase4_input(d$x, d$env, d$dist)
  spec <- internal(".limit_effect_spec")(d$env, p$basis, names(d$env)[1L])
  for (k in sort(unique(c(a$low, a$mid, a$high, ceiling(.20 * a$sites))))) {
    sets <- internal(".limit_draw")(p, k, 10000, "conventional", .5, 50000L, 10000L)
    for (component in c("total", "relative")) {
      z <- internal(".limit_fit")(p, sets, component, "constant", FALSE, spec, TRUE)
      rows[[length(rows) + 1L]] <- cbind(data.frame(community = a$community, n = a$sites,
        order = k, component, environmental_variables = ncol(d$env), M = 10000),
        as.data.frame(as.list(z$diagnostics)))
    }
  }
}
phase4_write(do.call(rbind, rows), "screen", "full")
