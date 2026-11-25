# Rscript docs/manuscript_validation/support.R OUT
source("docs/manuscript_validation/common.R")
out <- commandArgs(TRUE)[1L]
meta <- read_manifest()
rows <- list()
for (i in seq_len(nrow(meta))) {
  a <- load_dd2(meta[i, ], Sys.getenv("DD2_DATA"))
  P <- (a$x > 0) * 1L
  for (alpha in c(0, .5)) for (k in as.integer(meta[i, c("low", "mid", "high")])) {
    for (M in c(10000L, 50000L)) {
      s <- draw(P, k, M, alpha, seed = 530000L + k)
      d <- suppressWarnings(ugdm_design(a$env, s$idx))
      b <- contrast_basis(a$env, d)
      spread_spans <- matrix(apply(d$Xdir, 2L, function(z) diff(range(z))), ncol = 3L, byrow = TRUE)
      spread <- apply(b[, 1:3, drop = FALSE] <= spread_spans + 1e-10, 1L, all)
      position <- b[, 4L] <= apply(d$Xpos, 2L, function(z) diff(range(z))) + 1e-10
      rows[[length(rows) + 1L]] <- data.frame(case = meta$id[i],
        sampler = if (alpha == 0) "conventional" else "species",
        order = k, budget = M, term = c(paste0("spread.", d$vars), paste0("position.", d$vars)),
        supported = unname(c(spread, position)))
    }
  }
}
write_study(do.call(rbind, rows), out, "contrast_support")
