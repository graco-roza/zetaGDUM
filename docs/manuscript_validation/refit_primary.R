# Rscript docs/manuscript_validation/refit_primary.R OUT CASE
# Recheck near-zero breakpoint-minus-pair intervals at the larger MC budget.
source("docs/manuscript_validation/common.R")
args <- commandArgs(TRUE)
out <- args[1L]
case <- args[2L]
meta <- read_manifest()
meta <- meta[match(case, meta$id), , drop = FALSE]
a <- load_dd2(meta, Sys.getenv("DD2_DATA"))
P <- (a$x > 0) * 1L
n <- nrow(P)
range <- internal(".pair_distance_range")(a$dist)
rows <- list()
low <- NULL
for (alpha in c(0, .5)) {
  fit_order <- function(k) {
    message(case, " 50000 ", alpha, " ", k)
    s <- draw(P, k, 50000, alpha, seed = 675000L + k)
    d <- suppressWarnings(ugdm_design(a$env, s$idx, a$dist, dist_range = range))
    y <- cbind(total = ugdm_response(a$x, s$idx, "total", working_weights = "constant")$y,
               relative = ugdm_response(a$x, s$idx, "relative", working_weights = "constant")$y)
    om <- omission(s$idx, s$sampling_weights, n)
    base_fit <- internal(".fit_logit_nnls")
    Xspread <- -cbind(d$Xdir, d$Xdist)
    Xposition <- cbind(d$Xpos, -d$Xpos)
    Xfull <- cbind(Xspread, Xposition)
    value <- function(keep) {
      vapply(seq_len(ncol(y)), function(j) {
        yy <- y[keep, j]
        w <- s$sampling_weights[keep]
        full <- base_fit(Xfull[keep, , drop = FALSE], yy, w)
        spread <- base_fit(Xspread[keep, , drop = FALSE], yy, w)
        position <- base_fit(Xposition[keep, , drop = FALSE], yy, w)
        if (!full$converged || !spread$converged || !position$converged ||
            !is.finite(full$D2) || full$D2 <= 1e-10 ||
            max(spread$D2, position$D2) > full$D2 + 1e-6) return(NA_real_)
        (full$D2 + position$D2 - spread$D2) / (2 * full$D2)
      }, numeric(1L))
    }
    center <- setNames(value(rep(TRUE, nrow(s$idx))), colnames(y))
    draws <- t(vapply(om$keep, value, numeric(ncol(y))))
    list(center = center, delete = draws, min_ess = min(om$ess))
  }
  if (is.null(low) || choose(n, meta$low) > 50000L) low <- fit_order(meta$low)
  high <- fit_order(meta$high)
  i <- interval(high$center - low$center, high$delete - low$delete, "jackknife")
  rows[[length(rows) + 1L]] <- data.frame(case, sampler = if (alpha == 0) "conventional" else "species",
    order = meta$high, from = meta$low, low_min_ess = low$min_ess, high_min_ess = high$min_ess, i)
}
write_study(do.call(rbind, rows), out, paste(case, "refit50000", sep = "_"))
