# Rscript docs/manuscript_validation/profiles.R OUT CASE [B=100]
source("docs/manuscript_validation/common.R")
Rcpp::sourceCpp("docs/manuscript_validation/profile_helper.cpp",
  cacheDir = Sys.getenv("PROFILE_CPP_CACHE", tempdir()))
args <- commandArgs(TRUE)
out <- args[1L]
case <- args[2L]
B <- if (length(args) > 2L) as.integer(args[3L]) else 100L
meta <- read_manifest()
meta <- meta[match(case, meta$id), , drop = FALSE]
d <- load_dd2(meta, Sys.getenv("DD2_DATA"))
n <- nrow(d$x)
p <- (d$x > 0) * 1L
q <- d$x / rowSums(d$x)
traits <- cbind(occupancy = colMeans(p), local_share = colSums(q) / colSums(p))
representations <- list(total = p * 1, relative = q)
orders <- as.integer(meta[c("low", "mid", "high")])
W <- site_weights(n, B, seed = 42006)
tt <- ncol(traits)
profiles <- function(a) {
  sim <- a[, 2:(tt + 1L), drop = FALSE] / a[, 1L]
  dis <- a[, (tt + 3L):(2L * tt + 2L), drop = FALSE] / a[, tt + 2L]
  out <- cbind(sim, dis, sim - dis)
  colnames(out) <- paste0("profile.", rep(c("similarity", "dissimilarity", "gap"), each = tt),
                          ".", rep(colnames(traits), 3))
  out
}
idx_check <- sample_sets(n, 3L, M = 20L, seed = 93)
for (name in names(representations)) {
  p0 <- profile_totals(representations[[name]], idx_check, traits)
  pp <- species_contributions(d$x, idx_check, if (name == "total") "presence" else "relative")
  expected <- cbind(rowSums(pp$similarity), pp$similarity %*% traits,
                    rowSums(pp$dissimilarity), pp$dissimilarity %*% traits)
  stopifnot(max(abs(p0 - expected)) < 1e-12)
}
rows <- list()
for (alpha in c(0, .5)) {
  sampler <- if (alpha == 0) "conventional" else "species"
  fits <- list()
  for (k in orders) {
    message(case, " profile ", sampler, " ", k)
    p0 <- draw(p, k, 10000L, alpha, seed = 530000L + k)
    mats <- lapply(representations, function(x) profile_totals(x, p0$idx, traits))
    boot <- lapply(seq_len(B), function(b) {
      if (p0$exact) {
        w <- internal(".product_sampling_weights")(W[, b], p0$idx, p0$sampling_weights)
        lapply(mats, function(x) profiles(matrix(colSums(x * w), 1L))[1L, ])
      } else {
        s <- draw(p, k, 10000L, alpha, omega = W[, b], seed = 635000L + k + 101L * b)
        lapply(representations, function(x) profiles(matrix(colSums(
          profile_totals(x, s$idx, traits) * s$sampling_weights), 1L))[1L, ])
      }
    })
    fits[[as.character(k)]] <- lapply(names(mats), function(response) {
      mat <- mats[[response]] * p0$sampling_weights
      u <- colSums(mat)
      delete <- profiles(matrix(u, n, length(u), byrow = TRUE) - inclusion(p0$idx, mat, n))
      list(center = profiles(matrix(u, 1L))[1L, ], jackknife = delete,
           bayesian = do.call(rbind, lapply(boot, function(x) x[[response]])))
    })
    names(fits[[as.character(k)]]) <- names(mats)
  }
  for (k in orders) for (response in names(representations)) {
    fit <- fits[[as.character(k)]][[response]]
    null <- fits[[as.character(orders[1L])]][[response]]
    for (method in c("bayesian", "jackknife")) {
      rows[[length(rows) + 1L]] <- data.frame(case, sampler, order = k, from = NA_integer_,
        response, method, interval(fit$center, fit[[method]], method))
      if (k != orders[1L]) {
        rows[[length(rows) + 1L]] <- data.frame(case, sampler, order = k, from = orders[1L],
          response, method, interval(fit$center - null$center, fit[[method]] - null[[method]], method))
      }
      if (response == "relative") {
        incidence <- fits[[as.character(k)]]$total
        rows[[length(rows) + 1L]] <- data.frame(case, sampler, order = k, from = NA_integer_,
          response = "relative_minus_total", method,
          interval(fit$center - incidence$center, fit[[method]] - incidence[[method]], method))
      }
    }
  }
}
write_study(do.call(rbind, rows), out, paste(case, "profiles", sep = "_"))
