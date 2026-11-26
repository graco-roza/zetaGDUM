library(zetaGDUM)
within_dir <- "docs/within_order/results"
within_grid <- expand.grid(n = c(70L, 150L), order = c(2L, 4L, 8L, 12L),
                           response = c("incidence", "relative"))

within_reference <- function() {
  env <- data.frame(strong = seq(-1, 1, length.out = 61L))
  env$weak <- env$strong[(seq_len(61L) * 11L) %% 61L + 1L]
  env$irrelevant <- env$strong[(seq_len(61L) * 23L) %% 61L + 1L]
  ugdm_design(env, cbind(1:60, 2:61))
}

within_sites <- function(n, scenario = "baseline") {
  latent <- matrix(rnorm(n * 3L), n)
  if (scenario == "spatial") {
    for (i in 2:n) latent[i, ] <- .65 * latent[i-1L, ] + sqrt(1 - .65^2) * latent[i, ]
  }
  if (scenario == "correlated") latent[, 3L] <- .75 * latent[, 1L] + sqrt(1 - .75^2) * latent[, 3L]
  env <- as.data.frame(2 * pnorm(latent) - 1)
  names(env) <- c("strong", "weak", "irrelevant")
  optima <- expand.grid(strong = seq(-1, 1, length.out = 10L), weak = seq(-1, 1, length.out = 4L))
  first <- if (scenario == "broad") .7 else .30
  second <- if (scenario == "equal") first else .75
  intensity <- exp(-.5 * (outer(env$strong, optima$strong, "-") / first)^2 -
                     .5 * (outer(env$weak, optima$weak, "-") / second)^2)
  if (scenario == "position") intensity <- intensity * exp(1.5 * env$strong^2)
  if (scenario == "equal") {
    optima <- expand.grid(strong = seq(-1, 1, length.out = 7L), weak = seq(-1, 1, length.out = 7L))
    intensity <- exp(-.5 * (outer(env$strong, optima$strong, "-") / first)^2 -
                       .5 * (outer(env$weak, optima$weak, "-") / second)^2)
  }
  if (scenario == "sparse") intensity <- intensity * .35
  counts <- cbind(1L, matrix(rpois(length(intensity), 3 * intensity + .015), n))
  list(env = env, counts = counts, blocks = cut(seq_len(n), min(12L, n %/% 5L)))
}

within_data <- function(sites, k, M, response, ref = within_reference(), idx = NULL) {
  if (is.null(idx)) idx <- sample_sets(nrow(sites$env), k, M)
  basis <- ugdm_design(sites$env, idx, ref = ref)
  d <- ugdm_interaction_design(basis, sites$env, position = "nonlinear", spread = "constant")
  r <- ugdm_response(sites$counts, idx, if (response == "incidence") "total" else "relative",
                     working_weights = "constant")
  list(design = d, y = r$y, trials = r$trials, idx = idx, base = basis)
}

within_fit <- function(a) ugdm_fit_interactions(a$design, a$y, a$trials,
  ridge = 1e-4, spread_smooth = 1e-4, position_smooth = 1e-4)

within_write <- function(a, name) {
  dir.create(within_dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(within_dir, name)
  tmp <- tempfile(tmpdir = within_dir)
  if (endsWith(name, ".gz")) {
    con <- gzfile(tmp, "wt")
    write.csv(a, con, row.names = FALSE)
    close(con)
  } else write.csv(a, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not write result.")
}
