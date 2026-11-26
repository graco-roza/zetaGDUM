test_that("global score controls preserve exact uniform and weighted targets", {
  set.seed(896)
  n <- 9L
  x <- cbind(2L, matrix(rpois(n * 5L, 1), n))
  env <- data.frame(temperature = seq(-2, 3, length.out = n))
  orders <- c(2L, 4L)
  exact <- function(k) {
    sets <- t(combn(n, k))
    list(idx = sets, sampling_weights = rep(1, nrow(sets)), exact = TRUE)
  }
  exact_paths <- lapply(orders, exact)
  unequal <- lapply(exact_paths, function(a) {
    counts <- rep(c(1L, 2L, 4L), length.out = nrow(a$idx))
    list(idx = a$idx[rep(seq_along(counts), counts), , drop = FALSE],
         sampling_weights = rep(sum(counts) / (length(counts) * counts), counts))
  })
  pilot <- sample_sets_coupled(x, orders, M = 120L, sampler = "species", exact_limit = 0L)
  theta <- matrix(seq(-.25, .25, length.out = 10L), nrow = 5L)
  n_coef <- nrow(theta)
  for (type in c("relative", "incidence")) {
    for (weighting in if (type == "relative") "constant" else c("constant", "richness")) {
      plan <- ugdm_global_score_plan(x, env, orders, pilot, theta = c(theta), type = type,
                                     weighting = weighting)
      raw <- ugdm_global_score_plan(x, env, orders, type = type, weighting = weighting,
                                    controls = character())
      a <- ugdm_global_site_scores(plan, exact_paths, c(theta))
      b <- ugdm_global_site_scores(plan, unequal, c(theta))
      c <- ugdm_global_site_scores(raw, exact_paths, c(theta))
      expect_equal(a$pseudovalues, b$pseudovalues, tolerance = 1e-12)
      expect_equal(a$pseudovalues, c$pseudovalues, tolerance = 1e-12)
      expect_equal(a$score, colMeans(a$pseudovalues), tolerance = 1e-12)
      order_data <- ugdm_path_data(x, env, exact_paths)
      for (j in seq_along(orders)) {
        d <- order_data[[j]]
        trial <- if (weighting == "richness") d$mean_richness else NULL
        ref <- ugdm_score_pseudovalues(d$design, d[[type]], theta[, j], trials = trial,
                                      method = "inclusion")
        expect_equal(unname(a$pseudovalues[, n_coef * (j - 1L) + seq_len(n_coef)]),
                     unname(ref$pseudovalues), tolerance = 1e-12)
      }
    }
  }
  intercept <- function(a) cbind(intercept = rep(1, nrow(a$design$idx)))
  simple <- ugdm_global_score_plan(x, env, orders, pilot, c(0, .2),
    controls = "response", design_matrix = intercept)
  simple_exact <- ugdm_global_site_scores(simple, exact_paths, c(0, .2))
  expect_equal(dim(simple_exact$pseudovalues), c(n, 2L))
  expect_equal(simple_exact$pseudovalues,
               ugdm_global_site_scores(simple, unequal, c(0, .2))$pseudovalues)
  expect_equal(simple_exact$pilot_evaluations, 240L)
})
