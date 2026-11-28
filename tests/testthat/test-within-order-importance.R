test_that("within-order grouped attribution has the declared full total", {
  set.seed(12)
  n <- 25L
  env <- data.frame(a = rnorm(n), b = rnorm(n))
  idx <- sample_sets(n, 3L, M = 250L)
  base <- ugdm_design(env, idx)
  d <- ugdm_position_design(base, env, "nonlinear")
  y <- plogis(-.7 + .25 * d$dir[, 1L] + .2 * d$pos[, 1L] +
                rnorm(nrow(idx), sd = .13))
  fit <- ugdm_fit_position(d, y, distance = FALSE, ridge = .001)
  a <- ugdm_within_order_importance(fit, d, y)
  expect_equal(sum(a$predictors$absolute), a$predictors$total, tolerance = 1e-7)
  expect_equal(sum(a$components$absolute), a$components$total, tolerance = 1e-7)
  expect_setequal(names(a$predictors$absolute), c("a", "b"))
  expect_setequal(names(a$components$absolute), c("spread", "position"))
  b <- ugdm_within_order_importance(fit, d, y, evaluation = d, evaluation_y = y)
  expect_equal(a$predictors$absolute, b$predictors$absolute, tolerance = 1e-7)
})

test_that("held-out importance validates basis and penalty coordinates", {
  set.seed(91)
  n <- 32L
  env <- data.frame(a = rnorm(n), b = rnorm(n))
  train <- seq_len(21L)
  idx <- t(combn(train, 2L))
  ref <- ugdm_design(env, idx, basis_sites = train)
  d <- ugdm_position_design(ref, env, "nonlinear")
  y <- plogis(-.7 + d$pos[, 1L] + .1 * d$pos[, 4L])
  fit <- ugdm_fit_position(d, y, distance = FALSE, ridge = .001)
  ii <- t(combn(setdiff(seq_len(n), train), 2L))
  held <- ugdm_position_design(ugdm_design(env, ii, ref = ref), env, "nonlinear")
  validation <- ugdm_within_order_importance(fit, d, y, evaluation = held,
    evaluation_y = plogis(-.7 - held$pos[, 1L]))
  expect_lt(validation$predictors$total, 0)
  expect_true(all(is.na(validation$predictors$relative)))
  independent <- ugdm_position_design(ugdm_design(env, ii), env, "nonlinear")
  expect_error(ugdm_within_order_importance(fit, d, y, evaluation = independent,
    evaluation_y = rep(.5, nrow(ii))), "training reference")
  fit$penalty[2L, 2L] <- 1
  expect_error(ugdm_within_order_importance(fit, d, y), "coefficient coordinates")
})

test_that("Owen values match sampled blocked orderings and cache each coalition", {
  groups <- list(a = list(spread = 2L, position = 3L),
                 b = list(spread = 4L, position = 5L),
                 distance = list(distance = 6L))
  X <- matrix(c(1, .25, -.5, .3, -.2, .1), 1L,
              dimnames = list(NULL, c("intercept", "x", "x2", "z", "z2", "distance")))
  p <- list(X = X, y = .8, weights = 1, test_X = X, test_y = .57,
            test_weights = 1, groups = groups)
  counter <- new.env(parent = emptyenv())
  mock <- function(keep, start) {
    key <- paste(keep, collapse = ",")
    counter[[key]] <- if (is.null(counter[[key]])) 1L else counter[[key]] + 1L
    list(coef = c(-.27, rep(.8, length(keep) - 1L)), converged = TRUE)
  }
  exact <- zetaGDUM:::.within_importance_game(p, fit_coalition = mock)
  expect_equal(unlist(as.list(counter)) |> unname(), rep(1L, length(ls(counter))))
  expect_equal(sum(exact$hierarchy$absolute), exact$predictors$total)
  expect_equal(as.numeric(tapply(exact$hierarchy$absolute, exact$hierarchy$predictor, sum)),
    unname(exact$predictors$absolute))
  expect_equal(exact$diagnostics$null_deviance,
    zetaGDUM:::.binom_dev(.57, plogis(-.27), 1))
  set.seed(30)
  mc <- zetaGDUM:::.within_importance_game(p, permutations = 2500L,
    exact_max = 0L, fit_coalition = mock)
  expect_true(all(abs(mc$hierarchy$absolute - exact$hierarchy$absolute) <=
    4 * mc$hierarchy$mc_se))
  expect_equal(sum(mc$components$absolute), exact$predictors$total)
  bad <- function(keep, start) {
    if (identical(keep, c(1L, 2L))) return(list(converged = FALSE))
    mock(keep, start)
  }
  expect_error(zetaGDUM:::.within_importance_game(p, fit_coalition = bad), "Coalition 1 failed")
  expect_warning(result <- zetaGDUM:::.within_importance_game(p,
    fit_coalition = bad, failure = "warn"), "coalition fits failed")
  expect_equal(result$diagnostics$failures[["1"]]$columns, c("intercept", "x"))
  expect_true(anyNA(result$hierarchy$absolute))
})
