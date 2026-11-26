source("docs/manuscript_validation/common.R")
write5 <- function(x, name) write_study(x, "docs/phase5/runs", name)
qweight <- function(x, w, probs) {
  use <- is.finite(x) & is.finite(w) & w > 0
  if (!any(use)) return(rep(NA_real_, length(probs)))
  x <- x[use]
  w <- w[use]
  i <- order(x)
  vapply(probs, function(p) x[i[which(cumsum(w[i]) / sum(w) >= p - 1e-12)[1L]]], numeric(1))
}
model_data <- function(a, component, model, keep = rep(TRUE, length(a$relative))) {
  y <- a[[component]]
  keep <- keep & is.finite(y) & a$sampling_weights > 0
  d <- ugdm_model_design(internal(".subset_order_design")(a$design, keep), model)
  list(design = d, y = y[keep], trials = rep(1, sum(keep)), sampling_weights = a$sampling_weights[keep])
}
fit5 <- function(a, model, free_sign = FALSE, smooth = .001, ridge = .0001, start = NULL) {
  ugdm_fit_variant(a$design, a$y, a$trials, model = model,
    sampling_weights = a$sampling_weights, free_sign = free_sign, smooth = smooth, ridge = ridge, start = start)
}
effect5 <- function(a, f, env, local = NULL, diagnose = TRUE) {
  d <- a$design
  if (!f$converged || !is.finite(f$D2)) return(rep(NA_real_, 14L))
  b <- d$environment_basis
  s <- (quantile(env[[1]], c(.25, .75), names = FALSE) - b$center[1]) / b$scale[1]
  spread <- if (identical(f$model, "B")) diff(s) * f$beta[1] else
    sum(f$beta[1:3] * colSums(diff(predict(b$splines[[1]], s))))
  if (diagnose) {
    dg <- ugdm_variant_diagnostics(f, d, a$y, a$trials, a$sampling_weights)
    imp <- ugdm_importance(f, d, sampling_weights = a$sampling_weights)
  } else {
    dg <- list(shapley = rep(NA_real_, 2), relative_shapley = rep(NA_real_, 2),
               condition = NA_real_, active = sum(f$beta < 1e-8))
    imp <- NA_real_
  }
  c(spread = unname(.1 * spread), position = unname(diff(s) * f$lambda[1]),
    importance = unname(imp[1]), shapley_spread = unname(dg$shapley[1]),
    shapley_position = unname(dg$shapley[2]), relative_spread = unname(dg$relative_shapley[1]),
    relative_position = unname(dg$relative_shapley[2]), condition = dg$condition, active = dg$active,
    D2 = f$D2, local_spread = if (is.null(local)) NA_real_ else local$spread[1],
    local_position = if (is.null(local)) NA_real_ else local$position[1],
    local_residual = if (is.null(local)) NA_real_ else local$residual[1],
    local_predicted = if (is.null(local)) NA_real_ else local$predicted[1])
}
effect_names <- c("spread", "position", "importance", "shapley_spread", "shapley_position",
                  "relative_spread", "relative_position", "condition", "active", "D2",
                  "local_spread", "local_position", "local_residual", "local_predicted")
