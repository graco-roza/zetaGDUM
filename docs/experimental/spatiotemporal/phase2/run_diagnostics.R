root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase2_load.R"))
st2_load(root)
for (f in c("factorial", "prediction_validation")) source(file.path(root, "phase2", paste0(f, ".R")))
out <- Sys.getenv("ST2_OUTPUT", file.path(root, "phase2/results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save <- function(x, name) utils::write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "")
sim <- st2_factorial_world(88114, "additive", region = TRUE)
ref <- st2_basis(sim$enabled, sim$env)
p <- sim$enabled

zeta <- list()
for (mode in c("all", "consecutive")) {
  z <- st2_zeta(p, 5, 4, mode = mode)
  zeta[[length(zeta) + 1L]] <- z$surface
  if (mode == "consecutive") {
    index <- expand.grid(spatial_order = 1:5, temporal_order = 1:4)
    save(cbind(index, occurrence_mean = rowMeans(z$contributions$presence),
               conditional_minimum = rowSums(z$contributions$abundance) /
                 rowSums(z$contributions$presence),
               zeta_abundance = rowSums(z$contributions$abundance)),
         "decomposition_surface")
  }
}
save(do.call(rbind, zeta), "zeta_surfaces")
d <- st2_order_diagnostic(p, sim$env, 5, 4, budget = 64, mc_replicates = 5, basis = ref)
save(d$orders, "order_diagnostic")
save(d$weights, "order_weights")

# Change every site by the same standardized amount, keeping initial position
# and spread fixed. Means and local deviations are recomputed by construction.
set.seed(2017)
a <- list()
for (k in 1:3) {
  design <- st2_transition_design(p, sim$env, k, sets = st2_draw_sets(60, k, 80), basis = ref)
  model <- st2_fit_transition(p, design)
  grid <- expand.grid(delta = seq(-2, 2, length.out = 25), r0 = 0:k)
  select <- design
  select$X <- matrix(0, nrow(grid) * k, ncol(design$X), dimnames = list(NULL, colnames(design$X)))
  select$X[, "position_change.climate"] <- rep(k * grid$delta, each = k)
  select$row_id <- rep(seq_len(nrow(grid)), each = k)
  select$sets <- matrix(rep(seq_len(k), times = nrow(grid)), nrow(grid), k, byrow = TRUE)
  states <- list(z = matrix(rep(seq_len(k), nrow(grid)) <= rep(grid$r0, each = k),
                            ncol = 1, dimnames = list(NULL, model$species[1L])))
  states$z <- matrix(rep(states$z, length(model$species)), ncol = length(model$species),
                     dimnames = list(NULL, model$species))
  if (k > 1L) states$other <- matrix(rep(grid$r0, each = k), nrow(states$z), ncol(states$z)) - states$z
  else states$other <- NULL
  pr <- st2_predict_transition(model, select, states)
  q <- st2_transition_summary(pr, select, states)
  average <- function(x) {
    r <- rowsum(x, select$row_id, reorder = FALSE)[, 1L]
    r
  }
  a[[k]] <- data.frame(spatial_order = k, grid, gain_probability = .st_ratio(average(pr$gain * !states$z), k-grid$r0),
    loss_probability = .st_ratio(average(pr$loss * states$z), grid$r0),
    expansion = q$expansion[, 1L], contraction = q$contraction[, 1L],
    redistribution = q$redistribution[, 1L], pure_redistribution = q$pure_redistribution[, 1L])
}
save(do.call(rbind, a), "environmental_response")

# Native species-informed proposal on one predeclared initial window, benchmarked
# against independent uniform draws. HT summary targets the uniform population.
sampling <- list()
for (k in c(2L, 4L, 6L)) for (r in 1:16) {
  seed <- 84101L + r
  for (method in c("uniform", "species_informed")) {
    alpha <- if (method == "uniform") 0 else .5
    count <- .st2_minimum(p$x, 1:3)
    draw <- zetaGDUM::sample_sets_species(count, k, 100L, alpha = alpha, seed = seed)
    q <- exp(-lchoose(nrow(count), k)) / draw$sampling_weights
    basis_test <- st2_transition_design(p, sim$env, k, 3L, sets = draw$idx, proposal = q, basis = ref)
    fit <- st2_fit_transition(p, basis_test)
    exact <- sum(.st_choose_prob(colSums(count > 0), k, nrow(count)))
    sampling[[length(sampling) + 1L]] <- data.frame(spatial_order = k, replicate = r, method, truth = exact,
      HT = mean(draw$sampling_weights * draw$shared), shared_fraction = mean(draw$shared > 0),
      effective_set_size = sum(draw$sampling_weights)^2 / sum(draw$sampling_weights^2),
      climate_gain = fit$fit$gain$coefficients[["position_change.climate"]], status = fit$fit$gain$status)
  }
}
save(do.call(rbind, sampling), "sampling_precision")

# Published pairwise accounting on fixed paired samples, with all signs retained.
if (requireNamespace("ecopart", quietly = TRUE)) {
  a <- .st_slice(p, 1); b <- .st_slice(p, 12)
  idx <- st2_draw_sets(nrow(a), 2, 50)
  shapley <- st2_beta_accounting(a, b, idx)
  save(shapley, "native_beta_change")
  parts <- lapply(c("sorensen", "bray-curtis"), function(index) {
    values <- t(vapply(seq_len(nrow(idx)), function(j) st2_tatsumi(a[idx[j, ], ], b[idx[j, ], ], index), numeric(4)))
    colnames(values) <- c("loss_homogenising", "loss_differentiating",
                           "gain_homogenising", "gain_differentiating")
    data.frame(index, set = seq_len(nrow(idx)), values, row.names = NULL)
  })
  save(do.call(rbind, parts), "tatsumi_pair")
}
# Detection experiment with the same latent communities, one visit and no latent model.
detection <- list()
for (r in 1:16) {
  sim <- st2_factorial_world(10000L + r, "null")
  a <- sim$enabled
  observed <- a
  pdet <- plogis(.85 + 1.35 * sim$env[, , 1L])
  detected <- array(runif(length(a$x)), dim(a$x)) < rep(as.vector(pdet), dim(a$x)[3L])
  observed$x <- a$x * detected
  for (panel in c("latent", "observed")) {
    x <- if (panel == "latent") a else observed
    design <- st2_transition_design(x, sim$env, 1)
    f <- st2_fit_transition(x, design)
    for (event in c("gain","loss")) detection[[length(detection)+1L]] <- data.frame(
      replicate = r, panel, event, coefficient = f$fit[[event]]$coefficients[["position_change.climate"]],
      events = f$fit[[event]]$events, opportunities = f$fit[[event]]$opportunities)
  }
}
save(do.call(rbind, detection), "detection_sensitivity")

# Fixed-neighbourhood sensitivity and count outcome family comparison on held-out histories.
extras <- list()
for (r in 1:6) {
  sim <- st2_factorial_world(10550L + r, "additive", region = TRUE, count = TRUE)
  a <- sim$enabled
  n <- dim(a$x)[1]
  test <- which(sim$blocks == 6L)
  training <- st2_subset_panel(a, setdiff(seq_len(n), test))
  testing <- st2_subset_panel(a, test)
  ref <- st2_basis(training, sim$env)
  d0 <- st2_transition_design(training, sim$env, 1L, basis = ref)
  d1 <- st2_transition_design(testing, sim$env, 1L, basis = ref)
  z <- st2_states(testing, d1)
  blocks <- sim$blocks
  ids <- dimnames(a$x)[[1L]]
  neighbourhood <- setNames(lapply(seq_len(n), function(i) ids[setdiff(which(blocks == blocks[i]), i)]), ids)
  m <- st2_fit_dynamic(training, d0)
  mn <- st2_fit_dynamic(training, d0, sim$env, neighbourhood)
  for (method in c("M2", "M2_neighbourhood")) {
    risk <- st2_predict_dynamic(if (method == "M2") m else mn, d1, sim$env)
    extras[[length(extras)+1L]] <- data.frame(replicate = r, method, response = "site_gain",
      error = st2_scores(risk$gain, z$y, !z$z)[["log_loss"]], metric = "log_loss", status = "scored")
  }
  for (method in c("lognormal", "truncated_nbinom")) {
    if (method == "truncated_nbinom" && !suppressWarnings(requireNamespace("glmmTMB", quietly = TRUE))) next
    model <- suppressWarnings(st2_fit_positive(training, d0, family = method))
    forecast <- st2_predict_positive_parameters(model, d1, z)
    pred <- if (method == "lognormal") exp(forecast$location + forecast$scale^2 / 2) else {
      mu <- exp(forecast$location)
      mu / (1 - dnbinom(0, mu = mu, size = forecast$scale))
    }
    extras[[length(extras)+1L]] <- data.frame(replicate = r, method, response = "positive_count",
      error = mean(abs(pred[z$y] - z$after[z$y])), metric = "MAE",
      status = paste(model$fit$new$status, model$fit$persistent$status, sep = "/"))
  }
}
save(do.call(rbind, extras), "model_sensitivities")
