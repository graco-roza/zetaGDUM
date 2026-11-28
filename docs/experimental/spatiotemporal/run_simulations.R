# Run from repository root. Optional environment settings support quick pilots.
root <- "docs/experimental/spatiotemporal"
source(file.path(root, "load.R"))
st_load(root)
reps <- as.integer(Sys.getenv("ST_REPLICATES", "32"))
if (is.na(reps) || reps < 2L) stop("Use at least two replicates.")
output <- Sys.getenv("ST_OUTPUT", root)
dir.create(file.path(output, "results"), recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, name) {
  utils::write.csv(x, file.path(output, "results", paste0(name, ".csv")), row.names = FALSE, na = "")
}
set.seed(911)
terms <- c(initial = "position_initial.gradient", change = "position_change.gradient",
           interaction = "position_by_change.gradient")

scenario_summary <- function(p, scenario, replicate, observation = "latent") {
  n <- dim(p$x)[1]
  o <- st_occupancy(p, orders = seq_len(min(6L, n)))
  zeta <- aggregate(cbind(before, after, gains, losses, persistent, change) ~ spatial_order,
                    o$contributions, sum)
  species <- o$species
  b <- st_beta_change(p)
  finite_mean <- function(v) if (any(is.finite(v))) mean(v[is.finite(v)]) else NA_real_
  beta <- do.call(rbind, lapply(split(b, b$component), function(z) {
    valid <- is.finite(z$gain_contribution + z$loss_contribution)
    data.frame(component = z$component[1], before = finite_mean(z$before),
               after = finite_mean(z$after), change = finite_mean(z$change),
               gain_contribution = finite_mean(z$gain_contribution[valid]),
               loss_contribution = finite_mean(z$loss_contribution[valid]),
               change_on_counterfactual_support = finite_mean(z$change[valid]),
               endpoint_fraction = mean(is.finite(z$change)), counterfactual_fraction = mean(valid))
  }))
  rows <- data.frame(scenario, replicate, observation, sites = n)
  data <- data.frame(n_expansion = sum(species$expansion), n_contraction = sum(species$contraction),
                     n_redistribution = sum(species$redistribution),
                     n_redistribution_only = sum(species$redistribution_only),
                     n_expansion_and_redistribution = sum(species$expansion & species$redistribution),
                     n_contraction_and_redistribution = sum(species$contraction & species$redistribution),
                     appearances = sum(species$regional_appearance),
                     disappearances = sum(species$regional_disappearance),
                     changed_total = sum(abs(species$total_change) > 1e-10),
                     changed_positive_mean = sum(abs(species$conditional_change) > 1e-10, na.rm = TRUE),
                     changed_persistent_mean = sum(abs(species$persistent_conditional_change) > 1e-10, na.rm = TRUE),
                     changed_dominance = sum(abs(species$dominance_change) > 1e-10, na.rm = TRUE),
                     shared_site_gains = sum(species$gained), shared_site_losses = sum(species$lost),
                     pair_sorensen_change = beta$change[beta$component == "sorensen"],
                     multisite_sorensen_change = unname(st_beta_multisite(p, time = length(p$times))["sorensen"] -
                                                        st_beta_multisite(p)["sorensen"]))
  list(summary = cbind(rows, data), zeta = cbind(rows, zeta), beta = cbind(rows, beta), species = species)
}

scenarios <- c(as.list(st_scenarios),
               list(combined = c("expansion", "redistribution", "abundance_only")))
names(scenarios) <- c(st_scenarios, "combined")
cases <- list()
for (i in seq_along(scenarios)) {
  for (r in seq_len(reps)) {
    s <- simulate_st_gdum(scenarios[[i]], seed = 10000L + r)
    case <- scenario_summary(s$latent, names(scenarios)[i], r)
    cases[[length(cases) + 1L]] <- case[1:3]
    if (r == 1L && names(scenarios)[i] %in% c("expansion", "contraction", "abundance_only", "redistribution", "combined")) {
      write_result(case$species[setdiff(names(case$species), c("new_sites", "lost_sites", "persistent_sites"))],
                   paste0("species_", names(scenarios)[i]))
    }
  }
}
for (r in seq_len(reps)) {
  s <- simulate_st_gdum("spatial_filtering", n_sites = 48, abundance = "count",
                        detection = 0.7, variable_effort = TRUE, missing = 0.1, seed = 20000L + r)
  p <- st_panel(s$observed, s$latent$times, "complete_sites")
  matched <- s$latent
  matched$x <- matched$x[dimnames(p$x)[[1]], , , drop = FALSE]
  cases[[length(cases) + 1L]] <- scenario_summary(p, "detection_failure", r, "observed")[1:3]
  cases[[length(cases) + 1L]] <- scenario_summary(matched, "detection_failure", r, "matched_latent")[1:3]
}
for (group in c("summary", "zeta", "beta")) {
  write_result(do.call(rbind, lapply(cases, `[[`, group)), paste0("scenario_", group))
}

recover <- function(fit, event, term, truth, scenario, r, order = 1L, pool = "all") {
  value <- fit$fit[[event]]
  b <- if (term %in% names(value$coefficients)) value$coefficients[[term]] else NA_real_
  se <- if (term %in% names(value$jackknife_se)) value$jackknife_se[[term]] else NA_real_
  data.frame(scenario, replicate = r, spatial_order = order, pool, event, term,
             status = value$status, truth, estimate = b, se)
}

recovery <- list()
entries <- 0L
for (r in seq_len(reps)) {
  for (scenario in c("gains_only", "losses_only", "interaction", "null",
                     "uncommon_effect", "widespread_effect", "persistent_effect", "abundance_only",
                     "detection_failure")) {
    imperfect <- scenario == "detection_failure"
    s <- simulate_st_gdum(if (imperfect) "spatial_filtering" else scenario, n_sites = 44, n_species = 64,
                         times = c(0, 1, 2, 4, 7), abundance = if (scenario == "abundance_only") "gamma" else "constant",
                         detection = if (imperfect) 0.65 else 1,
                         variable_effort = imperfect, missing = if (imperfect) 0.08 else 0,
                         seed = 30000L + r)
    p <- if (imperfect) st_panel(s$observed, s$latent$times, "complete_sites") else s$latent
    d <- st_design(p, s$env)
    scale <- d$ref$scale[[1]]
    center <- d$ref$center[[1]]
    block <- seq_len(dim(p$x)[1])
    if (scenario == "abundance_only") {
      f <- st_abundance_gdum(p, d, terms["change"], blocks = block)
      entries <- entries + 1L
      recovery[[entries]] <- data.frame(scenario, replicate = r, spatial_order = 1, pool = "all",
                                        event = "positive_abundance", term = unname(terms["change"]),
                                        status = if (f$fit$converged) "estimated" else "nonconverged",
                                        truth = 0.65 * scale, estimate = unname(f$fit$coefficients[2]),
                                        se = unname(f$jackknife_se[2]))
      next
    }
    target <- switch(scenario, uncommon_effect = !s$species$widespread,
                     widespread_effect = s$species$widespread, persistent_effect = s$species$persistent,
                     NULL)
    pools <- if (is.null(target)) list(all = s$species$species) else
      list(target = s$species$species[target], complement = s$species$species[!target])
    if (scenario == "losses_only") pools <- list(low_persistence = s$species$species[!s$species$persistent])
    for (pool in names(pools)) {
      used <- if (scenario == "interaction") unname(terms) else unname(terms[c("initial", "change")])
      f <- st_gdum(p, d, used, pool = pools[[pool]], blocks = block)
      gain_truth <- if (scenario == "interaction") 1.25 * center * scale else
        if (scenario == "gains_only") scale else if (pool == "target") 1.2 * scale else 0
      loss_truth <- if (scenario == "losses_only") -scale else 0
      for (event in c("gain", "loss")) {
        entries <- entries + 1L
        recovery[[entries]] <- recover(f, event, terms["change"],
                                       if (event == "gain") gain_truth else loss_truth, scenario, r, pool = pool)
      }
      if (scenario == "interaction") {
        entries <- entries + 1L
        recovery[[entries]] <- recover(f, "gain", terms["interaction"], 1.25 * scale^2, scenario, r)
      }
    }
  }
}
recovery <- do.call(rbind, recovery)
recovery$error <- recovery$estimate - recovery$truth
recovery$covered_95 <- abs(recovery$error) <= 1.96 * recovery$se
recovery$reject_zero <- abs(recovery$estimate) > 1.96 * recovery$se
write_result(recovery, "parameter_recovery")

# Held-site-out predictions. Whole site histories are withheld, and for k=2
# evaluation uses only pairs entirely within the train or held-out panel.
holdout <- list()
for (r in seq_len(reps)) {
  s <- simulate_st_gdum("interaction", n_sites = 36, n_species = 64, seed = 40000 + r)
  folds <- sample(rep(1:3, each = 12L))
  for (k in 1:2) {
    for (fold in 1:3) {
      train <- s$latent
      test <- s$latent
      train$x <- train$x[folds != fold, , , drop = FALSE]
      test$x <- test$x[folds == fold, , , drop = FALSE]
      dt <- st_design(train, s$env, k)
      ds <- st_design(test, s$env, k, ref = dt$ref)
      actual <- st_gdum(test, ds)$counts
      for (architecture in c("additive", "interaction")) {
        gain_terms <- unname(terms[if (architecture == "additive") 1:2 else 1:3])
        model <- st_gdum(train, dt, gain_terms, loss_terms = character())
        pred <- st_predict(model, ds, actual$before)
        probability <- pmax(pmin(pred$gain_probability, 1 - 1e-10), 1e-10)
        log_score <- -sum(actual$gain * log(probability) +
                            (actual$gain_opportunities - actual$gain) * log1p(-probability)) /
          sum(actual$gain_opportunities)
        holdout[[length(holdout) + 1L]] <- data.frame(replicate = r, spatial_order = k,
                                                       fold, architecture, gain_log_loss = log_score)
      }
    }
  }
}
write_result(do.call(rbind, holdout), "site_holdout")

# Quantify the effect of overlapping sets on the uncertainty calculation.
dependence <- list()
spatial <- list()
for (r in seq_len(reps)) {
  s <- simulate_st_gdum("null", n_sites = 24, n_species = 48, seed = 50000L + r)
  for (k in 1:2) {
    d <- st_design(s$latent, s$env, k)
    f <- st_gdum(s$latent, d, unname(terms[1:2]), blocks = seq_len(24))
    for (event in c("gain", "loss")) {
      fit <- f$fit[[event]]
      j <- terms["change"]
      dependence[[length(dependence) + 1L]] <- data.frame(replicate = r, spatial_order = k, event,
        estimate = unname(fit$coefficients[j]), independent_se = unname(fit$working_independent_se[j]),
        site_se = unname(fit$jackknife_se[j]))
    }
  }
  for (scenario in c("spatial_filtering", "null")) {
    s <- simulate_st_gdum(scenario, n_sites = 36, n_species = 48, seed = 60000L + r)
    for (component in c("turnover", "total", "nestedness", "relative")) {
      p <- s$latent$x[, 1, , drop = FALSE]
      e <- s$env[, 1, , drop = FALSE]
      f <- st_spatial_gdum(p, e, component = component)[[1]]$fit
      spatial[[length(spatial) + 1L]] <- data.frame(replicate = r, scenario, component,
        spread_height = f$height[[1]], position_slope = f$lambda[[1]], D2 = f$D2)
    }
  }
}
write_result(do.call(rbind, dependence), "overlap_dependence")
write_result(do.call(rbind, spatial), "spatial_recovery")

demo <- simulate_st_gdum(c("interaction", "losses_only"), n_sites = 28, seed = 71)
surface <- st_order_diagnostic(demo$latent, 7, 5, demo$env[, , 1], budget = 1500)
write_result(surface, "order_surface")

scaling <- lapply(c(40L, 80L, 160L, 320L), function(n) {
  s <- simulate_st_gdum("null", n_sites = n, seed = 301)
  dt <- system.time(d <- st_zeta_decline(s$latent, 7, 5))
  data.frame(n_sites = n, spatial_max = 7L, temporal_max = 5L,
             rectangles_at_top_order = tail(d$rectangles, 1), elapsed_seconds = unname(dt["elapsed"]))
})
write_result(do.call(rbind, scaling), "scaling")
write_result(data.frame(r = as.character(getRversion()), zetaGDUM = as.character(packageVersion("zetaGDUM")),
                        replicates = reps), "run_metadata")
