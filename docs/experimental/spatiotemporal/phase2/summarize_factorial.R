root <- "docs/experimental/spatiotemporal"
out <- Sys.getenv("ST2_OUTPUT", file.path(root, "phase2/results"))
d <- utils::read.csv(file.path(out, "factorial_contrasts.csv"))
save <- function(x, name) utils::write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "")
a <- subset(d, target > 0 & ((scenario %in% c("gain", "loss") & response == "presence") |
                            (scenario == "abundance" & response == "abundance")))
a$target_cell <- a$group == a$target
a$squared_error <- a$error^2
keep <- c("architecture", "estimand", "scenario","target_cell")
moment <- lapply(split(subset(a, group > 0), a[a$group > 0, keep], drop = TRUE), function(z) {
  rmse <- sqrt(mean(z$squared_error, na.rm = TRUE))
  zero <- sqrt(mean(z$truth^2, na.rm = TRUE))
  valid <- is.finite(z$truth) & is.finite(z$estimate)
  correlation <- if (sum(valid) >= 2L && sd(z$truth[valid]) > 0 && sd(z$estimate[valid]) > 0)
    cor(z$truth[valid], z$estimate[valid]) else NA_real_
  cbind(z[1, keep], bias = mean(z$error, na.rm = TRUE), rmse = rmse, zero_contrast_rmse = zero,
        rmse_ratio = if (zero > 0) rmse / zero else NA_real_, correlation = correlation,
        scored = sum(is.finite(z$estimate)), row.names = NULL)
})
save(do.call(rbind, moment), "factorial_summary")

group_energy <- lapply(split(subset(a, group > 0), a[a$group > 0, c("architecture", "estimand", "scenario", "target","replicate")],
                             drop = TRUE), function(z) {
  rms <- aggregate(cbind(truth = truth^2, estimate = estimate^2) ~ group, z, mean)
  t <- z$target[1L]
  if (!nrow(rms)) return(NULL)
  data.frame(z[1L, c("architecture", "estimand", "scenario", "target", "replicate")],
             truth_top = rms$group[which.max(rms$truth)],
             estimated_top = rms$group[which.max(rms$estimate)],
             truth_fraction = sum(rms$truth[rms$group == t]) / sum(rms$truth),
             estimated_fraction = sum(rms$estimate[rms$group == t]) / sum(rms$estimate),
             row.names = NULL)
})
save(do.call(rbind, group_energy), "factorial_group_discrimination")

# Evaluate finite-panel distinguishability from group-unlabelled aggregate
# truth surfaces. Each pair's distance is scaled by typical within-label noise.
pooled <- subset(a, group == 0 & architecture == "M2" & estimand == "trajectory_surface")
pooled$order <- paste(pooled$spatial_order, pooled$temporal_order)
distances <- lapply(split(pooled, pooled$scenario), function(z) {
  groups <- lapply(split(z, z$target), function(g) {
    x <- xtabs(truth ~ replicate + order, g)
    center <- colMeans(x)
    list(center = center, noise = mean(sweep(x, 2L, center)^2))
  })
  pair <- t(combn(seq_along(groups), 2L))
  rows <- lapply(seq_len(nrow(pair)), function(i) {
    a <- groups[[pair[i,1L]]]; b <- groups[[pair[i,2L]]]
    between <- sqrt(mean((a$center - b$center)^2))
    noise <- sqrt((a$noise + b$noise)/2)
    data.frame(scenario=z$scenario[1L], group_a=pair[i,1L], group_b=pair[i,2L],
               between=between, within=noise, distance_noise_ratio=between/noise)
  })
  do.call(rbind, rows)
})
save(do.call(rbind, distances), "factorial_surface_similarity")
