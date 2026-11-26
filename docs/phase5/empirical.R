source("docs/phase5/common.R")
args <- commandArgs(TRUE)
meta <- read_manifest()[as.integer(args[2L]), ]
data <- load_dd2(meta, args[1L])
env <- data$env
x <- data$x
dm <- data$dist
orders <- unlist(meta[c("low", "mid", "high")], use.names = FALSE)
model_names <- c("A","B","C")
responses <- c("incidence","relative")
n <- nrow(x)
base <- suppressWarnings(ugdm_design(env, cbind(seq_len(n-1L), 2:n), dm, dist_range = range(as.dist(dm))))
key <- as.character(meta$id)
stamp <- function(z) cbind(community = meta$community, n = n, z)
point <- local <- intervals <- support <- list()
message(meta$community, " full-model point grid")
for (M in c(1500L,5000L)) for (coupled in c(FALSE,TRUE)) {
  paths <- sample_sets_coupled(x, orders, M, coupled = coupled, seed = 372 + M)
  a <- ugdm_path_data(x, env, paths, dm, base)
  for (j in seq_along(orders)) for (component in responses) for (m in model_names) {
    a_m <- model_data(a[[j]], component, m)
    if (!length(a_m$y)) next
    f <- fit5(a_m, m)
    e <- setNames(effect5(a_m, f, env), effect_names)
    point[[length(point) + 1L]] <- data.frame(order = orders[j], M, coupled, component, model = m,
      converged = f$converged, exact = a[[j]]$exact, as.list(e))
    if (M == 5000 && coupled && f$converged && is.finite(f$D2)) {
      lc <- ugdm_local_distinctiveness(x, f, a_m$design,
        if (component == "relative") "relative" else "normalised", a_m$sampling_weights)
      local[[length(local) + 1L]] <- data.frame(order = orders[j], component, model = m,
        site_number = seq_len(nrow(lc)), lc[c("observed","predicted","residual","baseline","spread","position")])
      # Predeclared deterministic exemplars from observed pairs near a common anchor.
      d0 <- a_m$design
      if (orders[j] == 2L) {
        ref <- d0$idx[1L, ]
        mat <- matrix(d0$idx %in% ref, nrow(d0$idx))
        overlap <- rowSums(mat)
        diff_pos <- rowSums((sweep(d0$Xpos,2,d0$Xpos[1L,]) / 2)^2)
        diff_sp <- abs(d0$Xdir[,1] - d0$Xdir[1L,1])
        member <- list(one_site = which(overlap == 1L),
          whole_set = which(overlap == 0L),
          near_fixed_position = which(diff_sp >= quantile(diff_sp,.60)))
        to <- vapply(member, function(i) i[which.min(diff_pos[i])], integer(1))
        u <- ugdm_configuration_effects(f, d0, rep(1L,3), to)
        support[[length(support) + 1L]] <- data.frame(component, model = m, name = names(member),
          u$effects, position_shift_norm = sqrt(rowSums(u$mean_position_change^2)))
      }
    }
  }
}
write5(stamp(do.call(rbind, point)), paste0("empirical_points_", key))
write5(stamp(do.call(rbind, local)), paste0("empirical_local_", key))
write5(stamp(do.call(rbind, support)), paste0("empirical_config_", key))

message(meta$community, " matched uncertainty")
orders <- orders[1:2]
M <- 2000L
a <- ugdm_path_data(x, env, sample_sets_coupled(x, orders, M, exact_limit = 0, seed = 571), dm, base)
weights <- site_weights(n, 799, seed = 661)
keep_effects <- c(1L,2L,4L,5L,6L,7L)
names6 <- effect_names[keep_effects]
for (component in responses) {
  cp <- vector("list", length(orders))
  for (j in seq_along(orders)) {
    cp[[j]] <- t(vapply(seq_len(39L), function(r) {
      idx <- sample_sets_weighted(weights[,r], orders[j], M)
      tmp <- ugdm_path_data(x, env, list(list(idx = idx, sampling_weights = rep(1,M))), dm, base)[[1L]]
      zz <- model_data(tmp, component, "B")
      effect5(zz, fit5(zz, "B"), env)[keep_effects]
    }, numeric(length(names6))))
  }
  for (m in model_names) {
    q <- lapply(a, model_data, component = component, model = m)
    f <- lapply(q, fit5, model = m)
    eq <- lapply(seq_along(q), function(j) effect5(q[[j]], f[[j]], env)[keep_effects])
    D <- lapply(seq_along(q), function(j) {
      z <- q[[j]]
      t(vapply(seq_len(n), function(i) {
        keep <- rowSums(z$design$idx == i) == 0
        zz <- model_data(a[[j]], component, m, keep = rowSums(a[[j]]$design$idx == i) == 0)
        start <- f[[j]]
        if (!any(keep) || !start$converged || !length(zz$y)) return(rep(NA_real_, length(names6)))
        fz <- fit5(zz, m, start = start)
        effect5(zz, fz, env)[keep_effects]
      }, numeric(length(names6))))
    })
    for (j in 1:3) {
      point <- if (j == 3L) eq[[2]] - eq[[1]] else eq[[j]]
      delete <- if (j == 3L) D[[2]] - D[[1]] else D[[j]]
      cp_ <- if (j == 3L) cp[[2]] - cp[[1]] else cp[[j]]
      a_ok <- which(colSums(is.finite(delete)) == n & is.finite(point))
      if (!length(a_ok)) next
      js <- internal(".jackknife_sd")(delete)
      pe <- ugdm_pseudovalue_bootstrap(point[a_ok], delete[,a_ok,drop=FALSE], weights = weights)
      for (ii in seq_along(a_ok)) {
        t <- a_ok[ii]
        dr <- pe$draws[,ii]
        bb <- cp_[,t]
        comp <- if (t == 4 && j == 3L) dr else if (t == 4 && 3L %in% a_ok)
          dr - pe$draws[,match(3L,a_ok)] else if (t <= 2L) dr else NULL
        intervals[[length(intervals) + 1L]] <- data.frame(component, model = m,
          order = if (j == 3L) orders[2L] else orders[j], previous = if (j == 3L) orders[1L] else NA,
          effect = names6[t], point = point[t],
          jack_se = js[t], jack_lo = point[t] - qt(.975,n-1)*js[t],
          jack_hi = point[t] + qt(.975,n-1)*js[t],
          pseudo_lo = quantile(dr,.025,names=FALSE), pseudo_hi = quantile(dr,.975,names=FALSE),
          p_positive = mean(dr > 0), p_practical = mean(dr > .01), n_valid_delete = n,
          p_hypothesis = if (is.null(comp)) NA_real_ else mean(comp > 0),
          p_hypothesis_practical = if (is.null(comp)) NA_real_ else mean(comp > .01),
          cp_lo = if (m=="B") qweight(bb, rep(1,length(bb)), .025) else NA_real_,
          cp_hi = if (m=="B") qweight(bb, rep(1,length(bb)), .975) else NA_real_)
      }
    }
  }
}
write5(stamp(do.call(rbind, intervals)), paste0("empirical_uncertainty_", key))
