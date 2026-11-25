library(zetaGDUM)
internal <- function(name) get(name, asNamespace("zetaGDUM"))
draw <- internal(".inference_sets")
subset_design <- internal(".subset_order_design")
inclusion <- internal(".inclusion_sums")
kish <- internal(".kish")
jk_sd <- internal(".jackknife_sd")
info <- internal(".order_information")
loss <- internal(".binom_dev")
read_manifest <- function() read.csv("docs/manuscript_validation/communities.csv")
write_study <- function(value, out, name) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = out)
  write.csv(value, tmp, row.names = FALSE)
  if (!file.rename(tmp, file.path(out, paste0(name, ".csv")))) stop("Failed to write ", name)
}
interval <- function(point, draws, type, min_valid = .95) {
  if (is.null(dim(draws))) stop("Draws must be a matrix.")
  finite <- colSums(is.finite(draws))
  if (type == "jackknife") {
    se <- jk_sd(draws)
    span <- qt(.975, nrow(draws) - 1L) * se
    lim <- rbind(point - span, point + span)
  } else {
    lim <- vapply(seq_along(point), function(j) {
      x <- draws[, j]
      if (finite[j] < min_valid * length(x)) return(rep(NA_real_, 2L))
      quantile(x, c(.025, .975), na.rm = TRUE, names = FALSE)
    }, numeric(2L))
    se <- apply(draws, 2L, sd, na.rm = TRUE)
    se[finite < min_valid * nrow(draws)] <- NA_real_
  }
  data.frame(term = names(point), estimate = unname(point), low = lim[1L, ], high = lim[2L, ],
             se = unname(se), valid = finite, attempted = nrow(draws), row.names = NULL)
}

load_dd2 <- function(meta, root) {
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  table <- function(tag) {
    matches <- all_files[basename(all_files) == paste0(meta$id, "-", tag, ".csv")]
    if (length(matches) != 1L) stop("Expected exactly one ", tag, " table for ", meta$id)
    read.csv(matches, row.names = 1L, check.names = FALSE)
  }
  x <- table("S")
  env <- table("E")
  xy <- table("C")
  if (meta$id == "N63TBIRD") {
    stopifnot(identical(rownames(x), rownames(env)), nrow(x) == nrow(xy),
              identical(rownames(xy), as.character(seq_len(nrow(xy)))))
    rownames(xy) <- rownames(x)
    env <- env[, c("AprMay", "DJF", "Urb", "Br", "Co", "Op")]
    env[c("Urb", "Op")] <- lapply(env[c("Urb", "Op")], function(v) log1p(100 * v))
    stopifnot(cor(xy$Lat, env$AprMay, use = "complete.obs") < -.9)
  } else {
    stopifnot(setequal(rownames(x), rownames(env)), setequal(rownames(x), rownames(xy)))
  }
  env <- env[rownames(x), , drop = FALSE]
  coord <- match(c("long", "lat"), sub("^lon$", "long", tolower(names(xy))))
  xy <- xy[rownames(x), coord, drop = FALSE]
  keep <- complete.cases(x, env, xy) & rowSums(x) > 0
  x <- as.matrix(x[keep, , drop = FALSE])
  x <- x[, colSums(x) > 0, drop = FALSE]
  env <- env[keep, , drop = FALSE]
  xy <- xy[keep, , drop = FALSE]
  if (meta$id != "N63TBIRD") env[] <- lapply(env, function(v) {
    if (sd(v) > 0 && mean(as.numeric(scale(v))^3) > 2) log1p(v - min(v)) else v
  })
  stopifnot(nrow(x) == meta$sites, ncol(x) == meta$species, all(x >= 0),
            any(x > 0 & x != 1), all(abs(xy[[1]]) <= 180), all(abs(xy[[2]]) <= 90))
  list(x = x, env = env, dist = coords_dist(xy, lonlat = TRUE))
}

contrast_basis <- function(env, design) {
  tr <- design$environment_basis
  qs <- apply(env, 2L, quantile, probs = c(.25, .75))
  spread <- lapply(seq_along(design$vars), function(j) {
    x <- (qs[, j] - tr$center[j]) / tr$scale[j]
    B <- predict(tr$splines[[j]], x)
    c(.1 * (B[2L, ] - B[1L, ]), diff(x))
  })
  do.call(rbind, spread)
}

empirical_estimate <- function(d, y, w, contrast, diagnose = FALSE) {
  # Both manuscript responses are continuous and every defined set has working weight one.
  f <- ugdm_fit(d, y, rep(1, length(y)), sampling_weights = w)
  wt <- w / mean(w)
  spread <- internal(".fit_logit_nnls")(-cbind(d$Xdir, d$Xdist), y, wt)
  position <- internal(".fit_logit_nnls")(cbind(d$Xpos, -d$Xpos), y, wt)
  pp <- (f$D2 + position$D2 - spread$D2) / 2
  env_names <- d$vars
  names_effect <- function(values, type) setNames(as.numeric(values), paste0(type, ".", env_names))
  beta <- matrix(f$beta[seq_len(3L * length(env_names))], ncol = 3L, byrow = TRUE)
  imp <- ugdm_importance(f, d, sampling_weights = w)
  v <- c(names_effect(rowSums(beta * contrast[, 1:3, drop = FALSE]), "spread"),
         names_effect(contrast[, 4L] * f$lambda, "position"),
         setNames(imp, paste0("importance.", names(imp))),
         shapley.spread = f$D2 - pp, shapley.position = pp,
         share.position = if (is.finite(f$D2) && f$D2 > 1e-10) pp / f$D2 else NA_real_,
         fit.D2 = f$D2)
  if (!f$converged) v[grepl("^(spread|position|importance)\\.", names(v))] <- NA_real_
  if (!f$converged || !spread$converged || !position$converged ||
      any(c(spread$D2, position$D2) > f$D2 + 1e-6, na.rm = TRUE)) {
    v[grepl("^(shapley|share)\\.", names(v))] <- NA_real_
  }
  diagnostics <- if (diagnose) {
    h <- info(d, f, w)
    data.frame(converged = f$converged, converged_spread = spread$converged,
      converged_position = position$converged, rank = h$rank, parameters = h$parameters,
      lambda_min = h$min_eigenvalue, condition = h$condition, active = h$active)
  }
  list(estimate = v, diagnostics = diagnostics, fit = f)
}

omission <- function(idx, w, n) {
  missed <- split(rep(seq_len(nrow(idx)), ncol(idx)), as.vector(idx))
  out <- lapply(seq_len(n), function(i) {
    keep <- rep(TRUE, nrow(idx))
    keep[missed[[as.character(i)]]] <- FALSE
    keep
  })
  list(keep = out, ess = vapply(out, function(i) kish(w[i]), numeric(1L)))
}
