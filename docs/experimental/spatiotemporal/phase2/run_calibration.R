root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase2_load.R"))
st2_load(root)
source(file.path(root, "phase2/factorial.R"))
out <- Sys.getenv("ST2_OUTPUT", file.path(root, "phase2/results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
replicates <- as.integer(Sys.getenv("ST2_NULL_REPS", "200"))
draws <- as.integer(Sys.getenv("ST2_BOOT_DRAWS", "40"))
if (is.na(replicates) || replicates < 2 || is.na(draws) || draws < 2) stop("Increase panel or bootstrap replicates.")

calibrate <- function(scenario, r) {
  x <- st2_factorial_world(54000L + r, scenario, region = TRUE)
  base <- x$enabled
  basis <- st2_basis(base, x$env)
  rows <- list()
  for (k in 1:2) {
    sets <- st2_draw_sets(60L, k, 60L)
    design <- st2_transition_design(base, x$env, k, sets = sets, basis = basis)
    fit <- st2_fit_transition(base, design)
    for (event in c("gain", "loss")) {
      f <- fit$fit[[event]]
      term <- "position_change.climate"
      if (f$status == "estimated") {
        u <- st2_transition_uncertainty(fit, event, x$blocks, draws,
                                        temporal_length = 3L, exact_delete = r <= 4L,
                                        terms = term)
        se <- u$se[1L, ]
        estimate <- f$coefficients[[term]]
        a <- f$coefficients[["intercept"]]
        change <- .5 * k
        contrast <- plogis(a + change * estimate) - plogis(a)
        contrast_se <- change * plogis(a) * (1 - plogis(a)) * se
      } else {
        se <- setNames(rep(NA_real_, 5L), c("naive", "delete_site", "delete_spatial_block", "moving_temporal", "crossed"))
        estimate <- contrast <- contrast_se <- NA_real_
      }
      rows[[length(rows) + 1L]] <- data.frame(scenario, replicate = r, spatial_order = k, event,
        status = f$status, rank = f$rank, method = names(se), estimate, se = unname(se),
        contrast, contrast_se = unname(contrast_se), exact_delete = r <= 4L, row.names = NULL)
    }
  }
  do.call(rbind, rows)
}
jobs <- expand.grid(scenario = c("null", "shock"), r = seq_len(replicates), stringsAsFactors = FALSE)
workers <- as.integer(Sys.getenv("ST2_WORKERS", "3"))
args <- seq_len(nrow(jobs))
work <- function(i) calibrate(jobs$scenario[i], jobs$r[i])
runs <- if (.Platform$OS.type == "windows" || workers == 1L) lapply(args, work) else
  parallel::mclapply(args, work, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
if (any(vapply(runs, inherits, logical(1), "try-error"))) stop("At least one null calibration worker failed.")
raw <- do.call(rbind, runs)
utils::write.csv(raw, file.path(out, "null_panels.csv"), row.names = FALSE, na = "")
groups <- split(raw, raw[c("scenario", "spatial_order", "event", "method")], drop = TRUE)
summary <- lapply(groups, function(z) {
  moment <- function(x, se) {
    good <- is.finite(x) & is.finite(se)
    value <- x[good]
    coverage <- mean(abs(value) <= 1.96 * se[good])
    c(n = sum(good), bias = mean(value), rmse = sqrt(mean(value^2)), sampling_sd = sd(value),
      mean_se = mean(se[good]), coverage = coverage, type1 = 1 - coverage,
      mc_se = sqrt(coverage * (1 - coverage) / sum(good)))
  }
  coef <- moment(z$estimate, z$se)
  contr <- moment(z$contrast, z$contrast_se)
  cbind(z[1L, c("scenario", "spatial_order", "event", "method")], as.list(coef),
        as.list(setNames(contr, paste0("contrast_", names(contr)))),
        boundary = mean(z$status == "boundary"), aliased = mean(z$status == "aliased"),
        nonconverged = mean(z$status == "nonconverged"))
})
utils::write.csv(do.call(rbind, summary), file.path(out, "null_summary.csv"), row.names = FALSE, na = "")
