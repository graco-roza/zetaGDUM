source("docs/within_order/simulation_design.R")
args <- commandArgs(TRUE)
B <- if (length(args)) as.integer(args[1L]) else 24L
set.seed(83107)
cells <- expand.grid(order = c(2L,4L,8L,12L), response = c("incidence","relative"),
  M = c(700L,2800L,11200L), stringsAsFactors = FALSE)
ref <- within_reference()
rows <- lapply(seq_len(nrow(cells)), function(i) {
  z <- cells[i, ]
  set.seed(1533 + z$order)
  site <- within_sites(70L)
  assessment <- within_data(within_sites(1700L), z$order, 5000L, z$response, ref)
  compute <- function(M) {
    a <- within_data(site, z$order, M, z$response, ref)
    f <- within_fit(a)
    v <- suppressWarnings(ugdm_within_order_importance(f, a$design, a$y, a$trials,
      evaluation = assessment$design, evaluation_y = assessment$y,
      evaluation_trials = assessment$trials, failure = "warn"))
    v$predictors$absolute
  }
  x <- t(vapply(seq_len(B), function(j) compute(z$M), numeric(3L)))
  cbind(z[rep(1,nrow(x)), ], replicate = seq_len(nrow(x)), as.data.frame(x))
})
within_write(do.call(rbind, rows), "integration_study.csv")
raw <- do.call(rbind, rows)
precision <- do.call(rbind, lapply(split(raw, raw[c("order","response","M")]), function(a) {
  b <- a[order(a$replicate), ]
  finite <- is.finite(b$strong) & is.finite(b$weak) & is.finite(b$irrelevant)
  split_at <- floor(nrow(b) / 2)
  z <- cbind(b[c("strong","weak","irrelevant")],
             strong_weak = b$strong - b$weak)
  sd_pool <- function(i) {
    if (sum(finite[i]) < 2L) return(NA_real_)
    max(apply(z[i[finite[i]], , drop = FALSE], 2L, sd))
  }
  data.frame(b[1L, c("order","response","M")], valid = sum(finite),
    pilot_se = sd_pool(seq_len(split_at)),
    confirmation_se = sd_pool(seq.int(split_at + 1L, nrow(b))))
}))
precision$pilot_selected <- FALSE
for (k in unique(precision$order)) for (response in unique(precision$response)) {
  can <- which(precision$order == k & precision$response == response &
                 is.finite(precision$pilot_se) & precision$pilot_se <= .01)
  if (length(can)) precision$pilot_selected[can[which.min(precision$M[can])]] <- TRUE
}
within_write(precision, "integration_budget.csv")
