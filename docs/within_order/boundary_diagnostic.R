source("docs/within_order/simulation_design.R")

cells <- data.frame(n = c(70L, 150L, 70L, 150L), order = c(2L, 4L, 8L, 12L),
                    response = c("incidence", "incidence", "relative", "relative"))
fixed <- within_reference()
results <- lapply(seq_len(nrow(cells)), function(j) {
  cell <- cells[j, ]
  do.call(rbind, lapply(seq_len(200L), function(i) {
    set.seed(14018L + 1871L * i + 45L * cell$n + cell$order)
    a <- within_data(within_sites(cell$n), cell$order, 700L, cell$response, fixed)
    fit <- within_fit(a)
    layout <- zetaGDUM:::.interaction_layout(a$design, fit$has_distance)
    active <- fit$coef[layout$dir] <= 1e-8
    data.frame(cell, replicate = i, converged = fit$converged,
      active_strong = sum(active[1:3]), active_weak = sum(active[4:6]),
      active_irrelevant = sum(active[7:9]))
  }))
})
within_write(do.call(rbind, results), "boundary_diagnostic.csv")
