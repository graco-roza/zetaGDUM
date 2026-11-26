source("docs/phase4/common.R")
# Single almost ubiquitous but very low abundance species. Site totals are identical.
n <- 250L
s <- support_case(n)
s$comm[229:n, 1L] <- 0
s$comm[229:n, -1] <- s$comm[229:n, -1] * (1 + 1e-7)
k <- c(3, 4, 5, 13, 30, 50, 100, 150, 200, n)
h <- zeta_order_support(s$comm, k)
phase4_write(cbind(h, exact_relative_response_sd =
  (1 / (1 + 1e7)) * sqrt(h$richness_variance) * ifelse(k >= 5, 1, NA_real_)),
  "retention", "rare_core")
